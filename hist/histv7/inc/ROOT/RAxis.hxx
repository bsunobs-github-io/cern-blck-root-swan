/// \file ROOT/RAxis.hxx
/// \ingroup HistV7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-03-23
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2022, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAxis
#define ROOT7_RAxis

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

#include "ROOT/RAxisIterator.hxx"
#include "ROOT/RLogger.hxx"

namespace ROOT {
namespace Experimental {

namespace Internal {
namespace Hist {
/// The axis argument type for coordinates; corresponds to the alternatives of CoordArgVariant_t.
enum class ECoord {
   kDouble, ///< The axis takes `double` as coordinate.
   kFloat,  ///< The axis takes `float` as coordinate.
   kInt64,  ///< The axis takes `int64_t` as coordinate.
   kUInt64, ///< The axis takes `uint64_t` as coordinate.
   kLabel,  ///< The axis takes `std::string` as coordinate.

   kNumAxisCoords ///< Marker for determining the last element.
};

/// Find the value of ECoord that matches CoordType.
template <class CoordType, bool IsInt = std::is_integral_v<CoordType>>
static constexpr int kCoordArgVariant = -1;

template <>
constexpr int kCoordArgVariant<double> = static_cast<int>(ECoord::kDouble);

template <>
constexpr int kCoordArgVariant<float> = static_cast<int>(ECoord::kFloat);

template <class CoordType>
constexpr int
   kCoordArgVariant<CoordType, true> = static_cast<int>(std::is_signed_v<CoordType> ? ECoord::kInt64
                                                                                    : ECoord::kUInt64);

template <>
constexpr int kCoordArgVariant<std::string> = static_cast<int>(ECoord::kLabel);

template <>
constexpr int kCoordArgVariant<const char *> = static_cast<int>(ECoord::kLabel);
} // namespace Hist
} // namespace Internal

/**
 \class RAxisBase
 Histogram axis base class. Keeps track of the number of bins and overflow
 handling. Offers bin iteration.
 */
class RAxisBase {
public:
   using const_iterator = RAxisIterator;

   /// const_iterator representing an invalid bin.
   static constexpr const_iterator kInvalidBin{std::numeric_limits<const_iterator::value_type>::min()};

   /// Underlying axis kind
   enum class EKind {
      kEquidistant, ///< The axis consists of bins with computed size
      kIrregular,   ///< The axis contains bins with consecutive, user-provided bin boorders.

      kNumAxisKinds ///< Marker for determining the last element.
   };

   /// Handling of entries outside the (for growing axis: initial) axis range.
   enum class EOverflow {
      kNone,          ///< The axis does not account for under-/overflow and cannot groow
      kWrap,          ///< The axis has no over-/underflow and instead wraps around (modulo axis size).
      kGrow,          ///< The axis extends the range by adding additional bins outside the initial range.
      kMerge,         ///< The axis extends the range by merging neighboring bins.
      kUnder = 0b100, ///< The axis counts entries below the axis range; test bit-wise.
      kOver = 0b1000, ///< The axis counts entries above the axis range; test bit-wise.

      kNumAxisOverflows, ///< Marker for determining the last element.

      kBoth = kUnder | kOver, ///< The axis has both under- and overflow.
   };

   /** The bin coordinate is given as one of these to find the bin. This needs
    to match the axis's coordinate type; FindBin() will throw otherwise.
    */
   using CoordArgVariant_t = std::variant<double, float, int64_t, uint64_t, std::string>;

protected:
   ///\name Inaccessible copy, assignment
   /// The copy and move constructors and assignment operators are protected to
   /// prevent slicing.
   ///\{
   RAxisBase(const RAxisBase &) = default;
   RAxisBase(RAxisBase &&) = default;
   RAxisBase &operator=(const RAxisBase &) = default;
   RAxisBase &operator=(RAxisBase &&) = default;
   ///\}

   /// Default construct a RAxisBase (for use by derived classes for I/O)
   RAxisBase() noexcept(noexcept(std::string())) = default;

   /// Virtual destructor needed in this inheritance-based design
   virtual ~RAxisBase() = default;

   /// Construct a RAxisBase.
   ///
   ///\param[in] title - axis title used for graphics and text representation.
   ///\param[in] overflowHandling - how to treat coordinates beyond the axis range.
   RAxisBase(std::string_view title, EKind kind, Internal::Hist::ECoord coordType, EOverflow overflowHandling) noexcept
      : fTitle(title), fKind(kind), fCoordType(coordType), fOverflowHandling(overflowHandling)
   {
   }

   /// Check if two axis have the same bin borders
   ///
   /// Default implementation should work for any RAxis type, but is quite
   /// inefficient as it does virtual GetBinFrom calls in a loop. RAxis
   /// implementations are encouraged to provide optimized overrides for common
   /// axis binning comparison scenarios.
   virtual bool HasSameBinBordersAs(const RAxisBase &other) const
   {
      // Axis underlying value type and overflow handling must match
      if (fCoordType != other.fCoordType || fOverflowHandling != other.fOverflowHandling)
         return false;

      // Number of normal bins must match
      if (GetNBinsNoOver() != other.GetNBinsNoOver())
         return false;

      // Left borders of normal bins must match
      for (auto i = begin(); i != end(); ++i)
         if (GetBinFromV(i) != other.GetBinFromV(i))
            return false;

      // Right border of the last normal bin (aka maximum) must also match
      if (GetMaximumV() != other.GetMaximumV())
         return false;

      // If all of these checks passed, the two axes have the same bin borders
      return true;
   }

public:
   /// Get the axis's title
   const std::string &GetTitle() const { return fTitle; }

   /// Whether this axis can grow (and thus has no overflow bins).
   bool CanGrow() const noexcept { return fOverflowHandling == EOverflow::kGrow; };

   /// Whether this axis wraps its coordinates (and thus has no overflow bins).
   /// Wrapping means that values above (below) the maximum are reduced (increased)
   /// by multiples of the axis length, until the value fits into the axis. Example:
   ///
   ///   x = 2.3 * pi
   /// For a wrapping axis -pi .. pi, the bin will be filled at 0.3 * pi.
   bool DoesWrap() const noexcept { return fOverflowHandling == EOverflow::kWrap; };

   /// Whether the axis accounts for entries below the axis range.
   bool HasUnderflow() const noexcept
   {
      return static_cast<int>(fOverflowHandling) & static_cast<int>(EOverflow::kUnder);
   }

   /// Whether the axis accounts for entries above the axis range.
   bool HasOverflow() const noexcept
   {
      return static_cast<int>(fOverflowHandling) & static_cast<int>(EOverflow::kOver);
   }

   /// Get the number of bins, excluding under- and overflow.
   virtual int GetNBinsNoOver() const noexcept = 0;

   /// Get the number of bins, including under- and overflow.
   int GetNBins() const noexcept { return GetNBinsNoOver() + GetNOverflowBins(); }

   /// Get the number of over- and underflow bins.
   int GetNOverflowBins() const noexcept
   {
      int ret = 0;
      if (HasUnderflow())
         ++ret;
      if (HasUnderflow())
         ++ret;
      return ret;
   };

   ///\name Iterator interfaces
   ///\{

   /// Get a const_iterator pointing to the first regular bin.
   const_iterator begin() const noexcept { return const_iterator{1}; }

   /// Get a const_iterator pointing past the last regular bin.
   const_iterator end() const noexcept { return const_iterator{GetNBinsNoOver() + 1}; }

   /// Get a const_iterator pointing to the underflow bin if it exists, or the first bin.
   const_iterator begin_underflow() const
   {
      if (HasUnderflow())
         return const_iterator{0};
      return begin();
   }

   /// Get a const_iterator pointing past the overflow bin if it exists, or the end().
   const_iterator end_overflow() const
   {
      if (HasOverflow())
         return const_iterator{GetNBinsNoOver() + 2};
      return end();
   }
   ///\}

   /// Find the bin (returning `begin_underflow()` for underflow and `end_overflow()`
   /// for overflow) for the given coordinate.
   /// \note Passing the coordinate of a bin border can either return the bin above or
   /// below the bin border. I.e. don't do that for reliable results!
   virtual const_iterator FindBin(const CoordArgVariant_t &x) const noexcept = 0;

   ///\name Getters for bin coordinates.
   ///\{

   /// Get the bin center for the given bin index.
   /// For non-floating point axis coordinates, this is the bin's coordinate value.
   /// The result of this method on an overflow or underflow bin is unspecified.
   virtual double GetBinCenter(const_iterator bin) const = 0;

   /// Get the low bin border ("left edge") for the given bin index.
   /// For non-floating point axis coordinates, this is the bin's coordinate value.
   /// The result of this method on an underflow bin is unspecified.
   virtual CoordArgVariant_t GetBinFromV(const_iterator bin) const = 0;

   /// Get the high bin border ("right edge") for the given bin index.
   /// For non-floating point axis coordinates, this is the coordinate value of the subsequent bin.
   /// The result of this method on an overflow bin is unspecified.
   CoordArgVariant_t GetBinToV(const_iterator bin) const
   {
      if (bin < begin())
         return GetMinimumV();
      return GetBinFromV(bin + 1);
   }

   /// Get the low end of the axis range.
   CoordArgVariant_t GetMinimumV() const { return GetBinFromV(begin()); }

   /// Get the high end of the axis range.
   CoordArgVariant_t GetMaximumV() const { return GetBinToV(end() - 1); }

   /// If the coordinate `x` is within 10 ULPs of a bin low edge coordinate,
   /// return the bin for which this is a low edge. If it's not a bin edge,
   /// return `kInvalidBin`.
   virtual const_iterator GetBinIndexForLowEdge(double x) const;

   ///\}

   /// Check if two axes use the same binning convention, i.e.
   ///
   /// - Either they are both growable or neither of them is growable.
   /// - Minimum, maximum, and all bin borders in the middle are the same.
   /// - Bin labels must match (exactly including order, for now).
   bool HasSameBinningAs(const RAxisBase &other) const;

private:
   std::string fTitle; ///< Title of this axis, used for graphics / text.

   EKind fKind : 2; ///< Derived axis type.
   static_assert(1 >> 2 < static_cast<int>(EKind::kNumAxisKinds), "Too few bits for fKind.");

   Internal::Hist::ECoord fCoordType : 3; ///< Type of the axis coordinates.
   static_assert(1 >> 3 < static_cast<int>(Internal::Hist::ECoord::kNumAxisCoords), "Too few bits for fCoordType.");

   EOverflow fOverflowHandling : 4; ///< Overflow handling of the axis.
   static_assert(1 >> 4 < static_cast<int>(EOverflow::kNumAxisOverflows), "Too few bits for fOverflowHandling.");
};

namespace Internal {
template <class CoordType, bool IsFloat>
struct RAxisEquidistantBinCalc {
   using const_iterator = RAxisBase::const_iterator;

   CoordType fLow;         ///< The lower limit of the axis
   CoordType fInvBinWidth; ///< The inverse of the bin width
   int fNBinsNoOver;       ///< Number of bins excluding under- and overflow.

   RAxisEquidistantBinCalc(CoordType low, CoordType high, int nBinsNoOver)
      : fLow(low), fInvBinWidth(GetInvBinWidth(nBinsNoOver, low, high)), fNBinsNoOver(nBinsNoOver)
   {
   }

   /// Determine the inverse bin width.
   /// \param nbinsNoOver - number of bins without unter-/overflow
   /// \param low - lower axis boundary
   /// \param high - upper axis boundary
   static CoordType GetInvBinWidth(int nbinsNoOver, CoordType low, CoordType high)
   {
      if (high - low == 0)
         return std::numeric_limits<CoordType>::min();
      return nbinsNoOver / (high - low);
   }

   const_iterator begin() const { return const_iterator{1}; }

   /// Find the bin index for the given coordinate, ignoring potential overflows.
   /// \note Passing a bin border coordinate can either return the bin above or
   /// below the bin border. I.e. don't do that for reliable results!
   const_iterator::value_type FindBinRaw(CoordType x) const noexcept
   {
      return 1 + const_iterator::value_type(x - fLow) * fInvBinWidth;
   }

   /// Get the bin center for the given bin index.
   /// For the bin == 1 (the first bin) of 2 bins for an axis (0., 1.), this
   /// returns 0.25.
   /// The result of this method on an overflow or underflow bin is unspecified.
   double GetBinCenter(const_iterator bin) const { return fLow + (bin - begin() + 0.5) / fInvBinWidth; }

   /// Get the low bin border for the given bin index.
   /// For the bin == 1 (the first bin) of 2 bins for an axis (0., 1.), this
   /// returns 0.
   /// The result of this method on an overflow bin is unspecified.
   CoordType GetBinFrom(const_iterator bin) const
   {
      if (*bin > fNBinsNoOver) {
         // handle overflow bin gracefully
         return GetTo();
      }
      return fLow + (bin - begin()) / fInvBinWidth;
   }

   CoordType GetBinWidth() const noexcept { return 1. / fInvBinWidth; }

   /// Upper end of the axis.
   CoordType GetTo() const noexcept { return fLow + GetLength(); }

   /// Lower end of the axis.
   CoordType GetFrom() const noexcept { return fLow; }

   /// Length of the axis.
   CoordType GetLength() const noexcept { return fNBinsNoOver / fInvBinWidth; }
};

template <class CoordType>
struct RAxisEquidistantBinCalc<CoordType, /* IsFloat */ false> {
   using const_iterator = RAxisBase::const_iterator;

   static constexpr bool kIsSigned = std::is_signed_v<CoordType>;
   using SameSignedInt64_t = std::conditional_t<kIsSigned, int64_t, uint64_t>;

   SameSignedInt64_t fLow;    ///< The lower limit of the axis
   SameSignedInt64_t fLength; ///< The difference between higher and lower limit.
   int fNBinsNoOver;          ///< The number of bins.

   RAxisEquidistantBinCalc(CoordType low, CoordType high, int nBinsNoOver) noexcept
      : fLow(low), fLength(high - low), fNBinsNoOver(nBinsNoOver)
   {
   }

   const_iterator begin() const noexcept { return const_iterator{1}; }

   /// Find the bin index for the given coordinate, ignoring potential overflows.
   const_iterator::value_type FindBinRaw(CoordType x) const noexcept
   {
      return 1 + const_iterator::value_type{((x - fLow) * fNBinsNoOver) / fLength};
   }

   /// Get the bin center for the given bin index.
   /// For the bin == 1 (the first bin) of 2 bins for an axis (0., 1.), this
   /// returns 0.25.
   /// The result of this method on an overflow or underflow bin is unspecified.
   double GetBinCenter(const_iterator bin) const noexcept
   {
      return fLow + (bin - begin() + 0.5) * fLength / fNBinsNoOver;
   }

   /// Get the low bin border for the given bin index.
   /// For the bin == 1 (the first bin) of 2 bins for an axis (0., 1.), this
   /// returns 0.
   /// The result of this method on an overflow bin is unspecified.
   CoordType GetBinFrom(const_iterator bin) const noexcept
   {
      if (*bin > fNBinsNoOver) {
         // handle overflow bin gracefully
         return fLow + fLength;
      }
      return fLow + (bin - begin()) * fLength / fNBinsNoOver;
   }

   CoordType GetBinWidth() const noexcept { return fLength / fNBinsNoOver; }

   /// Upper end of the axis.
   CoordType GetTo() const noexcept { return fLow + GetLength(); }

   /// Lower end of the axis.
   CoordType GetFrom() const noexcept { return fLow; }

   /// Length of the axis.
   CoordType GetLength() const noexcept { return fLength; }
};

} // namespace Internal

/**
 Axis with equidistant bin borders. Defined by lower `L` and upper `U` limit and
 the number of bins `N`. All bins have the same width `(U-L)/N`.
 */
template <class CoordType>
class RAxisEquidistant : public RAxisBase {
public:
   using CoordType_t = CoordType;

   static_assert(std::is_same_v<CoordType, std::string> && "Please use RAxisLabels instead.");

protected:
   static constexpr bool kIsFloat = std::is_floating_point_v<CoordType>;
   using BinCalc_t = Internal::RAxisEquidistantBinCalc<CoordType_t, kIsFloat>;

   BinCalc_t fBinCalc;

   /// See RAxisBase::HasSameBinBordersAs
   bool HasSameBinBordersAs(const RAxisBase &other) const override;

   /// Given rawbin (`<1` for underflow, `>GetNBinsNoOver()` for overflow),
   /// determine the bin number taking into account how over/underflow
   /// should be handled.
   ///
   /// \param[in] rawbin for which to determine the bin number.
   /// \return Returns the bin number adjusted for potential over- and underflow
   /// bins. Returns `kInvalidBin` if the axis cannot handle the over- / underflow.
   ///
   const_iterator ClampBinIndex(const_iterator::value_type rawbin) const
   {
      // iterator is 1-based:
      const_iterator ret{rawbin + 1};

      // FIXME: Handle kGrow!

      // Underflow: Put in underflow bin if any, otherwise ignore
      if (ret < begin())
         return begin_underflow();

      // Overflow: Put in overflow bin if any, otherwise ignore
      // `rawbin` is not an integer, cannot compare `rawbin > GetLastBin()`.
      if (ret >= end())
         return end_overflow();

      // Bin index is in range and has been corrected for over/underflow
      return ret;
   }

public:
   RAxisEquidistant() = default;

   /// Initialize a RAxisEquidistant.
   /// \param[in] title - axis title used for graphics and text representation.
   /// \param nbinsNoOver - number of bins in the axis, excluding under- and overflow
   ///   bins.
   /// \param low - the low axis range. Any coordinate below that is considered
   ///   as underflow. The first bin's lower edge is at this value.
   /// \param high - the high axis range. Any coordinate above that is considered
   ///   as overflow. The last bin's higher edge is at this value.
   explicit RAxisEquidistant(std::string_view title, int nbinsNoOver, CoordType low, CoordType high, EOverflow overflow = EOverflow::kBoth) noexcept
      : RAxisBase(title, EKind::kEquidistant, Internal::Hist::kCoordArgVariant<CoordType_t>, overflow), fBinCalc(low, high, nbinsNoOver)
   {
   }

   /// Initialize a RAxisEquidistant.
   /// \param nbinsNoOver - number of bins in the axis, excluding under- and overflow
   ///   bins.
   /// \param low - the low axis range. Any coordinate below that is considered
   ///   as underflow. The first bin's lower edge is at this value.
   /// \param high - the high axis range. Any coordinate above that is considered
   ///   as overflow. The last bin's higher edge is at this value.
   explicit RAxisEquidistant(int nbinsNoOver, CoordType low, CoordType high) noexcept
      : RAxisEquidistant("", nbinsNoOver, low, high)
   {
   }

   /// Get the number of bins, excluding under- and overflow.
   int GetNBinsNoOver() const noexcept final { return fBinCalc.fNBinsNoOver; }

   /// Find the adjusted bin index (returning `kUnderflowBin` for underflow and
   /// `GetOverflowBin()` for overflow) for the given coordinate.
   /// \note Passing a bin border coordinate can either return the bin above or
   /// below the bin border. I.e. don't do that for reliable results!
   const_iterator FindBin(const CoordArgVariant_t &x) const noexcept final
   {
      auto coord = std::get<Internal::Hist::kCoordArgVariant<CoordType>>(x);
      if (DoesWrap())
         coord = std::fmod(coord, fBinCalc.GetLength());

      const_iterator::value_type unclamped = FindBinRaw(coord);
      return ClampBinIndex(unclamped);
   }

   /// Get the width of the bins.
   CoordType_t GetBinWidth() const noexcept { return fBinCalc.GetBinWidth(); }

   /// Get the bin center for the given bin index.
   /// For the bin == 1 (the first bin) of 2 bins for an axis (0., 1.), this
   /// returns 0.25.
   /// The result of this method on an overflow or underflow bin is unspecified.
   double GetBinCenter(const_iterator bin) const final { return fBinCalc.GetBinCenter(bin); }

   /// Get the low bin border for the given bin index.
   /// For the bin == 1 (the first bin) of 2 bins for an axis (0., 1.), this
   /// returns 0.
   /// The result of this method on an underflow bin is unspecified.
   CoordType_t GetBinFrom(const_iterator bin) const { return fBinCalc.GetBinFrom(bin); }

   /// Get the low bin border for the given bin index.
   /// For the bin == 1 (the first bin) of 2 bins for an axis (0., 1.), this
   /// returns 0.
   /// The result of this method on an underflow bin is unspecified.
   CoordArgVariant_t GetBinFromV(const_iterator bin) const final { return GetBinFrom(bin); }

   /// If the coordinate `x` is within 10 ULPs of a bin low edge coordinate,
   /// return the bin for which this is a low edge. If it's not a bin edge,
   /// return `kInvalidBin`.
   const_iterator GetBinIndexForLowEdge(double x) const final;
};


/**
 An axis with non-equidistant bins (also known as "variable binning"). It is
 defined by an array of bin borders - one more than the number of
 (non-overflow-) bins it has! As an example, an axis with two bin needs three
 bin borders:
   - lower edge of the first bin;
   - higher edge of the first bin, identical to the lower edge of the second
     bin;
   - higher edge of the second bin

 This axis cannot grow; the size of new bins would not be well defined.
 */
template <class CoordType>
class RAxisIrregular : public RAxisBase {
public:

   using CoordType_t = CoordType;
   static_assert(std::is_same_v<CoordType, std::string> && "Please use RAxisLabels instead.");

private:
   /// Bin borders, one more than the number of regular bins.
   std::vector<double> fBinBorders;

protected:
   /// See RAxisBase::HasSameBinBordersAs
   bool HasSameBinBordersAs(const RAxisBase &other) const override;

   /// Find the raw bin index (not adjusted) for the given coordinate `x`.
   /// The resulting raw bin is 1-based.
   /// \note Passing a bin border coordinate can either return the bin above or
   /// below the bin border. I.e. don't do that for reliable results!
   const_iterator FindBinRaw(const CoordArgVariant_t &x) const noexcept
   {
      const auto bBegin = fBinBorders.begin();
      const auto bEnd = fBinBorders.end();
      // lower_bound finds the first bin border that is >= x.
      auto iNotLess = std::lower_bound(bBegin, bEnd, x);
      return iNotLess - bBegin;
   }

public:
   RAxisIrregular() = default;

   /// Construct a RAxisIrregular from a vector of bin borders.
   /// \note The bin borders must be sorted in increasing order!
   explicit RAxisIrregular(const std::vector<double> &binborders) : RAxisBase(), fBinBorders(binborders)
   {
#ifdef R__DO_RANGE_CHECKS
      if (!std::is_sorted(fBinBorders.begin(), fBinBorders.end()))
         R__LOG_ERROR("HIST") << "Bin borders must be sorted!";
#endif // R__DO_RANGE_CHECKS
   }

   /// Construct a RAxisIrregular from a vector of bin borders.
   /// \note The bin borders must be sorted in increasing order!
   /// Faster, noexcept version taking an rvalue of binborders. The compiler will
   /// know when it can take this one.
   explicit RAxisIrregular(std::vector<double> &&binborders) noexcept : RAxisBase(), fBinBorders(std::move(binborders))
   {
#ifdef R__DO_RANGE_CHECKS
      if (!std::is_sorted(fBinBorders.begin(), fBinBorders.end()))
         R__LOG_ERROR("HIST") << "Bin borders must be sorted!";
#endif // R__DO_RANGE_CHECKS
   }

   /// Construct a RAxisIrregular from a vector of bin borders.
   /// \note The bin borders must be sorted in increasing order!
   explicit RAxisIrregular(std::string_view title, const std::vector<double> &binborders)
      : RAxisBase(title), fBinBorders(binborders)
   {
#ifdef R__DO_RANGE_CHECKS
      if (!std::is_sorted(fBinBorders.begin(), fBinBorders.end()))
         R__LOG_ERROR("HIST") << "Bin borders must be sorted!";
#endif // R__DO_RANGE_CHECKS
   }

   /// Construct a RAxisIrregular from a vector of bin borders.
   /// \note The bin borders must be sorted in increasing order!
   /// Faster, noexcept version taking an rvalue of binborders. The compiler will
   /// know when it can take this one.
   explicit RAxisIrregular(std::string_view title, std::vector<double> &&binborders) noexcept
      : RAxisBase(title), fBinBorders(std::move(binborders))
   {
#ifdef R__DO_RANGE_CHECKS
      if (!std::is_sorted(fBinBorders.begin(), fBinBorders.end()))
         R__LOG_ERROR("HIST") << "Bin borders must be sorted!";
#endif // R__DO_RANGE_CHECKS
   }

   /// Convert to RAxisConfig.
   operator RAxisConfig() const
   {
      return RAxisConfig(GetTitle(), GetBinBorders());
   }

   /// Get the number of bins, excluding under- and overflow.
   int GetNBinsNoOver() const noexcept final
   {
      return fBinBorders.size() - 1;
   }

   /// Find the bin index (adjusted with under- and overflow) for the given coordinate `x`.
   /// \note Passing a bin border coordinate can either return the bin above or
   /// below the bin border. I.e. don't do that for reliable results!
   int FindBin(double x) const noexcept final
   {
      int rawbin = FindBinRaw(x);
      // No need for ClampBinNumber(rawbin) here; lower_bound() is the
      // answer: e.g. for x < *bBegin, rawbin is -1.
      if (rawbin < GetFirstBin())
         return kUnderflowBin;
      if (rawbin >= GetLastBin() + 1)
         return GetOverflowBin();
      return rawbin;
   }

   /// Get the bin center of the bin with the given index.
   /// The result of this method on an overflow or underflow bin is unspecified.
   double GetBinCenter(int bin) const final
   {
      return 0.5 * (fBinBorders[bin - 1] + fBinBorders[bin]);
   }

   /// Get the lower bin border for a given bin index.
   /// The result of this method on an underflow bin is unspecified.
   double GetBinFrom(int bin) const final
   {
      if (bin == GetOverflowBin())
         return fBinBorders[GetLastBin()];
      return fBinBorders[bin - 1];
   }

   /// If the coordinate `x` is within 10 ULPs of a bin low edge coordinate,
   /// return the bin for which this is a low edge. If it's not a bin edge,
   /// return `kInvalidBin`.
   int GetBinIndexForLowEdge(double x) const noexcept final;

   /// This axis cannot be extended.
   bool CanGrow() const noexcept final
   {
      return false;
   }

   /// Access to the bin borders used by this axis.
   const std::vector<double> &GetBinBorders() const noexcept
   {
      return fBinBorders;
   }
};


namespace Internal {

template <>
struct AxisConfigToType<RAxisConfig::kIrregular> {
   using Axis_t = RAxisIrregular;

   Axis_t operator()(const RAxisConfig &cfg) { return RAxisIrregular(cfg.GetTitle(), cfg.GetBinBorders()); }
};

} // namespace Internal

/**
 \class RAxisLabels
 A RAxisGrow that has a label assigned to each bin and a bin width of 1.

 While filling still works through coordinates (i.e. arrays of doubles),
 RAxisLabels allows to convert a string to a bin number or the bin's coordinate
 center. The number of labels and the number of bins reported by RAxisGrow might
 differ: the RAxisGrow will only grow when seeing a Fill(), while the RAxisLabels
 will add a new label whenever `GetBinCenter()` is called.

 Implementation details:
 Filling happens often; `GetBinCenter()` needs to be fast. Thus the unordered_map.
 The painter needs the reverse: it wants the label for bin 0, bin 1 etc. The axis
 should only store the bin labels once; referencing them is (due to re-allocation,
 hashing etc) non-trivial. So instead, build a `vector<string_view>` for the few
 times the axis needs to be painted.
 */
class RAxisLabels : public RAxisGrow {
private:
   /// Map of label (view on `fLabels`'s elements) to bin index
   std::unordered_map<std::string, int /*bin number*/> fLabelsIndex;

public:
   /// Construct a RAxisLables from a `vector` of `string_view`s, with title.
   explicit RAxisLabels(std::string_view title, const std::vector<std::string_view> &labels)
      : RAxisGrow(title, labels.size(), 0., static_cast<double>(labels.size()))
   {
      for (size_t i = 0, n = labels.size(); i < n; ++i)
         fLabelsIndex[std::string(labels[i])] = i;
   }

   /// Construct a RAxisLables from a `vector` of `string`s, with title.
   explicit RAxisLabels(std::string_view title, const std::vector<std::string> &labels)
      : RAxisGrow(title, labels.size(), 0., static_cast<double>(labels.size()))
   {
      for (size_t i = 0, n = labels.size(); i < n; ++i)
         fLabelsIndex[labels[i]] = i;
   }

   /// Construct a RAxisLables from a `vector` of `string_view`s
   explicit RAxisLabels(const std::vector<std::string_view> &labels) : RAxisLabels("", labels) {}

   /// Construct a RAxisLables from a `vector` of `string`s
   explicit RAxisLabels(const std::vector<std::string> &labels) : RAxisLabels("", labels) {}

   /// Convert to RAxisConfig.
   operator RAxisConfig() const { return RAxisConfig(GetTitle(), GetBinLabels()); }

   /// Get the bin index with label.
   int FindBinByName(const std::string &label)
   {
      auto insertResult = fLabelsIndex.insert({label, -1});
      if (insertResult.second) {
         // we have created a new label
         int idx = fLabelsIndex.size() - 1;
         insertResult.first->second = idx;
         return idx;
      }
      return insertResult.first->second;
   }

   /// Get the center of the bin with label.
   double GetBinCenterByName(const std::string &label)
   {
      return FindBinByName(label) + 0.5; // bin *center*
   }

   /// Build a vector of labels. The position in the vector defines the label's bin.
   std::vector<std::string_view> GetBinLabels() const
   {
      std::vector<std::string_view> vec(fLabelsIndex.size());
      for (const auto &kv : fLabelsIndex)
         vec.at(kv.second) = kv.first;
      return vec;
   }

   /// Result of an RAxisLabels label set comparison
   enum LabelsCmpFlags {
      /// Both axes have the same labels, mapping to the same bins
      kLabelsCmpSame = 0,

      /// The other axis doesn't have some labels from this axis
      kLabelsCmpSubset = 0b1,

      /// The other axis has some labels which this axis doesn't have
      kLabelsCmpSuperset = 0b10,

      /// The labels shared by both axes do not map into the same bins
      kLabelsCmpDisordered = 0b100,
   };

   /// Compare the labels of this axis with those of another axis
   LabelsCmpFlags CompareBinLabels(const RAxisLabels &other) const noexcept
   {
      // This will eventually contain the results of the labels comparison
      LabelsCmpFlags result = kLabelsCmpSame;
      size_t missing_in_other = 0;

      // First, check how this axis' labels map into the other axis
      for (const auto &kv : fLabelsIndex) {
         auto iter = other.fLabelsIndex.find(kv.first);
         if (iter == other.fLabelsIndex.cend()) {
            ++missing_in_other;
         } else if (iter->second != kv.second) {
            result = LabelsCmpFlags(result | kLabelsCmpDisordered);
         }
      }
      if (missing_in_other > 0)
         result = LabelsCmpFlags(result | kLabelsCmpSubset);

      // If this covered all labels in the other axis, we're done
      if (fLabelsIndex.size() == other.fLabelsIndex.size() + missing_in_other)
         return result;

      // Otherwise, we must check the labels of the other axis too
      for (const auto &kv : other.fLabelsIndex)
         if (fLabelsIndex.find(kv.first) == fLabelsIndex.cend())
            return LabelsCmpFlags(result | kLabelsCmpSuperset);
      return result;
   }
};

namespace Internal {

template <>
struct AxisConfigToType<RAxisConfig::kLabels> {
   using Axis_t = RAxisLabels;

   Axis_t operator()(const RAxisConfig &cfg) { return RAxisLabels(cfg.GetTitle(), cfg.GetBinLabels()); }
};

} // namespace Internal

///\name Axis Compatibility
///\{
enum class EAxisCompatibility {
   kIdentical, ///< Source and target axes are identical

   kContains, ///< The source is a subset of bins of the target axis

   /// The bins of the source axis have finer granularity, but the bin borders
   /// are compatible. Example:
   /// source: 0., 1., 2., 3., 4., 5., 6.; target: 0., 2., 5., 6.
   /// Note that this is *not* a symmetrical property: only one of
   /// CanMerge(source, target), CanMap(target, source) can return kContains.
   kSampling,

   /// The source axis and target axis have different binning. Example:
   /// source: 0., 1., 2., 3., 4., target: 0., 0.1, 0.2, 0.3, 0.4
   kIncompatible
};

/// Whether (and how) the source axis can be merged into the target axis.
EAxisCompatibility CanMap(const RAxisEquidistant &target, const RAxisEquidistant &source) noexcept;
///\}

} // namespace Hist
} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RAxis header guard
