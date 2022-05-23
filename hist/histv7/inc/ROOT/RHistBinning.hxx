/// \file ROOT/RHistBinning.hxx
/// \ingroup Hist ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2020-01-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RHistBinning
#define ROOT7_RHistBinning

#include "ROOT/RAxis.hxx"

namespace ROOT {
namespace Experimental {

/**
 \class RHistBinningBase
 Histogram axis configuration, to convert a n-dimensional coordinate to a global bin index.
 Base class that hides the actual axis types.
*/
template <int DIMENSIONS>
class RHistBinningBase {
public:

   /**
    \class const_iterator
    Random const_iterator through bins. Represents the bin index, not a bin
    content: the axis has no notion of any content.
    */
   class const_iterator : public std::iterator<std::random_access_iterator_tag, int /*value*/, int /*distance*/,
                                               const int * /*pointer*/, const int & /*ref*/> {
      int fCursor = 0; ///< Current iteration position
   };

   /// Get the product of all axes' bin numbers.
   int GetNBins() const
   {
      using namespace std;
      // Product of all axes' bin numbers.
      return accumulate(fAxes.begin(), fAxes.end(), 1, [](int old, const RAxisBase &ax) { return old * ax.GetNBins(); })
   }

   const std::array<unique_ptr<RAxisBase>, DIMENSIONS> &GetAxes() const { return fAxes; }

   const_iterator begin() const;
   const_iterator end() const;
   const_iterator cbegin() const { return begin(); }
   const_iterator cend() const { return begin(); }

protected:
   ///\name Inaccessible copy, assignment
   /// The copy and move constructors and assignment operators are protected to
   /// prevent slicing.
   ///\{

   /// Default construct a RAxisBase (for use by derived classes for I/O)
   RHistBinningBase() = default;
   RHistBinningBase(const RAxisBase &) = default;
   RHistBinningBase(RHistBinningBase &&) = default;
   RHistBinningBase &operator=(const RHistBinningBase &) = default;
   RHistBinningBase &operator=(RHistBinningBase &&) = default;

   ///\}
private:
   std::array<unique_ptr<RAxisBase>, DIMENSIONS> fAxes;
};

} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RHistBinning
