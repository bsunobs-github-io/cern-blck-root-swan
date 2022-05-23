/// \file ROOT/RAxisIterator.hxx
/// \ingroup HistV7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2012-05-17
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2022, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAxisIterator
#define ROOT7_RAxisIterator

#include <iterator>

namespace ROOT {
namespace Experimental {

/**
 \class RAxisIterator
 Random-access const_iterator through bins. Represents the bin index of an axis,
 not a histogram bin or a bin content.

 The internal counter starts at `1`; for the last bin its value is `nbins`. The
 iterator can also point to the underflow bin in which case its value is `0`, or
 at the overflow bin in which case its value is `nbins + 1`.

 The `end()` for non-overflow is thus `nbins + 1`, for overflow it's `nbins + 2`.

 Dereferncing the iterator returns the underlying counter.
 */

class RAxisIterator {
public:
   // iterator protocol:
   using iterator_category = std::random_access_iterator_tag;
   using value_type = int;
   using difference_type = int;
   using pointer = const int *;
   using reference = const int &;

private:
   value_type fCursor = 0; ///< Current iteration position

public:
   RAxisIterator() = default;

   /// Initialize a RAxisIterator with its position
   constexpr explicit RAxisIterator(value_type cursor) noexcept : fCursor(cursor) {}

   /// ++i
   constexpr RAxisIterator &operator++() noexcept
   {
      // Could check whether fCursor < fEnd - but what for?
      ++fCursor;
      return *this;
   }

   /// --i
   constexpr RAxisIterator &operator--() noexcept
   {
      // Could check whether fCursor > fBegin - but what for?
      --fCursor;
      return *this;
   }

   /// i++
   constexpr RAxisIterator operator++(difference_type) noexcept
   {
      RAxisIterator old(*this);
      ++(*this);
      return old;
   }

   // i--
   constexpr RAxisIterator operator--(difference_type) noexcept
   {
      RAxisIterator old(*this);
      --(*this);
      return old;
   }

   // i += 2
   constexpr RAxisIterator &operator+=(difference_type d) noexcept
   {
      fCursor += d;
      return *this;
   }

   // i -= 2
   constexpr RAxisIterator &operator-=(difference_type d) noexcept
   {
      fCursor -= d;
      return *this;
   }

   // i + 2
   constexpr RAxisIterator operator+(difference_type d) noexcept
   {
      RAxisIterator ret(*this);
      ret += d;
      return ret;
   }
   friend constexpr RAxisIterator operator+(difference_type d, RAxisIterator rhs) noexcept;

   // i - 2
   constexpr RAxisIterator operator-(difference_type d) noexcept
   {
      RAxisIterator ret(*this);
      ret -= d;
      return ret;
   }

   // i - j
   constexpr value_type operator-(const RAxisIterator &j) noexcept { return fCursor - j.fCursor; }

   // i[2]
   constexpr value_type operator[](difference_type d) noexcept { return fCursor + d; }

   // *i
   constexpr value_type operator*() const noexcept { return fCursor; }

   // i->
   constexpr const value_type *operator->() const noexcept { return &fCursor; }

   friend constexpr bool operator<(RAxisIterator lhs, RAxisIterator rhs) noexcept;
   friend constexpr bool operator>(RAxisIterator lhs, RAxisIterator rhs) noexcept;
   friend constexpr bool operator<=(RAxisIterator lhs, RAxisIterator rhs) noexcept;
   friend constexpr bool operator>=(RAxisIterator lhs, RAxisIterator rhs) noexcept;
   friend constexpr bool operator==(RAxisIterator lhs, RAxisIterator rhs) noexcept;
   friend constexpr bool operator!=(RAxisIterator lhs, RAxisIterator rhs) noexcept;
};

///\name RAxisIterator external operators
///\{

/// 2 + i
inline constexpr RAxisIterator operator+(int d, RAxisIterator rhs) noexcept
{
   return rhs + d;
}

/// i < j
inline constexpr bool operator<(RAxisIterator lhs, RAxisIterator rhs) noexcept
{
   return lhs.fCursor < rhs.fCursor;
}

/// i > j
inline constexpr bool operator>(RAxisIterator lhs, RAxisIterator rhs) noexcept
{
   return lhs.fCursor > rhs.fCursor;
}

/// i <= j
inline constexpr bool operator<=(RAxisIterator lhs, RAxisIterator rhs) noexcept
{
   return lhs.fCursor <= rhs.fCursor;
}

/// i >= j
inline constexpr bool operator>=(RAxisIterator lhs, RAxisIterator rhs) noexcept
{
   return lhs.fCursor >= rhs.fCursor;
}

/// i == j
inline constexpr bool operator==(RAxisIterator lhs, RAxisIterator rhs) noexcept
{
   return lhs.fCursor == rhs.fCursor;
}

/// i != j
inline constexpr bool operator!=(RAxisIterator lhs, RAxisIterator rhs) noexcept
{
   return lhs.fCursor != rhs.fCursor;
}
///\}

} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RAxisIterator header guard
