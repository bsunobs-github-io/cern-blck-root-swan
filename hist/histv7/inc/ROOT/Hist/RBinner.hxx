/// \file ROOT/RHist.hxx
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

#ifndef ROOT7_RBinner
#define ROOT7_RBinner

#include "ROOT/RLogger.hxx"

#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace ROOT {
namespace Experimental {

namespace Detail {
namespace Hist {

/// Abstract interface for converting coordinate values into a bin index.
class RBinnerBase {
protected:
   RBinnerBase() = default;

public:
   virtual ssize_t GetBin(std::vector<void *> coord) const = 0;

   static std::unique_ptr<RBinnerBase> Create(const std::vector<std::string> &axisTypes);
};

/// Helper to convert coordinate values into a bin index.
/// This class knows about the concrete axis types to be used.
template <class... Axis>
class RBinner : public RBinnerBase {
   std::tuple<Axis...*> fAxes;

public:
   RBinner(Axis... *axes): fAxes(axes...) {}
   ssize_t GetBin(std::vector<void *> coord) const final;

   template <class... Coord_t>
   ssize_t GetBin(Coord_t... coord) const;
};

} // namespace Hist
} // namespace Detail

#endif // ROOT7_RBinner
