// Storage for a whole spline
// And it owns corresponding Observable class

// Should I also hire this class for initialization of
// Observable class?

// But anyhow, maybe I need some factory to select the proper
// Observable instance for me

#ifndef _OBSERVABLE_SPLINES_
#define _OBSERVABLE_SPLINES_
#include "Framework/EventGen/EventRecord.h"
#include "ProfSpline/ExtendedSpline.h"
#include "ProfSpline/ObservableI.h"
#include <memory>

namespace genie {
namespace rew {
class ObservableSplines {
public:
  ObservableI *get_observable() const;
  const ExtendedSpline &get_bin(size_t bin_id) const;
  // the interpolation function
  double get_value(const EventRecord *, const std::vector<double> &) const;

  // This function is to restore binning information
  virtual size_t GetObservablesBinID(const std::vector<double> &) const;

  virtual void
  InitializeBins(const std::vector<std::vector<double>> &bin_edges);

  // TODO: figure out how to initialize bins
  // TODO: figure out if we want to initialize Observable here

private:
  std::unique_ptr<ObservableI> observable; // corresponding Observable
  std::vector<ExtendedSpline> bins;        // each bin
  std::vector<std::vector<double>> bin_edges;
};
} // namespace rew
} // namespace genie

#endif