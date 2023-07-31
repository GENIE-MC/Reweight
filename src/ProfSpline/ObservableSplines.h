// Storage for a whole spline
// And it owns corresponding Observable class

// Should I also hire this class for initialization of
// Observable class?

// But anyhow, maybe I need some factory to select the proper
// Observable instance for me

#ifndef _OBSERVABLE_SPLINES_
#define _OBSERVABLE_SPLINES_
#include "Framework/EventGen/EventRecord.h"
#include "ProfSpline/Ipol.h"
#include "ProfSpline/ObservableI.h"
#include <memory>

namespace genie {
namespace rew {
class ObservableSplines {
public:
  ObservableI *GetObservable() const;
  const Professor::Ipol &GetBin(size_t bin_id) const;
  const Professor::Ipol &operator[](size_t bin_id) const;
  const Professor::Ipol &operator[](const std::vector<double> &vars) const;
  // the interpolation function
  /// \uml{note takes eventrecord and nuisance parameters as input
  /// and return the differential cross section}
  double GetDXsec(const EventRecord *, const std::vector<double> &) const;

  // This function is to restore binning information
  /// \uml{note Calculate the bin id for a given set of observables}
  size_t GetObservablesBinID(const std::vector<double> &) const;

  /// \uml{note[right] Bin edges for each dimension
  /// Bin edges should be read from configuration file
  /// of comparsion package?}
  void InitializeBins(const std::vector<std::vector<double>> &bin_edges);

  // figure out how to initialize bins
  void InitialIpols(const std::vector<std::string> &lines);

  // TODO: figure out if we want to initialize Observable here

private:
  /// \uml{note[right] Actual class to calculate observable}
  std::unique_ptr<ObservableI> observable;

  std::vector<Professor::Ipol> bins;

  std::vector<std::vector<double>> bin_edges;
};
} // namespace rew
} // namespace genie

#endif