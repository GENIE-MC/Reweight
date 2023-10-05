// Storage for a whole spline
// And it owns corresponding Observable class

// Should I also hire this class for initialization of
// Observable class?

// But anyhow, maybe I need some factory to select the proper
// Observable instance for me

#ifndef _OBSERVABLE_SPLINES_
#define _OBSERVABLE_SPLINES_
#include "Framework/EventGen/EventRecord.h"
#include "ProfSpline/ObservableBins.h"
#include "ProfSpline/ObservableDiscreteBins.h"
#include "ProfSpline/ObservableI.h"
#include "Professor/Ipol.h"
#include <memory>
#include <vector>

namespace genie {
namespace rew {

/// \uml{note the actual calculation of weight and hold the information needed}
class ObservableSplines {
public:
  const Professor::Ipol &GetBin(size_t bin_id) const;
  // the interpolation function

  double GetDXsec(const EventRecord &, const std::vector<double> &) const;
  double GetRatio(const EventRecord &, const std::vector<double> &,
                  const std::vector<double> &) const;

  std::vector<int> GetObservablesBinID(const EventRecord &) const;

  size_t GetObservablesBinIDLinearized(const EventRecord &) const;
  // size_t GetObservablesBinIDLinearized(const std::vector<double> &) const;

  void InitializeBins(const std::vector<std::vector<double>> &bin_edges);

  // figure out how to initialize bins
  void InitializeIpols(const std::vector<std::string> &lines);

  void InitializeObservable(const std::string name,
                            const std::string config = "NoConfig");

  void
  InitializeDiscreteBins(const std::vector<std::string> &enabled_bin_names);

  // TODO: figure out if we want to initialize Observable here
  size_t GetNChannel() const;

  double GetValueInterpolated(size_t channel_id,
                              const std::vector<double> &obvs,
                              const std::vector<double> &paras) const;
  size_t toBinID(size_t channel_id, std::vector<int> bin_ids) const;

  size_t GetChannelID(const EventRecord &) const;
  double GetCellSize(std::vector<int> bin_ids) const;

private:
  // size_t GetObservablesBinID(const std::vector<double> &) const;

  const ObservableI *observable;

  std::vector<Professor::Ipol> bins;

  ObservableBins binning;

  std::unique_ptr<ObservableDiscreteBins> discrete_bins{};
};
} // namespace rew
} // namespace genie

#endif