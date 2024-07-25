// Storage for a whole spline
// And it owns corresponding Observable class

// Should I also hire this class for initialization of
// Observable class?

// But anyhow, maybe I need some factory to select the proper
// Observable instance for me

#ifndef _OBSERVABLE_SPLINES_
#define _OBSERVABLE_SPLINES_
#include "Framework/Algorithm/AlgId.h"
#include "Framework/EventGen/EventRecord.h"
#include "ProfSpline/KinematicVariables.h"
#include "ProfSpline/RwgKineSpace.h"
#include "Professor/Ipol.h"
#include "libxml/tree.h"
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

  void InitializeBins(const std::vector<std::vector<double>> &bin_edges);

  // figure out how to initialize bins
  void InitializeIpols(const std::vector<std::string> &lines);

  void InitializeObservable(const std::string name, const std::string config);
  void InitializeObservable(const std::string AlgID);
  void InitializeObservable(AlgId id);

  double GetValueInterpolated(const std::vector<double> &obvs,
                              const std::vector<double> &paras) const;
  size_t lookupBinID(const std::vector<double> &obvs) const;

  bool IsHandled(const EventRecord &event) const {
    return observable->IsHandled(event);
  }

  KinematicVariables CalcKinematicVariables(const EventRecord &event) const {
    return observable->CalcKinematicVariables(event);
  }

  ChannelIDs ChannelID(const EventRecord &event) const {
    return observable->ChannelID(event);
  }

  template <class BinningT, class FirstNeighbors, class channelT>
  ObservableSplines(BinningT &&bin_in, FirstNeighbors &&first_neighbour_in,
                    channelT &&channel_in, int dimension)
      : bin_edges(std::forward<BinningT>(bin_in)),
        first_neighbour(std::forward<FirstNeighbors>(first_neighbour_in)),
        dimension(dimension),
        channel(std::forward<channelT>(channel_in)) {}

  ObservableSplines() = default;

private:
  // size_t GetObservablesBinID(const std::vector<double> &) const;

  const RwgKineSpace *observable;

  std::vector<Professor::Ipol> bins;

  std::vector</*different bins*/ std::vector<
      /*dimensions*/ std::pair</*xmin*/ double, /*xmax*/ double>>>
      bin_edges{};
  std::vector<std::set<size_t>> first_neighbour{};

  size_t dimension{};
  ChannelIDs channel{};
  // std::string name{};
};
} // namespace rew
} // namespace genie

#endif
