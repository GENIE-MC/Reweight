#include "ProfSpline/ObservableBins.h"
#include "Framework/Messenger/Messenger.h"
namespace genie {
namespace rew {
std::vector<int>
ObservableBins::GetObservablesBinID(const std::vector<double> &vars) const {
  if (bin_edges.size() != vars.size()) {
    LOG("ObservableBins", pERROR) << "Number of dimensions in bin edges and "
                                     "number of dimensions in vars do not "
                                     "match";
    exit(1);
  }
  std::vector<int> bin_ids(vars.size());
  for (size_t i = 0; i < vars.size(); ++i) {
    bin_ids[i] = bin_edges[i].FindBin(vars[i]);
  }
  return bin_ids;
}

int ObservableBins::LinearizeBinID(const std::vector<int> &bin_ids) const {
  if (bin_ids.size() != bin_edges.size()) {
    LOG("ObservableBins", pERROR) << "Number of dimensions in bin edges and "
                                     "number of dimensions in bin_ids do not "
                                     "match";
    exit(1);
  }
  int bin_id = 0;
  for (size_t i = 0; i < bin_ids.size(); ++i) {
    bin_id *= bin_edges[i].GetNbins();
    if (bin_ids[i] == 0) {
      bin_id += 0;
    } else if (bin_ids[i] == bin_edges[i].GetNbins() + 1) {
      bin_id += bin_edges[i].GetNbins();
    } else {
      bin_id += bin_ids[i] - 1;
    }
  }
  return bin_id;
}

int ObservableBins::GetObservablesBinIDLinearized(
    const std::vector<double> &vars) const {
  return LinearizeBinID(GetObservablesBinID(vars));
}

void ObservableBins::InitializeBins(
    const std::vector<std::vector<double>> &bin_edges_in) {
  for (const auto &bin_edge : bin_edges_in) {
    bin_edges.emplace_back(bin_edge.size() - 1, bin_edge.data());
  }
}

} // namespace rew
} // namespace genie