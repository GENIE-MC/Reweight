#ifndef _OBSERVABLE_BINS_
#define _OBSERVABLE_BINS_

#include "TAxis.h"
#include <TNamed.h>
#include <string>
#include <vector>
namespace genie {
namespace rew {
class ObservableBins {
public:
  ObservableBins() = default;
  template <typename T>
  ObservableBins(T &&m_bin_edges) : bin_edges(std::forward<T>(m_bin_edges)){};
  void InitializeBins(const std::vector<std::vector<double>> & bin_edges_in);

  int GetObservablesBinIDLinearized(const std::vector<double> &) const;
  std::vector<int> GetObservablesBinID(const std::vector<double> &) const;
  int LinearizeBinID(const std::vector<int> &) const;

private:
  std::vector<TAxis> bin_edges;
};
} // namespace rew
} // namespace genie
#endif