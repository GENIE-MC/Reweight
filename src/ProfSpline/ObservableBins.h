#ifndef _OBSERVABLE_BINS_
#define _OBSERVABLE_BINS_

#include "TAxis.h"
#include <TNamed.h>
#include <string>
#include <vector>
namespace genie {
namespace rew {
/// \uml{note conventional binning for ordinary observables
/// Actually hold a set of TAxis in multiple dimensions
/// And provide interface to act like a single dimension binning.}
class ObservableBins {
public:
  ObservableBins() = default;
  template <typename T>
  ObservableBins(T &&m_bin_edges) : bin_edges(std::forward<T>(m_bin_edges)){};
  void InitializeBins(const std::vector<std::vector<double>> &bin_edges_in);

  int GetObservablesBinIDLinearized(const std::vector<double> &) const;
  int GetObservablesBinIDLinearized(const std::vector<int> &) const;
  
  // This the ROOT convention, the bin id starts from 1
  // 0 is underflow, nbin+1 is overflow
  std::vector<int> GetObservablesBinID(const std::vector<double> &) const;
  std::vector<double> GetObservablesBinLoc(const std::vector<double> &) const;
  // notice that the linearized bin id starts from 0
  int LinearizeBinID(const std::vector<int> &) const;
  const TAxis & GetAxis(size_t i) const { return bin_edges[i]; }

private:
  std::vector<TAxis> bin_edges;
};
} // namespace rew
} // namespace genie
#endif