#ifndef _OBSERVABLE_DISCRETE_BINS_H_
#define _OBSERVABLE_DISCRETE_BINS_H_
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/EventGen/EventRecord.h"
namespace genie {
namespace rew {
/// \uml{note Made this an interface to handle the possibliy of
/// different coice of different channel division. }
class ObservableDiscreteBinI : public Algorithm {
public:
  virtual ~ObservableDiscreteBinI() = default;
  virtual size_t GetBinID(const EventRecord &) const = 0;
  virtual size_t GetNBinMax() const = 0;
};

class ObservableDiscreteBinCCNC : public ObservableDiscreteBinI {
public:
  virtual ~ObservableDiscreteBinCCNC() = default;
  virtual size_t GetBinID(const EventRecord &) const;
  virtual size_t GetNBinMax() const;
};

class ObservableDiscreteBinInteraction : public ObservableDiscreteBinI {
public:
  virtual ~ObservableDiscreteBinInteraction() = default;
  virtual size_t GetBinID(const EventRecord &) const;
  virtual size_t GetNBinMax() const;
};

/// \uml{note Hold a set of ObservableDiscreteBinI pointers 
/// to handld multiple dimensions of discrete bins
/// And act like a single dimension binning.}
class ObservableDiscreteBins {
public:
  ObservableDiscreteBins(std::vector<std::string> enabled_bin_names);
  size_t GetBinID(const EventRecord &) const;
  size_t GetNBinMax() const;

private:
  // genie::AlgFactory will take care of deleting these pointers
  std::vector<const ObservableDiscreteBinI *> bins;
  size_t nbin_max{1};
};
} // namespace rew
} // namespace genie

#endif