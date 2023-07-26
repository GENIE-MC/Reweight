#ifndef _OBSERVABLE_MUON_MOMENTUM_
#define _OBSERVABLE_MUON_MOMENTUM_

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepVirtualList.h"
#include "ProfSpline/ObservableI.h"

namespace genie {
namespace rew {
class ObservableMuonMomentum : public ObservableI {
public:
  virtual std::vector<double>
  HandleEventRecord(const EventRecord *event) override;
};

} // namespace rew
} // namespace genie

#endif