#ifndef _OBSERVABLE_MUON_MOMENTUM_
#define _OBSERVABLE_MUON_MOMENTUM_

#include "Framework/EventGen/EventRecord.h"
#include "ProfSpline/ObservableI.h"

namespace genie {
namespace rew {
class ObservableMuonMomentum : public RwgKineSpace {
public:
  ObservableMuonMomentum();
  ObservableMuonMomentum(std::string config);
  virtual ~ObservableMuonMomentum() = default;

  virtual std::vector<double>
  KinematicVariables(const EventRecord &event) const override;

  virtual bool IsHandled(const EventRecord &event) const override;

private:
  virtual void LoadConfig(void) override;
  std::string channelid{};
};

} // namespace rew
} // namespace genie

#endif