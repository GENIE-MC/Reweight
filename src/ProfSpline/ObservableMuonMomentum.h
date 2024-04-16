#ifndef _OBSERVABLE_MUON_MOMENTUM_
#define _OBSERVABLE_MUON_MOMENTUM_

#include "Framework/EventGen/EventRecord.h"
#include "ProfSpline/RwgKineSpace.h"

namespace genie {
namespace rew {
class ObservableMuonMomentum : public RwgKineSpace {
public:
  ObservableMuonMomentum();
  ObservableMuonMomentum(std::string config);

public:
  virtual ~ObservableMuonMomentum() = default;

  virtual KinematicVariables
  CalcKinematicVariables(const EventRecord &event) const override;

  virtual ChannelIDs ChannelID(const EventRecord &event) const override;

  virtual bool IsHandled(const EventRecord &event) const override;

private:
  virtual void LoadConfig(void) override;
  std::string channelid{};
};

} // namespace rew
} // namespace genie

#endif
