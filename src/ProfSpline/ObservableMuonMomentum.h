#ifndef _OBSERVABLE_MUON_MOMENTUM_
#define _OBSERVABLE_MUON_MOMENTUM_

#include "Framework/EventGen/EventRecord.h"
#include "ProfSpline/RwgKineSpace.h"

namespace genie {
namespace rew {
class ObservableMuonMomentum : public RwgKineSpace {
protected:
  ObservableMuonMomentum();
  ObservableMuonMomentum(std::string config);

public:
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