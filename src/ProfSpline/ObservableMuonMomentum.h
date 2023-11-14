#ifndef _OBSERVABLE_MUON_MOMENTUM_
#define _OBSERVABLE_MUON_MOMENTUM_

#include "Framework/EventGen/EventRecord.h"
#include "ProfSpline/ObservableI.h"

namespace genie {
namespace rew {
class ObservableMuonMomentum : public ObservableI {
public:
  ObservableMuonMomentum();
  virtual ~ObservableMuonMomentum() = default;

  virtual std::vector<double>
  GetKinematicVariables(const EventRecord &event) const override;

  virtual bool IsHandled(const EventRecord &event) const override;
  virtual bool IsCC() const override;

  virtual void Configure(const Registry &config) override;
  virtual void Configure(string param_set) override;

private:
  void LoadConfig(void);
  std::string channelid{};
  bool isCC{true};
};

} // namespace rew
} // namespace genie

#endif