#ifndef _ObservablePMuEnu_H_
#define _ObservablePMuEnu_H_
#include "ProfSpline/ObservableI.h"

namespace genie {
namespace rew {
class ObservablePMuEnu : public RwgKineSpace {
public:
  ObservablePMuEnu();
  ObservablePMuEnu(std::string config);
  virtual std::vector<double>
  KinematicVariables(const EventRecord &event) const override;

  virtual ~ObservablePMuEnu() = default;

  virtual bool IsHandled(const EventRecord &event) const override;

private:
  virtual void LoadConfig(void) override;
  std::string channelid{};
};
} // namespace rew
} // namespace genie

#endif