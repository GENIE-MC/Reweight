#ifndef OBSERVABLEHADRONIZATION_H
#define OBSERVABLEHADRONIZATION_H
#include "ProfSpline/RwgKineSpace.h"

namespace genie::rew {
class ObservableHadronization : public RwgKineSpace {
public:
  using RwgKineSpace::RwgKineSpace;

  virtual KinematicVariables
  CalcKinematicVariables(const EventRecord &event) const override;

  virtual ChannelIDs ChannelID(const EventRecord &event) const override;

  virtual ~ObservableHadronization() = default;

  virtual bool IsHandled(const EventRecord &event) const override;

private:
  virtual void LoadConfig(void) override;
};

class ObservableHadronizationLowW final : public ObservableHadronization {
public:
  ObservableHadronizationLowW();
  ObservableHadronizationLowW(std::string config);
  virtual KinematicVariables
  CalcKinematicVariables(const EventRecord &event) const override;
  virtual ~ObservableHadronizationLowW() = default;
  virtual bool IsHandled(const EventRecord &event) const override;
};

class ObservableHadronizationHighW final : public ObservableHadronization {
public:
  ObservableHadronizationHighW();
  ObservableHadronizationHighW(std::string config);
  virtual KinematicVariables
  CalcKinematicVariables(const EventRecord &event) const override;
  virtual ~ObservableHadronizationHighW() = default;
  virtual bool IsHandled(const EventRecord &event) const override;
};

} // namespace genie::rew

#endif
