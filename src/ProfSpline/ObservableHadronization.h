#pragma once
#ifndef OBSERVABLEHADRONIZATION_H
#define OBSERVABLEHADRONIZATION_H
#include "ProfSpline/RwgKineSpace.h"


namespace genie {
namespace rew {
class ObservableHadronization : public RwgKineSpace {
public:
  ObservableHadronization();
  ObservableHadronization(std::string config);

public:
  virtual KinematicVariables
  CalcKinematicVariables(const EventRecord &event) const override;

  virtual ChannelIDs ChannelID(const EventRecord &event) const override;

  virtual ~ObservableHadronization() = default;

  virtual bool IsHandled(const EventRecord &event) const override;

private:
  virtual void LoadConfig(void) override;
  std::string channelid{};
};
} // namespace rew
} // namespace genie

#endif
