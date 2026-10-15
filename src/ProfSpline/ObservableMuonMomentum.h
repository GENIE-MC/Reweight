//____________________________________________________________________________
/*!

\class    genie::rew::ObservableMuonMomentum

\brief    RwgKineSpace implementation that uses the final-state muon momentum
          as the single kinematic variable for GENIE/Professor reweighting.

\author   Qiyu Yan <yanqiyu17 \at mails.ucas.ac.cn>
          University of Chinese Academy of Sciences (UCAS) and University of Warwick

\created  July 2023

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

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
