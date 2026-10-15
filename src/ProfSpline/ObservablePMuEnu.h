//____________________________________________________________________________
/*!

\class    genie::rew::ObservablePMuEnu

\brief    RwgKineSpace implementation that uses the final-state muon momentum
          and incoming neutrino energy, evaluated in the hit nucleon/nucleus
          rest frame, as the kinematic variables for GENIE/Professor
          reweighting.

\author   Qiyu Yan <yanqiyu17 \at mails.ucas.ac.cn>
          University of Chinese Academy of Sciences (UCAS) and University of Warwick

\created  August 2023

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _ObservablePMuEnu_H_
#define _ObservablePMuEnu_H_
#include "ProfSpline/RwgKineSpace.h"

namespace genie {
namespace rew {
class ObservablePMuEnu : public RwgKineSpace {
public:
  ObservablePMuEnu();
  ObservablePMuEnu(std::string config);

public:
  virtual KinematicVariables
  CalcKinematicVariables(const EventRecord &event) const override;

  virtual ChannelIDs ChannelID(const EventRecord &event) const override;

  virtual ~ObservablePMuEnu() = default;

  virtual bool IsHandled(const EventRecord &event) const override;

private:
  virtual void LoadConfig(void) override;
  std::string channelid{};
};
} // namespace rew
} // namespace genie

#endif
