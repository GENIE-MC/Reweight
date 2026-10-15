//____________________________________________________________________________
/*!

\class    genie::rew::ObservablePMuEnuW

\brief    RwgKineSpace implementation that uses the final-state muon momentum,
          incoming neutrino energy, and hadronic invariant mass W as the
          kinematic variables for GENIE/Professor reweighting.

\author   Qiyu Yan <yanqiyu17 \at mails.ucas.ac.cn>
          University of Chinese Academy of Sciences (UCAS) and University of Warwick

\created  August 2023

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _ObservablePMuEnuW_H_
#define _ObservablePMuEnuW_H_
#include "ProfSpline/RwgKineSpace.h"

namespace genie {
namespace rew {
class ObservablePMuEnuW : public RwgKineSpace {
public:
  ObservablePMuEnuW();
  ObservablePMuEnuW(std::string config);

public:
  virtual KinematicVariables
  CalcKinematicVariables(const EventRecord &event) const override;

  virtual ChannelIDs ChannelID(const EventRecord &event) const override;

  virtual ~ObservablePMuEnuW() = default;
};
} // namespace rew
} // namespace genie

#endif
