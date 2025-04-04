//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightINukeKinematics

\brief    Reweighting GENIE INTRANUKE/hA hadron transport model.

\author   Gray Putnam <gputnam \at fnal.gov>
          Fermi National Accelerator Laboratory
          
\created  Feb 14, 2025

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_INUKEKINEMATICS_H_
#define _G_REWEIGHT_INUKEKINEMATICS_H_

//#define _G_REWEIGHT_INUKEKINEMATICS_DEBUG_NTP_

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"
#include "RwCalculators/GReWeightINukeKinematicsParams.h"

using namespace genie::rew;
using namespace genie;

class TFile;
class TNtuple;
class TLorentzVector;

namespace genie {

 class HAIntranuke2018;
 class GHepParticle;

namespace rew   {

 class GReWeightINukeKinematics : public GReWeightModel
 {
 public:
   GReWeightINukeKinematics(): GReWeightModel("IntraNukeKin") {}
  ~GReWeightINukeKinematics() {}

   // implement the GReWeightI interface
   bool   AppliesTo      (const EventRecord & event) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

 private:
   GReWeightINukeKinematicsParams fRwParams;

#ifdef _G_REWEIGHT_INUKEKINEMATICS_DEBUG_NTP_
   TFile *              fTestFile;
   TNtuple *            fTestNtp;
#endif

 };

} // rew
} // genie

#endif
