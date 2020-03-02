//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightDeltaradAngle

\brief    Reweighting resonance decays

\author   Gray Yarbrough <gyarbrou \at vols.utk.edu>
          University of Tennessee, Knoxvile

          Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Mar 2, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_DELTARAD_ANGLE_H_
#define _G_REWEIGHT_DELTARAD_ANGLE_H_

#include <map>

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"

class TH1D;
class TNtupleD;
class TFile;

using namespace genie::rew;
using namespace genie;

namespace genie {
namespace rew   {

 class GReWeightDeltaradAngle : public GReWeightModel
 {
 public:
   GReWeightDeltaradAngle();
  ~GReWeightDeltaradAngle();

   // implement the GReWeightI interface
   bool   AppliesTo      (ScatteringType_t type, bool is_cc) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

 private:

   void   Init(void);
   double RewThetaDelta2NRad (const EventRecord & event);

   double fThetaDelta2NRadTwkDial;
 };

} // rew
} // genie

#endif // _G_REWEIGHT_DELTARAD_ANGLE_H_
