//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightFZone

\brief    Reweighting the formation zone model

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Sep 20, 2009

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_FZONE_H_
#define _G_REWEIGHT_FZONE_H_

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"

using namespace genie::rew;
using namespace genie;

namespace genie {

 class HAIntranuke2018;

namespace rew   {

 class GReWeightFZone : public GReWeightModel
 {
 public:
   GReWeightFZone();
  ~GReWeightFZone();

   // implement the GReWeightI interface
   bool   AppliesTo      (ScatteringType_t type, bool is_cc) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

   // other config options
   // set to match values used at event generation
   void SetCT0Pion(double ct0)    { fct0pion    = ct0; }
   void SetCT0Nucleon(double ct0) { fct0nucleon = ct0; }
   void SetK  (double k  )        { fK          = k;   }

 private:

   void Init(void);

   double fFZoneTwkDial; ///< formation zone tweaking dial

   double fct0pion;    ///<
   double fct0nucleon; ///<
   double fK;          ///<

   HAIntranuke2018* fFSIModel;

 };

} // rew
} // genie

#endif
