//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightFGM

\brief    Reweighting the Fermi Gas nuclear model

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Apr 26, 2010

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_FGM_H_
#define _G_REWEIGHT_FGM_H_

//#define _G_REWEIGHT_FGM_DEBUG_

#include <map>

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"

using namespace genie::rew;
using namespace genie;

class TH1D;
class TNtupleD;
class TFile;

namespace genie {

class NuclearModelI;

namespace rew   {

 class GReWeightFGM : public GReWeightModel
 {
 public:
   GReWeightFGM();
  ~GReWeightFGM();

   // implement the GReWeightI interface
   bool   AppliesTo      (const EventRecord & event) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

 private:

   void Init(void);

   double RewCCQEPauliSupViaKF   (const EventRecord & event);
   double RewCCQEMomDistroFGtoSF (const EventRecord & event);

   double fKFTwkDial;
   double fMomDistroTwkDial;

   const NuclearModelI * fFG;
   const NuclearModelI * fSF;

   std::map<int, TH1D *> fMapFGn;
   std::map<int, TH1D *> fMapFGp;
   std::map<int, TH1D *> fMapSFn;
   std::map<int, TH1D *> fMapSFp;

#ifdef _G_REWEIGHT_FGM_DEBUG_
   TFile *    fTestFile;
   TNtupleD * fTestNtp;
#endif

 };

} // rew
} // genie

#endif
