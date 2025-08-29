//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightResonanceDecay

\brief    Reweighting resonance decays

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Apr 26, 2010

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_RESDEC_H_
#define _G_REWEIGHT_RESDEC_H_

//#define _G_REWEIGHT_RESDEC_DEBUG_

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

 class GReWeightResonanceDecay : public GReWeightModel
 {
 public:
   GReWeightResonanceDecay();
  ~GReWeightResonanceDecay();

   // implement the GReWeightI interface
   bool   AppliesTo      (const EventRecord & event) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

   // various config options
   void RewNue      (bool tf ) { fRewNue     = tf; }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf; }
   void RewNumu     (bool tf ) { fRewNumu    = tf; }
   void RewNumubar  (bool tf ) { fRewNumubar = tf; }
   void RewCC       (bool tf ) { fRewCC      = tf; }
   void RewNC       (bool tf ) { fRewNC      = tf; }

 private:

   void   Init              (void);
   double RewBR             (const EventRecord & event);
   double RewThetaDelta2Npi (const EventRecord & event);

   double fBR1gammaTwkDial;
   double fBR1etaTwkDial;
   double fThetaDelta2NpiTwkDial;

   bool   fRewNue;       ///< reweight nu_e?
   bool   fRewNuebar;    ///< reweight nu_e_bar?
   bool   fRewNumu;      ///< reweight nu_mu?
   bool   fRewNumubar;   ///< reweight nu_mu_bar?
   bool   fRewCC;        ///< reweight CC?
   bool   fRewNC;        ///< reweight NC?

   std::map<int, TH1D*> fMpBR1gammaDef; // resonance pdg -> X + 1gamma, default BR = f(W)
   std::map<int, TH1D*> fMpBR1etaDef;   // resonance pdg -> X + 1eta,   default BR = f(W)

#ifdef _G_REWEIGHT_RESDEC_DEBUG_
   TFile *    fTestFile;
   TNtupleD * fTestNtp;
#endif
 };

} // rew
} // genie

#endif
