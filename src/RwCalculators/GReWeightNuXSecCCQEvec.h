//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecCCQEvec

\brief    Reweighting vector form factors in GENIE CCQE neutrino cross
          section calculations.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  May 24, 2010

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_CCQE_VEC_H_
#define _G_REWEIGHT_NU_XSEC_CCQE_VEC_H_

//#define _G_REWEIGHT_CCQE_VEC_DEBUG_

#include <map>
#include <string>

// GENIE/Reweight includes
#include "Framework/EventGen/EventRecord.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "RwCalculators/GReWeightModel.h"


class TFile;
class TNtupleD;

namespace genie {

class XSecAlgorithmI;
class XSecIntegratorI;

namespace rew   {

 class GReWeightNuXSecCCQEvec : public GReWeightModel
 {
 public:
   GReWeightNuXSecCCQEvec();
  ~GReWeightNuXSecCCQEvec();

   // implement the GReWeightI interface
   bool   AppliesTo      (const EventRecord &event) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

   // various config options
   void RewNue      (bool tf ) { fRewNue     = tf;   }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf;   }
   void RewNumu     (bool tf ) { fRewNumu    = tf;   }
   void RewNumubar  (bool tf ) { fRewNumubar = tf;   }

 private:

   void Init (void);

   XSecAlgorithmI * fXSecModel_bba;  ///< CCQE model with BBA05  f/f (default)
   XSecAlgorithmI * fXSecModel_dpl;  ///< CCQE model with dipole f/f ("maximally" tweaked)

   XSecIntegratorI * fXSecIntegrator_bba;  ///< Integrator for BBA05 cross section
   XSecIntegratorI * fXSecIntegrator_dpl;  ///< Integrator for Dipole cross section

   double fFFTwkDial;    ///< tweaking dial (0: bba/default, +1: dipole)

   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?

#ifdef _G_REWEIGHT_CCQE_VEC_DEBUG_
   TFile *    fTestFile;
   TNtupleD * fTestNtp;
#endif
 };

} // rew   namespace
} // genie namespace

#endif
