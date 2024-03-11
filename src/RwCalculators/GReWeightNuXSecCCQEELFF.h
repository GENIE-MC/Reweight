//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecCCQEELFFELFF

\brief    Reweighting CCQE GENIE neutrino cross sections

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_CCQE_ELFF_H_
#define _G_REWEIGHT_NU_XSEC_CCQE_ELFF_H_

//#define _G_REWEIGHT_CCQE_DEBUG_

#include <string>

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"

class TFile;
class TNtupleD;

namespace genie {

class XSecAlgorithmI;
class XSecIntegratorI;
class Registry;

namespace rew   {

 class GReWeightNuXSecCCQEELFF : public GReWeightModel
 {
 public:
   static const int kModeZExp             = 0;

   GReWeightNuXSecCCQEELFF();
   GReWeightNuXSecCCQEELFF(std::string model, std::string type);
  ~GReWeightNuXSecCCQEELFF();

   // implement the GReWeightI interface
   bool   AppliesTo      (const EventRecord &event) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

   // various config options
   void SetMode     (int mode) { fMode       = mode; }
   void RewNue      (bool tf ) { fRewNue     = tf;   }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf;   }
   void RewNumu     (bool tf ) { fRewNumu    = tf;   }
   void RewNumubar  (bool tf ) { fRewNumubar = tf;   }
   // z-expansion specific options
   void SetZExpPath    (string p){ fZExpPath    = p;   }

 private:

   void   Init                (void);
   double CalcWeightZExp      (const EventRecord & event);
   XSecAlgorithmI * fXSecModelDef;    ///< default model
   XSecAlgorithmI * fXSecModel;       ///< tweaked model
   XSecIntegratorI * fXSecIntegratorDef; ///< integrator for default model
   XSecIntegratorI * fXSecIntegrator; ///< integrator for tweaked model
   Registry *       fXSecModelConfig; ///< config in tweaked model
   string fFFModel; ///< String name of form factor model
   bool fModelIsZExp;             ///< Using Zexp form factors?


   std::string fManualModelName; ///< If using a tweaked model that isn't the same as default, name
   std::string fManualModelType; ///< If using a tweaked model that isn't the same as default, type

   int    fMode;         ///< 0: Ma, 1: Norm and MaShape, 2: Z-Expansion
   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?

   // unused // int     fZExpCurrIdx; ///< current coefficient index
   int     fZExpMaxCoef; ///< max number of coefficients to use
   string  fZExpPath;    ///< algorithm path to get coefficients

   struct ZExpPara {
     bool   fQ4limit;
     int    fKmax;
     double fT0;
     double fTcut;
     double fGep0;
     double fGmp0;
     double fGen0;
     double fGmn0;
     //double fZ_An[11];
     double *fZ_APn;
     double *fZ_BPn;
     double *fZ_ANn;
     double *fZ_BNn;
   } fZExpParaDef, fZExpPara, fZExpParaTwkDial;

#ifdef _G_REWEIGHT_CCQE_DEBUG_
   TFile *    fTestFile;
   TNtupleD * fTestNtp;
#endif
 };

} // rew   namespace
} // genie namespace

#endif
