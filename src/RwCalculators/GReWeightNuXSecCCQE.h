//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecCCQE

\brief    Reweighting CCQE GENIE neutrino cross sections

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_CCQE_H_
#define _G_REWEIGHT_NU_XSEC_CCQE_H_

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

 class GReWeightNuXSecCCQE : public GReWeightModel
 {
 public:
   static const int kModeMa               = 0;
   static const int kModeNormAndMaShape   = 1;
   static const int kModeZExp             = 2;

   static const int fZExpMaxSyst          = 4; ///< maximum number of systematics

   GReWeightNuXSecCCQE();
   GReWeightNuXSecCCQE(std::string model, std::string type);
  ~GReWeightNuXSecCCQE();

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
   void SetMaPath   (string p) { fMaPath     = p;    }
   // z-expansion specific options
   void SetZExpPath    (string p){ fZExpPath    = p;   }
   // RunningMa specific options
   void SetE0Path    (string p){ fE0Path    = p;   }

 private:

   void   Init                (void);
   double CalcWeightNorm      (const EventRecord & event);
   double CalcWeightMaShape   (const EventRecord & event);
   double CalcWeightMa        (const EventRecord & event);
   double CalcWeightZExp      (const EventRecord & event);
   double CalcWeightRPA       (const EventRecord & event);
   double CalcWeightCoulomb   (const EventRecord & event);

   XSecAlgorithmI * fXSecModelDef;    ///< default model
   XSecAlgorithmI * fXSecModel;       ///< tweaked model

   XSecIntegratorI * fXSecIntegratorDef; ///< integrator for default model
   XSecIntegratorI * fXSecIntegrator; ///< integrator for tweaked model

   Registry *       fXSecModelConfig; ///< config in tweaked model
   string fFFModel; ///< String name of form factor model
   bool fModelIsDipole;           ///< Using dipole form factors?
   bool fModelIsZExp;             ///< Using Zexp form factors?
   bool fModelIsMArunAxial;       ///< Using  running axial mass form factors?
   std::string fManualModelName; ///< If using a tweaked model that isn't the same as default, name
   std::string fManualModelType; ///< If using a tweaked model that isn't the same as default, type

   int    fMode;         ///< 0: Ma, 1: Norm and MaShape, 2: Z-Expansion
   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?
   string fMaPath;       ///< M_{A} path in config Registry
   double fNormTwkDial;  ///<
   double fNormDef;      ///<
   double fNormCurr;     ///<
   double fMaTwkDial;    ///<
   double fMaDef;        ///<
   double fMaCurr;       ///<
   double fE0TwkDial;    ///<
   double fE0Def;        ///<
   double fE0Curr;       ///<
   string fE0Path;       ///< E_{0} for RunningMA path in config Registry



   // unused // int     fZExpCurrIdx; ///< current coefficient index
   int     fZExpMaxCoef; ///< max number of coefficients to use
   string  fZExpPath;    ///< algorithm path to get coefficients
   double  fZExpTwkDial[fZExpMaxSyst]; ///<
   double  fZExpDef    [fZExpMaxSyst]; ///<
   double  fZExpCurr   [fZExpMaxSyst]; ///< array of current parameter values

   double fRPATwkDial; ///< 0 = default, 1 = RPA off (changes Nieves CCQE only)

   /// Scales the EM potential used for Coulomb corrections (changes Nieves CCQE
   /// only)
   double fCoulombTwkDial;

   /// Copy of the default cross section model with RPA turned off
   XSecAlgorithmI* fXSecModelDefNoRPA;

   /// Copy of the default cross section model with a tweaked value of the EM
   /// potential used for Coulomb corrections
   XSecAlgorithmI* fXSecModelDefCoulomb;

#ifdef _G_REWEIGHT_CCQE_DEBUG_
   TFile *    fTestFile;
   TNtupleD * fTestNtp;
#endif
 };

} // rew   namespace
} // genie namespace

#endif
