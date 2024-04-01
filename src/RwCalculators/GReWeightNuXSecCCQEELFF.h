//____________________________________________________________________________
/*!

  \class    genie::rew::GReWeightNuXSecCCQEELFFELFF

  \brief    Reweighting CCQE GENIE neutrino cross sections
            Z expansion vector form factor model

  \author   Liang Liu <liangliu \at fnal.gov>
  Fermi National Accelerator Laboratory

  \created  March 25, 2024

  \cpright  Copyright (c) 2003-2024, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  */
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_CCQE_ELFF_H_
#define _G_REWEIGHT_NU_XSEC_CCQE_ELFF_H_

#include <string>

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"
#include "TRandom3.h"

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

        double GetOneSigma(const EventRecord & event);
        double UpdateXSec(const EventRecord & event);
        void XSecPartialDerivative(const EventRecord & event);

        XSecAlgorithmI * fXSecModelDef;    ///< default model
        XSecAlgorithmI * fXSecModel;       ///< tweaked model
        Registry *       fXSecModelConfig; ///< config in tweaked model
        string fFFModel; ///< String name of form factor model
        bool fModelIsZExp;             ///< Using Zexp form factors?


        std::string fManualModelName; ///< If using a tweaked model that isn't the same as default, name
        std::string fManualModelType; ///< If using a tweaked model that isn't the same as default, type

        int    fMode;         ///< 0: ZExp; TODO: currently there is only one model implemented. 
        bool   fRewNue;       ///< reweight nu_e CC?
        bool   fRewNuebar;    ///< reweight nu_e_bar CC?
        bool   fRewNumu;      ///< reweight nu_mu CC?
        bool   fRewNumubar;   ///< reweight nu_mu_bar CC?

        // unused // int     fZExpCurrIdx; ///< current coefficient index
        int     fZExpMaxCoef; ///< max number of coefficients to use
        string  fZExpPath;    ///< algorithm path to get coefficients

        // parameters in z expansion vector model
        // define the default, current and tweak dial
        struct ZExpPara {
          bool   fQ4limit;
          int    fKmax;
          double fT0;
          double fTcut;
          double fGep0;
          double fGmp0;
          double fGen0;
          double fGmn0;
          double *fZ_APn;
          double *fZ_BPn;
          double *fZ_ANn;
          double *fZ_BNn;
        } fZExpParaDef, fZExpPara, fZExpParaTwkDial;
        // tweek dial and scale factor in propagation method
        double fZExpTwkDial;
        double fZExp_Scale;

        bool fIsSinglePara;  
        bool fIsAllPara;        // flag of propagation method
        double A_f[16];

        // List of the uncertainties of parameters from Kaushik
        // ap1, ap2, ap3, ap4, 
        // bp1, bp2, bp3, bp4, 
        // an1, an2, an3, an4,
        // bn1, bn2, bn3, bn4
        double errors[16]; 
        double error_mat[16][16] = 
        {{ 9.37324e-05,0.000153491,-0.00104771,7.3357e-05,-0.000116564,0.00105353,-0.000441308,-0.00681675,     0,     0,     0,     0,     0,     0,     0,     0},
        {0.000153491,0.00272255,0.0020043,-0.0193067,-0.00104321,-0.00112326,0.02234,-0.0162318,     0,     0,     0,     0,     0,     0,     0,     0},
        {-0.00104771,0.0020043,0.0216328,-0.0375785,0.000483311,-0.0176086,0.0429679,0.0794517,     0,     0,     0,     0,     0,     0,     0,     0},
        {7.3357e-05,-0.0193067,-0.0375785,0.165852,0.00580325,0.0259044,-0.187424,0.0155864,     0,     0,     0,     0,     0,     0,     0,     0},
        {-0.000116564,-0.00104321,0.000483311,0.00580325,0.00089678,-0.00171593,-0.00890089,0.0248036,     0,     0,     0,     0,     0,     0,     0,     0},
        {0.00105353,-0.00112326,-0.0176086,0.0259044,-0.00171593,0.0282469,-0.0571237,-0.150721,     0,     0,     0,     0,     0,     0,     0,     0},
        {-0.000441308,0.02234,0.0429679,-0.187424,-0.00890089,-0.0571237,0.353319,0.0423899,     0,     0,     0,     0,     0,     0,     0,     0},
        {-0.00681675,-0.0162318,0.0794517,0.0155864,0.0248036,-0.150721,0.0423899,1.12658,     0,     0,     0,     0,     0,     0,     0,     0},
        {     0,     0,     0,     0,     0,     0,     0,     0,0.000317543,4.97647e-05,-0.00552168,0.00614938,     0,     0,     0,     0},
        {     0,     0,     0,     0,     0,     0,     0,     0,4.97647e-05,0.00401998,0.00255274,-0.0260365,     0,     0,     0,     0},
        {     0,     0,     0,     0,     0,     0,     0,     0,-0.00552168,0.00255274,0.102708,-0.136429,     0,     0,     0,     0},
        {     0,     0,     0,     0,     0,     0,     0,     0,0.00614938,-0.0260365,-0.136429,0.312044,     0,     0,     0,     0},
        {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,0.00558202,-0.0170784,-0.0436313,0.18835},
        {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,-0.0170784,0.107663,0.064924,-0.856043},
        {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,-0.0436313,0.064924,0.570347,-1.20359},
        {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,0.18835,-0.856043,-1.20359,7.84867 }};
    };

  } // rew   namespace
} // genie namespace

#endif
