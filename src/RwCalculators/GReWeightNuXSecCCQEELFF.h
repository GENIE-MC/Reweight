//____________________________________________________________________________
/*!

  \class    genie::rew::GReWeightNuXSecCCQEELFFELFF

  \brief    Reweighting CCQE GENIE neutrino cross sections
            Z expansion vector form factor model

  \author   Liang Liu <liangliu \at fnal.gov>
  Fermi National Accelerator Laboratory

  \created  March 25, 2024

  \cpright  Copyright (c) 2003-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  */
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_CCQE_ELFF_H_
#define _G_REWEIGHT_NU_XSEC_CCQE_ELFF_H_

#include <string>

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"
#include "TRandom3.h"
#include "TMatrixDSym.h"

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
          std::vector<double> fZ_APn;
          std::vector<double> fZ_BPn;
          std::vector<double> fZ_ANn;
          std::vector<double> fZ_BNn;
        } fZExpParaDef, fZExpPara, fZExpParaTwkDial;
        // tweek dial and scale factor in propagation method
        double fZExpTwkDial;
        double fZExp_Scale;

        // Two methods are provided to calculate the uncertainties of XSec
        // 1. propagation of errors: it is based on grwght1p
        // 2. Cholesky decomposition: it is based on grwghtnp
        bool fIsSinglePara;     // it will be used for Cholesky decomposition
        bool fIsAllPara;        // flag of propagation method
        std::vector<double> A_f;

        // List of the uncertainties of parameters from Kaushik
        // ap1, ap2, ap3, ap4,
        // bp1, bp2, bp3, bp4,
        // an1, an2, an3, an4,
        // bn1, bn2, bn3, bn4
        std::vector<double> errors;
        TMatrixDSym error_mat;
    };

  } // rew   namespace
} // genie namespace

#endif
