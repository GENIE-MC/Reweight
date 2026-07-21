//____________________________________________________________________________
/*!

  \class    genie::rew::GReWeightNuXSecCCQEZAFF

  \brief    Reweighting CCQE GENIE neutrino cross sections
            Z expansion axial form factor model

  \author   Liang Liu <liangliu \at fnal.gov>
  Fermi National Accelerator Laboratory

  \created  March 25, 2024

  \cpright  Copyright (c) 2003-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  */
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_CCQE_ZAFF_H_
#define _G_REWEIGHT_NU_XSEC_CCQE_ZAFF_H_

#include <string>

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"
#include "TRandom3.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"

class TFile;
class TNtupleD;

namespace genie {

  class XSecAlgorithmI;
  class XSecIntegratorI;
  class Registry;

  namespace rew   {

    class GReWeightNuXSecCCQEZAFF : public GReWeightModel
    {
      public:
        static const int kModeZExp             = 0;

        GReWeightNuXSecCCQEZAFF();
        GReWeightNuXSecCCQEZAFF(std::string model, std::string type);
        ~GReWeightNuXSecCCQEZAFF();

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
        // finite-difference step for XSecPartialDerivative, as a fraction of the
        // per-coefficient 1-sigma; built-in default 0.1, per-job override
        // (grwght1p --fd-delta). Must be > 0 (checked at first use).
        void SetFiniteDiffDelta (double d){ fFiniteDiffDelta = d; }

        // How the xsec 1-sigma is estimated from the coefficient covariance:
        //   kSigmaPropagation - analytic error propagation via finite-difference
        //                       derivatives (default; exact for the quadratic
        //                       coefficient dependence)
        //   kSigmaCholesky    - MC sampling: universes a' = a + L*z with L the
        //                       Cholesky factor of the covariance, sigma = sample
        //                       std dev of the recomputed xsec (grwghtnp-style)
        enum ESigmaEstimator { kSigmaPropagation = 0, kSigmaCholesky = 1 };
        void SetSigmaEstimator (ESigmaEstimator m){ fSigmaEstimator = m; }
        // number of universes for kSigmaCholesky (default 1000; must be >= 2,
        // checked at first use). Driver option: grwght1p --n-universes.
        void SetNUniverses     (int n){ fNUniverses = n; }

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
          std::vector<double> fZ_An;
        } fZExpParaDef, fZExpPara, fZExpParaTwkDial;
        // tweek dial and scale factor in propagation method
        double fZExpTwkDial;
        double fZExp_Scale;
        double fFiniteDiffDelta; ///< finite-difference step (fraction of 1-sigma) in XSecPartialDerivative

        // Two methods are provided to calculate the uncertainties of XSec
        // (see ESigmaEstimator): propagation of errors, or Cholesky-sampled
        // universes. fIsSinglePara/fIsAllPara placeholders retired in favour
        // of fSigmaEstimator.
        ESigmaEstimator fSigmaEstimator; ///< how GetOneSigma estimates sigma_xsec
        int             fNUniverses;     ///< universes for kSigmaCholesky
        TMatrixD        fLch;            ///< cached Cholesky factor L of error_mat
        double GetOneSigmaCholesky(const EventRecord & event);
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
