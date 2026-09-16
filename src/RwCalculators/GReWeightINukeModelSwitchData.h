//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightINukeModelSwitchData

\created  Apr 4, 2025

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_INUKE_DATA_H_
#define _G_REWEIGHT_INUKE_DATA_H_

//#define _G_REWEIGHT_INUKE_DATA_DEBUG_NTP_

// GENIE/Reweight includes
#include "Framework/Numerical/Spline.h"

using namespace genie::rew;
using namespace genie;

namespace genie {

namespace rew   {
 typedef enum EModelSwitchType {
    kNoSwitch = 0,
    kRwINukeG4,
    kRwINukeINCL,
 } ModelSwitch_t;

 class GReWeightINukeModelSwitchData
 {
 public:
   GReWeightINukeModelSwitchData();
  ~GReWeightINukeModelSwitchData();

   static const GReWeightINukeModelSwitchData *Instance();

   bool IsHandled(GSyst_t syst) const;
   double FateFraction(ModelSwitch_t model, GSyst_t syst, double KE) const;

   const Spline *G4PA_Tot() const {return fG4PA_Tot;}
   const Spline *G4FracPA_Tot() const {return fG4FracPA_Tot;}
   const Spline *G4FracPA_Inel() const {return fG4FracPA_Inel;}
   const Spline *G4FracPA_CEx() const {return fG4FracPA_CEx;}
   const Spline *G4FracPA_Abs() const {return fG4FracPA_Abs;}
   const Spline *G4FracPA_PiPro() const {return fG4FracPA_PiPro;}

   const Spline *G4NA_Tot() const {return fG4NA_Tot;}
   const Spline *G4FracNA_Tot() const {return fG4FracNA_Tot;}
   const Spline *G4FracNA_Inel() const {return fG4FracNA_Inel;}
   const Spline *G4FracNA_CEx() const {return fG4FracNA_CEx;}
   const Spline *G4FracNA_Abs() const {return fG4FracNA_Abs;}
   const Spline *G4FracNA_PiPro() const {return fG4FracNA_PiPro;}

   const Spline *INCLPA_Tot() const {return fINCLPA_Tot;}
   const Spline *INCLFracPA_Tot() const {return fINCLFracPA_Tot;}
   const Spline *INCLFracPA_Inel() const {return fINCLFracPA_Inel;}
   const Spline *INCLFracPA_CEx() const {return fINCLFracPA_CEx;}
   const Spline *INCLFracPA_Abs() const {return fINCLFracPA_Abs;}
   const Spline *INCLFracPA_PiPro() const {return fINCLFracPA_PiPro;}

   const Spline *INCLNA_Tot() const {return fINCLNA_Tot;}
   const Spline *INCLFracNA_Tot() const {return fINCLFracNA_Tot;}
   const Spline *INCLFracNA_Inel() const {return fINCLFracNA_Inel;}
   const Spline *INCLFracNA_CEx() const {return fINCLFracNA_CEx;}
   const Spline *INCLFracNA_Abs() const {return fINCLFracNA_Abs;}
   const Spline *INCLFracNA_PiPro() const {return fINCLFracNA_PiPro;}

 private:
   Spline *fG4PA_Tot;
   Spline *fG4FracPA_Tot;
   Spline *fG4FracPA_Inel;
   Spline *fG4FracPA_CEx;
   Spline *fG4FracPA_Abs;
   Spline *fG4FracPA_PiPro;

   Spline *fG4NA_Tot;
   Spline *fG4FracNA_Tot;
   Spline *fG4FracNA_Inel;
   Spline *fG4FracNA_CEx;
   Spline *fG4FracNA_Abs;
   Spline *fG4FracNA_PiPro;

   Spline *fINCLPA_Tot;
   Spline *fINCLFracPA_Tot;
   Spline *fINCLFracPA_Inel;
   Spline *fINCLFracPA_CEx;
   Spline *fINCLFracPA_Abs;
   Spline *fINCLFracPA_PiPro;

   Spline *fINCLNA_Tot;
   Spline *fINCLFracNA_Tot;
   Spline *fINCLFracNA_Inel;
   Spline *fINCLFracNA_CEx;
   Spline *fINCLFracNA_Abs;
   Spline *fINCLFracNA_PiPro;

   // constants
   static constexpr double fMinKinEnergy = 0.;
   static constexpr double fMaxKinEnergy = 1000.;

 };

} // rew
} // genie

#endif
