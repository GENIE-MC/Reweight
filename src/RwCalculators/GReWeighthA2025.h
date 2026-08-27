//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightINuke

\brief    Applies an approximate reweight from hA2018 to hA2025

\author   Mohamed Ismail <msi10 \at pitt.edu>
    			University of Pittsburgh

\created  Sep 2025

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef GREWEIGHTHA2025_H
#define GREWEIGHTHA2025_H

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"

namespace genie {

 class HAIntranuke2018;

namespace rew   {

 class GReWeighthA2025 : public GReWeightModel {

 public:

   GReWeighthA2025();
  ~GReWeighthA2025();

   // Implement the GReWeightI interface
   bool   AppliesTo      (const EventRecord& event) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
 };

} // rew
} // genie

#endif // GREWEIGHTHA2025_H
