//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightINukeExtra

\brief    Adds some additional functionality to reweighting hA2018.
          Originally developed for use by SBND and ICARUS.

\author   Gray Putnam <gputnam \at fnal.gov>
          Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Oct 1, 2025

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_INUKE_EXTRA_H_
#define _G_REWEIGHT_INUKE_EXTRA_H_

// GENIE/Reweight includes
#include "RwCalculators/GReWeightINuke.h"
#include "RwCalculators/GReWeightINukeParamsExtra.h"

namespace genie {

namespace rew   {

 class GReWeightINukeExtra : public GReWeightINuke
 {
 public:
   GReWeightINukeExtra();
  ~GReWeightINukeExtra();

   // Adjust some details of the base implementation to accommodate
   // new dials while keeping the originals untouched
   virtual bool IsHandled( GSyst_t syst ) const override final;
 };

} // rew
} // genie

#endif
