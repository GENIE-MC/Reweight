//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightDISNuclMod

\brief    Reweighting the DIS nuclear modification model

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Apr 26, 2010

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_DISNUCLMOD_H_
#define _G_REWEIGHT_DISNUCLMOD_H_

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"

using namespace genie::rew;
using namespace genie;

namespace genie {
namespace rew   {

 class GReWeightDISNuclMod : public GReWeightModel
 {
 public:
   GReWeightDISNuclMod();
  ~GReWeightDISNuclMod();

   // implement the GReWeightI interface
   bool   AppliesTo      (const EventRecord & event) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

 private:
   void Init(void);

   double fNuclModTwkDial;
 };

} // rew
} // genie

#endif
