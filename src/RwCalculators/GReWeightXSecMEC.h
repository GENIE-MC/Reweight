//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightXSecMEC

\brief    Model-independent tweak dials for reweighting MEC events

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Sep 11, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_MEC_H_
#define _G_REWEIGHT_NU_XSEC_MEC_H_

#include <map>
#include <string>

// GENIE includes
#include "Framework/Interaction/InteractionType.h"
#include "RwCalculators/GReWeightModel.h"
#include "RwFramework/GSyst.h"

namespace genie {

namespace rew   {

 class GReWeightXSecMEC : public GReWeightModel {
 public:

   GReWeightXSecMEC();
   GReWeightXSecMEC(std::string model, std::string type);
  ~GReWeightXSecMEC();

   // implement the GReWeightI interface
   bool   AppliesTo      (ScatteringType_t type, bool is_cc) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord& event);

 private:

   void Init(void);
   double CalcWeightNorm(const EventRecord& event);

   /// Simple struct containing tweak dial information for the
   /// normalization of one MEC interaction type (CC, NC, EM)
   struct NormMapEntry {
     NormMapEntry() : fTwkDial(0.), fNormDef(0.), fNormCurr(0.) {}
     NormMapEntry(double twk_dial, double def, double curr)
       : fTwkDial( twk_dial ), fNormDef( def ), fNormCurr( curr ) {}
     double fTwkDial;
     double fNormDef;
     double fNormCurr;
   };

   /// Map linking MEC interaction modes to normalization tweak dial values
   std::map<InteractionType_t, NormMapEntry> fNormMap;

   /// Lookup table linking MEC GSyst_t tweak dial enum labels to
   /// interaction modes
   static std::map<GSyst_t, InteractionType_t> fGSystToIntTypeMap;
 };

} // rew   namespace
} // genie namespace

#endif // _G_REWEIGHT_NU_XSEC_MEC_H_
