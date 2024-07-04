//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London
*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

// GENIE/Generator includes
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightNuXSecNC.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecNC::GReWeightNuXSecNC()  :
GReWeightModel("NC")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecNC::~GReWeightNuXSecNC()
{

}
//_______________________________________________________________________________________
bool GReWeightNuXSecNC::IsHandled(GSyst_t syst) const
{
   switch(syst) {
     case ( kXSecTwkDial_NC ) :
       return true;
       break;
     default:
       return false;
       break;
   }
   return false;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecNC::AppliesTo(const EventRecord & event) const
{
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if (!is_cc) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecNC::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

   switch(syst) {
     case ( kXSecTwkDial_NC ) :
       fNCTwkDial = twk_dial;
       break;
     default:
       return;
       break;
   }
}
//_______________________________________________________________________________________
void GReWeightNuXSecNC::Reset(void)
{
  fNCTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightNuXSecNC::Reconfigure(void)
{
  fNCTwkDial = TMath::Max(0., fNCTwkDial);
}
//_______________________________________________________________________________________
double GReWeightNuXSecNC::CalcWeight(const genie::EventRecord & event)
{
  bool tweaked = (TMath::Abs(fNCTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  Interaction * interaction = event.Summary();
  bool is_nc  = interaction->ProcInfo().IsWeakNC();
  if(!is_nc) return 1.;

  bool is_qel = event.Summary()->ProcInfo().IsQuasiElastic();
  if(is_qel && !fRewQE) return 1.;

  bool is_res = event.Summary()->ProcInfo().IsResonant();
  if(is_res && !fRewRES) return 1.;

  bool is_dis = event.Summary()->ProcInfo().IsDeepInelastic();
  if(is_dis && !fRewDIS) return 1.;

  int nupdg = interaction->InitState().ProbePdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  return fNCTwkDial;
}
//_______________________________________________________________________________________
void GReWeightNuXSecNC::Init(void)
{
  fNCTwkDial = 0.;

  this->RewNue       (true);
  this->RewNuebar    (true);
  this->RewNumu      (true);
  this->RewNumubar   (true);

  this->RewQE  (true );
  this->RewRES (false); // assume GReWeightNuXSecNCRES is going to be your 1st choice
  this->RewDIS (false); // assume GReWeightNuXSecDIS   is going to be your 1st choice

}
//_______________________________________________________________________________________
