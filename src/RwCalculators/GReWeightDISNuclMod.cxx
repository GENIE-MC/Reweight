//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

// GENIE/Generator includes
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightDISNuclMod.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightDISNuclMod::GReWeightDISNuclMod() :
GReWeightModel("DISNuclMod")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightDISNuclMod::~GReWeightDISNuclMod()
{

}
//_______________________________________________________________________________________
bool GReWeightDISNuclMod::IsHandled (GSyst_t syst) const
{
  switch(syst) {
    case(kXSecTwkDial_DISNuclMod) :
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}
//_______________________________________________________________________________________
bool GReWeightDISNuclMod::AppliesTo (ScatteringType_t type, bool /*is_cc*/) const
{
  if (type==kScDeepInelastic) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightDISNuclMod::SetSystematic(GSyst_t syst, double twk_dial)
{
  switch(syst) {
    case(kXSecTwkDial_DISNuclMod) :
      fNuclModTwkDial = twk_dial;
      break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightDISNuclMod::Reset(void)
{
  fNuclModTwkDial = 0.;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightDISNuclMod::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightDISNuclMod::CalcWeight(const EventRecord & /*event*/)
{
  if ( fNuclModTwkDial != 0.0 ) {
    // fail hard if this reweighting is actually attempted
    LOG("ReW",pFATAL) << "Not implemented.";
    exit(-1);
  }
  return 1.;
}
//_______________________________________________________________________________________
void GReWeightDISNuclMod::Init(void)
{
  fNuclModTwkDial = 0.;

}
//_______________________________________________________________________________________
