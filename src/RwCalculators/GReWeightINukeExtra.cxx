//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

// GENIE/Generator includes
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Messenger/Messenger.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightINukeExtra.h"
#include "RwCalculators/GReWeightINukeParamsExtra.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightINukeExtra::GReWeightINukeExtra() : GReWeightINuke()
{
  // Switch to a different name to avoid collisions with an instance
  // of the base class
  fName = "IntraNukeExtra";

  // Switch to using the derived "params" object
  fINukeRwParams = std::make_shared< GReWeightINukeParamsExtra >();
}
//_______________________________________________________________________________________
GReWeightINukeExtra::~GReWeightINukeExtra()
{
}
//_______________________________________________________________________________________
bool GReWeightINukeExtra::IsHandled(GSyst_t syst) const {
  bool handle;

  switch (syst) {
  // Nucleon MFP variations in kinetic energy bins
  case (kINukeTwkDial_MFPLoE_N):
  case (kINukeTwkDial_MFPM1E_N):
  case (kINukeTwkDial_MFPM2E_N):
  case (kINukeTwkDial_MFPHiE_N):

  // Map hA2018 nucleon fates to either G4 or INCL
  case (kINukeTwkDial_G4_N):
  case (kINukeTwkDial_INCL_N):

  // Map hA2018 nucleon fates to G4 in kinetic energy bins
  case (kINukeTwkDial_G4LoE_N):
  case (kINukeTwkDial_G4M1E_N):
  case (kINukeTwkDial_G4M2E_N):
  case (kINukeTwkDial_G4HiE_N):

  // Map hA2018 nucleon fates to INCL in kinetic energy bins
  case (kINukeTwkDial_INCLLoE_N):
  case (kINukeTwkDial_INCLM1E_N):
  case (kINukeTwkDial_INCLM2E_N):
  case (kINukeTwkDial_INCLHiE_N):
    handle = true;
    break;

  default:
    handle = false;
  }

  return handle;
}
