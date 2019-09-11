//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// GENIE/Generator includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/InteractionType.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Registry/Registry.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightXSecMEC.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

namespace {

  // Helper function to help us initialize the static std::map
  // owned by the GReWeightXSecMEC class
  std::map<GSyst_t, InteractionType_t> make_gsyst_to_inttype_map() {
    std::map<GSyst_t, InteractionType_t> temp_map;

    temp_map[ kXSecTwkDial_NormCCMEC ] = kIntWeakCC;
    temp_map[ kXSecTwkDial_NormNCMEC ] = kIntWeakNC;
    temp_map[ kXSecTwkDial_NormEMMEC ] = kIntEM;

    return temp_map;
  }

}

// Define the static std::map owned by the GReWeightXSecMEC class
// TODO: switch to something better for GENIE 4. C++11 makes initializing
// static std::map objects a lot easier. - S. Gardiner
std::map<GSyst_t, InteractionType_t> GReWeightXSecMEC::fGSystToIntTypeMap
  = make_gsyst_to_inttype_map();

//_______________________________________________________________________________________
GReWeightXSecMEC::GReWeightXSecMEC()
  : GReWeightModel("MEC")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightXSecMEC::GReWeightXSecMEC(std::string /*model*/, std::string /*type*/)
  : GReWeightModel("MEC")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightXSecMEC::~GReWeightXSecMEC()
{

}
//_______________________________________________________________________________________
bool GReWeightXSecMEC::IsHandled(GSyst_t syst) const
{
  // If we have an entry for this knob in the GSyst_t -> InteractionType_t
  // map, then this calculator can handle it. Otherwise, it can't.
  bool handle = fGSystToIntTypeMap.count( syst );
  return handle;
}
//_______________________________________________________________________________________
bool GReWeightXSecMEC::AppliesTo(ScatteringType_t type, bool /*is_cc*/) const
{
  // Weights can be calculated for CC, NC, and EM MEC events
  if ( type == kScMEC ) return true;
  return false;
}
//_______________________________________________________________________________________
void GReWeightXSecMEC::SetSystematic(GSyst_t syst, double twk_dial)
{
  if ( !this->IsHandled(syst) ) return;

  // We've already checked that there is an entry for the given knob in the map
  // during the previous call to IsHandled(). Therefore, just retrieve the
  // stored value this time.
  InteractionType_t type = fGSystToIntTypeMap.at( syst );

  // Store the new tweak dial value in the entry for the interaction type of
  // interest
  fNormMap.at(type).fTwkDial = twk_dial;
}
//_______________________________________________________________________________________
void GReWeightXSecMEC::Reset(void)
{
  // Reset all of the normalization tweak dials to their defaults
  std::map<InteractionType_t, NormMapEntry>::iterator it = fNormMap.begin();
  std::map<InteractionType_t, NormMapEntry>::iterator end = fNormMap.end();
  while ( it != end ) {
    it->second.fTwkDial = 0.;
    it->second.fNormCurr = it->second.fNormDef;
    ++it;
  }

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightXSecMEC::Reconfigure(void)
{
  GSystUncertainty* fracerr = GSystUncertainty::Instance();

  // Loop over all of the normalization tweak dials to update their current
  // values
  std::map<GSyst_t, InteractionType_t>::const_iterator it = fGSystToIntTypeMap.cbegin();
  std::map<GSyst_t, InteractionType_t>::const_iterator end = fGSystToIntTypeMap.cend();
  while ( it != end ) {

    GSyst_t syst = it->first;
    InteractionType_t type = it->second;

    // Note: this assumes that the error is symmetric.
    // TODO: consider changing this to handle asymmetric errors on the normalization
    double frac_err_norm = fracerr->OneSigmaErr( syst );

    NormMapEntry& entry = fNormMap.at( type );
    entry.fNormCurr = std::max(0.,
      entry.fNormDef * (1. + entry.fTwkDial * frac_err_norm));

    ++it;
  }

}
//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeight(const genie::EventRecord& event)
{
  bool is_mec = event.Summary()->ProcInfo().IsMEC();
  if ( !is_mec ) return 1.;

  double weight = this->CalcWeightNorm(event);
  return weight;
}
//_______________________________________________________________________________________
void GReWeightXSecMEC::Init(void) {
  // Set the default normalization for each interaction type (tweak dial = 0
  // corresponds to a normalization factor of 1)
  std::map<GSyst_t, InteractionType_t>::const_iterator it = fGSystToIntTypeMap.cbegin();
  std::map<GSyst_t, InteractionType_t>::const_iterator end = fGSystToIntTypeMap.cend();

  while ( it != end ) {
    fNormMap[ it->second ] = NormMapEntry(0., 1., 1.);
    ++it;
  }

}
//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeightNorm(const genie::EventRecord& event)
{
  InteractionType_t type = event.Summary()->ProcInfo().InteractionTypeId();

  // Find the tweak dial information for the current event's interaction type.
  // If a match isn't found, then just return a weight of unity.
  std::map<InteractionType_t, NormMapEntry>::const_iterator it = fNormMap.find( type );
  if ( it == fNormMap.cend() ) {
    LOG("ReW", pWARN) << "Unrecognized MEC event encountered in"
      << " GReWeightXSecMEC::CalcWeightNorm()";
    return 1.;
  }

  // If the tweak dial is set to zero (or is really small) then just return a
  // weight of unity
  double twk_dial = it->second.fTwkDial;
  bool tweaked = ( std::abs(twk_dial) > controls::kASmallNum );
  if ( !tweaked ) return 1.;

  double weight = it->second.fNormCurr;

  return weight;
}
//_______________________________________________________________________________________
