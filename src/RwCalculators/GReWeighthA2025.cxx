//____________________________________________________________________________
/*
 Copyright (c) 2003-2026, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Mohamed Ismail <msi10 \at pitt.edu>
          University of Pittsburgh

*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>

#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TVector.h>

// GENIE/Generator includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/HadronTransport/HAIntranuke2018.h"
#include "Physics/HadronTransport/Intranuke2018.h"
#include "Physics/HadronTransport/INukeHadroData2018.h"
#include "Physics/HadronTransport/INukeHadroFates2018.h"
#include "Physics/HadronTransport/INukeUtils2018.h"

#include "Physics/HadronTransport/INukeHadroData2025.h"
#include "Physics/HadronTransport/INukeHadroFates2025.h"
#include "Physics/HadronTransport/INukeUtils2025.h"
#include "Physics/HadronTransport/HAIntranuke2025.h"
#include "Physics/HadronTransport/Intranuke2025.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeighthA2025.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystUncertainty.h"

//______________________________________________________________________________
genie::rew::GReWeighthA2025::GReWeighthA2025()
  : GReWeightModel("IntraNuke2018to2025")
{
  // Look up the FSI model for the current tune. Also check whether FSIs are
  // actually enabled.
  AlgConfigPool* conf_pool = AlgConfigPool::Instance();
  Registry* gpl = conf_pool->GlobalParameterList();
  RgAlg fsi_alg = gpl->GetAlg( "HadronTransp-Model" );
  bool fsi_enabled = gpl->GetBool( "HadronTransp-Enable" );

  if ( !fsi_enabled ) {
    LOG( "ReW", pFATAL ) << "FSIs are not enabled for the current tune."
      << " Refusing to reweight FSIs.";
    std::exit( 1 );
  }

  if ( fsi_alg.name != "genie::HAIntranuke2018" ) {
    LOG( "ReW", pFATAL ) << "Reweighting events produced with the FSI model "
      << fsi_alg << " is not currently supported.";
    std::exit( 1 );
  }

}
//______________________________________________________________________________
genie::rew::GReWeighthA2025::~GReWeighthA2025()
{
}
//______________________________________________________________________________
bool genie::rew::GReWeighthA2025::IsHandled(GSyst_t syst) const
{
  if ( syst == kINukehA2025_cex ) return true;
  return false;
}
//______________________________________________________________________________
bool genie::rew::GReWeighthA2025::AppliesTo( const genie::EventRecord& evrec )
  const
{
  auto type = evrec.Summary()->ProcInfo().ScatteringTypeId();
  switch (type) {
    case kScCoherentProduction:
    case kScDiffractive:
    case kScNuElectronElastic:
    case kScAMNuGamma:
    case kScCoherentElastic:
      return false;
    default:
      return true;
  }
}
//______________________________________________________________________________
void genie::rew::GReWeighthA2025::SetSystematic(GSyst_t /*syst*/,
  double /*val*/)
{
 // if(this->IsHandled(syst)) {
 //    fINukeRwParams.SetTwkDial(syst, val);
 // }
}
//______________________________________________________________________________
void genie::rew::GReWeighthA2025::Reset(void)
{
 // fINukeRwParams.Reset();
 // this->Reconfigure();
}
//_______________________________________________________________________________________
void genie::rew::GReWeighthA2025::Reconfigure(void)
{
 // fINukeRwParams.Reconfigure();
}
//______________________________________________________________________________
double genie::rew::GReWeighthA2025::CalcWeight(const EventRecord & event)
{
  // Non-trivial weights can only be returned for a complex nuclear target
  GHepParticle* tgt = event.TargetNucleus();
  if ( !tgt ) return 1.0;
  double A = tgt->A();
  double Z = tgt->Z();
  if ( A <= 1 ) return 1.0;
  if ( Z <= 1 ) return 1.0;

  // Get the pre-FSI nuclear remnant. The A and Z values for this particle
  // (distinct from both the post-FSI remnant and the target nucleus)
  // are used to compute mean free paths in Intranuke2018.
  GHepParticle* pre_fsi_remnant = 0;
  TObjArrayIter piter( &event );
  // First loop over all particles in the event record and find the
  // final-state nuclear remnant (i.e., the post-FSI nuclear remnant)
  while( GHepParticle* p = dynamic_cast<GHepParticle*>(piter.Next()) ) {
    if ( p->Status() == genie::kIStFinalStateNuclearRemnant ) {
      // The pre-FSI nuclear remnant is set to be the first mother for
      // the post-FSI one by Intranuke2018
      int mother_idx = p->FirstMother();
      pre_fsi_remnant = event.Particle( mother_idx );
      break;
    }
  }

  // Skip this event if we couldn't find a pre-FSI nuclear remnant that
  // is an ion (something went wrong in the search above)
  bool pre_fsi_remnant_ok = false;
  if ( pre_fsi_remnant ) {
   pre_fsi_remnant_ok = genie::pdg::IsIon( pre_fsi_remnant->Pdg() );
  }
  if ( !pre_fsi_remnant_ok ) {
    LOG( "ReW", pWARN ) << "Could not find a suitable pre-FSI remnant"
      << " in GReWeighthA2025::CalcWeight()";
    return 1.;
  }

  // Store the nucleon and proton numbers for the pre-FSI remnant.
  int remnA = pre_fsi_remnant->A();
  int remnZ = pre_fsi_remnant->Z();
  LOG( "ReW", pDEBUG ) << "Found pre-FSI remnant with A = " << remnA
    << ", Z = " << remnZ << ". Target had A = " << A << ", Z = " << Z;

  INukeHadroData2018* hd2018 = INukeHadroData2018::Instance();
  INukeHadroData2025* hd2025 = INukeHadroData2025::Instance();

  double event_weight  = 1.0;

  // Loop over GHepParticle entries and calculate a per-particle weight for
  // hadrons in the nucleus. The product of these particle weights is the
  // overall event weight.
  int ip = -1;
  GHepParticle* p = 0;
  TIter event_iter( &event );
  while ( (p = dynamic_cast< GHepParticle* >(event_iter.Next())) ) {
    ++ip;

    // Skip particles other than pions (hA2025 has the same fate fractions
    // as hA2018 for nucleons)
    int pdgc = p->Pdg();
    bool is_pion = pdg::IsPion( pdgc );
    if ( !is_pion ) continue;

    // This weight calculator only accounts for differences in fate fractions,
    // so skip particles that did not interact during the cascade
    auto fsi_code = static_cast< INukeFateHA_t >( p->RescatterCode() );
    bool interacted = ( fsi_code != kIHAFtNoInteraction );
    if ( !interacted ) continue;

    // Skip particles with an unhandled fate (this shouldn't occur, it's a
    // fallback for robustness)
    if ( fsi_code != kIHAFtAbs && fsi_code != kIHAFtInelas
      && fsi_code != kIHAFtCEx && fsi_code != kIHAFtPiProd )
    {
	    LOG("ReW", pWARN) << "Unhandled hadron fate "
        << INukeHadroFates::AsString( static_cast< INukeFateHA_t >(fsi_code) )
        << " encountered in GReWeighthA2025::CalcWeight() for particle"
        << " with index " << ip << " and PDG code = " << pdgc;
      continue;
    }

    // Retrieve the kinetic energy
	  double ke = p->KinE();

    // Convert the kinetic energy to MeV (expected units for the calls to
    // FracADep below)
	  double ke_in_MeV = ke / units::MeV;

    double fate_frac2018 = hd2018->FracADep( pdgc, fsi_code, ke_in_MeV, remnA );
    double fate_frac2025 = hd2025->FracADep( pdgc, fsi_code, ke_in_MeV, remnA );

	  LOG("ReW", pDEBUG)
      << "GReWeighthA2025 reweighted hadron at position = " << ip
      << " with PDG code = " << pdgc
      << ", FSI code = "  << fsi_code
      << ", KE= "  << ke
      << ", A= "  << remnA
      << " (" << INukeHadroFates::AsString((INukeFateHA_t)fsi_code) << ") :"
      << " frac2018 = "  << fate_frac2018
      <<", frac2025 = " << fate_frac2025;

    double frac_ratio = 1.0;
    if ( fate_frac2018 != 0.0 ) {
	    frac_ratio = fate_frac2025 / fate_frac2018;
	  }

    // Update the current event weight
    event_weight *= frac_ratio;

  } // particle loop

  return event_weight;
}
