//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London
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

// GENIE/Reweight includes
#include "RwCalculators/GReWeightINuke.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightINuke::GReWeightINuke() :
GReWeightModel("IntraNuke")
{
#ifdef _G_REWEIGHT_INUKE_DEBUG_NTP_
  fTestFile = new TFile("./intranuke_reweight_test.root","recreate");
  fTestNtp  = new TNtuple("testntp","","pdg:E:mfp_twk_dial:d:d_mfp:fate:interact:w_mfp:w_fate");
#endif

  // Look up the FSI model for the current tune. Also check whether FSIs are
  // actually enabled.
  AlgConfigPool* conf_pool = AlgConfigPool::Instance();
  Registry* gpl = conf_pool->GlobalParameterList();
  RgAlg fsi_alg = gpl->GetAlg( "HadronTransp-Model" );
  bool fsi_enabled = gpl->GetBool( "HadronTransp-Enable" );

  if ( !fsi_enabled ) {
    LOG( "ReW", pERROR ) << "FSIs are not enabled for the current tune."
      << " Refusing to reweight FSIs.";
    std::exit( 1 );
  }

  AlgId id( fsi_alg );

  AlgFactory* algf = AlgFactory::Instance();

  Algorithm* alg = algf->AdoptAlgorithm( id );
  fFSIModel = dynamic_cast< HAIntranuke2018* >( alg );

  if ( !fFSIModel ) {
    LOG( "ReW", pERROR ) << "Reweighting events produced with the FSI model "
      << fsi_alg << " is not currently supported.";
    std::exit( 1 );
  }

  fFSIModel->AdoptSubstructure();

  fINukeRwParams = std::make_shared< GReWeightINukeParams >();
}
//_______________________________________________________________________________________
GReWeightINuke::~GReWeightINuke()
{
#ifdef _G_REWEIGHT_INUKE_DEBUG_NTP_
  assert(fTestFile);
  assert(fTestNtp);
  fTestFile->cd();
  fTestNtp->Write();
  fTestFile->Close();
  delete fTestFile;
  //delete fTestNtp;
#endif
}
//_______________________________________________________________________________________
bool GReWeightINuke::IsHandled(GSyst_t syst) const
{
   bool handle;

   switch(syst) {
     case ( kINukeTwkDial_MFP_pi      ) :
     case ( kINukeTwkDial_MFP_N       ) :
     case ( kINukeTwkDial_FrCEx_pi    ) :
       //     case ( kINukeTwkDial_FrElas_pi   ) :
     case ( kINukeTwkDial_FrInel_pi   ) :
     case ( kINukeTwkDial_FrAbs_pi    ) :
     case ( kINukeTwkDial_FrPiProd_pi ) :
     case ( kINukeTwkDial_FrCEx_N     ) :
       //     case ( kINukeTwkDial_FrElas_N    ) :
     case ( kINukeTwkDial_FrInel_N    ) :
     case ( kINukeTwkDial_FrAbs_N     ) :
     case ( kINukeTwkDial_FrPiProd_N  ) :
          handle = true;
          break;

     default:
          handle = false;
   }

   return handle;
}
//_______________________________________________________________________________________
bool GReWeightINuke::AppliesTo(const EventRecord & event) const
{
  auto type = event.Summary()->ProcInfo().ScatteringTypeId();
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
//_______________________________________________________________________________________
void GReWeightINuke::SetSystematic(GSyst_t syst, double val)
{
  if(this->IsHandled(syst)) {
     fINukeRwParams->SetTwkDial(syst, val);
  }
}
//_______________________________________________________________________________________
void GReWeightINuke::Reset(void)
{
  fINukeRwParams->Reset();
  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightINuke::Reconfigure(void)
{
  fINukeRwParams->Reconfigure();
}
//_______________________________________________________________________________________
double GReWeightINuke::CalcWeight(const EventRecord & event)
{
  // get the atomic mass number for the hit nucleus
  GHepParticle * tgt = event.TargetNucleus();
  if (!tgt) return 1.0;
  double A = tgt->A();
  double Z = tgt->Z();
  if (A<=1) return 1.0;
  if (Z<=1) return 1.0;

  fINukeRwParams->SetTargetA( A );

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
      << " in GReWeightINuke::CalcWeight()";
    return 1.;
  }

  // Store the nucleon and proton numbers for the pre-FSI remnant.
  // These will be used for mean free path calculations in
  // genie::utils::rew::MeanFreePathWeight
  int remnA = pre_fsi_remnant->A();
  int remnZ = pre_fsi_remnant->Z();
  LOG( "ReW", pDEBUG ) << "Found pre-FSI remnant with A = " << remnA
    << ", Z = " << remnZ << ". Target had A = " << A << ", Z = " << Z;

  // Tell the FSI model about the remnant A and Z. Normally these values are
  // set at the start of Intranuke2018::TransportHadrons.
  fFSIModel->SetRemnA( remnA );
  fFSIModel->SetRemnZ( remnZ );

  double event_weight  = 1.0;

  // Loop over stdhep entries and only calculate weights for particles.
  // All particles that are not hadrons generated inside the nucleus are given weights of 1.0
  int ip=-1;
  GHepParticle * p = 0;
  TIter event_iter(&event);
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
     ip++;

     // Skip particles not rescattered by the actual hadron transport code
     int  pdgc       = p->Pdg();
     bool is_pion    = pdg::IsPion   (pdgc);
     bool is_nucleon = pdg::IsNucleon(pdgc);
     bool is_kaon = pdg::IsKaon( pdgc );
     if(!is_pion && !is_nucleon && !is_kaon)
     {
        continue;
     }

     // Skip particles with code other than 'hadron in the nucleus'
     GHepStatus_t ist  = p->Status();
     if(ist != kIStHadronInTheNucleus)
     {
        continue;
     }

     // Calculate how the nuclear (A, Z) values changed due to FSIs experienced
     // by the current particle
     int deltaA, deltaZ;
     this->CalcDeltaAZ( event, *p, deltaA, deltaZ );

     // Kaon FSIs can't currently be reweighted. Just update (A, Z) based on
     // the particle's daughters and move on.
     if ( is_kaon ) {
       this->UpdateRemnantAZ( deltaA, deltaZ );
       continue;
     }

     // Determine the interaction type for current hadron in nucleus, if any
     int fsi_code = p->RescatterCode();
     LOG("ReW", pDEBUG)
        << "Attempting to reweight hadron at position = " << ip
        << " with PDG code = " << pdgc
        << " and FSI code = "  << fsi_code
        << " (" << INukeHadroFates::AsString((INukeFateHA_t)fsi_code) << ")";
     if(fsi_code == -1 || fsi_code == (int)kIHAFtUndefined) {
       LOG("ReW", pFATAL) << "INTRANUKE didn't set a valid rescattering code for event in position: " << ip;
       LOG("ReW", pFATAL) << "Here is the problematic event:";
       LOG("ReW", pFATAL) << event;
       exit(1);
     }
     bool escaped    = (fsi_code == (int)kIHAFtNoInteraction);
     bool interacted = !escaped;

     // Get 4-momentum and 4-position
     TLorentzVector x4 (p->Vx(), p->Vy(), p->Vz(), 0.    );
     TLorentzVector p4 (p->Px(), p->Py(), p->Pz(), p->E());

     // Init current hadron weights
     double w_mfp  = 1.0;
     double w_fate = 1.0;

     // Check which weights need to be calculated (only if relevant params were tweaked)
     bool calc_w_mfp  = fINukeRwParams->MeanFreePathParams(pdgc)->IsTweaked();
     bool calc_w_fate = fINukeRwParams->FateParams(pdgc)->IsTweaked();

     // Compute weight to account for changes in the total rescattering probability
     double mfp_scale_factor = 1.;
     if(calc_w_mfp)
     {
        mfp_scale_factor = fINukeRwParams->MeanFreePathParams(pdgc)->ScaleFactor(p4);
        w_mfp = utils::rew::MeanFreePathWeight( pdgc, x4, p4, A, Z,
          mfp_scale_factor, interacted, *fFSIModel );
     } // calculate mfp weight?

     // Compute weight to account for changes in relative fractions of reaction channels
     if(calc_w_fate && interacted)
     {
        double fate_fraction_scale_factor =
             fINukeRwParams->FateParams(pdgc)->ScaleFactor(
                  GSyst::INukeFate2GSyst((INukeFateHA_t)fsi_code,pdgc), p4);
        w_fate = fate_fraction_scale_factor;
     }

     // Calculate the current hadron weight
     double hadron_weight = w_mfp * w_fate;

     LOG("ReW", pNOTICE)
        << "Reweighted hadron at position = " << ip
        << " with PDG code = " << pdgc
        << ", FSI code = "  << fsi_code
        << " (" << INukeHadroFates::AsString((INukeFateHA_t)fsi_code) << ") :"
        << " w_mfp = "  << w_mfp
        <<", w_fate = " << w_fate;

     // Debug info
#ifdef _G_REWEIGHT_INUKE_DEBUG_NTP_
     // TODO: Fix the debugging functions here for the hA2018 updates
     double d        = utils::intranuke2018::Dist2Exit(x4,p4,A);
     double d_mfp    = utils::intranuke2018::Dist2ExitMFP(pdgc,x4,p4,A,Z);
     double Eh       = p->E();
     double iflag    = (interacted) ? 1 : 0;
     fTestNtp->Fill(pdgc, Eh, mfp_scale_factor, d, d_mfp, fsi_code, iflag, w_mfp, w_fate);
#endif

     // Update the current event weight
     event_weight *= hadron_weight;

     // Now that we've handled the current hadron, update the nuclear (A, Z)
     // values to what they were when FSIs were simulated for the next hadron.
     this->UpdateRemnantAZ( deltaA, deltaZ );

  }//particle loop

  return event_weight;
}
//_______________________________________________________________________________________
void GReWeightINuke::CalcDeltaAZ( const EventRecord& event,
  const GHepParticle& p, int& deltaA, int& deltaZ )
{
  // Compute the total nucleon number and electric charge (in units of the up
  // quark charge) for all "stable final state" daughters of the current
  // GHepParticle
  int myA = 0;
  int myQ = 0;
  genie::utils::rew::TallyAQ( event, p, myA, myQ );

  deltaA = genie::utils::rew::GetParticleA( p.Pdg() ) - myA;
  // Convert to units of the elementary charge
  deltaZ = ( p.Charge() - myQ ) / 3;

  // Deal with apparent charge conservation issues in the absorption fate. Note
  // that there are also some failure modes for "too few nucleon" cases which
  // appear to yield irrecoverable changes to (A, Z). I don't attempt to
  // address those here since the information is simply lost. The good news is
  // that it shouldn't be a big deal for FSI reweighting.
  // -- S. Gardiner, 19 June 2021
  if ( p.RescatterCode() == genie::kIHAFtAbs ) {

    if ( p.Pdg() == genie::kPdgPiM ) {
      int daught = p.FirstDaughter();
      int d_pdg = event.Particle( daught )->Pdg();
      // If the first daughter is a nucleon cluster, then a multinucleon
      // absorption reaction was simulated which doesn't suffer from the
      // charge conservation issue. If the first daughter is a pi minus,
      // then the absorption simulation failed and terminated early.
      if ( d_pdg != genie::kPdgCompNuclCluster
        && d_pdg != genie::kPdgPiM ) deltaZ--;
    } // Pi-

    else if ( p.Pdg() == genie::kPdgKP ) {
      int daught = p.FirstDaughter();
      if ( event.Particle(daught)->Pdg()
        == genie::kPdgCompNuclCluster )
      {
        deltaZ++;
      }
    } // K+
  } // absorption fate
}
//_______________________________________________________________________________________
void GReWeightINuke::UpdateRemnantAZ( int deltaA, int deltaZ ) {
  int remnA = fFSIModel->GetRemnA();
  int remnZ = fFSIModel->GetRemnZ();

  remnA += deltaA;
  remnZ += deltaZ;

  fFSIModel->SetRemnA( remnA );
  fFSIModel->SetRemnZ( remnZ );
}
