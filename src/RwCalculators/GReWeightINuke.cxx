//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

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
#include "Framework/Interaction/ScatteringType.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/PhysUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/DeepInelastic/EventGen/DISHadronicSystemGenerator.h"
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

  // Also get the formation zone settings from the DISHadronicSystemGenerator
  const genie::Algorithm* dis_hs_alg = algf->GetAlgorithm(
    "genie::DISHadronicSystemGenerator", "Default" );
  const genie::Registry& dis_hs_reg = dis_hs_alg->GetConfig();

  fCt0_pion = dis_hs_reg.GetDouble( "FZONE-ct0pion" );
  fCt0_nucleon = dis_hs_reg.GetDouble( "FZONE-ct0nucleon" );
  fKpt2 = dis_hs_reg.GetDouble( "FZONE-KPt2" );
  fNucl_R0 = dis_hs_reg.GetDouble( "NUCL-R0" );
  fNucl_NR = dis_hs_reg.GetDouble( "NUCL-NR" );
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
bool GReWeightINuke::AppliesTo(ScatteringType_t type, bool /*is_cc*/) const
{
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
     fINukeRwParams.SetTwkDial(syst, val);
  }
}
//_______________________________________________________________________________________
void GReWeightINuke::Reset(void)
{
  fINukeRwParams.Reset();
  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightINuke::Reconfigure(void)
{
  fINukeRwParams.Reconfigure();
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

  fINukeRwParams.SetTargetA( A );

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
     if(!is_pion && !is_nucleon)
     {
        continue;
     }

     // Skip particles with code other than 'hadron in the nucleus'
     GHepStatus_t ist  = p->Status();
     if(ist != kIStHadronInTheNucleus)
     {
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
     TLorentzVector p4 (p->Px(), p->Py(), p->Pz(), p->E());

     // NOTE: Unlike Intranuke, the new Intranuke2018 shifts the
     // position of the original particle with the kIStHadronInTheNucleus
     // status code. This screws up the mean free path reweighting.
     // See line 315 of Generator/src/Physics/HadronTransport/Intranuke2018.cxx.
     // As a stopgap solution, reconstruct the starting position of the
     // particle using other information from the event. For events other than
     // DIS, this is easy: just use the 4-position of the first mother of the
     // current particle. For DIS, apply the appropriate formation zone step
     // from the location of the initial struck nucleon (see the source code
     // for DISHadronicSystemGenerator). Plan on fixing Generator for
     // GENIE v3.2 to make this section of the reweighting code obsolete.
     // -- S. Gardiner, 11 June 2021
     TLorentzVector x4;
     genie::ScatteringType_t scatter_type = event.Summary()
      ->ProcInfo().ScatteringTypeId();
     if ( scatter_type == genie::kScDeepInelastic ) {
       // Compute the usual tracking radius for FSIs
       // TODO: add a check that this is consistent with the one used
       // by Intranuke2018
       double tracking_R = fNucl_NR * fNucl_R0 * std::pow( A, 1./3. );

       // Get the 3-momentum of the pre-fragmentation hadronic system
       genie::GHepParticle* had_syst = event.FinalStateHadronicSystem();
       TVector3 p3_had_syst = had_syst->P4()->Vect();

       double mass = p->Mass(); // On-shell mass of the current particle
       int pdg = p->Pdg();
       double ct0 = genie::pdg::IsNucleon( pdg ) ? fCt0_nucleon : fCt0_pion;

       // Compute the formation zone that was used for this particle
       double fzone = genie::utils::phys::FormationZone( mass, p4,
         p3_had_syst, ct0, fKpt2 );

       // Step one formation zone length away from the vertex in the direction
       // of the particle's 3-momentum.
       TVector3 step3 = p4.Vect().Unit();
       step3.SetMag( fzone );
       TLorentzVector step4( step3, 0. );

       const TLorentzVector& vtx4 = *event.HitNucleon()->X4();

       x4 = vtx4 + step4;

       // See hard-coded stuff in DISHadronicSystemGenerator
       double r_max = tracking_R + 2.; // fm
       double r = x4.P();
       if ( r > r_max ) {
         double scale = r_max / r;
         x4 *= scale;
       }
     }
     else {
       GHepParticle* temp_mother = event.Particle( p->FirstMother() );
       x4 = *temp_mother->X4();
     }

     // Init current hadron weights
     double w_mfp  = 1.0;
     double w_fate = 1.0;

     // Check which weights need to be calculated (only if relevant params were tweaked)
     bool calc_w_mfp  = fINukeRwParams.MeanFreePathParams(pdgc)->IsTweaked();
     bool calc_w_fate = fINukeRwParams.FateParams(pdgc)->IsTweaked();

     // Compute weight to account for changes in the total rescattering probability
     double mfp_scale_factor = 1.;
     if(calc_w_mfp)
     {
        mfp_scale_factor = fINukeRwParams.MeanFreePathParams(pdgc)->ScaleFactor();
        w_mfp = utils::rew::MeanFreePathWeight( pdgc, x4, p4, A, Z,
          mfp_scale_factor, interacted, *fFSIModel );
     } // calculate mfp weight?

     // Compute weight to account for changes in relative fractions of reaction channels
     if(calc_w_fate && interacted)
     {
        double fate_fraction_scale_factor =
             fINukeRwParams.FateParams(pdgc)->ScaleFactor(
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

  }//particle loop

  return event_weight;
}
//_______________________________________________________________________________________
