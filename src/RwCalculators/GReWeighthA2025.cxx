//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Mohamed Ismail <msi10@pitt.edu>
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

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeighthA2025::GReWeighthA2025() :
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
    LOG( "ReW", pFATAL ) << "FSIs are not enabled for the current tune."
      << " Refusing to reweight FSIs.";
    std::exit( 1 );
  }

  AlgId id( fsi_alg );

  AlgFactory* algf = AlgFactory::Instance();

  Algorithm* alg = algf->AdoptAlgorithm( id );
  fFSIModel = dynamic_cast< HAIntranuke2018* >( alg );

  if ( !fFSIModel ) {
    LOG( "ReW", pFATAL ) << "Reweighting events produced with the FSI model "
      << fsi_alg << " is not currently supported.";
    std::exit( 1 );
  }

  fFSIModel->AdoptSubstructure();
}
//_______________________________________________________________________________________
GReWeighthA2025::~GReWeighthA2025()
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
bool GReWeighthA2025::IsHandled(GSyst_t syst) const
{
   bool handle;

   switch(syst) {
     case ( kINukehA2025_cex ) :
          handle = true;
          break;

     default:
          handle = false;
   }

   return handle;
}
//_______________________________________________________________________________________
bool GReWeighthA2025::AppliesTo(ScatteringType_t type, bool /*is_cc*/) const
{
  if (type != kScCoherentProduction ) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________

void GReWeighthA2025::SetSystematic(GSyst_t syst, double val)
{
 // if(this->IsHandled(syst)) {
 //    fINukeRwParams.SetTwkDial(syst, val);
 // }
}
//_______________________________________________________________________________________
void GReWeighthA2025::Reset(void)
{
 // fINukeRwParams.Reset();
 // this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeighthA2025::Reconfigure(void)
{
 // fINukeRwParams.Reconfigure();
}

//_______________________________________________________________________________________
double GReWeighthA2025::CalcWeight(const EventRecord & event)
{
  // get the atomic mass number for the hit nucleus
  GHepParticle * tgt = event.TargetNucleus();
  if (!tgt) return 1.0;
  double A = tgt->A();
  double Z = tgt->Z();
  if (A<=1) return 1.0;
  if (Z<=1) return 1.0;

  //fINukeRwParams.SetTargetA( A );

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
  // These will be used for mean free path calculations in
  // genie::utils::rew::MeanFreePathWeight
  int remnA = pre_fsi_remnant->A();
  int remnZ = pre_fsi_remnant->Z();
  LOG( "ReW", pDEBUG ) << "Found pre-FSI remnant with A = " << remnA
    << ", Z = " << remnZ << ". Target had A = " << A << ", Z = " << Z;


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

     
     // Get 4-momentum and 4-position
     TLorentzVector x4 (p->Vx(), p->Vy(), p->Vz(), 0.    );
     TLorentzVector p4 (p->Px(), p->Py(), p->Pz(), p->E());
     
     int fsi_code = p->RescatterCode();
     bool escaped    = (fsi_code == (int)kIHAFtNoInteraction);
     bool interacted = !escaped;
     
     // new 
    double fate_frac2018 = 0.0;
    double fate_frac2025 = 0.0;
    if (is_pion  and fsi_code > 1 ){ 
      INukeHadroData2018 * hd2018 = INukeHadroData2018::Instance();
      INukeHadroData2025 * hd2025 = INukeHadroData2025::Instance();

	  // convert to MeV and
	  
	  double KE = p4.Energy() - p4.M(); // kinetic energy
	  double ke = KE / units::MeV;
	  
	  // get particle fate
	  auto fate_rescatter = p->RescatterCode();
	  
	  // get fate frac

	  if (fate_rescatter == kIHAFtAbs){
	   fate_frac2018 = hd2018->FracADep(pdgc, kIHAFtAbs, ke, remnA);
	   fate_frac2025 = hd2025->FracADep(pdgc, kIHAFtAbs, ke, remnA);	  
	  }
	  
	   else if (fate_rescatter == kIHAFtInelas){
	   fate_frac2018 = hd2018->FracADep(pdgc, kIHAFtInelas, ke, remnA);
	   fate_frac2025 = hd2025->FracADep(pdgc, kIHAFtInelas, ke, remnA);	  
	  }
	  
	 else if (fate_rescatter == kIHAFtCEx){
	   fate_frac2018 = hd2018->FracADep(pdgc, kIHAFtCEx, ke, remnA);
	   fate_frac2025 = hd2025->FracADep(pdgc, kIHAFtCEx, ke, remnA);	  
	  }
	  
	 else if (fate_rescatter == kIHAFtPiProd){
	   fate_frac2018 = hd2018->FracADep(pdgc, kIHAFtPiProd, ke, remnA);
	   fate_frac2025 = hd2025->FracADep(pdgc, kIHAFtPiProd, ke, remnA);	  
	  }
	  
	LOG("ReW", pNOTICE)
        << "New:::::Reweighted hadron at position = " << ip
        << " with PDG code = " << pdgc
        << ", FSI code = "  << fsi_code
        << ", KE= "  << ke
        << ", A= "  << remnA
        << " (" << INukeHadroFates::AsString((INukeFateHA_t)fsi_code) << ") :"
        << " frac2018 = "  << fate_frac2018
        <<", frac2025 = " << fate_frac2025;
	  


    }
    double frac_ratio = 0.0;
    if (fate_frac2018 != 0.0) {
	    frac_ratio = fate_frac2025 / fate_frac2018;
	} else {
	frac_ratio =1.0 ; 
	}
        

 



     // Update the current event weight
     event_weight *= frac_ratio;



  }//particle loop

  return event_weight;
}

