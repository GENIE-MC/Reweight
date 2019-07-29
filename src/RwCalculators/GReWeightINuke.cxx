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
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/HadronTransport/INukeHadroData.h"
#include "Physics/HadronTransport/INukeHadroFates.h"
#include "Physics/HadronTransport/INukeUtils.h"

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
    case kScCoherent:
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
     TLorentzVector x4 (p->Vx(), p->Vy(), p->Vz(), 0.    );
     TLorentzVector p4 (p->Px(), p->Py(), p->Pz(), p->E());

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
        w_mfp = utils::rew::MeanFreePathWeight(pdgc,x4,p4,A,Z,mfp_scale_factor,interacted);
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
     double d        = utils::intranuke::Dist2Exit(x4,p4,A);
     double d_mfp    = utils::intranuke::Dist2ExitMFP(pdgc,x4,p4,A,Z);
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
