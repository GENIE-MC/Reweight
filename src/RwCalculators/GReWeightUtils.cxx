//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author:  Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

// GENIE/Generator includes
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/HadronTransport/Intranuke2018.h"
#include "Physics/HadronTransport/INukeHadroData2018.h"
#include "Physics/HadronTransport/INukeHadroFates2018.h"
#include "Physics/HadronTransport/INukeUtils2018.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightUtils.h"

using namespace genie;
using namespace genie::rew;
using namespace genie::controls;

//____________________________________________________________________________
double genie::utils::rew::MeanFreePathWeight(
  int pdgc, const TLorentzVector & x4, const TLorentzVector & p4,
  double A, double Z,
  double mfp_scale_factor, bool interacted, const Intranuke2018& fsi_model )
{
   LOG("ReW", pINFO)
     << "Calculating mean free path weight: "
     << "target A = " << A << ", target Z = " << Z
     << ", mfp_scale = " << mfp_scale_factor
     << ", interacted = " << interacted;

   // Get the nominal survival probability
   double pdef = utils::intranuke2018::ProbSurvival(
      pdgc, x4, p4, A, Z, 1., fsi_model );
   LOG("ReW", pINFO)  << "Probability(default mfp) = " << pdef;
   if(pdef<=0) return 1.;

   // Get the survival probability for the tweaked mean free path
   double ptwk = utils::intranuke2018::ProbSurvival(
      pdgc, x4, p4, A, Z, mfp_scale_factor, fsi_model );
   LOG("ReW", pINFO)  << "Probability(tweaked mfp) = " << ptwk;
   if(ptwk<=0) return 1.;

   // Calculate weight
   double w_mfp = utils::rew::MeanFreePathWeight(pdef, ptwk, interacted);
   LOG("ReW", pINFO)  << "Mean free path weight = " << w_mfp;
   return w_mfp;
}
//____________________________________________________________________________
double genie::utils::rew::FZoneWeight(
  int pdgc, const TLorentzVector & vtx, const TLorentzVector & x4,
  const TLorentzVector & p4, double A, double Z,
  double fz_scale_factor, bool interacted, const Intranuke2018& fsi_model )
{
   // Calculate hadron start assuming tweaked formation zone
   TLorentzVector fz    = x4 - vtx;
   TLorentzVector fztwk = fz_scale_factor*fz;
   TLorentzVector x4twk = x4 + fztwk - fz;

   LOG("ReW", pDEBUG)  << "Formation zone = "<< fz.Vect().Mag() << " fm";

   // Get nominal survival probability.
   double pdef = utils::intranuke2018::ProbSurvival(
      pdgc, x4, p4, A, Z, 1., fsi_model );
   LOG("ReW", pDEBUG)  << "Survival probability (nominal) = "<< pdef;
   if(pdef<=0) return 1.;
   if(pdef>=1.){
     LOG("ReW", pERROR)
         <<  "Default formation zone takes hadron outside "
         <<  "nucleus so cannot reweight!" ;
     return 1.;
   }

   // Get tweaked survival probability.
   double ptwk = utils::intranuke2018::ProbSurvival(
      pdgc, x4twk, p4, A, Z, 1., fsi_model );
   if(ptwk<=0) return 1.;
   LOG("ReW", pDEBUG)  << "Survival probability (tweaked) = "<< ptwk;

   // Calculate weight
   double w_fz = utils::rew::MeanFreePathWeight(pdef, ptwk, interacted);
   LOG("ReW", pDEBUG)
       << "Particle weight for formation zone tweak = "<< ptwk;
   return w_fz;
}
//____________________________________________________________________________
double genie::utils::rew::MeanFreePathWeight(
       double pdef, double ptwk, bool interacted)
{
// Returns a weight to account for a change in hadron mean free path inside
// insidea nuclear medium.
//
// Inputs:
//   pdef : nominal survival probability
//   ptwk : survival probability for the tweaked mean free path
//   interacted : flag indicating whether the hadron interacted or escaped
//
// See utils::intranuke2018::ProbSurvival() for the calculation of probabilities.
//
  double w_mfp = 1.;

  if(interacted) {
     w_mfp = (1-pdef>0) ?  (1-ptwk)  /   (1-pdef)  : 1;
  } else {
     w_mfp = (pdef>0) ? ptwk / pdef : 1;
  }
  w_mfp = TMath::Max(0.,w_mfp);

  return w_mfp;
}
//____________________________________________________________________________
double genie::utils::rew::FateFraction(genie::rew::GSyst_t syst, double kinE,
  int target_A, double frac_scale_factor)
{
  if ( target_A < 1 ) LOG("ReW", pERROR) << "Invalid mass number A = "
    << target_A << " passed to genie::utils::rew::FateFraction";

  double fate_frac = 0.0;

  INukeHadroData2018 * hd = INukeHadroData2018::Instance();

  // convert to MeV and
  double ke = kinE / units::MeV;
  ke = TMath::Max(INukeHadroData2018::fMinKinEnergy,   ke);
  ke = TMath::Min(INukeHadroData2018::fMaxKinEnergyHA, ke);

  switch (syst) {

    //
    // pions
    //

    case (genie::rew::kINukeTwkDial_FrCEx_pi) :
    {
      fate_frac = hd->FracADep(kPdgPiP, kIHAFtCEx, ke, target_A);
    }
    break;

    //    case (genie::rew::kINukeTwkDial_FrElas_pi) :
    //    {
    //      fate_frac = hd->FracADep(kPdgPiP, kIHAFtElas, ke, target_A);
    //    }
    //    break;

    case (genie::rew::kINukeTwkDial_FrInel_pi) :
    {
      fate_frac = hd->FracADep(kPdgPiP, kIHAFtInelas, ke, target_A);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrAbs_pi) :
    {
      fate_frac = hd->FracADep(kPdgPiP, kIHAFtAbs, ke, target_A);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrPiProd_pi) :
    {
      fate_frac = hd->FracADep(kPdgPiP, kIHAFtPiProd,  ke, target_A);
    }
    break;

    //
    // nucleons
    //

    case (genie::rew::kINukeTwkDial_FrCEx_N) :
    {
      fate_frac = hd->FracAIndep(kPdgProton, kIHAFtCEx, ke);
    }
    break;

    //    case (genie::rew::kINukeTwkDial_FrElas_N) :
    //    {
    //      fate_frac = hd->Frac(kPdgProton, kIHAFtElas, ke);
    //    }
    //    break;

    case (genie::rew::kINukeTwkDial_FrInel_N) :
    {
      fate_frac = hd->FracAIndep(kPdgProton, kIHAFtInelas, ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrAbs_N) :
    {
      fate_frac = hd->FracAIndep(kPdgProton, kIHAFtAbs,    ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrPiProd_N) :
    {
      fate_frac = hd->FracAIndep(kPdgProton, kIHAFtPiProd,  ke);
    }
    break;

    default:
    {
      LOG("ReW", pDEBUG)
        << "Have reached default case and assigning fraction{fate} = 0";
      fate_frac = 0;
    }
    break;

  } // hadron_fate?

  fate_frac *= frac_scale_factor;

  return fate_frac;
}
//____________________________________________________________________________
double genie::utils::rew::WhichFateFractionScaleFactor(
    genie::rew::GSyst_t syst, double kinE, int target_A, double fate_frac)
{
  double fate_frac_nominal = FateFraction(syst, kinE, target_A, 1.0);

  // Avoid NaNs if both the nominal value and the tweaked value are zero
  if ( fate_frac_nominal == 0. && fate_frac == 0. ) return 1.;

  if ( fate_frac_nominal <= 0. ) {
    // We're having some sort of problem with the fate fraction calculation
    LOG("ReW", pERROR) << "Nonpositive nominal fate fraction"
      << " encountered in genie::utils::rew::WhichFateFractionScaleFactor()";
    return -99999.;
  }

  double scale = std::max(0., fate_frac) / fate_frac_nominal;

  return scale;
}
//____________________________________________________________________________
bool genie::utils::rew::HadronizedByAGKY(const EventRecord & event)
{
  Interaction * interaction = event.Summary();
  assert(interaction);

  bool is_dis  = interaction->ProcInfo().IsDeepInelastic();
  bool charm   = interaction->ExclTag().IsCharmEvent();
  bool by_agky = is_dis && !charm;

  return by_agky;
}
//____________________________________________________________________________
bool genie::utils::rew::HadronizedByAGKYPythia(const EventRecord & event)
{
// Check whether the event was hadronized by AGKY/KNO or AGKY/PYTHIA

  GHepStatus_t prefragm = kIStDISPreFragmHadronicState;
  bool found_string  = (event.FindParticle(kPdgString,  prefragm, 0) != 0);
  bool found_cluster = (event.FindParticle(kPdgCluster, prefragm, 0) != 0);
  bool handled_by_pythia = found_string || found_cluster;

  return handled_by_pythia;
}
//____________________________________________________________________________
TLorentzVector genie::utils::rew::Hadronic4pLAB(const EventRecord & event)
{
  GHepParticle * nu = event.Probe();                    // incoming v
  GHepParticle * N  = event.HitNucleon();               // struck nucleon
  GHepParticle * l  = event.FinalStatePrimaryLepton();  // f/s primary lepton

  assert(nu);
  assert(N);
  assert(l);

  // Compute the Final State Hadronic System 4p (PX = Pv + PN - Pl)

  const TLorentzVector & p4nu = *(nu->P4());
  const TLorentzVector & p4N  = *(N ->P4());
  const TLorentzVector & p4l  = *(l ->P4());

  TLorentzVector pX4 = p4nu + p4N - p4l;

  return pX4;
}
//____________________________________________________________________________
double genie::utils::rew::AGKYWeight(int /*pdgc*/, double /*xF*/, double /*pT2*/)
{
  return 1.0;
}
//____________________________________________________________________________
int genie::utils::rew::Sign(double twkdial)
{
  if(twkdial < 0.) return -1;
  if(twkdial > 0.) return +1;
  return 0;
}
//____________________________________________________________________________
