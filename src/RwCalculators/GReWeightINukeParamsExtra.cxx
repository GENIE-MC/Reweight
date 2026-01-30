//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Gray Putnam <gputnam \at fnal.gov>
         Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>

#include <TMath.h>
#include <TLorentzVector.h>

// GENIE/Generator includes
#include "Framework/Conventions/Controls.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/HadronTransport/INukeHadroData2018.h"
#include "Physics/HadronTransport/INukeHadroFates2018.h"
#include "Physics/NuclearState/NuclearUtils.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightINukeParamsExtra.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//___________________________________________________________________________
GReWeightINukeParamsExtra::GReWeightINukeParamsExtra(void)
{
  fParmPionFates = new GReWeightINukeParamsExtra::Fates(kRwINukePion);
  fParmNuclFates = new GReWeightINukeParamsExtra::Fates(kRwINukeNucl);
  fParmPionMFP   = new GReWeightINukeParamsExtra::MFP(kRwINukePion);
  fParmNuclMFP   = new GReWeightINukeParamsExtra::MFP(kRwINukeNucl);
}
//___________________________________________________________________________
GReWeightINukeParamsExtra::~GReWeightINukeParamsExtra(void)
{

}
//___________________________________________________________________________
//
// Fates nested class
//
//___________________________________________________________________________
GReWeightINukeParamsExtra::Fates::Fates(GReWeightINukeParams::HadronType_t ht)
  : GReWeightINukeParams::Fates::Fates( ht )
{
  fModelSwitch = kNoSwitch;
}
//___________________________________________________________________________
GReWeightINukeParamsExtra::Fates::~Fates(void)
{
  this->Reset();
}
//___________________________________________________________________________
void GReWeightINukeParamsExtra::Fates::Reset()
{
  GReWeightINukeParams::Fates::Reset();

  fSystKELow = -1.;
  fSystKEHigh = -1.;
}
//___________________________________________________________________________
void GReWeightINukeParamsExtra::Fates::SetSystKERange(GSyst_t syst) {
  if (syst == kINukeTwkDial_INCLLoE_N || syst == kINukeTwkDial_G4LoE_N) {
    fSystKELow = 0;
    fSystKEHigh = 0.15 / units::GeV;
  }
  else if (syst == kINukeTwkDial_INCLM1E_N || syst == kINukeTwkDial_G4M1E_N) {
    fSystKELow = 0.15 / units::GeV;
    fSystKEHigh = 0.3 / units::GeV;
  }
  else if (syst == kINukeTwkDial_INCLM2E_N || syst == kINukeTwkDial_G4M2E_N) {
    fSystKELow = 0.3 / units::GeV;
    fSystKEHigh = 0.6 / units::GeV;
  }
  else if (syst == kINukeTwkDial_INCLHiE_N || syst == kINukeTwkDial_G4HiE_N) {
    fSystKELow = 0.6 / units::GeV;
    fSystKEHigh = -1;
  }
}
//___________________________________________________________________________
void GReWeightINukeParamsExtra::Fates::SetTwkDial(GSyst_t syst, double val)
{
  // For dials that do not involve model transformations, the base class
  // implementation is already sufficient
  if (!this->IsModelTransform(syst)) {
    return GReWeightINukeParams::Fates::SetTwkDial( syst, val );
  }

  // check type of systematic
  if(!this->IsHandled(syst)) return;

  // can not explicitly set params designated as cushion terms
  if(this->IsIncluded(syst)) {
    if(this->IsCushionTerm(syst)) {
       LOG("ReW", pWARN)
         << "You may not set the value of cushion term " << GSyst::AsString(syst)
         << ". It is set automatically to maintain unitarity";
       return;
    }
  }

  if (this->IsG4ModelTransform(syst)) fModelSwitch = kRwINukeG4;
  else if (this->IsINCLModelTransform(syst)) fModelSwitch = kRwINukeINCL;

  // set the kinematic energy range for this systematic
  this->SetSystKERange(syst);

  // set all the syst values to 1 -- these will be scaled by the model difference
  int i=0;
  while( (syst = (fHadType == kRwINukePion) ?
          GSyst::NextPionFateSystematic(i++) :
          GSyst::NextNuclFateSystematic(i++)
       ) != kNullSystematic
     )
  {
    // Leave inelastic as cushion term to absorb the residual fraction in the systematic.
    // Since the fractions sum to 1, this will set the value correctly.
    if (syst != kINukeTwkDial_FrInel_N) {
      fSystValuesUser[syst]  = val;
      fIsCushion[syst] = false;
    }
  }
}
//___________________________________________________________________________
bool GReWeightINukeParamsExtra::Fates::IsInSystKERange(double KE) const {
  return (fSystKELow < 0 || KE >= fSystKELow) && (fSystKEHigh < 0. || KE < fSystKEHigh);
}
//___________________________________________________________________________
double GReWeightINukeParamsExtra::Fates::ScaleFactor(
      GSyst_t syst, double KE) const
{
  // check KE
  if (!IsInSystKERange(KE)) return 1.;

  // delegate to base class
  return GReWeightINukeParams::Fates::ScaleFactor(syst, KE);
}
//___________________________________________________________________________
double GReWeightINukeParamsExtra::Fates::OneSigmaErr(GSyst_t syst, double KE) const {
  // If we're not handling a model switch dial, then the base class
  // implementation is sufficient and we can delegate to it
  if ( fModelSwitch == kNoSwitch ) {
    return GReWeightINukeParams::Fates::OneSigmaErr( syst, KE );
  }

  GSystUncertainty * uncert = GSystUncertainty::Instance();
  const GReWeightINukeModelSwitchData* inuked
    = GReWeightINukeModelSwitchData::Instance();

  LOG("ReW", pNOTICE) << "OneSigmaErr FOR SYST: " << GSyst::AsString(syst) << " is switch: " << fModelSwitch << " is handled: " << inuked->IsHandled(syst);

  // Get the model error, and scale by the one-sigma
  double one_sigma = uncert->OneSigmaErr((fModelSwitch == kRwINukeG4) ? kINukeTwkDial_G4_N : kINukeTwkDial_INCL_N);

  // If we're scaling to G4 or INCL, fix the scale factor
  if (inuked->IsHandled(syst)) {
    double nom_frac = genie::utils::rew::FateFraction(syst, KE, fTargetA, 1.);
    double var_frac = inuked->FateFraction(fModelSwitch, syst, KE);

    double var_scale = 0;
    // protect against bad weight
    if (nom_frac > 1e-6) var_scale = (var_frac - nom_frac) / nom_frac;

    LOG("ReW", pNOTICE) << "OneSigmaErr Model Switch from " << nom_frac << " to " << var_frac << ", scale: " << var_scale;

    return var_scale * one_sigma;
  }

  return 0.;
}
//___________________________________________________________________________
bool GReWeightINukeParamsExtra::Fates::IsG4ModelTransform(GSyst_t syst) const
{
  return (syst == kINukeTwkDial_G4_N) || \
    (syst == kINukeTwkDial_G4LoE_N) || \
    (syst == kINukeTwkDial_G4M1E_N) || \
    (syst == kINukeTwkDial_G4M2E_N) || \
    (syst == kINukeTwkDial_G4HiE_N);
}
//___________________________________________________________________________
bool GReWeightINukeParamsExtra::Fates::IsINCLModelTransform(GSyst_t syst) const
{
  return (syst == kINukeTwkDial_INCL_N) || \
    (syst == kINukeTwkDial_INCLLoE_N) || \
    (syst == kINukeTwkDial_INCLM1E_N) || \
    (syst == kINukeTwkDial_INCLM2E_N) || \
    (syst == kINukeTwkDial_INCLHiE_N);
}
//___________________________________________________________________________
bool GReWeightINukeParamsExtra::Fates::IsModelTransform(GSyst_t syst) const
{
  return IsINCLModelTransform(syst) || IsG4ModelTransform(syst);
}
//___________________________________________________________________________
//
// MFP nested class
//
//___________________________________________________________________________
GReWeightINukeParamsExtra::MFP::MFP(GReWeightINukeParams::HadronType_t ht)
  : GReWeightINukeParams::MFP::MFP( ht )
{
}
//___________________________________________________________________________
GReWeightINukeParamsExtra::MFP::~MFP()
{

}
//___________________________________________________________________________
double GReWeightINukeParamsExtra::MFP::ScaleFactor(double KE) const
{
  // check KE
  if (!this->IsInSystKERange(KE)) return 1.;

  // delegate to base class
  return GReWeightINukeParams::MFP::ScaleFactor( KE );
}
//___________________________________________________________________________
void GReWeightINukeParamsExtra::MFP::Reset(void)
{
  GReWeightINukeParams::MFP::Reset();

  fSystKELow = -1.;
  fSystKEHigh = -1.;
}
//___________________________________________________________________________
bool GReWeightINukeParamsExtra::MFP::IsInSystKERange(double KE) const {
  return (fSystKELow < 0 || KE >= fSystKELow) && (fSystKEHigh < 0. || KE < fSystKEHigh);
}
//___________________________________________________________________________
void GReWeightINukeParamsExtra::MFP::SetSystKERange(GSyst_t syst) {
  if (syst == kINukeTwkDial_MFPLoE_N) {
    fSystKELow = 0;
    fSystKEHigh = 0.15 / units::GeV;
  }
  else if (syst == kINukeTwkDial_MFPM1E_N) {
    fSystKELow = 0.15 / units::GeV;
    fSystKEHigh = 0.3 / units::GeV;
  }
  else if (syst == kINukeTwkDial_MFPM2E_N) {
    fSystKELow = 0.3 / units::GeV;
    fSystKEHigh = 0.6 / units::GeV;
  }
  else if (syst == kINukeTwkDial_MFPHiE_N) {
    fSystKELow = 0.6 / units::GeV;
    fSystKEHigh = -1;
  }
}
//___________________________________________________________________________
void GReWeightINukeParamsExtra::MFP::SetTwkDial(GSyst_t syst, double val)
{
  this->SetSystKERange(syst);
  fTwkDial    = val;
  fIsIncluded = true;
}
