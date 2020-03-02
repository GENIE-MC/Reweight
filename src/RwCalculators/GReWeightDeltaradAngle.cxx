//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Gray Yarbrough <gyarbrou \at vols.utk.edu>
 University of Tennessee, Knoxvile

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "TH1D.h"
#include "TParticlePDG.h"
#include "TDecayChannel.h"

// GENIE/Generator includes
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightDeltaradAngle.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightDeltaradAngle::GReWeightDeltaradAngle() :
GReWeightModel("GReWeightDeltaradAngle")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightDeltaradAngle::~GReWeightDeltaradAngle()
{

}
//_________________
bool GReWeightDeltaradAngle::IsHandled(GSyst_t syst) const
{
  switch( syst ) {
    case( kRDcyTwkDial_Theta_Delta2NRad ):
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}
//_______________________________________________________________________________________
bool GReWeightDeltaradAngle::AppliesTo(ScatteringType_t type, bool /*is_cc*/) const
{
  if ( type == kScResonant) return true;
  return false;
}
//_______________________________________________________________________________________
void GReWeightDeltaradAngle::SetSystematic(GSyst_t syst, double twk_dial)
{
  switch( syst ) {
    case ( kRDcyTwkDial_Theta_Delta2NRad  ):
      fThetaDelta2NRadTwkDial = twk_dial;
      break;
    default:
      return;
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightDeltaradAngle::Reset(void)
{
  fThetaDelta2NRadTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightDeltaradAngle::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightDeltaradAngle::CalcWeight(const EventRecord & event)
{
  Interaction* interaction = event.Summary();

  bool is_res = interaction->ProcInfo().IsResonant();
  if ( !is_res ) return 1.;

  double wght = this->RewThetaDelta2NRad( event );

  return wght;
}
//_______________________________________________________________________________________
double GReWeightDeltaradAngle::RewThetaDelta2NRad(const EventRecord & event)
{
  bool tweaked = TMath::Abs( fThetaDelta2NRadTwkDial ) > controls::kASmallNum;
  if ( !tweaked ) return 1.;

  bool is_Delta_rad = false;
  int ir = -1; // resonance position
  int ig = -1; // gamma position
  int i  =  0;
  GHepParticle* p = 0;
  TIter iter( &event );
  while ( (p = dynamic_cast< GHepParticle* >(iter.Next())) ) {
    bool is_Deltap = ( p->Pdg() == kPdgP33m1232_DeltaP );
    if ( is_Deltap ) {
      ir = i;
      int fd = p->FirstDaughter();
      int ld = p->LastDaughter();
      int nd = 1 + ld - fd;
      if ( nd==2 ) {
        int fpdg = event.Particle( fd )->Pdg();
        int lpdg = event.Particle( ld )->Pdg();
        if ( fpdg == kPdgGamma && lpdg == kPdgProton) {
          is_Delta_rad = true;
          ig = fd;
        }
        if ( fpdg==kPdgProton && lpdg==kPdgGamma ) {
          is_Delta_rad = true;
          ig = ld;
        }
      }
    }

    if ( is_Delta_rad ) break;
    ++i;
  }

  if ( !is_Delta_rad ) return 1.;

  LOG("ReW", pDEBUG) << "A Delta+ -> p photon event:";
  LOG("ReW", pDEBUG) << "Resonance is at position: " << ir;
  LOG("ReW", pDEBUG) << "Gamma is at position: " << ig;

  // Get Delta and gamma 4-momenta
  TLorentzVector p4Res( *event.Particle(ir)->P4() );
  TLorentzVector p4Gamma( *event.Particle(ig)->P4() );

  // Boost gamma to the Delta rest frame
  TVector3 bv = -1. * p4Res.BoostVector();
  p4Gamma.Boost( bv );

  // Calculate the weight in the Delta rest frame
  // To conserve total cross section, normalize P1 and
  // P2 to integrate to unity on costheta in [-1, 1]
  // Radiative decays are isotropic by default
  double P1 = 1./2.;
  double costheta = p4Gamma.Vect().CosTheta();
  // Alternate distribution: cos^2(theta)
  double P2 = 3./2. * std::pow( costheta, 2 );
  double dial = fThetaDelta2NRadTwkDial;

  double wght = ( dial*P2 + (1. - dial)*P1 ) / P1;

  LOG("ReW", pDEBUG) << "Gamma Cos(ThetaCM) = "
    << costheta << ", weight = " << wght;

  return wght;
}
//_______________________________________________________________________________________
void GReWeightDeltaradAngle::Init(void)
{
  fThetaDelta2NRadTwkDial = 0.;
}
