//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool
*/
//____________________________________________________________________________

// GENIE/Generator includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/PhysUtils.h"
#include "Physics/HadronTransport/HAIntranuke2018.h"
#include "Physics/DeepInelastic/EventGen/DISHadronicSystemGenerator.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightFZone.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;
using namespace genie::utils;

//_______________________________________________________________________________________
GReWeightFZone::GReWeightFZone() :
GReWeightModel("FormationZone")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightFZone::~GReWeightFZone()
{

}
//_______________________________________________________________________________________
bool GReWeightFZone::IsHandled(GSyst_t syst) const
{
  switch(syst) {
    case (kHadrNuclTwkDial_FormZone) :
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}
//_______________________________________________________________________________________
bool GReWeightFZone::AppliesTo (const EventRecord & event) const
{
  auto type = event.Summary()->ProcInfo().ScatteringTypeId();
  if (type==kScDeepInelastic) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightFZone::SetSystematic(GSyst_t syst, double twk_dial)
{
  switch(syst) {
    case (kHadrNuclTwkDial_FormZone) :
      fFZoneTwkDial = twk_dial;
      break;
    default:
      return;
  }
}
//_______________________________________________________________________________________
void GReWeightFZone::Reset(void)
{
  fFZoneTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightFZone::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightFZone::CalcWeight(const EventRecord & event)
{
  // Physics parameter tweaked?
  bool tweaked = (TMath::Abs(fFZoneTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  // Skip events not involving nuclear targets.
  GHepParticle * tgt = event.TargetNucleus();
  if (!tgt) return 1.;
  double A = tgt->A();
  double Z = tgt->Z();
  if (A<=1) return 1.;
  if (Z<=1) return 1.0;


  // Skip event if was not DIS scattering.
  bool is_dis = event.Summary()->ProcInfo().IsDeepInelastic();
  if(!is_dis) {
    LOG("ReW", pDEBUG) << "Not a DIS event";
    return 1.;
  }

  //
  // Calculate event weight.
  //

  double event_weight = 1.;

  // hadronic system 4-/3-momentum
  TLorentzVector p4hadr = utils::rew::Hadronic4pLAB(event);
  TVector3 p3hadr = p4hadr.Vect();

  // vertex
  assert(event.HitNucleon());
  const TLorentzVector & vtx = *(event.HitNucleon()->X4());

  // formation zone: fractional 1sigma err
  GSystUncertainty * uncertainty = GSystUncertainty::Instance();
  double fracerr = uncertainty->OneSigmaErr(kHadrNuclTwkDial_FormZone);

  // Loop over particles calculate weights for all primary hadrons inside the nucleus.
  int ip=-1;
  GHepParticle * p = 0;
  TIter event_iter(&event);
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
      ip++;

     // Skip particles with code other than 'hadron in the nucleus'
     GHepStatus_t ist  = p->Status();
     if(ist != kIStHadronInTheNucleus)
     {
        continue;
     }
     // JIMTODO - Skip if is not a hadron
     // I am not sure why the above conditional would not catch this out - need to investigate.
     int pdgc = p->Pdg();   // hadron pdg code
     if (pdg::IsHadron(pdgc) == false)
     {
        continue;
     }

     // Determine whether it interacted or not
     int fsi_code = p->RescatterCode();
     if(fsi_code == -1 || fsi_code == (int)kIHAFtUndefined) {
       LOG("ReW", pFATAL) << "INTRANUKE didn't set a valid rescattering code for event in position: " << ip;
       LOG("ReW", pFATAL) << "Here is the problematic event:";
       LOG("ReW", pFATAL) << event;
//       exit(1);
      fsi_code = kIHAFtNoInteraction;
     }
     bool escaped    = (fsi_code == (int)kIHAFtNoInteraction);
     bool interacted = !escaped;

     LOG("ReW", pDEBUG)
        << "Attempting to reweight hadron at position = " << ip
        << " with PDG code = " << pdgc
        << ". The hadron "
        << ((interacted) ? "re-interacted" : "did not re-interact");

     // Default formation zone
     double m = p->Mass();
     const TLorentzVector * p4  = p->P4();

     double ct0=0.;
     pdg::IsNucleon(pdgc) ? ct0=fct0nucleon : ct0=fct0pion;
     double fz_def = phys::FormationZone(m,*p4,p3hadr,ct0,fK);

     double fz_scale_factor  = (1 + fFZoneTwkDial * fracerr);
     double fz_twk  = fz_def * fz_scale_factor;
     fz_twk = TMath::Max(0.,fz_twk);

     LOG("ReW", pDEBUG)
        << "Formation zone = " << fz_def << " fm (nominal), "
        << fz_twk << " fm (tweaked) - scale factor = " << fz_scale_factor;

     // Calculate hadron's position at the end of the formation zone step
     TVector3 step3v = p4->Vect();
     step3v.SetMag(fz_def);
     TLorentzVector step4v(step3v, 0.);

     TLorentzVector x4 = vtx + step4v;
     LOG("ReW", pDEBUG)  << "Vtx position: "<< print::X4AsString(&vtx);
     LOG("ReW", pDEBUG)  << "Hadron position: "<< print::X4AsString(&x4);

     // Calculate particle weight
     double hadron_weight = genie::utils::rew::FZoneWeight(
        pdgc, vtx, x4, *p4, A, Z, fz_scale_factor, interacted, *fFSIModel );

     // Update event weight
     event_weight *= hadron_weight;
  }

  return event_weight;
}
//_______________________________________________________________________________________
void GReWeightFZone::Init(void)
{
  fFZoneTwkDial = 0.;

  // Obtain the formation zone parameters needed for reweighting from the
  // DISHadronicSystemGenerator configuration.
  AlgFactory* algf = AlgFactory::Instance();
  const genie::Algorithm* dis_hs_alg = algf->GetAlgorithm(
    "genie::DISHadronicSystemGenerator", "Default" );
  const genie::Registry& dis_hs_reg = dis_hs_alg->GetConfig();

  fct0pion = dis_hs_reg.GetDouble( "FZONE-ct0pion" ); // fm
  fct0nucleon = dis_hs_reg.GetDouble( "FZONE-ct0nucleon" ); // fm
  fK = dis_hs_reg.GetDouble( "FZONE-KPt2" );

  // TODO: reduce code duplication from the constructor for GReWeightINuke
  // Look up the FSI model for the current tune. Also check whether FSIs are
  // actually enabled.
  AlgConfigPool* conf_pool = AlgConfigPool::Instance();
  Registry* gpl = conf_pool->GlobalParameterList();
  RgAlg fsi_alg = gpl->GetAlg( "HadronTransp-Model" );
  bool fsi_enabled = gpl->GetBool( "HadronTransp-Enable" );

  if ( !fsi_enabled ) {
    LOG( "ReW", pERROR ) << "FSIs are not enabled for the current tune."
      << " Refusing to reweight the formation zone.";
    std::exit( 1 );
  }

  AlgId id( fsi_alg );

  Algorithm* alg = algf->AdoptAlgorithm( id );
  fFSIModel = dynamic_cast< HAIntranuke2018* >( alg );

  if ( !fFSIModel ) {
    LOG( "ReW", pERROR ) << "Reweighting events produced with the FSI model "
      << fsi_alg << " is not currently supported.";
    std::exit( 1 );
  }

  fFSIModel->AdoptSubstructure();
}
//_______________________________________________________________________________________
