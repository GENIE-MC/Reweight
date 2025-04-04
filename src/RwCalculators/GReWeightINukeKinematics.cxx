#include <cassert>
#include <cstdlib>

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

// GENIE/Reweight includes
#include "RwCalculators/GReWeightINukeKinematics.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
bool GReWeightINukeKinematics::AppliesTo(const EventRecord & event) const
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
bool GReWeightINukeKinematics::IsHandled(GSyst_t syst) const
{
  switch (syst) {
    case ( kINukeKinematicsTwkDial_NP_N ):
    default:
    case ( kINukeKinematicsTwkDial_PP_N ):
      return true;
  }

  return false;

}

//_______________________________________________________________________________________
void GReWeightINukeKinematics::SetSystematic(GSyst_t syst, double val)
{
  if (IsHandled(syst)) {
     fRwParams.SetTwkDial(syst, val);
  }
}

//_______________________________________________________________________________________
void GReWeightINukeKinematics::Reset() {
  fRwParams.Reset();
}

//_______________________________________________________________________________________
void GReWeightINukeKinematics::Reconfigure() {}

// "Lab" energy of particle p in p + t -> f1 + f2 scattering process.
// That is -- the kinetic energy of p in the frame where t is at rest
double TLab(const TLorentzVector &p, const TLorentzVector &f1, const TLorentzVector &f2) {
  TLorentzVector t = f1 + f2 - p; // E.M. conservation
  double E_p = ((p + t).Mag2() - t.M2() - p.M2())/(2.0*t.M());

  return E_p - p.M();
}

// Scattering angle in p + t -> f1 + f2 scattering process,
// in center-of-momentum (C.O.M.) frame
double costhcm(const TLorentzVector &p, const TLorentzVector &f1, const TLorentzVector &f2) {
  TLorentzVector t = f1 + f2 - p; // E.M. conservation

  double Ecm2 = (f1 + f2).Mag2();

  double E_p_cm = (Ecm2 + p.M2() - t.M2()) / (2*sqrt(Ecm2));
  double P_p_cm = sqrt(E_p_cm*E_p_cm - p.M2());
  double E_f1_cm = (Ecm2 + f1.M2() - f2.M2()) / (2*sqrt(Ecm2));
  double P_f1_cm = sqrt(E_f1_cm*E_f1_cm - f1.M2());

  return ((p + f1).M2() - p.M2() - f1.M2() - 2*E_p_cm*E_f1_cm) / (2*P_f1_cm*P_p_cm);
}

double GReWeightINukeKinematics::CalcWeight(const EventRecord &event) {
  // to return
  double event_weight = 1.;

  // get the atomic mass number for the hit nucleus
  GHepParticle * tgt = event.TargetNucleus();
  if (!tgt) return 1.0;
  double A = tgt->A();
  double Z = tgt->Z();
  if (A<=1) return 1.0;
  if (Z<=1) return 1.0;

  fRwParams.SetTargetA( A );

  // Loop over stdhep entries and only calculate weights for particles.
  // All particles that are not hadrons generated inside the nucleus are given weights of 1.0
  int ip=-1;
  GHepParticle * p = 0;
  TIter event_iter(&event);
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
     ip++;

     // Skip particles not rescattered by the actual hadron transport code
     int  pdgc       = p->Pdg();
     // bool is_pion    = pdg::IsPion   (pdgc);
     bool is_nucleon = pdg::IsNucleon(pdgc);
     // bool is_kaon = pdg::IsKaon( pdgc );
     if (!is_nucleon)
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
     LOG("RWINukeKin", pDEBUG)
        << "Attempting to reweight hadron at position = " << ip
        << " with PDG code = " << pdgc
        << " and FSI code = "  << fsi_code
        << " (" << INukeHadroFates::AsString((INukeFateHA_t)fsi_code) << ")\n";

     // Inelastic, Charge Exchange 
     if (fsi_code == (int)kIHAFtInelas || fsi_code == (int)kIHAFtCEx) {
       const GHepParticle &f1 = *event.Particle(p->FirstDaughter());
       const GHepParticle &f2 = *event.Particle(p->LastDaughter());

       int np_pp = 0; // == 0 if nucleons in scatter are different, == 1 if they are the same
       if (fsi_code == (int)kIHAFtInelas) np_pp = (f1.Pdg() != f2.Pdg()) ? 0 : 1;
       else if (fsi_code == (int)kIHAFtCEx) np_pp =  (f1.Pdg() == f2.Pdg()) ? 0 : 1;

       double tlab_v = TLab(*p->P4(), *f1.P4(), *f2.P4());
       double costhcm_v = costhcm(*p->P4(), *f1.P4(), *f2.P4());

       double weight = fRwParams.CalcWeight(np_pp, tlab_v, costhcm_v);

       // std::cerr << "TLab: " << tlab_v << " costhcm: " << costhcm_v << " np/pp " << np_pp << " weight: " << weight << std::endl;

       event_weight *= weight;
     }
  }

  return event_weight;
}

