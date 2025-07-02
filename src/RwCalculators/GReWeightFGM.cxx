//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtupleD.h>

// GENIE/Generator includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/PauliBlocker.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/ParticleData/PDGCodes.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightFGM.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;
using namespace genie::utils;

//_______________________________________________________________________________________
GReWeightFGM::GReWeightFGM() :
GReWeightModel("FermiGasModel")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightFGM::~GReWeightFGM()
{
#ifdef _G_REWEIGHT_FGM_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightFGM::IsHandled(GSyst_t syst) const
{
  switch(syst) {
    case ( kSystNucl_CCQEPauliSupViaKF   ) :
    case ( kSystNucl_CCQEMomDistroFGtoSF ) :
        return true;
        break;
    default:
        return false;
        break;
  }
  return false;
}
//_______________________________________________________________________________________
bool GReWeightFGM::AppliesTo (const EventRecord & event) const
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
void GReWeightFGM::SetSystematic(GSyst_t syst, double val)
{
  switch(syst) {
    case ( kSystNucl_CCQEPauliSupViaKF   ) :
        fKFTwkDial = val;
        break;
    case ( kSystNucl_CCQEMomDistroFGtoSF ) :
        fMomDistroTwkDial = val;
        break;
    default:
        return;
        break;
  }
}
//_______________________________________________________________________________________
void GReWeightFGM::Reset(void)
{
  fKFTwkDial        = 0.;
  fMomDistroTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightFGM::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightFGM::CalcWeight(const EventRecord & event)
{
  double wght =
    this->RewCCQEPauliSupViaKF   (event) *
    this->RewCCQEMomDistroFGtoSF (event);

  return wght;
}
//_______________________________________________________________________________________
double GReWeightFGM::RewCCQEPauliSupViaKF(const EventRecord & event)
{
  bool kF_tweaked = (TMath::Abs(fKFTwkDial) > controls::kASmallNum);
  if ( !kF_tweaked ) return 1.;

  bool is_qe = event.Summary()->ProcInfo().IsQuasiElastic();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if ( !is_qe || !is_cc ) return 1.;

  const Target & target = event.Summary()->InitState().Tgt();
  if ( !target.IsNucleus() ) {
     return 1.;
  }

  // Don't bother to compute weights for very light nuclear targets
  double A = target.A();
  if ( A <= 4 ) {
    return 1.;
  }

  // Skip CCQE channels that do not produce a final-state nucleon
  // (e.g., Charm-CCQE)
  int final_nucleon_pdgc  = event.Summary()->RecoilNucleonPdg();
  if ( !pdg::IsNucleon(final_nucleon_pdgc) ) return 1.;

  // Use the PauliBlocker to look up the Fermi momentum. This is a bit
  // of a hack, but it will allow us to get the correct value regardless
  // of whether the current tune is using the local Fermi gas model
  // or not.
  AlgFactory* algf = AlgFactory::Instance();

  // TODO: Make this more robust. Right now it assumes that the
  // default PauliBlocker configuration is always being used.
  // This may not be the case forever.
  const PauliBlocker* pb = dynamic_cast<const PauliBlocker*>(
    algf->GetAlgorithm("genie::PauliBlocker", "Default"));

  assert( pb );

  int struck_nucleon_pdgc = target.HitNucPdg();
  double hit_nuc_radius = target.HitNucPosition();

  // Fermi momentum for the final nucleon
  double kFf = pb->GetFermiMomentum(target, final_nucleon_pdgc,
    hit_nuc_radius);

  // Get the fractional error for tweaking kF
  GSystUncertainty * uncertainty = GSystUncertainty::Instance();
  double kF_fracerr = uncertainty->OneSigmaErr(kSystNucl_CCQEPauliSupViaKF);

  // Tweaked kF for the final nucleon
  double kFf_twk = kFf * (1 + fKFTwkDial * kF_fracerr);

  // For QELEventGenerator, Pauli blocking is now handled by the cross
  // section models themselves. The calculation doesn't rely on the
  // analytic Pauli suppression factor used previously (which was only
  // really valid for the FGMBodekRitchie nuclear model). Here we check
  // whether the phase space for the CCQE event matches that used by
  // QELEventGenerator. Admittedly, this is also a bit of a hack.
  if ( event.DiffXSecVars() == kPSQELEvGen ) {

    // Final nucleon momentum
    double pNf = event.Summary()->Kine().HadSystP4().P();

    // This event was generated using QELEventGenerator, so
    // just check whether the final nucleon momentum is below the
    // tweaked Fermi momentum. If it is, set the event weight to
    // zero. Otherwise, set it to unity.
    double weight = 1.;
    if ( pNf < kFf_twk ) weight = 0.;

    return weight;
  }

  // ************
  // This code preserves the old behavior of this tweak dial for
  // cases where QELEventGenerator was not used
  // ************

  // Fermi momentum for the initial nucleon
  double kFi = pb->GetFermiMomentum(target, struck_nucleon_pdgc,
    hit_nuc_radius);

  // Tweaked kF for the initial nucleon
  double kFi_twk = kFi * (1 + fKFTwkDial * kF_fracerr);

  // hit nucleon mass
  double Mn = target.HitNucP4Ptr()->M(); // can be off m/shell

  // momentum transfer
  const Kinematics & kine = event.Summary()->Kine();
  bool selected = true;
  double q2 = kine.q2(selected);

  const double pmax = 0.5; // GeV

  // default nuclear suppression
  double R = utils::nuclear::RQEFG_generic(q2, Mn, kFi, kFf, pmax);

  // calculate tweaked nuclear suppression factor
  double Rtwk = nuclear::RQEFG_generic(q2, Mn, kFi_twk, kFf_twk, pmax);

  // calculate weight (ratio of suppression factors)

  double wght = 1.0;

  if ( R > 0. && Rtwk > 0. ) {
    wght = Rtwk/R;
  }

#ifdef _G_REWEIGHT_FGM_DEBUG_
  fTestNtp->Fill(-q2,wght);
#endif

  return wght;
}
//_______________________________________________________________________________________
double GReWeightFGM::RewCCQEMomDistroFGtoSF(const EventRecord & event)
{
  bool momdistro_tweaked = (TMath::Abs(fMomDistroTwkDial) > controls::kASmallNum);
  if(!momdistro_tweaked) return 1.;

  bool is_qe = event.Summary()->ProcInfo().IsQuasiElastic();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if(!is_qe || !is_cc) return 1.;

  GHepParticle * tgtnucleus = event.TargetNucleus();
  if(!tgtnucleus) return 1.; // scattering off free-nucleon

  GHepParticle * hitnucleon = event.HitNucleon();
  if(!hitnucleon) return 1.;

  const double kPmax = 0.5;
  double p = hitnucleon->P4()->Vect().Mag();
  if(p > kPmax) return 1.;

  TH1D * hfg = 0;
  TH1D * hsf = 0;

  int tgtpdg = tgtnucleus -> Pdg();
  int nucpdg = hitnucleon -> Pdg();

  map<int, TH1D*> & mapfg = pdg::IsNeutron(nucpdg) ? fMapFGn : fMapFGp;
  map<int, TH1D*> & mapsf = pdg::IsNeutron(nucpdg) ? fMapSFn : fMapSFp;

  map<int, TH1D*>::iterator it;
  it = mapfg.find(tgtpdg);
  if(it != mapfg.end()) { hfg = it->second; }
  it = mapsf.find(tgtpdg);
  if(it != mapsf.end()) { hsf = it->second; }

  bool have_weight_func = (hfg!=0) && (hsf!=0);
  if(!have_weight_func) {
     const int kNEv  = 20000;
     const int kNP   = 500;
     hfg = new TH1D("","",kNP,0.,kPmax);
     hsf = new TH1D("","",kNP,0.,kPmax);
     hfg -> SetDirectory(0);
     hsf -> SetDirectory(0);
     const Target & tgt = event.Summary()->InitState().Tgt();
     bool ok = true;
     for(int iev=0; iev<kNEv; iev++) {
       ok = fFG->GenerateNucleon(tgt);
       if(!ok) {
         delete hfg;
         delete hsf;
         return 1.;
       }
       hfg->Fill(fFG->Momentum());
     }//fg
     for(int iev=0; iev<kNEv; iev++) {
       ok = fSF->GenerateNucleon(tgt);
       if(!ok) {
         delete hfg;
         delete hsf;
         return 1.;
       }
       hsf->Fill(fSF->Momentum());
     }//sf
     hfg->Scale(1. / hfg->Integral("width"));
     hsf->Scale(1. / hsf->Integral("width"));
     mapfg.insert(map<int,TH1D*>::value_type(tgtpdg,hfg));
     mapsf.insert(map<int,TH1D*>::value_type(tgtpdg,hsf));
  }//create & store momentum distributions


  double f_fg = hfg->GetBinContent( hfg->FindBin(p) );
  double f_sf = hsf->GetBinContent( hsf->FindBin(p) );
  double dial = fMomDistroTwkDial;
  double wght = (f_sf * dial + f_fg * (1-dial)) / f_fg;

  return wght;
}
//_______________________________________________________________________________________
void GReWeightFGM::Init(void)
{
  fKFTwkDial        = 0.;
  fMomDistroTwkDial = 0.;

  AlgFactory * algf = AlgFactory::Instance();

  fFG = dynamic_cast<const NuclearModelI*> (
    algf->GetAlgorithm("genie::FGMBodekRitchie","Default"));
  fSF = dynamic_cast<const NuclearModelI*> (
    algf->GetAlgorithm("genie::SpectralFunc","Default"));

#ifdef _G_REWEIGHT_FGM_DEBUG_
  fTestFile = new TFile("./fgm_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","Q2:wght");
#endif
}
//_______________________________________________________________________________________
