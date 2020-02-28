//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TParticlePDG.h>
#include <TDecayChannel.h>

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
GReWeightModel("DeltaradAngle")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightDeltaradAngle::~GReWeightDeltaradAngle()
{
#ifdef _G_REWEIGHT_RESDEC_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_________________
bool GReWeightDeltaradAngle::IsHandled(GSyst_t syst) const
{
  switch(syst) {
  case( kRDcyTwkDial_Theta_Delta2NRad):
        return true;
        break;
     default:
        return false;
        break;
  }
  return false;
}
//_______________________________________________________________________________________
//_______________________________________________________________________________________
bool GReWeightDeltaradAngle::AppliesTo (ScatteringType_t type, bool /*is_cc*/) const
{
  if (type==kScResonant) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightDeltaradAngle::SetSystematic(GSyst_t syst, double twk_dial)
{
  switch(syst) {
       case (kRDcyTwkDial_Theta_Delta2NRad ) :
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
  Interaction * interaction = event.Summary();

  bool is_res = interaction->ProcInfo().IsResonant();
  if(!is_res) return 1.;

  bool is_cc  = interaction->ProcInfo().IsWeakCC();
  bool is_nc  = interaction->ProcInfo().IsWeakNC();
  if(is_cc && !fRewCC) return 1.;
  if(is_nc && !fRewNC) return 1.;

  int nupdg = interaction->InitState().ProbePdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  double wght =
    this->RewThetaDelta2NRad(event);

  return wght;
}
//_______________________________________________________________________________________
double GReWeightDeltaradAngle::RewThetaDelta2NRad(const EventRecord & event)
{
  bool tweaked = (TMath::Abs(fThetaDelta2NRadTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  bool is_Delta_rad = false;
  int ir = -1; // resonance position
  int ig = -1; // gamma position
  int i  =  0;
  GHepParticle * p = 0;
  TIter iter(&event);
  while((p=(GHepParticle*)iter.Next())) {
    bool is_Deltap = (p->Pdg()==kPdgP33m1232_DeltaP);
    if(is_Deltap) {
      ir = i;
      int fd = p->FirstDaughter();
      int ld = p->LastDaughter();
      int nd = 1 + ld - fd;
      if(nd==2) {
        int fpdg = event.Particle(fd)->Pdg();
        int lpdg = event.Particle(ld)->Pdg();
        if(fpdg==kPdgGamma && lpdg==kPdgProton) { is_Delta_rad = true; ig = fd; }
        if(lpdg==kPdgGamma && fpdg==kPdgProton) { is_Delta_rad = true; ig = ld; }
      }
    }
    
    if(is_Delta_rad) break;
    i++;
  }
   if(!is_Delta_rad) return 1.;
   LOG("ReW", pDEBUG) << "A Delta++ -> p photon event:";
   LOG("ReW", pDEBUG) << "Resonance is at position: " << ir;
   LOG("ReW", pDEBUG) << "Gamma is at position: " << ig;

    // Get Delta and pi+ 4-momentum vectors
  TLorentzVector p4res(*event.Particle(ir)->P4());
  TLorentzVector p4photonp(*event.Particle(ig)->P4());

  // Boost pi+ to the Delta CM
  TVector3 bv = -1*p4res.BoostVector();
  p4photonp.Boost(bv);

  //attempt at weight calculation
  double P1=1/3.14; //isotropic default normalization
  double costheta = p4photonp.Vect().CosTheta(); 
  double P2       = 2/3.14 * (costheta*costheta);//need a new function
  double dial     = fThetaDelta2NRadTwkDial;
  double wght = 1.;
  wght = (dial*P2+(1-dial)*P1)/P1;
   LOG("ReW", pDEBUG)
       << "Photon Cos(ThetaCM) = " << costheta << ", weight = " << wght;

  #ifdef _G_REWEIGHT_RESDEC_DEBUG_
  fTestNtp->Fill(costheta,wght);
#endif

  return wght;
    }
//_______________________________________________________________________________________
void GReWeightDeltaradAngle::Init(void)
{
  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);
  this->RewCC     (true);
  this->RewNC     (true);

  fThetaDelta2NRadTwkDial = 0.;

  const int respdgarray[] = {
    kPdgP33m1232_DeltaM, kPdgP33m1232_Delta0, kPdgP33m1232_DeltaP, kPdgP33m1232_DeltaPP,
    kPdgS31m1620_DeltaM, kPdgS31m1620_Delta0, kPdgS31m1620_DeltaP, kPdgS31m1620_DeltaPP,
    kPdgD33m1700_DeltaM, kPdgD33m1700_Delta0, kPdgD33m1700_DeltaP, kPdgD33m1700_DeltaPP,
    kPdgF35m1905_DeltaM, kPdgF35m1905_Delta0, kPdgF35m1905_DeltaP, kPdgF35m1905_DeltaPP,
    kPdgP31m1910_DeltaM, kPdgP31m1910_Delta0, kPdgP31m1910_DeltaP, kPdgP31m1910_DeltaPP,
    kPdgP33m1920_DeltaM, kPdgP33m1920_Delta0, kPdgP33m1920_DeltaP, kPdgP33m1920_DeltaPP,
    kPdgF37m1950_DeltaM, kPdgF37m1950_Delta0, kPdgF37m1950_DeltaP, kPdgF37m1950_DeltaPP,
    kPdgP11m1440_N0, kPdgP11m1440_NP,
    kPdgD13m1520_N0, kPdgD13m1520_NP,
    kPdgS11m1535_N0, kPdgS11m1535_NP,
    kPdgS11m1650_N0, kPdgS11m1650_NP,
    kPdgD15m1675_N0, kPdgD15m1675_NP,
    kPdgF15m1680_N0, kPdgF15m1680_NP,
    kPdgD13m1700_N0, kPdgD13m1700_NP,
    kPdgP11m1710_N0, kPdgP11m1710_NP,
    kPdgP13m1720_N0, kPdgP13m1720_NP,
    0
  };

  // init
  const int    kNW   = 30;
  const double kWmin = 0.9;
  const double kWmax = 5.0;
  unsigned int ires=0;
  int respdg = 0;
  while((respdg = respdgarray[ires++])) {
    fMpBR1gammaDef [respdg] = new TH1D("","",kNW,kWmin,kWmax);
    fMpBR1etaDef   [respdg] = new TH1D("","",kNW,kWmin,kWmax);
    fMpBR1gammaDef [respdg] -> SetDirectory(0);
    fMpBR1etaDef   [respdg] -> SetDirectory(0);
  }

  // find corresponding decay channels and store default BR
  PDGLibrary * pdglib = PDGLibrary::Instance();
  ires=0;
  while((respdg = respdgarray[ires++])) {
    TParticlePDG * res = pdglib->Find(respdg);
    if(!res) continue;
    for(int j=0; j<res->NDecayChannels(); j++) {
        TDecayChannel * dch = res->DecayChannel(j);
        int ngamma = 0;
        int neta   = 0;
        double br = dch->BranchingRatio();
        double mt = 0.;
        for(int k=0; k<dch->NDaughters(); k++) {
          int dpdg = dch->DaughterPdgCode(k);
          if(dpdg==kPdgGamma) ngamma++;
          if(dpdg==kPdgEta  ) neta++;
          mt += 0;
        }//decay channel f/s particles
        bool is_1gamma = (ngamma==1);
        bool is_1eta   = (neta  ==1);
        for(int ibin = 1; ibin <= fMpBR1gammaDef[respdg]->GetNbinsX(); ibin++) {
          double W = fMpBR1gammaDef[respdg]->GetBinLowEdge(ibin) +
                     fMpBR1gammaDef[respdg]->GetBinWidth(ibin);
          bool is_allowed = (W>mt);
          if(is_allowed && is_1gamma) { fMpBR1gammaDef[respdg]->Fill(W, br); }
          if(is_allowed && is_1eta  ) { fMpBR1etaDef  [respdg]->Fill(W, br); }
        }//W bins
    }//decay channels
  }//resonances

#ifdef _G_REWEIGHT_RESDEC_DEBUG_
  fTestFile = new TFile("./resdec_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","costheta:wght");
#endif
}
