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

#include <TMath.h>
#include <TFile.h>
#include <TNtupleD.h>

// GENIE/Generator includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightNuXSecCCQEvec.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecCCQEvec::GReWeightNuXSecCCQEvec() :
GReWeightModel("CCQEvec")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQEvec::~GReWeightNuXSecCCQEvec()
{
#ifdef _G_REWEIGHT_CCQE_VEC_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQEvec::IsHandled(GSyst_t syst) const
{
   switch(syst) {
    case ( kXSecTwkDial_VecFFCCQEshape ) :
       return true;
       break;
    default:
       return false;
       break;
   }
   return false;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQEvec::AppliesTo(ScatteringType_t type, bool is_cc) const
{
  if (type==kScQuasiElastic && is_cc) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEvec::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

  switch(syst) {
    case ( kXSecTwkDial_VecFFCCQEshape ) :
       fFFTwkDial = twk_dial;
       break;
    default:
       return;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEvec::Reset(void)
{
  fFFTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEvec::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEvec::CalcWeight(const genie::EventRecord & event)
{
  bool tweaked = (TMath::Abs(fFFTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  Interaction * interaction = event.Summary();

  bool is_qe = interaction->ProcInfo().IsQuasiElastic();
  bool is_cc = interaction->ProcInfo().IsWeakCC();
  if(!is_qe || !is_cc) return 1.;

  bool charm = interaction->ExclTag().IsCharmEvent(); // skip CCQE charm
  if(charm) return 1.;

  int nupdg = event.Probe()->Pdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  //
  // Calculate weight
  // Input tweaking dial changes elastic nucleon form factors
  // (twk dial: 0 -> default/BBA, twk dial: 1 -> dipole).
  // Calculated weight includes `shape only effect in dsigma/dQ2
  // (normalized to constant integrated cross section)
  //

  const KinePhaseSpace_t phase_space = event.DiffXSecVars();

  interaction->KinePtr()->UseSelectedKinematics();

  if (phase_space == kPSQ2fE) {
    interaction->SetBit(kIAssumeFreeNucleon);
  }

  double old_xsec   = event.DiffXSec();
  if (!fUseOldWeightFromFile || fNWeightChecksDone < fNWeightChecksToDo) {
    double calc_old_xsec = fXSecModel_bba->XSec(interaction, phase_space);
    if (fNWeightChecksDone < fNWeightChecksToDo) {
      if (std::abs(calc_old_xsec - old_xsec)/old_xsec > controls::kASmallNum) {
        LOG("ReW",pWARN) << "Warning - default dxsec does not match dxsec saved in tree. Does the config match?";
      }
      fNWeightChecksDone++;
    }
    if(!fUseOldWeightFromFile) {
      old_xsec = calc_old_xsec;
    }
  }

  double dial                = fFFTwkDial;
  double old_weight          = event.Weight();
  double dpl_xsec            = fXSecModel_dpl->XSec(interaction, phase_space);

  double def_integrated_xsec = fXSecIntegrator_bba->Integrate(fXSecModel_bba, interaction);
  double dpl_integrated_xsec = fXSecIntegrator_dpl->Integrate(fXSecModel_dpl, interaction);

  LOG("ReW", pERROR) << "def_integ = " << def_integrated_xsec;
  LOG("ReW", pERROR) << "dpl_integ = " << dpl_integrated_xsec;

  if ( def_integrated_xsec <= 0. || dpl_integrated_xsec <= 0. ) {
    LOG("ReW", pWARN) << "Non-positive total cross section encountered in"
      << " GReWeightNuXSecCCQEvec::CalcWeight()";
    return 1.;
  }

  double def_ratio = old_xsec / def_integrated_xsec;
  double dpl_ratio = dpl_xsec / dpl_integrated_xsec;

//  assert(def_ratio > 0.);
//  if(def_ratio <= 0) return 1.;

  double weight = old_weight * (dial * dpl_ratio + (1. - dial)*def_ratio) / def_ratio;

#ifdef _G_REWEIGHT_CCQE_VEC_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  fTestNtp->Fill(
    E,Q2,weight,def_integrated_xsec,dpl_integrated_xsec,def_xsec,dpl_xsec);
#endif

  interaction->KinePtr()->ClearRunningValues();

  if (phase_space == kPSQ2fE) {
    interaction->ResetBit(kIAssumeFreeNucleon);
  }

  return weight;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEvec::Init(void)
{
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * gpl = conf_pool->GlobalParameterList();
  RgAlg xsec_alg = gpl->GetAlg("XSecModel@genie::EventGenerator/QEL-CC");

  AlgId def_id(xsec_alg); // no "default" anymore
  AlgId elff_id(AlgId(xsec_alg).Name(),"Dipole");

  // I can't see why we'd want a non-default model name here, so this bit is unnecessary for now
  //~ if (fManualModelName.size()) {
    //~ def_id   = AlgId(fManualModelName,"Dipole");
    //~ elff_id  = AlgId(fManualModelName,"DipoleELFF");
  //~ }

  AlgFactory * algf = AlgFactory::Instance();

  Algorithm * alg0 = algf->AdoptAlgorithm(def_id);
  fXSecModel_bba = dynamic_cast<XSecAlgorithmI*>(alg0);
  fXSecModel_bba->AdoptSubstructure();

  Algorithm * alg1 = algf->AdoptAlgorithm(elff_id);
  fXSecModel_dpl = dynamic_cast<XSecAlgorithmI*>(alg1);
  fXSecModel_dpl->AdoptSubstructure();

  // Get the Algorithm objects that should be used to integrate the cross
  // sections. Use the "ReweightShape" configuration, which turns off averaging
  // over the initial state nucleon distribution.
  AlgId alg0_integ_ID( alg0->GetConfig().GetAlg("XSec-Integrator").name,
    "ReweightShape");

  fXSecIntegrator_bba = dynamic_cast<XSecIntegratorI*>(
    algf->AdoptAlgorithm(alg0_integ_ID));

  assert( fXSecIntegrator_bba );

  AlgId alg1_integ_ID( alg1->GetConfig().GetAlg("XSec-Integrator").name,
    "ReweightShape");

  fXSecIntegrator_dpl = dynamic_cast<XSecIntegratorI*>(
    algf->AdoptAlgorithm(alg1_integ_ID));

  assert( fXSecIntegrator_dpl );

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  fFFTwkDial = 0.;

#ifdef _G_REWEIGHT_CCQE_VEC_DEBUG_
  fTestFile = new TFile("./ccqevec_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:Q2:wght:sig0:sig:dsig0:dsig");
#endif

}
//_______________________________________________________________________________________
