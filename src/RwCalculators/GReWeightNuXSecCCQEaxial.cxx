//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Aaron Meyer <asmeyer \at uchicago.edu>
          University of Chicago/Fermilab
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
#include "RwCalculators/GReWeightNuXSecCCQEaxial.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecCCQEaxial::GReWeightNuXSecCCQEaxial() :
GReWeightModel("CCQEaxial")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQEaxial::~GReWeightNuXSecCCQEaxial()
{
#ifdef _G_REWEIGHT_CCQE_AXFF_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQEaxial::IsHandled(GSyst_t syst) const
{
   switch(syst) {
    case ( kXSecTwkDial_AxFFCCQEshape ) :
       return true;
       break;
    default:
       return false;
       break;
   }
   return false;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQEaxial::AppliesTo(ScatteringType_t type, bool is_cc) const
{
  if (type==kScQuasiElastic && is_cc) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEaxial::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

  switch(syst) {
    case ( kXSecTwkDial_AxFFCCQEshape ) :
       fFFTwkDial = twk_dial;
       break;
    default:
       return;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEaxial::Reset(void)
{
  fFFTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEaxial::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEaxial::CalcWeight(const genie::EventRecord & event)
{
  bool tweaked = (TMath::Abs(fFFTwkDial) > controls::kASmallNum);
  if ( !tweaked ) return 1.;

  Interaction * interaction = event.Summary();

  bool is_qe = interaction->ProcInfo().IsQuasiElastic();
  bool is_cc = interaction->ProcInfo().IsWeakCC();
  if ( !is_qe || !is_cc ) return 1.;

  bool charm = interaction->ExclTag().IsCharmEvent(); // skip CCQE charm
  if ( charm ) return 1.;

  int nupdg = event.Probe()->Pdg();
  if ( nupdg == kPdgNuMu     && !fRewNumu   ) return 1.;
  if ( nupdg == kPdgAntiNuMu && !fRewNumubar) return 1.;
  if ( nupdg == kPdgNuE      && !fRewNue    ) return 1.;
  if ( nupdg == kPdgAntiNuE  && !fRewNuebar ) return 1.;

  //
  // Calculate weight
  // Input tweaking dial changes elastic nucleon form factors
  // (twk dial: 0 -> default/dipole, twk dial: 1 -> zexp).
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
    double calc_old_xsec = fXSecModelDef->XSec(interaction, phase_space);
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

  double old_weight           = event.Weight();
  double dial                 = fFFTwkDial;
  double zexp_xsec            = fXSecModel_zexp->XSec(interaction, phase_space);

  double def_integrated_xsec  = fXSecIntegratorDef->Integrate(fXSecModelDef, interaction);
  double zexp_integrated_xsec = fXSecIntegrator_zexp->Integrate(fXSecModel_zexp, interaction);

  //assert(def_integrated_xsec > 0.);
  //assert(zexp_integrated_xsec > 0.);
  if ( def_integrated_xsec <= 0. || zexp_integrated_xsec <= 0. ) {
    LOG("ReW", pWARN) << "Non-positive total cross section encountered in"
      << " GReWeightNuXSecCCQEaxial::CalcWeight()";
    return 1.;
  }

  double def_ratio  = old_xsec  / def_integrated_xsec;
  double zexp_ratio = zexp_xsec / zexp_integrated_xsec;

  double weight = old_weight * (dial * zexp_ratio + (1. - dial)*def_ratio) / def_ratio;

#ifdef _G_REWEIGHT_CCQE_AXFF_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  fTestNtp->Fill(
    E,Q2,weight,def_integrated_xsec,zexp_integrated_xsec,old_xsec,zexp_xsec);
#endif

  interaction->KinePtr()->ClearRunningValues();

  if (phase_space == kPSQ2fE) {
    interaction->ResetBit(kIAssumeFreeNucleon);
  }

  return weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEaxial::CalcChisq()
{
  double chisq = TMath::Power(fFFTwkDial, 2.);
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEaxial::Init(void)
{
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * gpl = conf_pool->GlobalParameterList();
  RgAlg xsec_alg = gpl->GetAlg("XSecModel@genie::EventGenerator/QEL-CC");

  AlgId dipole_id (AlgId(xsec_alg).Name(),"Dipole");
  AlgId zexp_id   (AlgId(xsec_alg).Name(),"ZExp");

  // I can't see why we'd want a non-default model name here, so this bit is unnecessary for now
  //~ if (fManualModelName.size()) {
    //~ def_id   = AlgId(fManualModelName,"Dipole");
    //~ elff_id  = AlgId(fManualModelName,"DipoleELFF");
  //~ }

  AlgFactory * algf = AlgFactory::Instance();

  Algorithm * alg0 = algf->AdoptAlgorithm(dipole_id);
  fXSecModelDef = dynamic_cast<XSecAlgorithmI*>(alg0);
  fXSecModelDef->AdoptSubstructure();

  Algorithm * alg1 = algf->AdoptAlgorithm(zexp_id);
  fXSecModel_zexp = dynamic_cast<XSecAlgorithmI*>(alg1);
  fXSecModel_zexp->AdoptSubstructure();

  // Get the Algorithm objects that should be used to integrate the cross
  // sections. Use the "ReweightShape" configuration, which turns off averaging
  // over the initial state nucleon distribution.
  AlgId alg0_integ_ID( alg0->GetConfig().GetAlg("XSec-Integrator").name,
    "ReweightShape");

  fXSecIntegratorDef = dynamic_cast<XSecIntegratorI*>(
    algf->AdoptAlgorithm(alg0_integ_ID));

  assert( fXSecIntegratorDef );

  AlgId alg1_integ_ID( alg1->GetConfig().GetAlg("XSec-Integrator").name,
    "ReweightShape");

  fXSecIntegrator_zexp = dynamic_cast<XSecIntegratorI*>(
    algf->AdoptAlgorithm(alg1_integ_ID));

  assert( fXSecIntegrator_zexp );

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  fFFTwkDial = 0.;

#ifdef _G_REWEIGHT_CCQE_AXFF_DEBUG_
  fTestFile = new TFile("./ccqeaxil_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:Q2:wght:sig0:sig:dsig0:dsig");
#endif

}
//_______________________________________________________________________________________
