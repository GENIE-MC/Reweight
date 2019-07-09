//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

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
#include "Framework/Registry/Registry.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightNuXSecNCEL.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecNCEL::GReWeightNuXSecNCEL() :
GReWeightModel("NCEl"),
fManualModelName(),
fManualModelType()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecNCEL::GReWeightNuXSecNCEL(std::string model, std::string type) :
GReWeightModel("NCEl"),
fManualModelName(model),
fManualModelType(type)
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecNCEL::~GReWeightNuXSecNCEL()
{
#ifdef _G_REWEIGHT_NCEL_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightNuXSecNCEL::IsHandled(GSyst_t syst) const
{
   bool handle;
   switch(syst) {
     case ( kXSecTwkDial_MaNCEL  ) : handle = true; break;
     case ( kXSecTwkDial_EtaNCEL ) : handle = true; break;
     default:
          handle = false;
          break;
   }
   return handle;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecNCEL::AppliesTo(ScatteringType_t type, bool is_cc) const
{
  if (type==kScQuasiElastic && !is_cc) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecNCEL::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

  switch(syst) {
    case ( kXSecTwkDial_MaNCEL ) :
      fMaTwkDial = twk_dial;
      break;
    case ( kXSecTwkDial_EtaNCEL ) :
      fEtaTwkDial = twk_dial;
      break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecNCEL::Reset(void)
{
  fMaTwkDial   = 0.;
  fMaCurr      = fMaDef;
  fEtaTwkDial  = 0.;
  fEtaCurr     = fEtaDef;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecNCEL::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();

  int sign_matwk  = utils::rew::Sign(fMaTwkDial );
  int sign_etatwk = utils::rew::Sign(fEtaTwkDial);

  double fracerr_ma  = fracerr->OneSigmaErr(kXSecTwkDial_MaNCEL,  sign_matwk );
  double fracerr_eta = fracerr->OneSigmaErr(kXSecTwkDial_EtaNCEL, sign_etatwk);

  fMaCurr  = fMaDef  * (1. + fMaTwkDial  * fracerr_ma);
  fEtaCurr = fEtaDef * (1. + fEtaTwkDial * fracerr_eta);

  fMaCurr  = TMath::Max(0., fMaCurr  );
  fEtaCurr = TMath::Max(0., fEtaCurr );

  Registry & r = *fXSecModelConfig;

  r.Set(fMaPath,  fMaCurr );
  r.Set(fEtaPath, fEtaCurr);
  fXSecModel->Configure(r);

//LOG("ReW, pDEBUG) << *fXSecModel;
}
//_______________________________________________________________________________________
double GReWeightNuXSecNCEL::CalcWeight(const genie::EventRecord & event)
{
  Interaction * interaction = event.Summary();

  bool is_qe = interaction->ProcInfo().IsQuasiElastic();
  bool is_nc = interaction->ProcInfo().IsWeakNC();
  if(!is_qe || !is_nc) return 1.;

  int nupdg = event.Probe()->Pdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  bool tweaked_ma  = (TMath::Abs(fMaTwkDial ) > controls::kASmallNum);
  bool tweaked_eta = (TMath::Abs(fEtaTwkDial) > controls::kASmallNum);
  bool tweaked     = tweaked_ma || tweaked_eta;
  if(!tweaked) return 1.0;

  interaction->KinePtr()->UseSelectedKinematics();

  const KinePhaseSpace_t phase_space = event.DiffXSecVars();

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

  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, phase_space );
  double new_weight = old_weight * (new_xsec/old_xsec);

//LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
//LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
//LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
//LOG("ReW", pDEBUG) << "new weight = " << new_weight;

#ifdef _G_REWEIGHT_NCEL_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  fTestNtp->Fill(E,Q2,new_weight);
#endif

  interaction->KinePtr()->ClearRunningValues();

  if (phase_space == kPSQ2fE) {
    interaction->ResetBit(kIAssumeFreeNucleon);
  }

  return new_weight;
}
//_______________________________________________________________________________________
void GReWeightNuXSecNCEL::Init(void)
{
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * gpl = conf_pool->GlobalParameterList();
  RgAlg xsec_alg = gpl->GetAlg("XSecModel@genie::EventGenerator/QEL-NC");

  AlgId id(xsec_alg);

  AlgId twk_id(id);
  if (fManualModelName.size()) {
    twk_id = AlgId(fManualModelName,fManualModelType);
  }

  AlgFactory * algf = AlgFactory::Instance();

  Algorithm * algdef = algf->AdoptAlgorithm(id);
  fXSecModelDef = dynamic_cast<XSecAlgorithmI*>(algdef);
  fXSecModelDef->AdoptSubstructure();

  Algorithm * alg = algf->AdoptAlgorithm(twk_id);
  fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg);
  fXSecModel->AdoptSubstructure();

  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
//LOG("ReW", pDEBUG) << *fXSecModelConfig;

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

/*
  this->SetMaPath ("Ma");
  this->SetEtaPath("Eta");
*/

  this->SetMaPath( "QEL-Ma");
  this->SetEtaPath( "EL-Axial-Eta" );

  fMaTwkDial   = 0.;
  fMaDef       = fXSecModelConfig->GetDouble(fMaPath);
  fMaCurr      = fMaDef;
  fEtaTwkDial  = 0.;
  fEtaDef      = fXSecModelConfig->GetDouble(fEtaPath);
  fEtaCurr     = fEtaDef;

#ifdef _G_REWEIGHT_NCEL_DEBUG_
  fTestFile = new TFile("./ncel_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:Q2:wght");
#endif
}
//_______________________________________________________________________________________
