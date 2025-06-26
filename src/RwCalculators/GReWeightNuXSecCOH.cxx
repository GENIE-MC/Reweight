//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London
*/
//____________________________________________________________________________

#include <TMath.h>

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
#include "RwCalculators/GReWeightNuXSecCOH.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecCOH::GReWeightNuXSecCOH() :
GReWeightModel("CCCoh"),
fManualModelName(),
fManualModelType()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCOH::GReWeightNuXSecCOH(std::string model, std::string type) :
GReWeightModel("CCCoh"),
fManualModelName(model),
fManualModelType(type)
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCOH::~GReWeightNuXSecCOH()
{

}
//_______________________________________________________________________________________
bool GReWeightNuXSecCOH::IsHandled(GSyst_t syst) const
{
   bool handle;

   switch(syst) {
     case ( kXSecTwkDial_MaCOHpi ) :
     case ( kXSecTwkDial_R0COHpi ) :
     case ( kXSecTwkDial_NormCCCOHpi ) :
     case ( kXSecTwkDial_NormNCCOHpi ) :
       handle = true;
       break;
     default:
       handle = false;
       break;
   }
   return handle;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCOH::AppliesTo(const EventRecord & event) const
{
  auto type = event.Summary()->ProcInfo().ScatteringTypeId();
  if (type==kScCoherentProduction) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOH::SetSystematic(GSyst_t syst, double twk_dial)
{
   switch(syst) {
     case ( kXSecTwkDial_MaCOHpi ) :
       fMaTwkDial = twk_dial;
       break;
     case ( kXSecTwkDial_R0COHpi ) :
       fR0TwkDial = twk_dial;
       break;
    case ( kXSecTwkDial_NormCCCOHpi ) :
      fCCCOHNormTwkDial = twk_dial;
      break;
    case ( kXSecTwkDial_NormNCCOHpi ) :
      fNCCOHNormTwkDial = twk_dial;
      break;
     default:
       break;
   }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOH::Reset(void)
{
  fMaTwkDial   = 0.;
  fMaCurr      = fMaDef;
  fR0TwkDial   = 0.;
  fR0Curr      = fR0Def;

  fCCCOHNormTwkDial = 0.;
  fNCCOHNormTwkDial = 0.;

  fCurCOHNormCC = 1.;
  fCurCOHNormNC = 1.;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOH::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();

  double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaCOHpi);
  double fracerr_r0 = fracerr->OneSigmaErr(kXSecTwkDial_R0COHpi);

  fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
  fR0Curr = fR0Def * (1. + fR0TwkDial * fracerr_r0);

  fMaCurr = TMath::Max(0., fMaCurr  );
  fR0Curr = TMath::Max(0., fR0Curr  );

  Registry & r = *fXSecModelConfig;

  r.Set(fMaPath, fMaCurr);
  r.Set(fR0Path, fR0Curr);

  fXSecModel->Configure(r);

  // Note: this assumes that the error is symmetric.
  // TODO: consider changing this to handle asymmetric errors on the
  // normalization
  double frac_err_cc_norm = fracerr->OneSigmaErr( kXSecTwkDial_NormCCCOHpi );
  double frac_err_nc_norm = fracerr->OneSigmaErr( kXSecTwkDial_NormNCCOHpi );

  fCurCOHNormCC = std::max( 0., 1. + fCCCOHNormTwkDial * frac_err_cc_norm );
  fCurCOHNormNC = std::max( 0., 1. + fNCCOHNormTwkDial * frac_err_nc_norm );

//LOG("ReW", pDEBUG) << *fXSecModel;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCOH::CalcWeight(const genie::EventRecord & event)
{
  Interaction * interaction = event.Summary();

  bool is_coh = interaction->ProcInfo().IsCoherentProduction();
  if(!is_coh) return 1.;

  bool xsec_tweaked =
      (TMath::Abs(fMaTwkDial) > controls::kASmallNum) ||
      (TMath::Abs(fR0TwkDial) > controls::kASmallNum);

  bool norm_tweaked = (TMath::Abs(fCCCOHNormTwkDial) > controls::kASmallNum) ||
      (TMath::Abs(fNCCOHNormTwkDial) > controls::kASmallNum);

  bool tweaked = xsec_tweaked || norm_tweaked;

  if(!tweaked) return 1.0;

  bool is_cc  = interaction->ProcInfo().IsWeakCC();
  bool is_nc  = interaction->ProcInfo().IsWeakNC();
  if(is_cc && !fRewCC) return 1.;
  if(is_nc && !fRewNC) return 1.;

  int nupdg = interaction->InitState().ProbePdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  double old_weight = event.Weight();
  double new_weight = old_weight;
  if ( norm_tweaked ) {
    if ( is_cc ) new_weight *= fCurCOHNormCC;
    else if ( is_nc ) new_weight *= fCurCOHNormNC;
  }

  if ( !xsec_tweaked ) return new_weight;

  interaction->KinePtr()->UseSelectedKinematics();

  const KinePhaseSpace_t phase_space = event.DiffXSecVars();

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

  double new_xsec   = fXSecModel->XSec(interaction, phase_space);
  new_weight *= (new_xsec/old_xsec);

  interaction->KinePtr()->ClearRunningValues();

  return new_weight;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOH::Init(void)
{
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * gpl = conf_pool->GlobalParameterList();
// -->  RgAlg xsec_alg = gpl->GetAlg("XSecModel@genie::EventGenerator/COH-CC");
  RgAlg xsec_alg = gpl->GetAlg("XSecModel@genie::EventGenerator/COH-CC-PION");

  AlgId id(xsec_alg);

  AlgId twk_id(id);
  if (fManualModelName.size()) {
    twk_id = AlgId(fManualModelName,fManualModelType);
  }

  AlgFactory * algf = AlgFactory::Instance();
  Algorithm * alg_def = algf->AdoptAlgorithm(id);
  Algorithm * alg_twk = algf->AdoptAlgorithm(twk_id);

  fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg_twk);
  fXSecModel->AdoptSubstructure();

  fXSecModelDef = dynamic_cast<XSecAlgorithmI*>(alg_def);
  fXSecModelDef->AdoptSubstructure();

  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
//LOG("ReW", pNOTICE) << *fXSecModelConfig;

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);
  this->RewCC     (true);
  this->RewNC     (true);

  this->SetMaPath("COH-Ma");
  this->SetR0Path("COH-Ro");

  fMaTwkDial   = 0.;
  fMaDef       = fXSecModelConfig->GetDouble(fMaPath);
  fMaCurr      = fMaDef;
  fR0TwkDial   = 0.;
  fR0Def       = fXSecModelConfig->GetDouble(fR0Path);
  fR0Curr      = fR0Def;

  fCCCOHNormTwkDial = 0.;
  fNCCOHNormTwkDial = 0.;

  fCurCOHNormCC = 1.;
  fCurCOHNormNC = 1.;
}
//_______________________________________________________________________________________
