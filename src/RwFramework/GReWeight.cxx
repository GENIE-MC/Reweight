//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool 
*/
//____________________________________________________________________________

#include <vector>
#include <algorithm>

#include <TMath.h>
#include <TString.h>

// GENIE/Generator includes
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/RunOpt.h"
// GENIE/Reweight includes
#include "RwFramework/GReWeight.h"

using std::vector;

using namespace genie;
using namespace genie::rew;

//____________________________________________________________________________
GReWeight::GReWeight()
{
  // Disable cacheing that interferes with event reweighting
  RunOpt::Instance()->EnableBareXSecPreCalc(false);
}
//____________________________________________________________________________
GReWeight::~GReWeight()
{
  this->CleanUp();
}
//____________________________________________________________________________
void GReWeight::AdoptWghtCalc(string name, GReWeightI* wcalc)
{
  if(!wcalc) return;

  fWghtCalc.insert(map<string, GReWeightI*>::value_type(name,wcalc));
  
  if (std::find(fWghtCalcNames.begin(),fWghtCalcNames.end(),name) == fWghtCalcNames.end()) {
    fWghtCalcNames.push_back(name);
  }
}
//____________________________________________________________________________
GReWeightI* GReWeight::WghtCalc(string name)
{ 
  map<string, GReWeightI*>::iterator iter = fWghtCalc.find(name);
  if(iter != fWghtCalc.end()) return iter->second;
  
  return 0;
}
//____________________________________________________________________________
GSystSet & GReWeight::Systematics(void)
{ 
  return fSystSet; 
}
//____________________________________________________________________________
void GReWeight::Reconfigure(void)
{
  LOG("ReW", pNOTICE) << "Reconfiguring ...";

  vector<genie::rew::GSyst_t> svec = fSystSet.AllIncluded();

  map<string, GReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {

      GReWeightI * wcalc = it->second;

      vector<genie::rew::GSyst_t>::const_iterator parm_iter = svec.begin();
      for( ; parm_iter != svec.end(); ++parm_iter) {
          GSyst_t syst = *parm_iter;
          double val = fSystSet.Info(syst)->CurValue;
          wcalc->SetSystematic(syst, val);
      }//params

      wcalc->Reconfigure();

  }//weight calculators

  LOG("ReW", pDEBUG) << "Done reconfiguring";
}
//____________________________________________________________________________
double GReWeight::CalcWeight(const genie::EventRecord & event) 
{
// calculate weight for all tweaked physics parameters
//
  double weight = 1.0;
  map<string, GReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {
    GReWeightI * wcalc = it->second;
    double w = wcalc->CalcWeight(event); 
    LOG("ReW", pNOTICE) 
       << "Calculator: " << it->first << " => wght = " << w;	
    weight *= w;
  }
  return weight;
}
//____________________________________________________________________________
void GReWeight::CleanUp(void)
{
  map<string, GReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {
    GReWeightI * rw = it->second;
    if(rw) {
      delete rw;
      rw=0;
    }
  }
  fWghtCalc.clear();
}
//____________________________________________________________________________
void GReWeight::Print()
{
  vector<genie::rew::GSyst_t> syst_vec = this->Systematics().AllIncluded();
  int vec_size = syst_vec.size();

  LOG("ReW", pNOTICE) << "Current set of systematic params:";	
  for(int i = 0 ; i < vec_size ; i ++){
     LOG("ReW", pNOTICE) 
        << " --o "  << GSyst::AsString(syst_vec[i])
        << " is set at " << this->Systematics().Info(syst_vec[i])->CurValue;
  }		       	        
}
//____________________________________________________________________________
const std::vector<std::string> & GReWeight::WghtCalcNames() const
{
  return fWghtCalcNames;
}
//____________________________________________________________________________
