//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <c.andreopoulos \at cern.ch>
         University of Liverpool 
*/
//____________________________________________________________________________

// GENIE/Generator includes
#include "Framework/Messenger/Messenger.h"

// GENIE/Reweight includes
#include "RwFramework/GSystSet.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GSystSet::GSystSet()
{

}
//_______________________________________________________________________________________
GSystSet::GSystSet(const GSystSet & syst)
{
  this->Copy(syst);
}
//_______________________________________________________________________________________
GSystSet::~GSystSet()
{
  fSystematics.clear();
}
//_______________________________________________________________________________________
void GSystSet::Init(GSyst_t syst, double init, double min, double max, double step)
{
  if(syst == kNullSystematic) return;

  if(this->Added(syst)) {
    this->Remove(syst);    
  }

  GSystInfo * syst_info = new GSystInfo(init,min,max,step);
  fSystematics.insert( map<GSyst_t, GSystInfo*>::value_type(syst, syst_info) );     
}
//_______________________________________________________________________________________
void GSystSet::Remove(GSyst_t syst)
{
  fSystematics.erase(syst);
}
//_______________________________________________________________________________________
int GSystSet::Size(void) const
{
  return fSystematics.size();
}
//_______________________________________________________________________________________
bool GSystSet::Added(GSyst_t syst) const
{
  return (fSystematics.find(syst) != fSystematics.end());
}
//_______________________________________________________________________________________
vector<genie::rew::GSyst_t> GSystSet::AllIncluded(void)
{
  vector<GSyst_t> svec;

  map<GSyst_t, GSystInfo*>::const_iterator it = fSystematics.begin();
  for( ; it != fSystematics.end(); ++it) {
    GSyst_t syst = it->first;
    svec.push_back(syst);
  }
  return svec;
}
//_______________________________________________________________________________________
const GSystInfo * GSystSet::Info(GSyst_t syst) const
{
  if ( this->Added(syst) ) {
    return fSystematics.find(syst)->second;
  }
  return 0;
}
//_______________________________________________________________________________________
void GSystSet::Set(GSyst_t syst, double val)
{
  if ( this->Added(syst) ) {
    fSystematics[syst]->CurValue = val;
  }
  else {
    this->Init(syst);
    this->Set(syst,val);
  }
}
//_______________________________________________________________________________________
void GSystSet::Print(void)
{
  LOG("ReW", pNOTICE) 
     << "Considering " << this->Size() << " systematics";
				    
  vector<genie::rew::GSyst_t> svec = this->AllIncluded();

  unsigned int i=0;
  vector<genie::rew::GSyst_t>::const_iterator it = svec.begin();
  for( ; it != svec.end(); ++it) {
    GSyst_t syst = *it;
    LOG("ReW", pNOTICE) << "(" << i++ << ") : " << GSyst::AsString(syst);
  }
}
//_______________________________________________________________________________________
void GSystSet::Copy(const GSystSet & syst_set)
{
  return fSystematics.clear();

  map<GSyst_t, GSystInfo*>::const_iterator it = syst_set.fSystematics.begin();
  for( ; it != syst_set.fSystematics.end(); ++it) {
    GSyst_t     syst       = it->first;
    GSystInfo * syst_info  = it->second;

    double cur  = syst_info->CurValue;
    double init = syst_info->InitValue;
    double min  = syst_info->MinValue;
    double max  = syst_info->MaxValue;
    double step = syst_info->Step;

    this->Init(syst,init,min,max,step);
    this->Set(syst,cur);
  }
}
//_______________________________________________________________________________________
