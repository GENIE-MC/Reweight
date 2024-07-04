//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

// GENIE/Generator includes
#include "Framework/Messenger/Messenger.h"
// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightModel::GReWeightModel(std::string name) :
GReWeightI(),
fUseOldWeightFromFile(true),
fNWeightChecksToDo(20),
fNWeightChecksDone(0),
fFailedWeightCheck(false),
fName(name)
{

}
//_______________________________________________________________________________________
GReWeightModel::~GReWeightModel()
{
  if (fFailedWeightCheck && fUseOldWeightFromFile) {
    LOG("ReW",pWARN) << fName<< ": You used the weights from the files but the"
      <<" check against the calculated weights failed. Your weights are probably wrong!";
  }
}
//_______________________________________________________________________________________
void GReWeightModel::SetNWeightChecks(int n)
{
  fNWeightChecksToDo = n;
}
//_______________________________________________________________________________________
void GReWeightModel::UseOldWeightFromFile(bool should_we)
{
  fUseOldWeightFromFile = should_we;
}
//_______________________________________________________________________________________
