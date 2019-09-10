//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Messenger/Messenger.h"

#include "RwFramework/GSystUncertaintyTable.h"

using namespace genie;
using namespace genie::rew;

//____________________________________________________________________________
GSystUncertaintyTable::GSystUncertaintyTable() :
Algorithm("genie::rew::GSystUncertaintyTable")
{

}
//____________________________________________________________________________
GSystUncertaintyTable::GSystUncertaintyTable(std::string config) :
Algorithm("genie::rew::GSystUncertaintyTable", config)
{

}
//____________________________________________________________________________
GSystUncertaintyTable::~GSystUncertaintyTable()
{

}
//____________________________________________________________________________
void GSystUncertaintyTable::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GSystUncertaintyTable::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GSystUncertaintyTable::LoadConfig(void)
{
  // Load one-sigma uncertainties from the Registry for every defined
  // tweak dial. If there's a missing entry, set both the plus and minus
  // one-sigma values to zero.
  double plus_err, minus_err;
  int syst_knob_index = 1; // Index zero corresponds to kNullSystematic
  int last_index = static_cast<int>( kNTwkDials );

  while ( syst_knob_index < last_index ) {

    GSyst_t syst_knob = static_cast<GSyst_t>( syst_knob_index );

    std::string knob_name = GSyst::AsString( syst_knob );
    if ( !knob_name.empty() && knob_name != "-") {
      GetParamDef(knob_name + "@PlusOneSigma", plus_err, 0.);
      GetParamDef(knob_name + "@MinusOneSigma", minus_err, 0.);

      fErrorsMap[ syst_knob ] = MapEntry( plus_err, minus_err );

      LOG("ReW", pNOTICE) << "Reweight parameter "
       << knob_name << " has one-sigma uncertainties +" << plus_err
       << ", -" << minus_err;
    }
    else if ( knob_name == "-" ) {
      // This string is used for the now-invalid FrElas_N and FrElas_Pi
      // knobs, so just ignore it
      LOG("ReW", pNOTICE) << "Skipping deprecated tweak dial knob for "
       << "GSyst_t enum value " << syst_knob;
    }
    else {
      LOG("ReW", pWARN) << "Unrecognized"
        << " GSyst_t value " << syst_knob << " encountered.";
    }

    ++syst_knob_index;
  }

}
