//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include <cstdlib>

// GENIE/Generator includes
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Messenger/Messenger.h"

// GENIE/Reweight includes
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

GSystUncertainty * GSystUncertainty::fInstance = 0;
//____________________________________________________________________________
GSystUncertainty::GSystUncertainty()
{
//  fInstance = 0;
}
//____________________________________________________________________________
GSystUncertainty::~GSystUncertainty()
{
  fInstance = 0;
}
//____________________________________________________________________________
GSystUncertainty * GSystUncertainty::Instance()
{
  if(fInstance == 0) {
    LOG("ReW", pINFO) << "GSystUncertainty late initialization";
    static GSystUncertainty::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new GSystUncertainty;
    fInstance->SetDefaults();
  }
  return fInstance;
}
//____________________________________________________________________________
double GSystUncertainty::OneSigmaErr(GSyst_t s, int sign) const
{
  const std::map<GSyst_t, GSystUncertaintyTable::MapEntry>&
    error_map = fTable->GetErrorsMap();

  // Search the map for the uncertainties associated with this
  // tweak dial. If we couldn't find an entry, then just
  // return zero.
  std::map<GSyst_t, GSystUncertaintyTable::MapEntry>::const_iterator
    it = error_map.find( s );
  if ( it == error_map.end() ) return 0.;

  // Otherwise, get the appropriate one-sigma error
  // (plus, minus, or the mean of the two)
  if (sign > 0) return it->second.PlusOneSigmaErr();
  else if (sign < 0) return it->second.MinusOneSigmaErr();

  // Handle default argument (sign=0)
  // Case added for compatibility purposes since most existing weight
  // calcutators call GSystUncertainty::OneSigmaErr(GSyst_t) and the error
  // on most GSyst_t params is symmetric.
  double plus_err = it->second.PlusOneSigmaErr();
  double minus_err = it->second.MinusOneSigmaErr();
  double err = 0.5 * (plus_err + minus_err);
  return err;
}
//____________________________________________________________________________
void GSystUncertainty::SetUncertainty(
   GSyst_t s, double plus_err, double minus_err)
{
  std::map<GSyst_t, GSystUncertaintyTable::MapEntry>&
    error_map = (*fTable->GetErrorsMapPtr());

  std::string knob_name = GSyst::AsString( s );

  LOG("ReW", pINFO) << "Setting uncertainties for"
    << " Reweight tweak dial " << knob_name << " to +"
    << plus_err << ", -" << minus_err;

  error_map[s] = GSystUncertaintyTable::MapEntry(plus_err, minus_err);
}
//____________________________________________________________________________
void GSystUncertainty::SetDefaults(void)
{
  AlgFactory* alg_f = AlgFactory::Instance();
  fTable = dynamic_cast< GSystUncertaintyTable* >(
    alg_f->AdoptAlgorithm("genie::rew::GSystUncertaintyTable", "Default"));
  assert( fTable );
}
//____________________________________________________________________________
