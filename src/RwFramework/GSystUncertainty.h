//____________________________________________________________________________
/*!

\class    genie::rew::GSystUncertainty

\brief    Singleton class for looking up reweight tweak dial uncertainties

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Sep 1, 2009

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_SYST_UNCERTAINTY_H_
#define _G_SYST_UNCERTAINTY_H_

#include <map>

// GENIE/Reweight includes
#include "RwFramework/GSyst.h"
#include "RwFramework/GSystUncertaintyTable.h"

using std::map;

namespace genie {
namespace rew   {

class GSystUncertainty {

public:

  static GSystUncertainty * Instance (void);

  double OneSigmaErr    (GSyst_t syst, int sign=0) const;
  void   SetUncertainty (GSyst_t syst, double plus_err, double minus_err);

private:

  GSystUncertaintyTable* fTable;

  void SetDefaults(void);

  GSystUncertainty();
  GSystUncertainty(const GSystUncertainty& err);
 ~GSystUncertainty();

  static GSystUncertainty* fInstance;

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (GSystUncertainty::fInstance !=0) {
            delete GSystUncertainty::fInstance;
            GSystUncertainty::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

} // rew   namespace
} // genie namespace

#endif
