//____________________________________________________________________________
/*!

\class    genie::GSystUncertaintyTable

\brief    Algorithm that enables bi-directional 1-sigma parameter errors
          to be configured for Reweight via XML

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  September 10, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_SYST_UNCERTAINTY_TABLE_H_
#define _G_SYST_UNCERTAINTY_TABLE_H_

// Generator includes
#include "Framework/Algorithm/Algorithm.h"

// Reweight includes
#include "RwFramework/GSyst.h"

namespace genie {
namespace rew {

class GSystUncertaintyTable : public Algorithm {

public:
  GSystUncertaintyTable();
  GSystUncertaintyTable(string config);
  virtual ~GSystUncertaintyTable();

  // Override the Algorithm::Configure methods to load configuration
  // data in private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

  // Represents a pair of plus and minus one-sigma errors
  // for a model parameter accessible via a Reweight tweak dial
  class MapEntry {
    public:
      MapEntry() : fErrPlusOneSigma(0.), fErrMinusOneSigma(0.) {}
      MapEntry(double plus_one_sig_err, double minus_one_sig_err)
        : fErrPlusOneSigma( plus_one_sig_err ),
        fErrMinusOneSigma( minus_one_sig_err ) {}

      inline double PlusOneSigmaErr() const { return fErrPlusOneSigma; }
      inline double MinusOneSigmaErr() const { return fErrMinusOneSigma; }

      inline void SetPlusOneSigmaErr(double plus_err)
        { fErrPlusOneSigma = plus_err; }

      inline void SetMinusOneSigmaErr(double minus_err)
        { fErrMinusOneSigma = minus_err; }

    protected:
      double fErrPlusOneSigma;
      double fErrMinusOneSigma;
  };

  inline const std::map<GSyst_t, MapEntry>& GetErrorsMap() const
    { return fErrorsMap; }

  inline std::map<GSyst_t, MapEntry>* GetErrorsMapPtr() { return &fErrorsMap; }

private:

  void LoadConfig (void);

  std::map<GSyst_t, MapEntry> fErrorsMap;

};

} // rew namespace
} // genie namespace

#endif //_G_SYST_UNCERTAINTY_TABLE_H_
