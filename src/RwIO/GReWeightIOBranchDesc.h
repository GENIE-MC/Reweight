//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightIOBranchDesc

\brief

\author   Julia Yarba <yarba_j@fnal.gov>
          Fermilab

\created  Aug 1, 2016

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________


#ifndef GRWIOBRANCHDESC_H
#define GRWIOBRANCHDESC_H

#include <string>
#include <TObject.h>

class TRootIOCtor;

namespace genie
{
namespace rew
{

class GReWeightIOBranchDesc : public TObject
{

   public:
   
      GReWeightIOBranchDesc() : TObject(), 
                                fParameterName(), fParameterMean(0.), 
                                fParameterSigmaPlus(0.), fParameterSigmaMinus(0.) {} 
      GReWeightIOBranchDesc( const std::string&, const double, const double, const double );
      GReWeightIOBranchDesc( const GReWeightIOBranchDesc& );
      GReWeightIOBranchDesc( TRootIOCtor* );
      
      ~GReWeightIOBranchDesc() {}
      
      const std::string& GetParameterName()       const { return fParameterName; }
      double             GetParameterMean()       const { return fParameterMean; }
      double             GetParameterSigmaPlus()  const { return fParameterSigmaPlus; }
      double             GetParameterSigmaMinus() const { return fParameterSigmaMinus; }
      
      void               SetParameter( const std::string&, const double, const double, const double );
         
   private:
   
      std::string fParameterName;
      double      fParameterMean;
      double      fParameterSigmaPlus;
      double      fParameterSigmaMinus;


ClassDef(GReWeightIOBranchDesc,2)

};

} // end namespace rew
} // end namespace genie

#endif
