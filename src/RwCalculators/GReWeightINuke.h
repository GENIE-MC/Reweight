//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightINuke

\brief    Reweighting GENIE INTRANUKE/hA hadron transport model.

	  The reweighting code considers two sets of physics changes:
          * Change in the hadron mean free path \lambda, i.e. change in
            the total rescattering probability P_{rescat}.
          * Changes in probabilty for rescattering mode X, given a fixed
            total rescattering probability P(X | \lambda).
            X = {elastic, inelastic, charge exchange, pion production,
                 absorption + multi-nucleon knockout}.

          Physics changes are considered separately for pions and nucleons.
          Unitarity is explicitly conserved.

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Sep 10, 2009

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_INUKE_H_
#define _G_REWEIGHT_INUKE_H_

//#define _G_REWEIGHT_INUKE_DEBUG_NTP_

// Standard library includes
#include <memory>

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"
#include "RwCalculators/GReWeightINukeParams.h"

using namespace genie::rew;
using namespace genie;

class TFile;
class TNtuple;
class TLorentzVector;

namespace genie {

 class HAIntranuke2018;
 class GHepParticle;

namespace rew   {

 class GReWeightINuke : public GReWeightModel
 {
 public:
   GReWeightINuke();
   virtual ~GReWeightINuke();

   // implement the GReWeightI interface
   bool   AppliesTo      (const EventRecord & event) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

 protected:

   void CalcDeltaAZ( const EventRecord& event, const GHepParticle& p,
     int& deltaA, int& deltaZ );

   void UpdateRemnantAZ( int deltaA, int deltaZ );

   std::shared_ptr< GReWeightINukeParams > fINukeRwParams;

   HAIntranuke2018* fFSIModel;

#ifdef _G_REWEIGHT_INUKE_DEBUG_NTP_
   TFile *              fTestFile;
   TNtuple *            fTestNtp;
#endif

 };

} // rew
} // genie

#endif
