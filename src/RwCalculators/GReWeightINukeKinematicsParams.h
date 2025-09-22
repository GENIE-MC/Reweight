#ifndef _G_REWEIGHT_INTRANUKEKINEMATICS_PARAMS_H_
#define _G_REWEIGHT_INTRANUKEKINEMATICS_PARAMS_H_

#include <map>
#include <TLorentzVector.h>

// GENIE/Reweight includes
#include "RwFramework/GSyst.h"
#include "Framework/Numerical/BLI2D.h"

namespace genie {
namespace rew   {

 class GReWeightINukeKinematicsParams {

 public:
   GReWeightINukeKinematicsParams();
  ~GReWeightINukeKinematicsParams() {}

   void    Reset              (void);                   ///<
   void    SetTwkDial         (GSyst_t s, double val);  ///<
   void    SetTargetA         (int target_A); ///< Set the mass number of the hit nucleus
   double  CalcWeight         (int np_pp, float tlab, float costhcm); 
   double  CalcPionWeight     (TLorentzVector p, TLorentzVector f1, TLorentzVector f2, TLorentzVector f3);


   struct ReBounce {
     void SetUniverse(int u);
     void LoadUniverse(int u);
     void Reset();
     double CalcWeight(float tlab, float costhcm);

     ReBounce(std::string n);
     ~ReBounce();

     int fUniverse;
     std::string fName;
     std::shared_ptr<BLI2DNonUnifGrid> fRewt;
     std::map<int, std::shared_ptr<BLI2DNonUnifGrid>> fCache;
   };

 private:
    ReBounce fNPwgt;
    ReBounce fPPwgt;
    bool fFixPiPro;
    double fBiasPiPro;

 }; //GReWeightINukeKinematicsParams

}      // rew   namespace
}      // genie namespace

#endif // _G_REWEIGHT_INTRANUKEKINEMATICS_PARAMS_H_
