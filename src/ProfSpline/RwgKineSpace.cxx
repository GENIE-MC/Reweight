//____________________________________________________________________________
/*
 Copyright (c) 2003-2026, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Qiyu Yan <yanqiyu17 \at mails.ucas.ac.cn>
          University of Chinese Academy of Sciences (UCAS) and University of Warwick

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include "RwgKineSpace.h"

namespace genie {
namespace rew {
void RwgKineSpace::Configure(const Registry &config) {
  Algorithm::Configure(config);
  LoadConfig();
}
//____________________________________________________________________________
void RwgKineSpace::Configure(string param_set) {
  Algorithm::Configure(param_set);
  LoadConfig();
}
} // namespace rew

} // namespace genie