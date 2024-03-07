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