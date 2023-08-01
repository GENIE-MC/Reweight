#include "ProfSpline/ObservableF.h"

namespace genie {
namespace rew {
ObservableF &ObservableF::Instance() {
  static ObservableF instance;
  return instance;
}

bool ObservableF::register_creator(std::string name,
                                   std::function<ObservableI *()> creator) {
  creators[name] = creator;
  return true;
}

ObservableI *ObservableF::create(std::string name) const {
  auto it = creators.find(name);
  if (it == creators.end()) {
    return nullptr;
  }
  return it->second();
}
} // namespace rew
} // namespace genie