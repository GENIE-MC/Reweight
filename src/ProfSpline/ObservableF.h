#ifndef _OBSERVABLE_F_
#define _OBSERVABLE_F_
#include "ProfSpline/ObservableI.h"
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>
namespace genie {
namespace rew {
class ObservableF {
public:
  static ObservableF &Instance();

  bool register_creator(std::string name,
                        std::function<ObservableI *()> creator);

  ObservableI *create(std::string name) const;

  ObservableF(const ObservableF &) = delete;
  ObservableF(ObservableF &&) = delete;

private:
  ObservableF() = default;
  std::unordered_map<std::string, std::function<ObservableI *()>> creators;
};

#define REGISTER_OBSERVABLE(name, type)                                        \
  namespace {                                                                  \
  bool registered_##type = ObservableF::Instance().register_creator(           \
      name, [] { return new type(); });                                        \
  }

} // namespace rew
} // namespace genie

#endif