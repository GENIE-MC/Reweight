#ifndef _G_REWEIGHT_PROFESSOR_
#define _G_REWEIGHT_PROFESSOR_

#include "ProfSpline/ObservableSplines.h"
#include "RwCalculators/GReWeightModel.h"
#include <memory>
#include <string>

namespace genie {
namespace rew {
class GReWeightProfessor : public GReWeightModel {
public:
  GReWeightProfessor(std::string name);
  ~GReWeightProfessor();

  // inherited from GReWeightModel
  virtual bool AppliesTo(ScatteringType_t type, bool is_cc) const override;
  virtual bool IsHandled(GSyst_t syst) const override;
  virtual void SetSystematic(GSyst_t syst, double val) override;
  virtual void Reset(void) override;
  virtual void Reconfigure(void) override;
  virtual double CalcWeight(const genie::EventRecord &event) override;

  // Some GReWeightProfessor specific functions
  // TODO: we need to do
  //  - Get binning information and pass it to ObservableSplines
  //  - Get Observable information and pass it to ObservableSplines
  //  - Get a list of nuisance parameters and maintain a map and the vector

  void Initialize(std::string);


private:
  void ReadProf2Spline(std::string filepath);
  void ReadComparionConf(std::string filepath);

  // the observable splines
  std::unique_ptr<ObservableSplines> observable_splines;
  std::vector<double> systematics_values;
  std::vector<std::string> spline_vars;
  std::vector<double> var_min, var_max;
  size_t dimension;
};
} // namespace rew
} // namespace genie

#endif