#ifndef _G_REWEIGHT_PROFESSOR_
#define _G_REWEIGHT_PROFESSOR_

#include "ProfSpline/ObservableSplines.h"
#include "RwCalculators/GReWeightModel.h"
#include <cstddef>
#include <memory>
#include <string>

namespace genie {
namespace rew {
class GReWeightProfessor : public GReWeightModel {
public:
  GReWeightProfessor(std::string name);
  ~GReWeightProfessor(){};

  // inherited from GReWeightModel
  virtual bool AppliesTo (const EventRecord & event) const override;
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

  void ReadProf2Spline(std::string filepath);

  // void InitializeObservable(std::string name) {
  //   if (!observable_splines)
  //     observable_splines = std::make_unique<ObservableSplines>();
  //   observable_splines->InitializeObservable(name);
  // }
  // void InitializeBins(std::vector<std::vector<double>> binning) {
  //   if (!observable_splines)
  //     observable_splines = std::make_unique<ObservableSplines>();
  //   observable_splines->InitializeBins(binning);
  // }

  void SetSystematic(const std::vector<double> &m_systematics_values,
                     const std::vector<double> &m_orig_value) {
    systematics_values = m_systematics_values;
    orig_value = m_orig_value;
  }

  void ReadComparionXML(std::string filepath);

private:
  // the observable splines
  // std::unique_ptr<ObservableSplines> observable_splines;
  // std::vector<std::shared_ptr<ObservableSplines>> observable_splines;
  // std::map<std::string, std::shared_ptr<ObservableSplines>> observable_map_from_name;
  std::map<std::tuple<std::string, int, int>, std::shared_ptr<ObservableSplines>> observable_map_from_id;
  std::vector<double> systematics_values, orig_value;
  std::vector<std::string> spline_vars;
  std::vector<double> var_min, var_max;
};
} // namespace rew
} // namespace genie

#endif