#ifndef _G_REWEIGHT_PROFESSOR_
#define _G_REWEIGHT_PROFESSOR_

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_PROFESSOR2_ENABLED__
#include "ProfSpline/KinematicVariables.h"
#include "ProfSpline/ObservableSplines.h"
#endif
#include "RwCalculators/GReWeightModel.h"
#include <string>

namespace genie {
namespace rew {
class GReWeightProfessor : public GReWeightModel {
public:
  GReWeightProfessor(std::string name);
  ~GReWeightProfessor(){};

  // inherited from GReWeightModel
  virtual bool AppliesTo(const EventRecord &event) const override;
  virtual bool IsHandled(GSyst_t syst) const override;
  virtual void SetSystematic(GSyst_t syst, double val) override;
  virtual void Reset(void) override;
  virtual void Reconfigure(void) override;
  virtual double CalcWeight(const genie::EventRecord &event) override;


#ifdef __GENIE_PROFESSOR2_ENABLED__
  std::map<std::tuple<std::string /*configuration full id*/,
                      ChannelIDs /*Channel selection ID*/>,
           std::vector<std::string> /*the vars lines from prof2*/>
  ReadProf2Spline(std::string filepath);
#endif

  void SetSystematic(const std::vector<double> &m_systematics_values,
                     const std::vector<double> &m_orig_value) {
    systematics_values = m_systematics_values;
    orig_value = m_orig_value;
  }

  void ReadComparisonXML(std::string filepath, std::string spline_path);

 private:
#ifdef __GENIE_PROFESSOR2_ENABLED__
  std::vector<ObservableSplines> observables{};
#endif
  std::vector<double> systematics_values, orig_value;
  std::vector<std::string> spline_vars;
  std::vector<double> var_min, var_max;
};
} // namespace rew
} // namespace genie

#endif
