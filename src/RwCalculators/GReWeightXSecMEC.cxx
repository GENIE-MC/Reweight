//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// GENIE/Generator includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/InteractionType.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Physics/Multinucleon/XSection/SuSAv2MECPXSec.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightXSecMEC.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

#include <iomanip>

using namespace genie;
using namespace genie::rew;

namespace {

  // Helper function to help us initialize the static std::map
  // owned by the GReWeightXSecMEC class
  std::map<GSyst_t, InteractionType_t> make_gsyst_to_inttype_map() {
    std::map<GSyst_t, InteractionType_t> temp_map;

    temp_map[ kXSecTwkDial_NormCCMEC ] = kIntWeakCC;
    temp_map[ kXSecTwkDial_NormNCMEC ] = kIntWeakNC;
    temp_map[ kXSecTwkDial_NormEMMEC ] = kIntEM;

    return temp_map;
  }

// Neutrino Energies from 0.0 GeV to 3.0 GeV for
// Energy_Dependence_CCMEC parameter. 
// Using: GENIE v3.4.0, G18_10a (2p-2h Valencia) and modified 
// tune G18_10s (2p-2h SuSAv2)
static double const nu_energies[552] = {
0.000, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 
0.28, 0.285, 0.29, 0.295, 0.3, 0.305, 0.31, 0.315, 
0.32, 0.325, 0.33, 0.335, 0.34, 0.345, 0.35, 0.355, 
0.36, 0.365, 0.37, 0.375, 0.38, 0.385, 0.39, 0.395, 
0.4, 0.405, 0.41, 0.415, 0.42, 0.425, 0.43, 0.435, 
0.44, 0.445, 0.45, 0.455, 0.46, 0.465, 0.47, 0.475, 
0.48, 0.485, 0.49, 0.495, 0.5, 0.505, 0.51, 0.515, 
0.52, 0.525, 0.53, 0.535, 0.54, 0.545, 0.55, 0.555, 
0.56, 0.565, 0.57, 0.575, 0.58, 0.585, 0.59, 0.595, 
0.6, 0.605, 0.61, 0.615, 0.62, 0.625, 0.63, 0.635, 
0.64, 0.645, 0.65, 0.655, 0.66, 0.665, 0.67, 0.675, 
0.68, 0.685, 0.69, 0.695, 0.7, 0.705, 0.71, 0.715, 
0.72, 0.725, 0.73, 0.735, 0.74, 0.745, 0.75, 0.755, 
0.76, 0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795, 
0.8, 0.805, 0.81, 0.815, 0.82, 0.825, 0.83, 0.835, 
0.84, 0.845, 0.85, 0.855, 0.86, 0.865, 0.87, 0.875, 
0.88, 0.885, 0.89, 0.895, 0.9, 0.905, 0.91, 0.915, 
0.92, 0.925, 0.93, 0.935, 0.94, 0.945, 0.95, 0.955, 
0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 
1, 1.005, 1.01, 1.015, 1.02, 1.025, 1.03, 1.035, 
1.04, 1.045, 1.05, 1.055, 1.06, 1.065, 1.07, 1.075, 
1.08, 1.085, 1.09, 1.095, 1.1, 1.105, 1.11, 1.115, 
1.12, 1.125, 1.13, 1.135, 1.14, 1.145, 1.15, 1.155, 
1.16, 1.165, 1.17, 1.175, 1.18, 1.185, 1.19, 1.195, 
1.2, 1.205, 1.21, 1.215, 1.22, 1.225, 1.23, 1.235, 
1.24, 1.245, 1.25, 1.255, 1.26, 1.265, 1.27, 1.275, 
1.28, 1.285, 1.29, 1.295, 1.3, 1.305, 1.31, 1.315, 
1.32, 1.325, 1.33, 1.335, 1.34, 1.345, 1.35, 1.355, 
1.36, 1.365, 1.37, 1.375, 1.38, 1.385, 1.39, 1.395, 
1.4, 1.405, 1.41, 1.415, 1.42, 1.425, 1.43, 1.435, 
1.44, 1.445, 1.45, 1.455, 1.46, 1.465, 1.47, 1.475, 
1.48, 1.485, 1.49, 1.495, 1.5, 1.505, 1.51, 1.515, 
1.52, 1.525, 1.53, 1.535, 1.54, 1.545, 1.55, 1.555, 
1.56, 1.565, 1.57, 1.575, 1.58, 1.585, 1.59, 1.595, 
1.6, 1.605, 1.61, 1.615, 1.62, 1.625, 1.63, 1.635, 
1.64, 1.645, 1.65, 1.655, 1.66, 1.665, 1.67, 1.675, 
1.68, 1.685, 1.69, 1.695, 1.7, 1.705, 1.71, 1.715, 
1.72, 1.725, 1.73, 1.735, 1.74, 1.745, 1.75, 1.755, 
1.76, 1.765, 1.77, 1.775, 1.78, 1.785, 1.79, 1.795, 
1.8, 1.805, 1.81, 1.815, 1.82, 1.825, 1.83, 1.835, 
1.84, 1.845, 1.85, 1.855, 1.86, 1.865, 1.87, 1.875, 
1.88, 1.885, 1.89, 1.895, 1.9, 1.905, 1.91, 1.915, 
1.92, 1.925, 1.93, 1.935, 1.94, 1.945, 1.95, 1.955, 
1.96, 1.965, 1.97, 1.975, 1.98, 1.985, 1.99, 1.995, 
2, 2.005, 2.01, 2.015, 2.02, 2.025, 2.03, 2.035, 
2.04, 2.045, 2.05, 2.055, 2.06, 2.065, 2.07, 2.075, 
2.08, 2.085, 2.09, 2.095, 2.1, 2.105, 2.11, 2.115, 
2.12, 2.125, 2.13, 2.135, 2.14, 2.145, 2.15, 2.155, 
2.16, 2.165, 2.17, 2.175, 2.18, 2.185, 2.19, 2.195, 
2.2, 2.205, 2.21, 2.215, 2.22, 2.225, 2.23, 2.235, 
2.24, 2.245, 2.25, 2.255, 2.26, 2.265, 2.27, 2.275, 
2.28, 2.285, 2.29, 2.295, 2.3, 2.305, 2.31, 2.315, 
2.32, 2.325, 2.33, 2.335, 2.34, 2.345, 2.35, 2.355, 
2.36, 2.365, 2.37, 2.375, 2.38, 2.385, 2.39, 2.395, 
2.4, 2.405, 2.41, 2.415, 2.42, 2.425, 2.43, 2.435, 
2.44, 2.445, 2.45, 2.455, 2.46, 2.465, 2.47, 2.475, 
2.48, 2.485, 2.49, 2.495, 2.5, 2.505, 2.51, 2.515, 
2.52, 2.525, 2.53, 2.535, 2.54, 2.545, 2.55, 2.555, 
2.56, 2.565, 2.57, 2.575, 2.58, 2.585, 2.59, 2.595, 
2.6, 2.605, 2.61, 2.615, 2.62, 2.625, 2.63, 2.635, 
2.64, 2.645, 2.65, 2.655, 2.66, 2.665, 2.67, 2.675, 
2.68, 2.685, 2.69, 2.695, 2.7, 2.705, 2.71, 2.715, 
2.72, 2.725, 2.73, 2.735, 2.74, 2.745, 2.75, 2.755, 
2.76, 2.765, 2.77, 2.775, 2.78, 2.785, 2.79, 2.795, 
2.8, 2.805, 2.81, 2.815, 2.82, 2.825, 2.83, 2.835, 
2.84, 2.845, 2.85, 2.855, 2.86, 2.865, 2.87, 2.875, 
2.88, 2.885, 2.89, 2.895, 2.9, 2.905, 2.91, 2.915, 
2.92, 2.925, 2.93, 2.935, 2.94, 2.945, 2.95, 2.955, 
2.96, 2.965, 2.97, 2.975, 2.98, 2.985, 2.99, 2.995
};

// Cross-section ratios of the 2p-2h SuSAv2 to Valencia cross sections
// after normalising both distributions to be unity at 1.2 GeV.
// Using: GENIE v3.4.0, G18_10a (2p-2h Valencia) and modified 
// tune G18_10s (2p-2h SuSAv2)
/*
static double const xsec_ratios[552] = {
0.00000, 11.1108, 9.85901, 8.83019, 7.97319, 7.24971, 6.63189, 6.09899, 
5.63523, 5.22779, 4.86655, 4.54389, 4.25311, 3.98952, 3.75033, 3.56951, 
3.40238, 3.24497, 3.09793, 2.95988, 2.8295, 2.70762, 2.59203, 2.48365, 
2.38133, 2.28551, 2.1955, 2.11153, 2.03262, 1.95923, 1.89012, 1.82568, 
1.76537, 1.70866, 1.65564, 1.61708, 1.58105, 1.54681, 1.51424, 1.48303, 
1.45329, 1.42497, 1.39798, 1.37219, 1.34758, 1.32416, 1.3019, 1.28075, 
1.26067, 1.24159, 1.22347, 1.2063, 1.19004, 1.17464, 1.16008, 1.14888, 
1.13854, 1.12879, 1.1196, 1.11097, 1.10287, 1.09531, 1.08826, 1.08171, 
1.07566, 1.07009, 1.065, 1.06037, 1.05617, 1.05238, 1.04896, 1.04589, 
1.04315, 1.04072, 1.03859, 1.03581, 1.0332, 1.03085, 1.02876, 1.02691, 
1.02528, 1.02387, 1.02269, 1.02169, 1.02088, 1.02023, 1.01974, 1.01944, 
1.01929, 1.01928, 1.0194, 1.01964, 1.02004, 1.02056, 1.02119, 1.02085, 
1.0205, 1.02029, 1.02019, 1.02019, 1.02033, 1.02058, 1.02093, 1.02137, 
1.02194, 1.0226, 1.02333, 1.02416, 1.0251, 1.0261, 1.02717, 1.02836, 
1.0296, 1.0309, 1.03232, 1.03235, 1.0322, 1.03218, 1.03223, 1.03235, 
1.03258, 1.03287, 1.03324, 1.03369, 1.0342, 1.0348, 1.03546, 1.03617, 
1.03698, 1.03783, 1.03874, 1.03973, 1.04075, 1.04186, 1.04301, 1.04247, 
1.0417, 1.04096, 1.04033, 1.03974, 1.03922, 1.03877, 1.03836, 1.03805, 
1.03776, 1.03756, 1.0374, 1.03729, 1.03725, 1.03724, 1.0373, 1.03739, 
1.03756, 1.03774, 1.038, 1.03695, 1.03569, 1.03446, 1.0333, 1.03218, 
1.03111, 1.03009, 1.02912, 1.02819, 1.0273, 1.02647, 1.02567, 1.02492, 
1.0242, 1.02353, 1.02288, 1.02229, 1.02172, 1.0212, 1.0207, 1.01946, 
1.01806, 1.0167, 1.01537, 1.01408, 1.01282, 1.01159, 1.0104, 1.00923, 
1.0081, 1.00699, 1.00592, 1.00486, 1.00384, 1.00284, 1.00187, 1.00092, 
1, 1.0009, 1.00178, 1.00304, 1.00442, 1.00575, 1.00708, 1.00837, 
1.00965, 1.01091, 1.01214, 1.01337, 1.01456, 1.01575, 1.01689, 1.01803, 
1.01915, 1.02025, 1.02133, 1.02239, 1.02344, 1.02446, 1.02547, 1.02669, 
1.02797, 1.02925, 1.03048, 1.03172, 1.03293, 1.03413, 1.03532, 1.03647, 
1.03763, 1.03876, 1.03988, 1.04099, 1.04207, 1.04314, 1.0442, 1.04524, 
1.04628, 1.04728, 1.04828, 1.04935, 1.05041, 1.05148, 1.05251, 1.05355, 
1.05457, 1.05557, 1.05656, 1.05754, 1.05851, 1.05947, 1.06041, 1.06134, 
1.06227, 1.06316, 1.06406, 1.06495, 1.06582, 1.06669, 1.06754, 1.06834, 
1.06912, 1.06988, 1.07063, 1.07138, 1.0721, 1.07283, 1.07355, 1.07424, 
1.07493, 1.07562, 1.07628, 1.07694, 1.0776, 1.07823, 1.07887, 1.07949, 
1.0801, 1.0807, 1.0813, 1.08184, 1.08236, 1.08288, 1.08338, 1.08387, 
1.08436, 1.08483, 1.0853, 1.08576, 1.08621, 1.08665, 1.08709, 1.0875, 
1.08792, 1.08833, 1.08873, 1.08912, 1.08951, 1.08988, 1.09024, 1.09059, 
1.09091, 1.09123, 1.09155, 1.09185, 1.09215, 1.09244, 1.09273, 1.093, 
1.09327, 1.09353, 1.09378, 1.09403, 1.09428, 1.0945, 1.09473, 1.09496, 
1.09517, 1.09537, 1.09558, 1.09577, 1.09597, 1.09616, 1.09634, 1.09651, 
1.09669, 1.09686, 1.09701, 1.09716, 1.09731, 1.09746, 1.09759, 1.09772, 
1.09785, 1.09797, 1.09808, 1.09819, 1.0983, 1.09839, 1.09849, 1.0986, 
1.09871, 1.09882, 1.09893, 1.09903, 1.09912, 1.09922, 1.09931, 1.09939, 
1.09946, 1.09954, 1.09961, 1.09967, 1.09973, 1.09979, 1.09984, 1.09989, 
1.09993, 1.09998, 1.10001, 1.10006, 1.10012, 1.10017, 1.10022, 1.10026, 
1.10031, 1.10035, 1.10038, 1.10041, 1.10044, 1.10047, 1.10049, 1.10051, 
1.10053, 1.10054, 1.10054, 1.10055, 1.10056, 1.10056, 1.10055, 1.10056, 
1.10058, 1.10059, 1.10059, 1.1006, 1.10061, 1.1006, 1.1006, 1.1006, 
1.1006, 1.10058, 1.10057, 1.10056, 1.10055, 1.10053, 1.10051, 1.10049, 
1.10047, 1.10044, 1.10042, 1.1004, 1.10038, 1.10036, 1.10034, 1.10032, 
1.1003, 1.10028, 1.10026, 1.10023, 1.1002, 1.10018, 1.10015, 1.10012, 
1.10009, 1.10005, 1.10002, 1.09999, 1.09995, 1.09992, 1.09988, 1.09985, 
1.09982, 1.09979, 1.09976, 1.09972, 1.09969, 1.09966, 1.09963, 1.09959, 
1.09956, 1.09953, 1.09949, 1.09946, 1.09943, 1.0994, 1.09936, 1.09933, 
1.0993, 1.09927, 1.09923, 1.0992, 1.09917, 1.09914, 1.09911, 1.09909, 
1.09906, 1.09903, 1.099, 1.09898, 1.09895, 1.09893, 1.0989, 1.09888, 
1.09885, 1.09883, 1.09881, 1.09879, 1.09876, 1.09874, 1.09872, 1.0987, 
1.09868, 1.09866, 1.09864, 1.09862, 1.0986, 1.09858, 1.09857, 1.09855, 
1.09853, 1.09852, 1.0985, 1.09848, 1.09847, 1.09845, 1.09844, 1.09842, 
1.09841, 1.09839, 1.09838, 1.09837, 1.09836, 1.09834, 1.09833, 1.09832, 
1.09831, 1.0983, 1.09829, 1.09828, 1.09827, 1.09826, 1.09825, 1.09824, 
1.09823, 1.09822, 1.09821, 1.0982, 1.09819, 1.09818, 1.09818, 1.09817, 
1.09816, 1.09816, 1.09815, 1.09814, 1.09814, 1.09813, 1.09813, 1.09812, 
1.09812, 1.09811, 1.09811, 1.0981, 1.0981, 1.09809, 1.09809, 1.09808, 
1.09808, 1.09807, 1.09807, 1.09806, 1.09806, 1.09806, 1.09806, 1.09806, 
1.09805, 1.09805, 1.09805, 1.09805, 1.09804, 1.09804, 1.09804, 1.09804, 
1.09803, 1.09803, 1.09803, 1.09803, 1.09803, 1.09802, 1.09802, 1.09802, 
1.09802, 1.09802, 1.09802, 1.09802, 1.09802, 1.09802, 1.09802, 1.09801, 
1.09801, 1.09801, 1.09801, 1.09801, 1.09801, 1.09801, 1.098, 1.098
};
*/

// Cross-section ratios of the 2p-2h SuSAv2 to Valencia cross sections
// after normalising both distributions to be unity at 10.019 GeV.
// Using: GENIE v3.4.0, G18_10a (2p-2h Valencia) and modified 
// tune G18_10s (2p-2h SuSAv2)
static double const xsec_ratios[552] = {
0.00000, 5.83079, 5.49906, 5.23727, 5.0254, 4.85042, 4.70346, 4.57829, 
4.4704, 4.37645, 4.29389, 4.22078, 4.15557, 4.09706, 4.04425, 3.64996, 
3.33725, 3.09439, 2.90032, 2.74168, 2.60958, 2.49787, 2.40218, 2.31928, 
2.24677, 2.18282, 2.12598, 2.07515, 2.02941, 1.98804, 1.95043, 1.9161, 
1.88464, 1.8557, 1.82899, 1.76818, 1.71209, 1.66243, 1.61816, 1.57845, 
1.54263, 1.51014, 1.48056, 1.4535, 1.42865, 1.40576, 1.3846, 1.36498, 
1.34674, 1.32975, 1.31386, 1.29899, 1.28504, 1.27192, 1.25956, 1.25008, 
1.24133, 1.23307, 1.22525, 1.21783, 1.2108, 1.20411, 1.19774, 1.19168, 
1.1859, 1.18037, 1.17509, 1.17004, 1.1652, 1.16056, 1.15611, 1.15184, 
1.14773, 1.14378, 1.13998, 1.13836, 1.13703, 1.13574, 1.13448, 1.13327, 
1.13208, 1.13094, 1.12982, 1.12873, 1.12767, 1.12664, 1.12564, 1.12466, 
1.12371, 1.12278, 1.12187, 1.12099, 1.12013, 1.11928, 1.11846, 1.11865, 
1.11899, 1.11931, 1.11963, 1.11994, 1.12025, 1.12055, 1.12085, 1.12114, 
1.12143, 1.12171, 1.12199, 1.12226, 1.12253, 1.1228, 1.12306, 1.12331, 
1.12357, 1.12381, 1.12406, 1.12484, 1.1257, 1.12654, 1.12738, 1.12821, 
1.12903, 1.12984, 1.13064, 1.13143, 1.13222, 1.13299, 1.13376, 1.13451, 
1.13526, 1.136, 1.13673, 1.13745, 1.13817, 1.13888, 1.13958, 1.13987, 
1.14008, 1.14028, 1.14049, 1.14069, 1.14089, 1.14109, 1.14129, 1.14149, 
1.14168, 1.14188, 1.14207, 1.14226, 1.14245, 1.14264, 1.14283, 1.14301, 
1.14319, 1.14338, 1.14356, 1.14259, 1.14137, 1.14016, 1.13895, 1.13776, 
1.13657, 1.13538, 1.13421, 1.13304, 1.13188, 1.13072, 1.12957, 1.12843, 
1.12729, 1.12616, 1.12503, 1.12392, 1.1228, 1.1217, 1.1206, 1.11912, 
1.11754, 1.11598, 1.11442, 1.11287, 1.11133, 1.10979, 1.10826, 1.10673, 
1.10521, 1.1037, 1.10219, 1.10069, 1.0992, 1.09771, 1.09623, 1.09475, 
1.09328, 1.09182, 1.09036, 1.08895, 1.08755, 1.08616, 1.08477, 1.08339, 
1.08201, 1.08064, 1.07927, 1.0779, 1.07654, 1.07518, 1.07383, 1.07248, 
1.07114, 1.0698, 1.06846, 1.06713, 1.0658, 1.06448, 1.06316, 1.062, 
1.0609, 1.05981, 1.05871, 1.05762, 1.05653, 1.05544, 1.05436, 1.05328, 
1.0522, 1.05112, 1.05005, 1.04898, 1.04791, 1.04684, 1.04578, 1.04472, 
1.04366, 1.04261, 1.04155, 1.04066, 1.03981, 1.03896, 1.03812, 1.03728, 
1.03644, 1.0356, 1.03477, 1.03393, 1.0331, 1.03227, 1.03144, 1.03061, 
1.02978, 1.02896, 1.02813, 1.02731, 1.02649, 1.02567, 1.02485, 1.02416, 
1.0235, 1.02285, 1.0222, 1.02155, 1.02091, 1.02026, 1.01961, 1.01897, 
1.01832, 1.01768, 1.01704, 1.0164, 1.01576, 1.01512, 1.01448, 1.01384, 
1.01321, 1.01257, 1.01194, 1.01143, 1.01098, 1.01052, 1.01007, 1.00961, 
1.00916, 1.00871, 1.00825, 1.0078, 1.00735, 1.0069, 1.00645, 1.006, 
1.00555, 1.0051, 1.00466, 1.00421, 1.00376, 1.00332, 1.00287, 1.00257, 
1.00234, 1.0021, 1.00187, 1.00163, 1.0014, 1.00116, 1.00093, 1.0007, 
1.00046, 1.00023, 1, 1.00024, 1.00047, 1.0007, 1.00094, 1.00117, 
1.0014, 1.00163, 1.00187, 1.00202, 1.00214, 1.00225, 1.00237, 1.00249, 
1.0026, 1.00272, 1.00283, 1.00295, 1.00307, 1.00318, 1.0033, 1.00342, 
1.00353, 1.00365, 1.00376, 1.00388, 1.00399, 1.00411, 1.00423, 1.0043, 
1.00434, 1.00439, 1.00443, 1.00448, 1.00453, 1.00457, 1.00462, 1.00466, 
1.00471, 1.00476, 1.0048, 1.00485, 1.00489, 1.00494, 1.00499, 1.00503, 
1.00508, 1.00512, 1.00517, 1.00519, 1.0052, 1.00521, 1.00521, 1.00522, 
1.00523, 1.00524, 1.00524, 1.00525, 1.00526, 1.00527, 1.00527, 1.00528, 
1.00529, 1.0053, 1.0053, 1.00531, 1.00532, 1.00533, 1.00533, 1.00533, 
1.00533, 1.00532, 1.00531, 1.00531, 1.0053, 1.0053, 1.00529, 1.00529, 
1.00528, 1.00527, 1.00527, 1.00526, 1.00526, 1.00525, 1.00525, 1.00524, 
1.00523, 1.00523, 1.00522, 1.00522, 1.00521, 1.0052, 1.00519, 1.00519, 
1.00518, 1.00517, 1.00516, 1.00516, 1.00515, 1.00514, 1.00513, 1.00513, 
1.00512, 1.00511, 1.00511, 1.0051, 1.00509, 1.00508, 1.00508, 1.00507, 
1.00506, 1.00505, 1.00504, 1.00503, 1.00502, 1.00501, 1.005, 1.00499, 
1.00499, 1.00498, 1.00497, 1.00496, 1.00495, 1.00494, 1.00493, 1.00492, 
1.00491, 1.0049, 1.00489, 1.00488, 1.00487, 1.00486, 1.00485, 1.00484, 
1.00483, 1.00482, 1.0048, 1.00479, 1.00478, 1.00477, 1.00476, 1.00475, 
1.00474, 1.00472, 1.00471, 1.0047, 1.00469, 1.00468, 1.00467, 1.00466, 
1.00465, 1.00464, 1.00463, 1.00462, 1.00461, 1.0046, 1.00459, 1.00458, 
1.00457, 1.00456, 1.00455, 1.00454, 1.00453, 1.00452, 1.00451, 1.0045, 
1.00449, 1.00448, 1.00447, 1.00447, 1.00446, 1.00445, 1.00444, 1.00444, 
1.00443, 1.00442, 1.00442, 1.00441, 1.0044, 1.0044, 1.00439, 1.00438, 
1.00437, 1.00437, 1.00436, 1.00435, 1.00435, 1.00434, 1.00433, 1.00432, 
1.00432, 1.00431, 1.00431, 1.0043, 1.00429, 1.00429, 1.00428, 1.00427, 
1.00427, 1.00426, 1.00425, 1.00425, 1.00424, 1.00423, 1.00423, 1.00422, 
1.00422, 1.00421, 1.0042, 1.0042, 1.00419, 1.00418, 1.00417, 1.00416, 
1.00416, 1.00415, 1.00414, 1.00413, 1.00412, 1.00412, 1.00411, 1.0041, 
1.00409, 1.00409, 1.00408, 1.00407, 1.00406, 1.00405, 1.00405, 1.00404, 
1.00403, 1.00402, 1.00401, 1.004, 1.00399, 1.00398, 1.00397, 1.00396, 
1.00395, 1.00394, 1.00393, 1.00393, 1.00392, 1.00391, 1.0039, 1.00389
};

std::unique_ptr<TGraph> ratioGraph_nu = std::unique_ptr<TGraph>(new TGraph(552, nu_energies, xsec_ratios));

  //// MECGenerator::SelectEmpiricalKinematics() uses bogus hard-coded
  //// limits which are copied below for consistency.
  //// TODO: Do something better in MECGenerator, then change this
  //// code for consistency
  //Range1D_t getQ2LimitsEmpiricalMEC() {
  //  const double Q2min =  0.01;
  //  const double Q2max =  8.00;
  //  Range1D_t rQ2( Q2min, Q2max );
  //  return rQ2;
  //}

  //Range1D_t getWLimitsEmpiricalMEC() {
  //  const double Wmin  =  1.88;
  //  const double Wmax  =  3.00;
  //  Range1D_t rW( Wmin, Wmax );
  //  return rW;
  //}

}

// Define the static std::map owned by the GReWeightXSecMEC class
// TODO: switch to something better for GENIE 4. C++11 makes initializing
// static std::map objects a lot easier. - S. Gardiner
std::map<GSyst_t, InteractionType_t> GReWeightXSecMEC::fGSystToIntTypeMap
  = make_gsyst_to_inttype_map();

//_______________________________________________________________________________________
GReWeightXSecMEC::GReWeightXSecMEC()
  : GReWeightModel("MEC")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightXSecMEC::GReWeightXSecMEC(std::string /*model*/, std::string /*type*/)
  : GReWeightModel("MEC")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightXSecMEC::~GReWeightXSecMEC()
{
  // Delete the adopted algorithms if needed
  if ( fXSecAlgCCDef ) delete fXSecAlgCCDef;
  if ( fXSecAlgCCAlt_Nieves ) delete fXSecAlgCCAlt_Nieves;
  if ( fXSecAlgCCAlt_SuSAv2 ) delete fXSecAlgCCAlt_SuSAv2;
  if ( fXSecAlgCCAlt_Empirical ) delete fXSecAlgCCAlt_Empirical;
  if ( fXSecAlgCCAlt_Martini ) delete fXSecAlgCCAlt_Martini;
}
//_______________________________________________________________________________________
bool GReWeightXSecMEC::IsHandled(GSyst_t syst) const
{
  // Some MEC tweak dials are independent of interaction type. Check
  // whether we can handle these first.
  if ( syst == kXSecTwkDial_DecayAngMEC ) return true;
  if ( syst == kXSecTwkDial_DecayAng2MEC ) return true;
  if ( syst == kXSecTwkDial_FracPN_CCMEC ) return true;
  if ( syst == kXSecTwkDial_FracDelta_CCMEC ) return true;
  if ( syst == kXSecTwkDial_XSecShape_CCMEC ) return true;
  if ( syst == kXSecTwkDial_XSecShape_CCMEC_Empirical ) return true;
  if ( syst == kXSecTwkDial_XSecShape_CCMEC_Martini ) return true;
  if ( syst == kXSecTwkDial_EnergyDependence_CCMEC ) return true;

  // If we have an entry for a knob that is interaction type dependent in the
  // GSyst_t -> InteractionType_t map, then this calculator can handle it.
  // Otherwise, it can't.
  bool handle = fGSystToIntTypeMap.count( syst );
  return handle;
}
//_______________________________________________________________________________________
//bool GReWeightXSecMEC::AppliesTo(ScatteringType_t type, bool /*is_cc*/) const
bool GReWeightXSecMEC::AppliesTo(const EventRecord & event) const
{
  auto type = event.Summary()->ProcInfo().ScatteringTypeId();
  // Weights can be calculated for CC, NC, and EM MEC events
  if ( type == kScMEC ) return true;
  return false;
}
//_______________________________________________________________________________________
void GReWeightXSecMEC::SetSystematic(GSyst_t syst, double twk_dial)
{
  if ( !this->IsHandled(syst) ) return;

  // Handle the knobs that are independent of interaction type first
  if ( syst == kXSecTwkDial_DecayAngMEC ) {
    fDecayAngTwkDial = twk_dial;
//    std::cout << "fDecayAngTwkDial " << fDecayAngTwkDial << std::endl;
    return;
  }
  if ( syst == kXSecTwkDial_DecayAng2MEC ) {
    fDecayAng2TwkDial = twk_dial;
 //   std::cout << "fDecayAng2TwkDial " << fDecayAng2TwkDial << std::endl;
    return;
  }
  else if ( syst == kXSecTwkDial_FracPN_CCMEC ) {
    fFracPN_CCTwkDial = twk_dial;
    return;
  }
  else if ( syst == kXSecTwkDial_FracDelta_CCMEC ) {
    fFracDelta_CCTwkDial = twk_dial;
    return;
  }
  else if ( syst == kXSecTwkDial_XSecShape_CCMEC ) {
    fCCXSecShapeTwkDial = twk_dial;
    return;
  }
  else if ( syst == kXSecTwkDial_XSecShape_CCMEC_Empirical ) {
    fCCXSecShapeEmpiricalTwkDial = twk_dial;
    return;
  }
  else if ( syst == kXSecTwkDial_XSecShape_CCMEC_Martini ) {
    fCCXSecShapeMartiniTwkDial = twk_dial;
    return;
  }
  else if ( syst == kXSecTwkDial_EnergyDependence_CCMEC ) {
    fEnergyDependenceTwkDial = twk_dial;
    return;
  }

  // We've already checked that there is an entry for the given knob in the map
  // during the previous call to IsHandled(). Therefore, just retrieve the
  // stored value this time.
  InteractionType_t type = fGSystToIntTypeMap.at( syst );

  // Store the new tweak dial value in the entry for the interaction type of
  // interest
  fNormMap.at(type).fTwkDial = twk_dial;
}
//_______________________________________________________________________________________
void GReWeightXSecMEC::Reset(void)
{
  // Reset all of the normalization tweak dials to their defaults
  std::map<InteractionType_t, NormMapEntry>::iterator it = fNormMap.begin();
  std::map<InteractionType_t, NormMapEntry>::iterator end = fNormMap.end();
  while ( it != end ) {
    it->second.fTwkDial = 0.;
    it->second.fNormCurr = it->second.fNormDef;
    ++it;
  }

  fDecayAngTwkDial = 0.;
  fDecayAng2TwkDial = 0.;
  fCCXSecShapeTwkDial = 0.;
  fCCXSecShapeEmpiricalTwkDial = 0.;
  fCCXSecShapeMartiniTwkDial = 0.;

  fFracPN_CCTwkDial = 0.;
  fFracDelta_CCTwkDial = 0.;

  fEnergyDependenceTwkDial = 0.;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightXSecMEC::Reconfigure(void)
{
  GSystUncertainty* fracerr = GSystUncertainty::Instance();

  // Loop over all of the normalization tweak dials to update their current
  // values
  std::map<GSyst_t, InteractionType_t>::const_iterator it = fGSystToIntTypeMap.cbegin();
  std::map<GSyst_t, InteractionType_t>::const_iterator end = fGSystToIntTypeMap.cend();
  while ( it != end ) {

    GSyst_t syst = it->first;
    InteractionType_t type = it->second;

    // Note: this assumes that the error is symmetric.
    // TODO: consider changing this to handle asymmetric errors on the normalization
    double frac_err_norm = fracerr->OneSigmaErr( syst );

    NormMapEntry& entry = fNormMap.at( type );
    entry.fNormCurr = std::max(0.,
      entry.fNormDef * (1. + entry.fTwkDial * frac_err_norm));

    ++it;
  }

}
//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeight(const genie::EventRecord& event)
{
  bool is_mec = event.Summary()->ProcInfo().IsMEC();
  if ( !is_mec ) return 1.;

  double weight = this->CalcWeightNorm( event );
    std::cout << "Weight norm: " << weight << std::endl;
  weight *= this->CalcWeightAngularDist( event );
    std::cout << "Weight ang: " << weight << std::endl;
  weight *= this->CalcWeightPNDelta( event );
    std::cout << "Weight pndel: " << weight << std::endl;
  weight *= this->CalcWeightXSecShape( event );
    std::cout << "Weight xsecshape: " << weight << std::endl;
  weight *= this->CalcWeightXSecShape_Empirical( event );
    std::cout << "Weight xsecshape_empirical: " << weight << std::endl;
  weight *= this->CalcWeightXSecShape_Martini( event );
    std::cout << "Weight xsecshape_martini: " << weight << std::endl;
  weight *= this->CalcWeightEnergyDependence( event );
    std::cout << "Weight edep: " << weight << std::endl;
  return weight;
}
//_______________________________________________________________________________________
// Could expand on ideas to further look into effects of reweighting in dependence on 
// initial state nucleon pair content
/*
double GReWeightXSecMEC::SetFlagispnevent(const genie::EventRecord& event)
{
  bool is_mec = event.Summary()->ProcInfo().IsMEC();
  if ( !is_mec ) return 1.;

  bool Flagispnevent  = this->CalcWeightNorm( event );
  Flagispnevent *= this->CalcWeightAngularDist( event );
  Flagispnevent *= this->CalcWeightPNDelta( event );
  Flagispnevent *= this->CalcWeightXSecShape( event );
  return Flagispnevent;
}
*/
//_______________________________________________________________________________________
void GReWeightXSecMEC::Init(void) {

  // Set the tweak dials to their default values
  fDecayAngTwkDial = 0.;
  fDecayAng2TwkDial = 0.;
  fFracPN_CCTwkDial = 0.;
  fFracDelta_CCTwkDial = 0.;
  fCCXSecShapeTwkDial = 0.;
  fCCXSecShapeEmpiricalTwkDial = 0.;
  fCCXSecShapeMartiniTwkDial = 0.;
  fEnergyDependenceTwkDial = 0.;

  // Set the default normalization for each interaction type (tweak dial = 0
  // corresponds to a normalization factor of 1)
  std::map<GSyst_t, InteractionType_t>::const_iterator it = fGSystToIntTypeMap.cbegin();
  std::map<GSyst_t, InteractionType_t>::const_iterator end = fGSystToIntTypeMap.cend();

  while ( it != end ) {
    fNormMap[ it->second ] = NormMapEntry(0., 1., 1.);
    ++it;
  }

  // Get the default CCMEC cross section model from the current tune
  // TODO: add retrieval of NCMEC, EMMEC
  genie::AlgFactory* algf = genie::AlgFactory::Instance();
  genie::AlgConfigPool* conf_pool = genie::AlgConfigPool::Instance();
  genie::Registry* gpl = conf_pool->GlobalParameterList();

  RgAlg cc_def_id = gpl->GetAlg( "XSecModel@genie::EventGenerator/MEC-CC" );
  fXSecAlgCCDef = dynamic_cast< XSecAlgorithmI* >( algf->AdoptAlgorithm(
    cc_def_id.name, cc_def_id.config) );
  assert( fXSecAlgCCDef );
  fXSecAlgCCDef->AdoptSubstructure();

  //// Pull out the integrator used for the default model so we
  //// can integrate the alternate one using the same kinematic limits.
  //// This ensures proper PDF normalization for the XSecShape weight
  //// calculation.
  //fXSecIntegrator = dynamic_cast<const XSecIntegratorI*>(
  //  fXSecAlgCCDef->SubAlg("NumericalIntegrationAlg") );
  //assert( fXSecIntegrator );

  // Get the "fast" configuration of MECXSec to use to integrate
  // the alternate MEC model in the CalcWeightShape member function
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI*>( algf->AdoptAlgorithm(
    "genie::MECXSec", "Fast") );
  assert( fXSecIntegrator );

  // Get an alternate (Valencia) CCMEC cross section model for reshaping 
  // the default SuSAv2 model (for XSecShape_CCMEC parameter)
  // TODO: change hard-coding here, or add different knobs for different
  // target models
  AlgId alt_id( "genie::NievesSimoVacasMECPXSec2016", "Default" );
  //  AlgId alt_id( "genie::SuSAv2MECPXSec", "Default" );
  //  AlgId alt_id( "genie::EmpiricalMECPXSec2015", "Reweight" );
  fXSecAlgCCAlt_Nieves = dynamic_cast< XSecAlgorithmI* >( algf->AdoptAlgorithm(alt_id) );
  assert( fXSecAlgCCAlt_Nieves );
  fXSecAlgCCAlt_Nieves->AdoptSubstructure();

  // Get an alternate (SuSAv2) CCMEC cross section model for reshaping 
  // the default Valencia model (for XSecShape_CCMEC parameter)
  //  TODO: This does not actually work yet (so I set the alternative model to 
  //  "genie::EmpiricalMECPXSec2015" instead of "genie::SuSAv2MECPXSec"). This may be 
  //  fixed by looking into the details of "AdoptAlgorithm".  
  //  AlgId alt_id2( "genie::NievesSimoVacasMECPXSec2016", "Default" );
  //  AlgId alt_id2( "genie::SuSAv2MECPXSec", "Default" );
  AlgId alt_id2( "genie::EmpiricalMECPXSec2015", "Reweight" );
  fXSecAlgCCAlt_SuSAv2 = dynamic_cast< XSecAlgorithmI* >( algf->AdoptAlgorithm(alt_id2) );
  assert( fXSecAlgCCAlt_SuSAv2 );
  fXSecAlgCCAlt_SuSAv2->AdoptSubstructure();

  // Get another alternate CCMEC cross section model (Empirical) for reshaping the
  // default (SuSAv2 or Valencia) (for the XSecShape_CCMEC_Empirical parameter)
  AlgId alt_id3( "genie::EmpiricalMECPXSec2015", "Reweight" );
  fXSecAlgCCAlt_Empirical = dynamic_cast< XSecAlgorithmI* >( algf->AdoptAlgorithm(alt_id3) );
  assert( fXSecAlgCCAlt_Empirical );
  fXSecAlgCCAlt_Empirical->AdoptSubstructure();  

  // Get another alternate CCMEC cross section model (Martini) for reshaping the
  // default (SuSAv2 or Valencia) (for the XSecShape_CCMEC_Martini parameter)
  // TODO: Change the following line once the Martini model is available in GENIE. Currently
  // the XSecShape_CCMEC_Martini parameter will set the alternative to the Empirical model.
  //  AlgId alt_id4( "genie::SuSAv2MECPXSec", "Default" );
  AlgId alt_id4( "genie::EmpiricalMECPXSec2015", "Reweight" );
  // AlgId alt_id4( "genie::MartiniMECPXSec2024", "Default" );
  fXSecAlgCCAlt_Martini = dynamic_cast< XSecAlgorithmI* >( algf->AdoptAlgorithm(alt_id4) );
  assert( fXSecAlgCCAlt_Martini );
  fXSecAlgCCAlt_Martini->AdoptSubstructure();
}
//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeightNorm(const genie::EventRecord& event)
{
  InteractionType_t type = event.Summary()->ProcInfo().InteractionTypeId();

  // Find the tweak dial information for the current event's interaction type.
  // If a match isn't found, then just return a weight of unity.
  std::map<InteractionType_t, NormMapEntry>::const_iterator it = fNormMap.find( type );
  if ( it == fNormMap.cend() ) {
    LOG("ReW", pWARN) << "Unrecognized MEC event encountered in"
      << " GReWeightXSecMEC::CalcWeightNorm()";
    return 1.;
  }

  // If the tweak dial is set to zero (or is really small) then just return a
  // weight of unity
  double twk_dial = it->second.fTwkDial;
  bool tweaked = ( std::abs(twk_dial) > controls::kASmallNum );
  if ( !tweaked ) return 1.;

  double weight = it->second.fNormCurr;

  return weight;
}
//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeightAngularDist(const genie::EventRecord& event)
{
  // Only tweak dial values on the interval [0, 1] make sense for the angular
  // distribution. Enforce this here regardless of what the user requested.
  // double twk_dial = std::max( std::min(1., fDecayAngTwkDial), 0. );
  // Since there are yet only ad-hoc assumption on how DecayAngMEC on the angular distribution 
  // dependence, I added more freedom in varying this tweak dial, changed it to [-1, 1] - L. Bathe-Peters
  double twk_dial = std::max( std::min(1., fDecayAngTwkDial), -1. );

  // Second tweak dial that chagnes the frequeny of the harmonic function that the 
  // DecayAngMEC parameter is reweighting to.
  double twk_dial2 = fDecayAng2TwkDial;

  // If the tweak dial is set to zero (or is really small) then just return a
  // weight of unity
  bool tweaked = ( ( std::abs(twk_dial) > controls::kASmallNum ) || ( std::abs(twk_dial2) > controls::kASmallNum ) );
  if ( !tweaked ) return 1.;

  // Get the daughters of the recoiled two-nucleon cluster
  // TODO: Consider using something less fragile here. Right now, this relies
  // on the observation that MECGenerator.cxx always places the recoiling
  // nucleon cluster at position 5 in the GENIE event record.
  const int recoil_nucleon_cluster_pos = 5;
  GHepParticle* nucleon_cluster = event.Particle( recoil_nucleon_cluster_pos );
  assert( nucleon_cluster );

  // Make sure that the retrieved particle is really a two-nucleon cluster. If
  // it isn't, just complain and return a unit weight.
  int cluster_pdg = nucleon_cluster->Pdg();
  if ( !pdg::Is2NucleonCluster(cluster_pdg) ) {
    LOG("ReW", pERROR) << "Invalid two-nucleon cluster PDG code " << cluster_pdg
      << " encountered in GReWeightXSecMEC::CalcWeightAngularDist()";
    return 1.;
  }

  // The two-nucleon cluster should have exactly two daughters (the two
  // final-state nucleons). If it doesn't complain and return a unit weight.
  int first = nucleon_cluster->FirstDaughter();
  int last = nucleon_cluster->LastDaughter();
  if ( !nucleon_cluster->HasDaughters() || (last - first) != 1 ) {
    LOG("ReW", pERROR) << "Invalid number of daughters for a two-nucleon"
      << " cluster encountered in GReWeightXSecMEC::CalcWeightAngularDist()";
    return 1.;
  }

  // Get the two final-state nucleons
  GHepParticle* N1 = event.Particle( first );
  GHepParticle* N2 = event.Particle( last );

  // If one of them isn't really a nucleon, complain and return a unit weight
  if ( !pdg::IsNucleon(N1->Pdg()) || !pdg::IsNucleon(N2->Pdg()) ) {
    LOG("ReW", pERROR) << "Non-nucleon daughter of a two-nucleon"
      << " cluster encountered in GReWeightXSecMEC::CalcWeightAngularDist()";
    return 1.;
  }

  // Get the 4-momenta of the two outgoing nucleons
  TLorentzVector p4N1 = *N1->P4();
  TLorentzVector p4N2 = *N2->P4();

  // Boost the 4-momenta of the two nucleons from the lab frame to their
  // CM frame (which is also the rest frame of the recoiling nucleon cluster)
  TLorentzVector p4Cluster = p4N1 + p4N2;
  TVector3 boostToCM = -p4Cluster.BoostVector();

  p4N1.Boost( boostToCM );
  p4N2.Boost( boostToCM );

  // Also get the 4-momenta of the initial and final leptons. These will be
  // used to compute the 4-momentum transfer
  TLorentzVector p4Probe = *event.Probe()->P4();
  TLorentzVector p4Lep = *event.FinalStatePrimaryLepton()->P4();

  TLorentzVector q4 = p4Probe - p4Lep;

  // Boost the 4-momentum transfer into the two-nucleon CM frame
  q4.Boost( boostToCM );

  // Use the 3-momentum transfer in the two-nucleon CM frame as the reference
  // z-axis for the altered angular distribution
  TVector3 q3 = q4.Vect().Unit();

  // Determine a rotation axis and angle that will cause the 3-momentum to
  // point along the +z direction
  TVector3 zvec(0., 0., 1.);
  TVector3 rot = ( q3.Cross(zvec) ).Unit();
  double angle = zvec.Angle( q3 );

  // Handle the edge case where q3 is along -z, so the
  // cross product above vanishes
  if ( q3.Perp() == 0. && q3.Z() < 0. ) {
    rot = TVector3(0., 1., 0.);
    angle = constants::kPi;
  }

  // If the rotation vector is non-null (within numerical precision) then
  // rotate the CM frame 3-momentum of nucleon #1 into a frame where q3 points along +z
  TVector3 p3N1 = p4N1.Vect();
  if ( rot.Mag() >= controls::kASmallNum ) {
    p3N1.Rotate(angle, rot);
  }

  // We now have what we need. Compute the emission angles for nucleon #1 relative to the
  // 3-momentum transfer in the rest frame of the recoiling nucleon cluster.
  double theta_N1 = p3N1.Theta();
  //double phi_N1 = p3N1.Phi();

  // Default model (used by all current GENIE MEC implementations) is to decay
  // the recoiling nucleon cluster isotropically. The alternate model is to
  // decay it according to (3/2)*cos^2(theta) in the CM frame, with the
  // 3-momentum transfer along the +z direction. The tweak dial linearly
  // interpolates between purely isotropic (0) and purely the alternate
  // distribution (1).
  // TODO: come up with something better for the alternate distribution
  // double weight = 3.*twk_dial*std::pow(std::cos(theta_N1), 2) + (1. - twk_dial);
  double weight = 3.*twk_dial*std::pow(std::cos(twk_dial2*theta_N1), 2) + (1. - twk_dial);

  return weight;
}
//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeightPNDelta(const genie::EventRecord& event)
{
  // Only handle CC events for now (and return unit weight for the others)
  // TODO: Add capability to tweak nucleon pair isospin for NC, EM
  InteractionType_t type = event.Summary()->ProcInfo().InteractionTypeId();
  if ( type != kIntWeakCC ) return 1.;

  // If the tweak dial is set to zero (or is really small) then just return a
  // weight of unity
  bool tweaked = ( std::abs(fFracPN_CCTwkDial) > controls::kASmallNum
    || std::abs(fFracDelta_CCTwkDial) > controls::kASmallNum );
  if ( !tweaked ) return 1.;

  // Enforce that the current event involves an initial nucleon cluster.
  // Determine whether the cluster is p+n
  GHepParticle* initial_nucleon_cluster = event.HitNucleon();
  assert( initial_nucleon_cluster );

  int two_nuc_pdg = initial_nucleon_cluster->Pdg();
  assert( pdg::Is2NucleonCluster(two_nuc_pdg) );

  // Whether the current event involved an initial pn nucleon cluster
  bool is_pn_event = ( two_nuc_pdg == kPdgClusterNP );

  // Whether the current event involved a virtual delta resonance
  // (currently only used by the Valencia model). If this is the
  // case, then a resonance will be set in the interaction.
  //  bool is_delta_event = event.Summary()->ExclTag().KnownResonance();

  // Calculate the model's default fraction of initial pn pairs. For empirical
  // MEC, this is just a fixed number from the configuration. For Valencia, it's
  // dependent on the kinematics, so we'll need to compute ratios of
  // differential cross sections.
  double pn_frac_def = 0.;

  // Also calculate the model's default fraction of internal deltas. This
  // is only used by the Valencia model for now.
  double delta_frac_def = 0.;

  std::string cc_def_alg_name = fXSecAlgCCDef->Id().Name();
  // std::cout << "GENIE model: " << cc_def_alg_name << std::endl;

  if ( cc_def_alg_name == "genie::EmpiricalMECPXSec2015" ) {
    // For GENIE's empirical MEC model, the pn fraction is not dependent on
    // kinematics. We can just retrieve the value from the model configuration.
    pn_frac_def = fXSecAlgCCDef->GetConfig().GetDouble( "EmpiricalMEC-FracPN_CC" );
    // The empirical MEC model doesn't account for internal delta resonances
    // explicitly. We've set the default delta fraction to zero above, but just in
    // case, let's repeat that here.
    delta_frac_def = 0.;
  }
  else if ( cc_def_alg_name == "genie::NievesSimoVacasMECPXSec2016" ) {
    // For the Valencia MEC model, the pn fraction can vary with q0 and q3. We
    // can get the pn fraction for this event's kinematics by computing the
    // differential cross section for each case. A similar thing is done in
    // MECGenerator in order to decide which kind of nucleon pair is hit. See
    // genie::MECGenerator::GenerateNSVInitialHadrons() for details.

    // Get the differential cross section for an initial pn pair and the total
    // for all pair types (this is how the Valencia calculation is organized in
    // GENIE). Note that the Valencia MEC model works in the kPSTlctl
    // phase space. Clone the input interaction so that we can modify
    // the PDG code of the initial nucleon cluster.
    Interaction* interaction = new Interaction( *event.Summary() );

    // TODO: When NC and/or EM interactions are added for Valencia, generalize
    // this for use those. Unlike CC, all three pair types can participate in
    // NC & EM.

    // Get the differential cross section for an initial pn pair. Clear any
    // set resonance so that we get the total differential cross section.
    interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( kPdgClusterNP );
    interaction->ExclTagPtr()->SetResonance( kNoResonance );
    double xsec_pn = fXSecAlgCCDef->XSec( interaction, kPSTlctl );

    // Get the total differential cross section (the resonance is still
    // cleared)
    interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( 0 );
    double xsec_tot = fXSecAlgCCDef->XSec( interaction, kPSTlctl );

    // Get the total differential cross section for an internal resonance
    interaction->ExclTagPtr()->SetResonance( kP33_1232 );
    double xsec_tot_delta = fXSecAlgCCDef->XSec( interaction, kPSTlctl );

    // Get the differential cross section for an internal resonance and an
    // initial pn pair (the resonance is still set)
    interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( kPdgClusterNP );
    double xsec_pn_delta = fXSecAlgCCDef->XSec( interaction, kPSTlctl );

    // We don't need the cloned interaction anymore, so delete it
    delete interaction;

    assert( xsec_tot > 0. );
    pn_frac_def = xsec_pn / xsec_tot;

    // Compute the delta fraction appropriate for the initial nucleon cluster
    // type sampled in this event. This allows us to maintain consistency when
    // potentially tweaking both the pn fraction and delta fraction.
    if ( is_pn_event ) {
      assert( xsec_pn > 0. );
      delta_frac_def = xsec_pn_delta / xsec_pn;
    }
    else {
      double xsec_nn = xsec_tot - xsec_pn;
      double xsec_nn_delta = xsec_tot_delta - xsec_pn_delta;
      assert( xsec_nn > 0. );
      delta_frac_def = xsec_nn_delta / xsec_nn;
    }

  }
  else if ( cc_def_alg_name == "genie::SuSAv2MECPXSec" ) {
    // Get the differential cross section for an initial pn pair and the total
    // for all pair types. Delta resonances are not used in the SuSAv2 model, so there is no extension
    // of the PNFrac_CCMEC to the DeltaNotDelta_parameter. Clone the input interaction.
    // TODO: Could think of a way to include delta resonances in the reweighting of the SuSAv2 model
    // predictions.
    Interaction* interaction = new Interaction( *event.Summary() );

    // Retrieve the fraction of initial pn pairs from GENIE. 
    pn_frac_def = dynamic_cast< const genie::SuSAv2MECPXSec* >( fXSecAlgCCDef )->PairRatio( interaction );

    // We don't need the cloned interaction anymore, so delete it
    delete interaction;
  } 
  else {
    LOG("ReW", pERROR) << "Unrecognized MEC model " << cc_def_alg_name
      << " encountered in genie::rew::GReWeightXSecMEC::CalcWeightPNDelta()";
    return 1.;
  }

  // Relax these sanity checks for now. Interpolation of the hadron tensors
  // can allow the pn or delta fraction to move a bit above one. Just
  // force it to be on [0, 1] instead of returning a unit weight.

  //// Check that the pn fraction computed above is sane. If not, complain and
  //// return a unit weight.
  //bool impossible_pp_or_nn_event = ( pn_frac_def == 1. && !is_pn_event );
  //if ( pn_frac_def < 0. || pn_frac_def > 1. || impossible_pp_or_nn_event ) {
  //  LOG("ReW", pERROR) << "Invalid pn fraction value " << pn_frac_def
  //    << " encountered in genie::rew::GReWeightXSecMEC::CalcWeightPNDelta()";
  //  return 1.;
  //}

  //// Do the same for the delta fraction.
  //bool impossible_delta_event = ( delta_frac_def == 1. && !is_delta_event );
  //if ( delta_frac_def < 0. || delta_frac_def > 1. || impossible_delta_event ) {
  //  LOG("ReW", pERROR) << "Invalid delta fraction value " << delta_frac_def
  //    << " encountered in genie::rew::GReWeightXSecMEC::CalcWeightPNDelta()";
  //  return 1.;
  //}

  // Force the default fractions to be on the interval [0, 1] for sanity's sake.
  // This will counteract interpolation problems.
  pn_frac_def = std::max( std::min(1., pn_frac_def), 0. );
  delta_frac_def = std::max( std::min(1., delta_frac_def), 0. );

  // TODO: add support for asymmetric errors here
  GSystUncertainty* gsu = GSystUncertainty::Instance();
  double frac_err_pn_cc = gsu->OneSigmaErr( kXSecTwkDial_FracPN_CCMEC );
  double frac_err_delta_cc = gsu->OneSigmaErr( kXSecTwkDial_FracDelta_CCMEC );

  // Compute the scaling factor for the pn fraction that corresponds to the
  // current tweak dial setting
  double pn_tweak_factor = ( 1. + fFracPN_CCTwkDial * frac_err_pn_cc );
  double delta_tweak_factor = ( 1. + fFracDelta_CCTwkDial * frac_err_delta_cc );

  // To conserve the total cross section (which is separately controlled by the
  // normalization tweak dials), enforce that the tweaked pn fraction lies on
  // the interval [0, 1]. Do the same for the delta fraction.
  double pn_frac_tweak = std::max( std::min(1., pn_frac_def * pn_tweak_factor), 0. );
  double delta_frac_tweak = std::max( std::min(1., delta_frac_def * delta_tweak_factor), 0. );

  // Assign the appropriate likelihood ratio as the weight. Note that
  // we've already checked that pn_frac_def lies in a reasonable range
  // above, so we can divide as shown without worrying about NaNs.
  double weight;
  if ( is_pn_event ) weight = pn_frac_tweak / pn_frac_def;
  else weight = ( 1. - pn_frac_tweak ) / ( 1. - pn_frac_def );

  // Also multiply by the tweaked delta fraction (skip these two lines if you want do not want to consider
  // delta resonances (FracPN parameter), otherwise delta resonances are considered ( DeltaNotDelta parameter).
  // This only works for an input tune with the 2p-2h Valencia model used for CCMEC.
  //  if ( is_delta_event ) weight *= delta_frac_tweak / delta_frac_def; // 1.; //delta_frac_tweak / delta_frac_def;
  // else weight *= ( 1. - delta_frac_tweak ) / ( 1. - delta_frac_def ); // 1.; //( 1. - delta_frac_tweak ) / ( 1. - delta_frac_def );
  LOG("ReW", pDEBUG) << "pn_twk_dial = " << fFracPN_CCTwkDial << ", frac_err = "
    << frac_err_pn_cc;
  LOG("ReW", pDEBUG) << "pn_frac_def = " << pn_frac_def << ", pn_frac_tweak = "
    << pn_frac_tweak;
  LOG("ReW", pDEBUG) << "delta_twk_dial = " << fFracDelta_CCTwkDial << ", frac_err = "
    << frac_err_delta_cc;
  LOG("ReW", pDEBUG) << "delta_frac_def = " << delta_frac_def << ", delta_frac_tweak = "
    << delta_frac_tweak << ", weight = " << weight;

  std::cout << "Weight: " << weight << std::endl;

  return weight;
}
//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeightXSecShape(const genie::EventRecord& event)
{
  // The XSecShape_CCMEC parameter reweights from the default SuSAv2 CCMEC
  // model to the Valencia CCMEC model (or vice versa in case a tune with
  // the Valencia CCMEC cross section model is used as input).
  // TODO: The reweighting from SuSAv2 to Valencia is very slow. Find ways speed
  // this up.
  // TODO: The reweighting from Valencia to SuSAv2 does not work yet. See comment
  // above for "AlgId alt_id2". Implement this.

  // Only handle CC events for now (and return unit weight for the others)
  // TODO: Add capability to tweak the shape for NC and EM
  InteractionType_t type = event.Summary()->ProcInfo().InteractionTypeId();
  if ( type != kIntWeakCC ) return 1.;

  // Only tweak dial values on the interval [0, 1] make sense for this
  // knob (0 is pure Valencia shape; 1 is pure GENIE empirical shape).
  // Enforce this here regardless of what the user requested.
  double twk_dial = std::max( std::min(1., fCCXSecShapeTwkDial), 0. );

  // If the tweak dial is set to zero (or is really small) then just return a
  // weight of unity
  bool tweaked = ( std::abs(twk_dial) > controls::kASmallNum );
  if ( !tweaked ) return 1.;

  // Compute the probability density for generating the selected kinematics
  // under the default cross section model
  // TODO: add check that the default model predictions match those stored
  // in the event record

  // Clone the input interaction so that we can clear the set nucleon cluster
  // PDG code and resonance flags. The total cross section stored in the
  // event record for the Valencia model includes all contributions.
  Interaction* interaction = new Interaction( *event.Summary() );

  // Double-check that the running value of the lepton kinetic energy
  // is set in the input interaction. If it isn't, set it manually
  // using the lepton 4-momentum.
  genie::Kinematics* kine_ptr = interaction->KinePtr();
  if ( !kine_ptr->KVSet(kKVTl) ) {

    // Final lepton mass
    double ml = interaction->FSPrimLepton()->Mass();
    // Final lepton 4-momentum
    const TLorentzVector& p4l = kine_ptr->FSLeptonP4();
    // Final lepton kinetic energy
    double Tl = p4l.E() - ml;
    // Final lepton scattering cosine
    double ctl = p4l.CosTheta();

    kine_ptr->SetKV( kKVTl, Tl );
    kine_ptr->SetKV( kKVctl, ctl );
  }

  // Set q0 and q3 in the interaction.
  if ( !kine_ptr->KVSet(kKVQ0) ) {

    // Final lepton 4-momentum
    const TLorentzVector& p4l = kine_ptr->FSLeptonP4();
    // Final lepton total energy
    double El = p4l.E();

    // Probe 4-momentum
    const InitialState& init_state = interaction->InitState();
    TLorentzVector* p4v = init_state.GetProbeP4( kRfLab );

    // Probe total energy
    double Ev = p4v->E();

    // Energy transfer
    double q0 = Ev - El;

    // Magnitude of the momentum transfer
    double q3 = ( (*p4v) - p4l ).Vect().Mag();

    kine_ptr->SetKV( kKVQ0, q0 );
    kine_ptr->SetKV( kKVQ3, q3 );
  }

  // Get the differential and total cross section for the default
  // MEC model (including contributions from both kinds of
  // initial nucleon clusters and both kinds of diagrams)
  // Save the hit nucleon cluster PDG code (we'll need to set it again
  // for the alternate modelNievesSimoVacasMECPXSec2016)
  int hit_nuc_pdg = interaction->InitState().Tgt().HitNucPdg();
  interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( 0 );
  interaction->ExclTagPtr()->SetResonance( kNoResonance );

  double diff_xsec_def = fXSecAlgCCDef->XSec( interaction,
    kPSTlctl);
  double tot_xsec_def = event.XSec();
  double prob_density_def = diff_xsec_def / tot_xsec_def;

  //double check_tot_xsec_def = fXSecAlgCCDef->Integral( interaction );
  //LOG("ReW", pDEBUG) << "check_tot_xsec_def = " << check_tot_xsec_def;

  // Get the names of the CCMEC cross section default and alternate models
  std::string cc_def_alg_name = fXSecAlgCCDef->Id().Name();
  std::string cc_alt_alg_name = fXSecAlgCCAlt_Nieves->Id().Name();
  std::string cc_def_alg_config = fXSecAlgCCDef->Id().Config();
  std::string cc_alt_alg_config = fXSecAlgCCAlt_Nieves->Id().Config();
  std::string cc_alt2_alg_name = fXSecAlgCCAlt_SuSAv2->Id().Name();

  // Set the hit nucleon cluster PDG code to its sampled value
  // (empirical MEC needs it)
  // Note: In order for this code to be more flexible, the following if-condition
  // makes sure that in case either fXSecAlgCCAlt_Nieves or fXSecAlgCCAlt_SuSAv2
  // is set to "genie::EmpiricalMECPXSec2015", one can still use the XSecShape_CCMEC
  // parameter. However, to reweight to the Empirical CCMEC model, the 
  // XSecShape_Empirical_CCMEC parameter was added, so that one can 
  // reweight to two models simultaneously and also in order to not
  // having to change the hard-coded part of the code.
  if ( ( cc_def_alg_name == "genie::SuSAv2MECPXSec" && cc_alt_alg_name == "genie::EmpiricalMECPXSec2015" ) || ( cc_def_alg_name == "genie::NievesSimoVacasMECPXSec2016" && cc_alt2_alg_name == "genie::EmpiricalMECPXSec2015" ) ) {
    interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( hit_nuc_pdg );
    // std::cout << "alternate model is empirical" << std::endl;
  }

  //// The empirical MEC differential cross section doesn't check the
  //// kinematic limits imposed by MECGenerator, so enforce them here.
  //// If we're outside, then don't bother calculating the alternative
  //// model differential cross section.
  //double diff_xsec_alt = 0.;
  //Range1D_t rW = getWLimitsEmpiricalMEC();
  //Range1D_t rQ2 = getQ2LimitsEmpiricalMEC();

  //double W = interaction->Kine().W();
  //double Q2 = interaction->Kine().Q2();

  //LOG("RwMEC", pERROR) << "NIEVES: W = " << W << ", Q2 = " << Q2;
  //if ( rW.min <= W && rW.max >= W && rQ2.min <= Q2 && rQ2.max >= Q2 ) {

  // If none of the following if statements are true, that means that the
  // input tune does not have the SuSAv2 or Valencia model as its default
  // CCMEC cross section model. In this case, just return a weight of 1 as 
  // there are currently no plans to have either model as the default model.
  double diff_xsec_alt = diff_xsec_def;
  double tot_xsec_alt = tot_xsec_def;

  if ( cc_def_alg_name == "genie::SuSAv2MECPXSec" ) {
    // If the SuSAv2 model is the default CCMEC cross section model from the 
    // input tune, compute the differential and total cross section of GENIE's 
    // Valencia MEC (alternative) model

    // LOG("ReW", pINFO) << "Input (default) CCMEC cross section model: SuSAv2";  
    // std::cout << "Input (default) CCMEC cross section model name: " << cc_def_alg_name << std::endl;
    // std::cout << "Alternative CCMEC cross section model name:     " << cc_alt_alg_name << std::endl;

    // std::cout << "Input (default) CCMEC cross section model config: " << cc_def_alg_config << std::endl;
    // std::cout << "Alternative CCMEC cross section model config:     " << cc_alt_alg_config << std::endl;

    diff_xsec_alt = fXSecAlgCCAlt_Nieves->XSec( interaction, kPSTlctl );

    //}

    tot_xsec_alt = this->GetXSecIntegral( fXSecAlgCCAlt_Nieves, interaction );

  }
  else if ( cc_def_alg_name == "genie::NievesSimoVacasMECPXSec2016" ) {
    // If the Valencia model is the default CCMEC cross section model from the 
    // input tune, compute the differential and total cross section of GENIE's 
    // SuSAv2 MEC (alternative) model

    // LOG("ReW", pINFO) << "Input (default) CCMEC cross section model: Valencia";
    // std::cout << "Input (default) CCMEC cross section model: " << cc_def_alg_name << std::endl;
    // std::cout << "Alternative CCMEC cross section model: " << cc_alt2_alg_name << std::endl; 
    // std::cout << "Alternative CCMEC cross section model: " << cc_alt_alg_name << std::endl;

    diff_xsec_alt = fXSecAlgCCAlt_SuSAv2->XSec( interaction, kPSTlctl ); 
    // diff_xsec_alt = fXSecAlgCCAlt_Nieves->XSec( interaction, kPSTlctl );
    //}

    tot_xsec_alt = this->GetXSecIntegral( fXSecAlgCCAlt_SuSAv2, interaction ); 
    // tot_xsec_alt = this->GetXSecIntegral( fXSecAlgCCAlt_Nieves, interaction );
    }
  else {
    // If the Empirical or Martini or any other model is the default CCMEC cross section model 
    // from the input tune, just return a weight of 1 as there are currently no 
    // plans to have either model as the default model.

    LOG("ReW", pWARN) << "MEC xsecshape reweighting for other CCMEC models but SuSAv2 or Valencia model not implemented";
    // std::cout << "Input (default) CCMEC cross section model: Other, i.e. " << cc_def_alg_name << std::endl;
    // std::cout << "MEC xsecshape reweighting for other CCMEC models but SuSAv2 or Valencia model not implemented" << std::endl;

    //  double diff_xsec_alt = diff_xsec_def; // already set above

    //}

    //  double tot_xsec_alt = tot_xsec_def; // already set above

  }


  //LOG("RwMEC", pERROR) << "diff_xsec_alt = " << diff_xsec_alt << ", tot_xsec_alt = " << tot_xsec_alt;

  //if ( tot_xsec_alt == 0. && diff_xsec_alt != 0. ) LOG("RwMEC", pERROR) << "OH NO!";

  // Protect against NaNs when the total cross section for the alternative
  // model is zero
  if ( tot_xsec_alt == 0. ) {
    diff_xsec_alt = 0.;
    tot_xsec_alt = 1.;
  }
  double prob_density_alt = diff_xsec_alt / tot_xsec_alt;

//  std::cout << "diff_xsec_def = " << std::setprecision(30) << diff_xsec_def << std::endl;
//  std::cout << "diff_xsec_alt = " << std::setprecision(30) << diff_xsec_alt << std::endl;
//  std::cout << "tot_xsec_def  = " << std::setprecision(30) << tot_xsec_def << std::endl;
//  std::cout << "tot_xsec_alt  = " << std::setprecision(30) << tot_xsec_alt << std::endl;

  // std::cout << "diff_xsec_def = " << diff_xsec_def << std::endl;
  // std::cout << "diff_xsec_alt = " << diff_xsec_alt << std::endl;
  // std::cout << "tot_xsec_def  = " << tot_xsec_def << std::endl;
  // std::cout << "tot_xsec_alt  = " << tot_xsec_alt << std::endl;

  // std::cout << "prob_density_def     = " << prob_density_def << std::endl;
  // std::cout << "prob_density_alt     = " << prob_density_alt << std::endl;

  // Compute a new probability density for this event by interpolating between
  // the two models while preserving the total cross section
  double tweaked_prob_density = (1. - twk_dial)*prob_density_def + twk_dial*prob_density_alt;

  // std::cout << "tweaked_prob_density = " << tweaked_prob_density << std::endl;

  // The weight is then the likelihood ratio
  double weight = tweaked_prob_density / prob_density_def;

  // std::cout << "weight = " << weight << std::endl;

  LOG("ReW", pDEBUG) << "xsec_def = " << diff_xsec_def << ", xsec_alt = " << diff_xsec_alt;
  LOG("ReW", pDEBUG) << "tot_xsec_def = " << tot_xsec_def << ", tot_xsec_alt = " << tot_xsec_alt;
  LOG("ReW", pDEBUG) << "twk_dial = " << fCCXSecShapeTwkDial << ", prob_density_def = "
    << prob_density_def << ", prob_density_alt = " << prob_density_alt << ", weight = " << weight;

  // We don't need the cloned interaction anymore, so delete it
  delete interaction;

  return weight;
}

//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeightXSecShape_Empirical(const genie::EventRecord& event)
{
  // Only handle CC events for now (and return unit weight for the others)
  // TODO: Add capability to tweak the shape for NC and EM
  InteractionType_t type = event.Summary()->ProcInfo().InteractionTypeId();
  if ( type != kIntWeakCC ) return 1.;

  // Only tweak dial values on the interval [0, 1] make sense for this
  // knob (0 is pure Valencia shape; 1 is pure GENIE empirical shape).
  // Enforce this here regardless of what the user requested.
  double twk_dial = std::max( std::min(1., fCCXSecShapeEmpiricalTwkDial), 0. );

  // If the tweak dial is set to zero (or is really small) then just return a
  // weight of unity
  bool tweaked = ( std::abs(twk_dial) > controls::kASmallNum );
  if ( !tweaked ) return 1.;

  // Compute the probability density for generating the selected kinematics
  // under the default cross section model
  // TODO: add check that the default model predictions match those stored
  // in the event record

  // Clone the input interaction so that we can clear the set nucleon cluster
  // PDG code and resonance flags. The total cross section stored in the
  // event record for the Valencia model includes all contributions.
  Interaction* interaction = new Interaction( *event.Summary() );

  // Double-check that the running value of the lepton kinetic energy
  // is set in the input interaction. If it isn't, set it manually
  // using the lepton 4-momentum.
  genie::Kinematics* kine_ptr = interaction->KinePtr();
  if ( !kine_ptr->KVSet(kKVTl) ) {

    // Final lepton mass
    double ml = interaction->FSPrimLepton()->Mass();
    // Final lepton 4-momentum
    const TLorentzVector& p4l = kine_ptr->FSLeptonP4();
    // Final lepton kinetic energy
    double Tl = p4l.E() - ml;
    // Final lepton scattering cosine
    double ctl = p4l.CosTheta();

    kine_ptr->SetKV( kKVTl, Tl );
    kine_ptr->SetKV( kKVctl, ctl );
  }

  // Set q0 and q3 in the interaction.
  if ( !kine_ptr->KVSet(kKVQ0) ) {

    // Final lepton 4-momentum
    const TLorentzVector& p4l = kine_ptr->FSLeptonP4();
    // Final lepton total energy
    double El = p4l.E();

    // Probe 4-momentum
    const InitialState& init_state = interaction->InitState();
    TLorentzVector* p4v = init_state.GetProbeP4( kRfLab );

    // Probe total energy
    double Ev = p4v->E();

    // Energy transfer
    double q0 = Ev - El;

    // Magnitude of the momentum transfer
    double q3 = ( (*p4v) - p4l ).Vect().Mag();

    kine_ptr->SetKV( kKVQ0, q0 );
    kine_ptr->SetKV( kKVQ3, q3 );
  }

  // Get the differential and total cross section for the default
  // MEC model (including contributions from both kinds of
  // initial nucleon clusters and both kinds of diagrams)
  // Save the hit nucleon cluster PDG code (we'll need to set it again
  // for the alternate modelNievesSimoVacasMECPXSec2016)
  int hit_nuc_pdg = interaction->InitState().Tgt().HitNucPdg();
  interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( 0 );
  interaction->ExclTagPtr()->SetResonance( kNoResonance );

  double diff_xsec_def = fXSecAlgCCDef->XSec( interaction,
    kPSTlctl);
  double tot_xsec_def = event.XSec();
  double prob_density_def = diff_xsec_def / tot_xsec_def;

  //double check_tot_xsec_def = fXSecAlgCCDef->Integral( interaction );
  //LOG("ReW", pDEBUG) << "check_tot_xsec_def = " << check_tot_xsec_def;

  // Get the names of the CCMEC cross section default and alternate models
  std::string cc_def_alg_name = fXSecAlgCCDef->Id().Name();
  // std::string cc_alt_alg_name = fXSecAlgCCAlt_Nieves->Id().Name();
  // std::string cc_alt2_alg_name = fXSecAlgCCAlt_SuSAv2->Id().Name();
  std::string cc_alt3_alg_name = fXSecAlgCCAlt_Empirical->Id().Name(); 
  // std::string cc_alt4_alg_name = fXSecAlgCCAlt_Martini->Id().Name();

  // Set the hit nucleon cluster PDG code to its sampled value
  // (empirical MEC needs it)
  if ( cc_alt3_alg_name == "genie::EmpiricalMECPXSec2015" ) {
    interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( hit_nuc_pdg );
    std::cout << "alternate model is empirical" << std::endl;
  }

  //// The empirical MEC differential cross section doesn't check the
  //// kinematic limits imposed by MECGenerator, so enforce them here.
  //// If we're outside, then don't bother calculating the alternative
  //// model differential cross section.
  //double diff_xsec_alt = 0.;
  //Range1D_t rW = getWLimitsEmpiricalMEC();
  //Range1D_t rQ2 = getQ2LimitsEmpiricalMEC();

  //double W = interaction->Kine().W();
  //double Q2 = interaction->Kine().Q2();

  //LOG("RwMEC", pERROR) << "NIEVES: W = " << W << ", Q2 = " << Q2;
  //if ( rW.min <= W && rW.max >= W && rQ2.min <= Q2 && rQ2.max >= Q2 ) {

  // std::cout << "Input (default) CCMEC cross section model: " << cc_def_alg_name << std::endl;
  // std::cout << "Alternative CCMEC cross section model: " << cc_alt3_alg_name << std::endl; 

  // Once the CCMEC Martini model is available, just change fXSecAlgCCAlt3 
  // (reweight from SuSAv2 or Valencia to the Empirical model) to 
  // fXSecAlgCCAlt4 (reweight from SuSAv2 or Valencia to the Martini model)
  // and this should just work
  double diff_xsec_alt = fXSecAlgCCAlt_Empirical->XSec( interaction, kPSTlctl ); 

  //}

  double tot_xsec_alt = this->GetXSecIntegral( fXSecAlgCCAlt_Empirical, interaction ); 

  //LOG("RwMEC", pERROR) << "diff_xsec_alt = " << diff_xsec_alt << ", tot_xsec_alt = " << tot_xsec_alt;

  //if ( tot_xsec_alt == 0. && diff_xsec_alt != 0. ) LOG("RwMEC", pERROR) << "OH NO!";

  // Protect against NaNs when the total cross section for the alternative
  // model is zero
  if ( tot_xsec_alt == 0. ) {
    diff_xsec_alt = 0.;
    tot_xsec_alt = 1.;
  }
  double prob_density_alt = diff_xsec_alt / tot_xsec_alt;

  // Compute a new probability density for this event by interpolating between
  // the two models while preserving the total cross section
  double tweaked_prob_density = (1. - twk_dial)*prob_density_def + twk_dial*prob_density_alt;

  // The weight is then the likelihood ratio
  double weight = tweaked_prob_density / prob_density_def;

  LOG("ReW", pDEBUG) << "xsec_def = " << diff_xsec_def << ", xsec_alt = " << diff_xsec_alt;
  LOG("ReW", pDEBUG) << "tot_xsec_def = " << tot_xsec_def << ", tot_xsec_alt = " << tot_xsec_alt;
  LOG("ReW", pDEBUG) << "twk_dial = " << fCCXSecShapeEmpiricalTwkDial << ", prob_density_def = "
    << prob_density_def << ", prob_density_alt = " << prob_density_alt << ", weight = " << weight;

  // We don't need the cloned interaction anymore, so delete it
  delete interaction;

  return weight;
}
//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeightXSecShape_Martini(const genie::EventRecord& event)
{
  // Only handle CC events for now (and return unit weight for the others)
  // TODO: Add capability to tweak the shape for NC and EM
  InteractionType_t type = event.Summary()->ProcInfo().InteractionTypeId();
  if ( type != kIntWeakCC ) return 1.;

  // Only tweak dial values on the interval [0, 1] make sense for this
  // knob (0 is pure Valencia shape; 1 is pure GENIE empirical shape).
  // Enforce this here regardless of what the user requested.
  double twk_dial = std::max( std::min(1., fCCXSecShapeMartiniTwkDial), 0. );

  // If the tweak dial is set to zero (or is really small) then just return a
  // weight of unity
  bool tweaked = ( std::abs(twk_dial) > controls::kASmallNum );
  if ( !tweaked ) return 1.;

  // Compute the probability density for generating the selected kinematics
  // under the default cross section model
  // TODO: add check that the default model predictions match those stored
  // in the event record

  // Clone the input interaction so that we can clear the set nucleon cluster
  // PDG code and resonance flags. The total cross section stored in the
  // event record for the Valencia model includes all contributions.
  Interaction* interaction = new Interaction( *event.Summary() );

  // Double-check that the running value of the lepton kinetic energy
  // is set in the input interaction. If it isn't, set it manually
  // using the lepton 4-momentum.
  genie::Kinematics* kine_ptr = interaction->KinePtr();
  if ( !kine_ptr->KVSet(kKVTl) ) {

    // Final lepton mass
    double ml = interaction->FSPrimLepton()->Mass();
    // Final lepton 4-momentum
    const TLorentzVector& p4l = kine_ptr->FSLeptonP4();
    // Final lepton kinetic energy
    double Tl = p4l.E() - ml;
    // Final lepton scattering cosine
    double ctl = p4l.CosTheta();

    kine_ptr->SetKV( kKVTl, Tl );
    kine_ptr->SetKV( kKVctl, ctl );
  }

  // Set q0 and q3 in the interaction.
  if ( !kine_ptr->KVSet(kKVQ0) ) {

    // Final lepton 4-momentum
    const TLorentzVector& p4l = kine_ptr->FSLeptonP4();
    // Final lepton total energy
    double El = p4l.E();

    // Probe 4-momentum
    const InitialState& init_state = interaction->InitState();
    TLorentzVector* p4v = init_state.GetProbeP4( kRfLab );

    // Probe total energy
    double Ev = p4v->E();

    // Energy transfer
    double q0 = Ev - El;

    // Magnitude of the momentum transfer
    double q3 = ( (*p4v) - p4l ).Vect().Mag();

    kine_ptr->SetKV( kKVQ0, q0 );
    kine_ptr->SetKV( kKVQ3, q3 );
  }

  // Get the differential and total cross section for the default
  // MEC model (including contributions from both kinds of
  // initial nucleon clusters and both kinds of diagrams)
  // Save the hit nucleon cluster PDG code (we'll need to set it again
  // for the alternate modelNievesSimoVacasMECPXSec2016)
  int hit_nuc_pdg = interaction->InitState().Tgt().HitNucPdg();
  interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( 0 );
  interaction->ExclTagPtr()->SetResonance( kNoResonance );

  double diff_xsec_def = fXSecAlgCCDef->XSec( interaction,
    kPSTlctl);
  double tot_xsec_def = event.XSec();
  double prob_density_def = diff_xsec_def / tot_xsec_def;

  //double check_tot_xsec_def = fXSecAlgCCDef->Integral( interaction );
  //LOG("ReW", pDEBUG) << "check_tot_xsec_def = " << check_tot_xsec_def;

  // Get the names of the CCMEC cross section default and alternate models
  std::string cc_def_alg_name = fXSecAlgCCDef->Id().Name();
  // std::string cc_alt_alg_name = fXSecAlgCCAlt_Nieves->Id().Name();
  // std::string cc_alt2_alg_name = fXSecAlgCCAlt_SuSAv2->Id().Name();
  // std::string cc_alt3_alg_name = fXSecAlgCCAlt_Empirical->Id().Name();
  std::string cc_alt4_alg_name = fXSecAlgCCAlt_Martini->Id().Name();

  // Set the hit nucleon cluster PDG code to its sampled value
  // (empirical MEC needs it)
  if ( cc_alt4_alg_name == "genie::EmpiricalMECPXSec2015" ) {
    interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( hit_nuc_pdg );
    std::cout << "alternate model is empirical" << std::endl;
  }

  //// The empirical MEC differential cross section doesn't check the
  //// kinematic limits imposed by MECGenerator, so enforce them here.
  //// If we're outside, then don't bother calculating the alternative
  //// model differential cross section.
  //double diff_xsec_alt = 0.;
  //Range1D_t rW = getWLimitsEmpiricalMEC();
  //Range1D_t rQ2 = getQ2LimitsEmpiricalMEC();

  //double W = interaction->Kine().W();
  //double Q2 = interaction->Kine().Q2();

  //LOG("RwMEC", pERROR) << "NIEVES: W = " << W << ", Q2 = " << Q2;
  //if ( rW.min <= W && rW.max >= W && rQ2.min <= Q2 && rQ2.max >= Q2 ) {

  // std::cout << "Input (default) CCMEC cross section model: " << cc_def_alg_name << std::endl;
  // std::cout << "Alternative CCMEC cross section model: " << cc_alt4_alg_name << std::endl;

  // Once the CCMEC Martini model is available, just change fXSecAlgCCAlt3 
  // (reweight from SuSAv2 or Valencia to the Empirical model) to 
  // fXSecAlgCCAlt4 (reweight from SuSAv2 or Valencia to the Martini model)
  // and this should just work
  double diff_xsec_alt = fXSecAlgCCAlt_Martini->XSec( interaction, kPSTlctl );

  //}

  double tot_xsec_alt = this->GetXSecIntegral( fXSecAlgCCAlt_Martini, interaction );

  //LOG("RwMEC", pERROR) << "diff_xsec_alt = " << diff_xsec_alt << ", tot_xsec_alt = " << tot_xsec_alt;

  //if ( tot_xsec_alt == 0. && diff_xsec_alt != 0. ) LOG("RwMEC", pERROR) << "OH NO!";

  // Protect against NaNs when the total cross section for the alternative
  // model is zero
  if ( tot_xsec_alt == 0. ) {
    diff_xsec_alt = 0.;
    tot_xsec_alt = 1.;
  }
  double prob_density_alt = diff_xsec_alt / tot_xsec_alt;

  // Compute a new probability density for this event by interpolating between
  // the two models while preserving the total cross section
  double tweaked_prob_density = (1. - twk_dial)*prob_density_def + twk_dial*prob_density_alt;

  // The weight is then the likelihood ratio
  double weight = tweaked_prob_density / prob_density_def;

  LOG("ReW", pDEBUG) << "xsec_def = " << diff_xsec_def << ", xsec_alt = " << diff_xsec_alt;
  LOG("ReW", pDEBUG) << "tot_xsec_def = " << tot_xsec_def << ", tot_xsec_alt = " << tot_xsec_alt;
  LOG("ReW", pDEBUG) << "twk_dial = " << fCCXSecShapeMartiniTwkDial << ", prob_density_def = "
    << prob_density_def << ", prob_density_alt = " << prob_density_alt << ", weight = " << weight;

  // We don't need the cloned interaction anymore, so delete it
  delete interaction;

  return weight;
}
//______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeightEnergyDependence(const genie::EventRecord& event)
{
  // Enforce tweak dial values on the interval [0, 1] for the energy
  // dependence dial, regardless of what the user requested.
  double twk_dial = std::max( std::min(1., fEnergyDependenceTwkDial), 0. );
  
  // If the tweak dial is set to zero (or close to zero) then just return a
  // weight of unity
  bool tweaked = ( std::abs(twk_dial) > controls::kASmallNum );
  if ( !tweaked ) return 1.;
    double weight = 1.;
    weight = CalcWeight2p2hEnergyDependence( event );
    return weight;
  }
//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeight2p2hEnergyDependence(const genie::EventRecord& event)
{
  double weight = 1.;

  // Get the neutrino
  GHepParticle* neutrino = event.Particle( 0 );
  double nu_pdg = neutrino->Pdg();

  // Get the neutrino energy
  double E_nu = neutrino->E();

  // Initialise cross section ratio
  double r = -999.;

  // Consider cross sections for Empirical, Valencia and SuSAv2 models in energy range [0.0,3.0].
  // Above 3 GeV the CC 2p-2h cross section predictions approaches more or less a constant value.

  // Check whether particle is neutrino ... (idea: can split energy dependence dial further into
  // neutrino and antineutrino cases, in addition to low and high energy split -> 4 cases)
  if( nu_pdg > 0 ){
    // If neutrino energy is in range, take ratio at that energy
    if( E_nu < 2.995 ){
      // r = 1.14523; // TODO: replace hard-coded with table value everywhere, TGraph
      r = ratioGraph_nu->Eval(E_nu);
      // std::cout << "E_nu is " << E_nu << " and r is " << ratioGraph_nu->Eval(E_nu) << std::endl;
    }
    else{ // else, take ratio at maximal neutrino energy
      // r = 1.3245;
      r = ratioGraph_nu->Eval(2.995);
      }
    if( E_nu < 0.245 ){ // if neutrino energy is below range, take ratio at minimal energy
      // r = 0.97832;
      r = ratioGraph_nu->Eval(0.245);
    }
    
    // Consider split into low and high energy dial (split is energy where r = 1); in this case change 
    // fEnergyDependenceTwkDial to fLowEnergyDependenceTwkDial ( for E_nu < 10.0 ) and 
    // fHighEnergyDependenceTwkDial (else)
    // if( E_nu < 1.2 ){
    if( E_nu < 10.0 ){
      // weight = fEnergyDependenceTwkDial + ( 1 - fEnergyDependenceTwkDial ) / r ; // For dial 1 being CV and dial 0 being tweaked
      weight = 1 - fEnergyDependenceTwkDial + fEnergyDependenceTwkDial / r ; // For dial 0 being CV and dial 1 being tweaked
    }
    else{
      // weight = fEnergyDependenceTwkDial + ( 1 - fEnergyDependenceTwkDial ) / r ;
      weight = 1 - fEnergyDependenceTwkDial + fEnergyDependenceTwkDial / r ;
    }
  // ... or antineutrino
  }
  else{
    if( E_nu < 2.995 ){ // if neutrino energy is in range, take ratio at that energy
      // r = 1.14523; // replace hard-coded with table value everywhere, TGraph
      r = ratioGraph_nu->Eval(E_nu);
    }
    else{ // else, take ratio at maximal neutrino energy
      // r = 1.3245;
      r = ratioGraph_nu->Eval(2.995);
    }
    if( E_nu < 0.245 ){ // if neutrino energy is below range, take ratio at minimal energy
      // r = 0.97832;
      r = ratioGraph_nu->Eval(0.245);
    }
    
    // Consider split into low and high energy dial (split is energy where r = 1); in this case change 
    // fEnergyDependenceTwkDial to fLowEnergyDependenceTwkDial ( for E_nu < 10.0 ) and 
    // fHighEnergyDependenceTwkDial (else)
    // if( E_nu < 1.2 ){
    if( E_nu < 10.0 ){
      // weight = fEnergyDependenceTwkDial + ( 1 - fEnergyDependenceTwkDial ) / r ;
      weight = 1 - fEnergyDependenceTwkDial + fEnergyDependenceTwkDial / r ;
    }
    else{
      // weight = fEnergyDependenceTwkDial + ( 1 - fEnergyDependenceTwkDial ) / r ;
      weight = 1 - fEnergyDependenceTwkDial + fEnergyDependenceTwkDial * r ;
    }
  }
  return weight;
}
//_______________________________________________________________________________________
double GReWeightXSecMEC::GetXSecIntegral(const XSecAlgorithmI* xsec_alg,
  const Interaction* interaction)
{
  double xsec = 0.;

  XSecSplineList* xssl = XSecSplineList::Instance();
  assert( xssl );

  // First check if a total cross section spline is already available
  // for the requested cross section model and interaction. If it is,
  // use it to get the integrated cross section.
  std::string curr_tune = xssl->CurrentTune();
  bool spline_computed = xssl->HasSplineFromTune( curr_tune )
    && xssl->SplineExists( xsec_alg, interaction );
  if ( spline_computed ) {
    const Spline* spl = xssl->GetSpline( xsec_alg, interaction );
    double Ev = interaction->InitState().ProbeE( kRfLab ); // kRfHitNucRest?
    if ( spl->ClosestKnotValueIsZero(Ev, "-") ) xsec = 0.;
    else xsec = spl->Evaluate( Ev );
  }
  // If not, fall back to doing the integration directly
  else {
    xsec = xsec_alg->Integral( interaction );
  }

  LOG("ReW", pDEBUG) << "MECshape spline check: "
    << ( spline_computed ? "found" : "not found" );

  return xsec;
}
