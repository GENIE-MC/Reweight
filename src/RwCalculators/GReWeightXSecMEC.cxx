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
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightXSecMEC.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

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

}
//_______________________________________________________________________________________
bool GReWeightXSecMEC::IsHandled(GSyst_t syst) const
{
  // Some MEC tweak dials are independent of interaction type. Check
  // whether we can handle these first.
  if ( syst == kXSecTwkDial_DecayAngMEC ) return true;
  if ( syst == kXSecTwkDial_FracPN_CCMEC ) return true;

  // If we have an entry for a knob that is interaction type dependent in the
  // GSyst_t -> InteractionType_t map, then this calculator can handle it.
  // Otherwise, it can't.
  bool handle = fGSystToIntTypeMap.count( syst );
  return handle;
}
//_______________________________________________________________________________________
bool GReWeightXSecMEC::AppliesTo(ScatteringType_t type, bool /*is_cc*/) const
{
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
    return;
  }
  else if ( syst == kXSecTwkDial_FracPN_CCMEC ) {
    fFracPN_CCTwkDial = twk_dial;
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

  fFracPN_CCTwkDial = 0.;

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
  weight *= this->CalcWeightAngularDist( event );
  weight *= this->CalcWeightPN( event );
  return weight;
}
//_______________________________________________________________________________________
void GReWeightXSecMEC::Init(void) {

  // Set the tweak dials to their default values
  fDecayAngTwkDial = 0.;
  fFracPN_CCTwkDial = 0.;

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
  double twk_dial = std::max( std::min(1., fDecayAngTwkDial), 0. );

  // If the tweak dial is set to zero (or is really small) then just return a
  // weight of unity
  bool tweaked = ( std::abs(twk_dial) > controls::kASmallNum );
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
  double weight = 3.*twk_dial*std::pow(std::cos(theta_N1), 2) + (1. - twk_dial);

  return weight;
}
//_______________________________________________________________________________________
double GReWeightXSecMEC::CalcWeightPN(const genie::EventRecord& event)
{
  // Only handle CC events for now (and return unit weight for the others)
  // TODO: Add capability to tweak nucleon pair isospin for NC, EM
  InteractionType_t type = event.Summary()->ProcInfo().InteractionTypeId();
  if ( type != kIntWeakCC ) return 1.;

  // If the tweak dial is set to zero (or is really small) then just return a
  // weight of unity
  bool tweaked = ( std::abs(fFracPN_CCTwkDial) > controls::kASmallNum );
  if ( !tweaked ) return 1.;

  // Enforce that the current event involves an initial nucleon cluster.
  // Determine whether the cluster is p+n
  GHepParticle* initial_nucleon_cluster = event.HitNucleon();
  assert( initial_nucleon_cluster );

  int two_nuc_pdg = initial_nucleon_cluster->Pdg();
  assert( pdg::Is2NucleonCluster(two_nuc_pdg) );

  bool is_pn_event = ( two_nuc_pdg == kPdgClusterNP );

  // Calculate the model's default fraction of initial pn pairs. For empirical
  // MEC, this is just a fixed number from the configuration. For Valencia, it's
  // dependent on the kinematics, so we'll need to compute ratios of
  // differential cross sections.
  double pn_frac_def = 0.;

  std::string cc_def_alg_name = fXSecAlgCCDef->Id().Name();
  if ( cc_def_alg_name == "genie::EmpiricalMECPXSec2015" ) {
    // For GENIE's empirical MEC model, the pn fraction is not dependent on
    // kinematics. We can just retrieve the value from the model configuration.
    pn_frac_def = fXSecAlgCCDef->GetConfig().GetDouble( "EmpiricalMEC-FracPN_CC" );
  }
  else if ( cc_def_alg_name == "genie::NievesSimoVacasMECPXSec2016" ) {
    // For the Valencia MEC model, the pn fraction can vary with q0 and q3. We
    // can get the pn fraction for this event's kinematics by computing the
    // differential cross section for each case. A similar thing is done in
    // MECGenerator in order to decide which kind of nucleon pair is hit. See
    // genie::MECGenerator::GenerateNSVInitialHadrons() for details.

    // Get the differential cross section for an initial pn pair and the total
    // for all pair types (this is how the Valencia calculation is organized in
    // GENIE). Note that the that the Valencia MEC model works in the kPSTlctl
    // phase space. Clone the input interaction so that we can modify
    // the PDG code of the initial nucleon cluster.
    Interaction* interaction = new Interaction( *event.Summary() );

    // TODO: When NC and/or EM interactions are added for Valencia, generalize
    // this for use those. Unlike CC, all three pair types can participate in
    // NC & EM.

    // Get the differential cross section for an initial pn pair
    interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( kPdgClusterNP );
    double xsec_pn = fXSecAlgCCDef->XSec( interaction, kPSTlctl );

    // Get the total differential cross section
    interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg( 0 );
    double xsec_tot = fXSecAlgCCDef->XSec( interaction, kPSTlctl );

    // We don't need the cloned interaction anymore, so delete it
    delete interaction;

    assert( xsec_tot > 0. );
    pn_frac_def = xsec_pn / xsec_tot;
  }
  else if ( cc_def_alg_name == "genie::SuSAv2MECPXSec" ) {
    // TODO: actually implement this
    LOG("ReW", pWARN) << "MEC pn reweighting for SuSAv2 not yet implemented";
    return 1.;
  }
  else {
    LOG("ReW", pERROR) << "Unrecognized MEC model " << cc_def_alg_name
      << " encountered in genie::rew::GReWeightXSecMEC::CalcWeightPN()";
    return 1.;
  }

  // Check that the pn fraction computed above is sane. If not, complain and
  // return a unit weight.
  bool impossible_pp_or_nn_event = ( pn_frac_def == 1. && !is_pn_event );
  if ( pn_frac_def <= 0. || pn_frac_def > 1. || impossible_pp_or_nn_event ) {
    LOG("ReW", pERROR) << "Invalid pn fraction value " << pn_frac_def
      << " encountered in genie::rew::GReWeightXSecMEC::CalcWeightPN()";
    return 1.;
  }

  // TODO: add support for asymmetric errors here
  GSystUncertainty* gsu = GSystUncertainty::Instance();
  double frac_err_pn_cc = gsu->OneSigmaErr( kXSecTwkDial_FracPN_CCMEC );

  // Compute the scaling factor for the pn fraction that corresponds to the
  // current tweak dial setting
  double pn_tweak_factor = ( 1. + fFracPN_CCTwkDial * frac_err_pn_cc );

  // To conserve the total cross section (which is separately controlled by the
  // normalization tweak dials), enforce that the tweaked pn fraction lies on
  // the interval [0, 1]
  double pn_frac_tweak = std::max( std::min(1., pn_frac_def * pn_tweak_factor), 0. );

  // Assign the appropriate likelihood ratio as the weight. Note that
  // we've already checked that pn_frac_def lies in a reasonable range
  // above, so we can divide as shown without worrying about NaNs.
  double weight;
  if ( is_pn_event ) weight = pn_frac_tweak / pn_frac_def;
  else weight = ( 1. - pn_frac_tweak ) / ( 1. - pn_frac_def );

  LOG("ReW", pERROR) << "twk_dial = " << fFracPN_CCTwkDial << ", frac_err = "
    << frac_err_pn_cc;
  LOG("ReW", pERROR) << "pn_frac_def = " << pn_frac_def << ", pn_frac_tweak = "
    << pn_frac_tweak << ", weight = " << weight;

  return weight;
}
