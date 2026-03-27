// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file hfFragmentationFunction.cxx
/// \brief charm hadron hadronization task
/// \author Christian Reckziegel <christian.reckziegel@cern.ch>, Federal University of ABC
/// \since 15.03.2024
///
/// The task store data relevant to the calculation of hadronization observables radial
/// profile and/or jet momentum fraction for charmed hadrons
#include <vector>
#include <string>

#include "TVector3.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
//#include "PWGHF/TableProducer/candidateSelectorD0.cxx"
//#include "EventFiltering/PWGHF/HFFilterHelpers.h" // where findBin() is defined

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetUtilities.h"


#include <typeinfo> // Include this header for typeid

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// calculate delta phi such that 0 < delta phi < pi
double deltaPhi(double phi1, double phi2)
{
  // Compute the absolute difference between phi1 and phi2
  double dphi = std::abs(phi1 - phi2);

  // Constrain angle between [min,min+2pi] = [-pi,-pi+2pi]
  dphi = RecoDecay::constrainAngle(dphi, -o2::constants::math::PI);

  // Return absolute value of distance
  return std::abs(dphi);
}

struct CutValues {
  double mass;
  double dca;
  double cosThetaStar;
  double ptK;
  double ptPi;
  double d0K;
  double d0pi;
  double d0d0;
  double cpa;
  double cpaXY;
  double minDecayLengthXY;
  double maxDecayLength;
  double maxDecayLengthXY;
  double minDecayLength;
  double normDauImpParXY;
};

int findBinPt(const std::vector<double>& bins, double value) {
  return std::upper_bound(bins.begin(), bins.end(), value) - bins.begin() - 1;
}

template <typename Candidate>
bool preSelection(const Candidate& candidate, const std::vector<double>& binsPt, const std::vector<CutValues>& cuts) {
  auto pTBin = findBinPt(binsPt, candidate.pt());
  const auto& cut = cuts[pTBin];

  return std::abs(candidate.m() - o2::constants::physics::MassD0) <= cut.mass &&
         candidate.cosThetaStar() >= cut.cosThetaStar &&
         candidate.ptProng0() >= cut.ptPi &&
         candidate.ptProng1() >= cut.ptK &&
         candidate.impactParameterProduct() <= cut.d0d0 &&
         candidate.cpa() >= cut.cpa &&
         candidate.cpaXY() >= cut.cpaXY &&
         candidate.decayLengthXYNormalised() >= cut.minDecayLengthXY &&
         candidate.decayLength() >= cut.minDecayLength &&
         candidate.decayLength() <= cut.maxDecayLength &&
         candidate.decayLengthXY() <= cut.maxDecayLengthXY &&
         std::abs(candidate.impactParameterNormalised0()) >= cut.normDauImpParXY &&
         std::abs(candidate.impactParameterNormalised1()) >= cut.normDauImpParXY;
}

// check if D0 passes pre-selection cuts
/*template <typename T1, typename T2, typename T3>
bool preSelection(const T1& candidate, T2 const& binsPt, T3 const& cuts) {
  auto candpT = candidate.pt();
  auto pTBin = findBin(binsPt, candpT);
  
  // candidate invariant mass
  if (std::abs(candidate.m() - o2::constants::physics::MassD0) > cuts->get(pTBin, "m")) {
    return false;
  }  
  // candidate DCA
  //if (std::abs(candidate.impactParameterXY()) > cuts->get(pTBin, "DCA")) {
  //  return false;
  //}
  // "cos theta*", 
  if (candidate.cosThetaStar() < cuts->get(pTBin, "cos theta*")) {
    return false;
  }
  // "pT K", 
  // "pT Pi", 
  if (candidate.ptProng0() < cuts->get(pTBin, "pT Pi") || candidate.ptProng1() < cuts->get(pTBin, "pT K")) {
    return false;
  }
  // cut on daughter DCA - need to add secondary vertex constraint here
  //if (std::abs(trackPion.dcaXY()) > cuts->get(pTBin, "d0pi") || std::abs(trackKaon.dcaXY()) > cuts->get(pTBin, "d0K")) {
  //  return false;
  //}
  // product of daughter impact parameters
  if (candidate.impactParameterProduct() > cuts->get(pTBin, "d0d0")) {
    return false;
  }
  // cosine of pointing angle
  if (candidate.cpa() < cuts->get(pTBin, "cos pointing angle")) {
    return false;
  }
  // cosine of pointing angle XY
  if (candidate.cpaXY() < cuts->get(pTBin, "cos pointing angle xy")) {
    return false;
  }
  // normalised decay length in XY plane
  if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "min norm decay length XY")) {
    return false;
  }
  // "max decay length", 
  if (candidate.decayLength() > cuts->get(pTBin, "max decay length")) {
    return false;
  }  
  // "max decay length XY", 
  if (candidate.decayLengthXY() > cuts->get(pTBin, "max decay length XY")) {
    return false;
  }
  // "min decay length", 
  if (candidate.decayLength() < cuts->get(pTBin, "min decay length")) {
    return false;
  }
  // impact parameter of prong 0/1 divided by its uncertainty
  if (std::abs(candidate.impactParameterNormalised0()) < cuts->get(pTBin, "norm dauImpPar XY") || std::abs(candidate.impactParameterNormalised1()) < cuts->get(pTBin, "norm dauImpPar XY")) {
    return false;
  }

  return true;
}*/

// creating table for storing distance data
namespace o2::aod
{
namespace jet_distance
{
DECLARE_SOA_COLUMN(JetHfDist, jetHfDist, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, int);
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);
DECLARE_SOA_COLUMN(HfPrompt, hfPrompt, bool);
DECLARE_SOA_COLUMN(HfMatch, hfMatch, bool);
DECLARE_SOA_COLUMN(HfMlScore0, hfMlScore0, float);
DECLARE_SOA_COLUMN(HfMlScore1, hfMlScore1, float);
DECLARE_SOA_COLUMN(HfMlScore2, hfMlScore2, float);
DECLARE_SOA_COLUMN(HfMatchedFrom, hfMatchedFrom, int);
DECLARE_SOA_COLUMN(HfSelectedAs, hfSelectedAs, int);
DECLARE_SOA_COLUMN(McJetHfDist, mcJetHfDist, float);
DECLARE_SOA_COLUMN(McJetPt, mcJetPt, float);
DECLARE_SOA_COLUMN(McJetEta, mcJetEta, float);
DECLARE_SOA_COLUMN(McJetPhi, mcJetPhi, float);
DECLARE_SOA_COLUMN(McJetNConst, mcJetNConst, int);
DECLARE_SOA_COLUMN(McHfPt, mcHfPt, float);
DECLARE_SOA_COLUMN(McHfEta, mcHfEta, float);
DECLARE_SOA_COLUMN(McHfPhi, mcHfPhi, float);
DECLARE_SOA_COLUMN(McHfY, mcHfY, float);
DECLARE_SOA_COLUMN(McHfPrompt, mcHfPrompt, bool);
DECLARE_SOA_COLUMN(McHfMatch, mcHfMatch, bool);
} // namespace jet_distance
DECLARE_SOA_TABLE(JetDistanceTable, "AOD", "JETDISTTABLE",
                  jet_distance::JetHfDist,
                  jet_distance::JetPt,
                  jet_distance::JetEta,
                  jet_distance::JetPhi,
                  jet_distance::JetNConst,
                  jet_distance::HfPt,
                  jet_distance::HfEta,
                  jet_distance::HfPhi,
                  jet_distance::HfMass,
                  jet_distance::HfY,
                  jet_distance::HfMlScore0,
                  jet_distance::HfMlScore1,
                  jet_distance::HfMlScore2);
DECLARE_SOA_TABLE(MCPJetDistanceTable, "AOD", "MCPJETDISTTABLE",
                  jet_distance::McJetHfDist,
                  jet_distance::McJetPt,
                  jet_distance::McJetEta,
                  jet_distance::McJetPhi,
                  jet_distance::McJetNConst,
                  jet_distance::McHfPt,
                  jet_distance::McHfEta,
                  jet_distance::McHfPhi,
                  jet_distance::McHfY,
                  jet_distance::McHfPrompt,
                  jet_distance::McHfMatch);
DECLARE_SOA_TABLE(MCDJetDistanceTable, "AOD", "MCDJETDISTTABLE",
                  jet_distance::JetHfDist,
                  jet_distance::JetPt,
                  jet_distance::JetEta,
                  jet_distance::JetPhi,
                  jet_distance::JetNConst,
                  jet_distance::HfPt,
                  jet_distance::HfEta,
                  jet_distance::HfPhi,
                  jet_distance::HfMass,
                  jet_distance::HfY,
                  jet_distance::HfPrompt,
                  jet_distance::HfMatch,
                  jet_distance::HfMlScore0,
                  jet_distance::HfMlScore1,
                  jet_distance::HfMlScore2,
                  jet_distance::HfMatchedFrom,
                  jet_distance::HfSelectedAs);
DECLARE_SOA_TABLE(MatchJetDistanceTable, "AOD", "MATCHTABLE",
                  jet_distance::McJetHfDist,
                  jet_distance::McJetPt,
                  jet_distance::McJetEta,
                  jet_distance::McJetPhi,
                  jet_distance::McJetNConst,
                  jet_distance::McHfPt,
                  jet_distance::McHfEta,
                  jet_distance::McHfPhi,
                  jet_distance::McHfY,
                  jet_distance::McHfPrompt,
                  jet_distance::JetHfDist,
                  jet_distance::JetPt,
                  jet_distance::JetEta,
                  jet_distance::JetPhi,
                  jet_distance::JetNConst,
                  jet_distance::HfPt,
                  jet_distance::HfEta,
                  jet_distance::HfPhi,
                  jet_distance::HfMass,
                  jet_distance::HfY,
                  jet_distance::HfPrompt,
                  jet_distance::HfMlScore0,
                  jet_distance::HfMlScore1,
                  jet_distance::HfMlScore2,
                  jet_distance::HfMatchedFrom,
                  jet_distance::HfSelectedAs);
} // namespace o2::aod

struct HfFragmentationFunction {
  // producing new table
  Produces<aod::JetDistanceTable> distJetTable;
  Produces<aod::MCPJetDistanceTable> mcpdistJetTable;
  Produces<aod::MCDJetDistanceTable> mcddistJetTable;
  Produces<aod::MatchJetDistanceTable> matchJetTable;

  // Tables for MC jet matching
  using JetMCDTable = soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>;
  using JetMCPTable = soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>;

  // slices for accessing proper HF mcdjets collision associated to mccollisions
  Preslice<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;
  Preslice<JetMCDTable> d0MCDJetsPerCollisionPreslice = aod::jet::collisionId;
  Preslice<JetMCPTable> d0MCPJetsPerMCCollisionPreslice = aod::jet::mcCollisionId;

  // Histogram registry: an object to hold your histograms
  HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  
  // pT array ranges
  // default values for the cuts
  /*static constexpr int NBinsPt = 25; // static constexpr int
  static constexpr int NCutVars = 15; // static constexpr int
  // default values for the pT bin edges (can be used to configure histogram axis)
  // offset by 1 from the bin numbers in cuts array
  double BinsPt[NBinsPt + 1] = {
    0,
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0,
    5.5,
    6.0,
    6.5,
    7.0,
    7.5,
    8.0,
    9.0,
    10.0,
    12.0,
    16.0,
    20.0,
    24.0,
    36.0,
    50.0,
    100.0};// constexpr double
  std::vector<double> vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};*/ // auto 
  std::vector<double> vecBinsPt = {0,
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0,
    5.5,
    6.0,
    6.5,
    7.0,
    7.5,
    8.0,
    9.0,
    10.0,
    12.0,
    16.0,
    20.0,
    24.0,
    36.0,
    50.0,
    100.0};
  
  /*double Cuts[NBinsPt][NCutVars] = {{0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},   // 0   < pT < 0.5
                                              {0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},   // 0.5 < pT < 1   
                                              {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},  // 1   < pT < 1.5 
                                              {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},  // 1.5 < pT < 2   
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0., 10., 10., 0.06, 0.5},  // 2   < pT < 2.5 
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0., 10., 10., 0.06, 0.5},  // 2.5 < pT < 3   
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  // 3   < pT < 3.5 
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  // 3.5 < pT < 4   
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 4   < pT < 4.5 
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 4.5 < pT < 5   
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 5   < pT < 5.5 
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 5.5 < pT < 6   
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 6   < pT < 6.5 
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 6.5 < pT < 7   
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 7   < pT < 7.5 
                                              {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 7.5 < pT < 8   
                                              {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 8   < pT < 9   
                                              {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 9   < pT < 10  
                                              {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 10  < pT < 12  
                                              {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 10000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   // 12  < pT < 16  
                                              {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  // 16  < pT < 20  
                                              {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  // 20  < pT < 24  
                                              {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  // 24  < pT < 36  
                                              {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  // 36  < pT < 50  
                                              {0.400, 300. * 1E-4, 1.0, 0.6, 0.6, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5}};*/ // 50  < pT < 100, constexpr double

  std::vector<CutValues> cutValues = {
    {0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},   /* 0   < pT < 0.5 */
    {0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},   /* 0.5 < pT < 1   */
    {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},  /* 1   < pT < 1.5 */
    {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5},  /* 1.5 < pT < 2   */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0., 10., 10., 0.06, 0.5},  /* 2   < pT < 2.5 */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0., 10., 10., 0.06, 0.5},  /* 2.5 < pT < 3   */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 3   < pT < 3.5 */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 3.5 < pT < 4   */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 4   < pT < 4.5 */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 4.5 < pT < 5   */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 5   < pT < 5.5 */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 5.5 < pT < 6   */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 6   < pT < 6.5 */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 6.5 < pT < 7   */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 7   < pT < 7.5 */
    {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 7.5 < pT < 8   */
    {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 8   < pT < 9   */
    {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 9   < pT < 10  */
    {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 10  < pT < 12  */
    {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 10000. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},   /* 12  < pT < 16  */
    {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 16  < pT < 20  */
    {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 20  < pT < 24  */
    {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 24  < pT < 36  */
    {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0., 10., 10., 0.06, 0.5},  /* 36  < pT < 50  */
    {0.400, 300. * 1E-4, 1.0, 0.6, 0.6, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.80, 0., 0., 10., 10., 0.06, 0.5}   /* 50  < pT < 100 */
  };
  
  // row labels
  /*std::vector<std::string> labelsPt = {
    "pT bin 0",
    "pT bin 1",
    "pT bin 2",
    "pT bin 3",
    "pT bin 4",
    "pT bin 5",
    "pT bin 6",
    "pT bin 7",
    "pT bin 8",
    "pT bin 9",
    "pT bin 10",
    "pT bin 11",
    "pT bin 12",
    "pT bin 13",
    "pT bin 14",
    "pT bin 15",
    "pT bin 16",
    "pT bin 17",
    "pT bin 18",
    "pT bin 19",
    "pT bin 20",
    "pT bin 21",
    "pT bin 22",
    "pT bin 23",
    "pT bin 24"}; // static const

  // column labels
  std::vector<std::string> labelsCutVar = {"m", "DCA", "cos theta*", "pT K", "pT Pi", "d0K", "d0pi", "d0d0", "cos pointing angle", "cos pointing angle xy", "min norm decay length XY", "max decay length", "max decay length XY", "min decay length", "norm dauImpPar XY"};*/ // static const

  // pre-selection HF cuts
  //Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{vecBinsPt}, "pT bin limits"};
  //Configurable<LabeledArray<double>> cuts{"cuts", {Cuts[0], NBinsPt, NCutVars, labelsPt, labelsCutVar}, "D0 candidate selection per pT bin"};
  // event selection cuts
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};

  std::vector<int> eventSelectionBits;

  void init(InitContext const&)
  {
    // initialise event selection:
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));

    // create histograms
    // collision system histograms
    std::vector<std::string> histLabels = {"mccollisions", "mc z_cut", "collisions", "sel8 + z_cut"};
    registry.add("h_collision_counter", ";# of collisions;", HistType::kTH1F, {{static_cast<int>(histLabels.size()), 0.0, static_cast<double>(histLabels.size())}});
    auto counter = registry.get<TH1>(HIST("h_collision_counter"));
    for (std::vector<std::string>::size_type iCounter = 0; iCounter < histLabels.size(); iCounter++) {
      counter->GetXaxis()->SetBinLabel(iCounter + 1, histLabels[iCounter].data());
    }
    registry.add("h_jet_counter", ";# of jets;", {HistType::kTH1F, {{6, 0., 3.0}}});
    auto jetCounter = registry.get<TH1>(HIST("h_jet_counter"));
    jetCounter->GetXaxis()->SetBinLabel(1, "particle level");
    jetCounter->GetXaxis()->SetBinLabel(2, "detector level");
    jetCounter->GetXaxis()->SetBinLabel(3, "particle matched jets");
    jetCounter->GetXaxis()->SetBinLabel(4, "detector matched jets");
    jetCounter->GetXaxis()->SetBinLabel(5, "mcd matched to mcp loop");
    jetCounter->GetXaxis()->SetBinLabel(6, "mcp matched to mcd loop");
    registry.add("h_d0_preselection",";# of candidates;", {HistType::kTH1F, {{3, 0., 3.}}});
    auto hfPreSelectionCounter = registry.get<TH1>(HIST("h_jet_counter"));
    hfPreSelectionCounter->GetXaxis()->SetBinLabel(1, "Free of cuts");
    hfPreSelectionCounter->GetXaxis()->SetBinLabel(2, "Pre-selected");
    hfPreSelectionCounter->GetXaxis()->SetBinLabel(3, "Failed pre-selection");
    // D0 candidate histograms from data
    registry.add("h_d0_jet_projection", ";z^{D^{0},jet}_{||};dN/dz^{D^{0},jet}_{||}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_d0_jet_distance_vs_projection", ";#DeltaR_{D^{0},jet};z^{D^{0},jet}_{||}", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
    registry.add("h_d0_jet_distance", ";#DeltaR_{D^{0},jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_d0_jet_pt", ";p_{T,D^{0} jet};dN/dp_{T,D^{0} jet}", {HistType::kTH1F, {{200, 0., 10.}}});
    registry.add("h_d0_jet_eta", ";#eta_{T,D^{0} jet};dN/d#eta_{D^{0} jet}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add("h_d0_jet_phi", ";#phi_{T,D^{0} jet};dN/d#phi_{D^{0} jet}", {HistType::kTH1F, {{250, -10., 10.}}});
    registry.add("h_d0_mass", ";m_{D^{0}} (GeV/c^{2});dN/dm_{D^{0}}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_d0_eta", ";#eta_{D^{0}} (GeV/c^{2});dN/d#eta_{D^{0}}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add("h_d0_phi", ";#phi_{D^{0}} (GeV/c^{2});dN/d#phi_{D^{0}}", {HistType::kTH1F, {{250, -10., 10.}}});
  }

  void processDummy(aod::TracksIU const&) {}
  PROCESS_SWITCH(HfFragmentationFunction, processDummy, "Dummy process function turned on by default", true);

  void processDataChargedSubstructure(aod::JetCollision const& collision,
                                      soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& jets,
                                      aod::CandidatesD0Data const&,
                                      aod::JetTracks const&)
  {
    // apply event selection and fill histograms for sanity check
    registry.fill(HIST("h_collision_counter"), 2.0);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST("h_collision_counter"), 3.0);

    for (const auto& jet : jets) {
      // fill jet counter histogram
      registry.fill(HIST("h_jet_counter"), 0.5);
      // obtaining jet 3-vector
      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      for (const auto& d0Candidate : jet.candidates_as<aod::CandidatesD0Data>()) {

        // obtaining jet 3-vector
        TVector3 d0Vector(d0Candidate.px(), d0Candidate.py(), d0Candidate.pz());

        // calculating fraction of the jet momentum carried by the D0 along the direction of the jet axis
        double zParallel = (jetVector * d0Vector) / (jetVector * jetVector);

        // calculating angular distance in eta-phi plane
        double axisDistance = jetutilities::deltaR(jet, d0Candidate);

        // filling histograms
        registry.fill(HIST("h_d0_jet_projection"), zParallel);
        registry.fill(HIST("h_d0_jet_distance_vs_projection"), axisDistance, zParallel);
        registry.fill(HIST("h_d0_jet_distance"), axisDistance);
        registry.fill(HIST("h_d0_jet_pt"), jet.pt());
        registry.fill(HIST("h_d0_jet_eta"), jet.eta());
        registry.fill(HIST("h_d0_jet_phi"), jet.phi());
        registry.fill(HIST("h_d0_mass"), d0Candidate.m());
        registry.fill(HIST("h_d0_eta"), d0Candidate.eta());
        registry.fill(HIST("h_d0_phi"), d0Candidate.phi());

        registry.fill(HIST("h_d0_preselection"), 0.0);

        // apply pre-selection cuts on D0 candidates
        if (preSelection(d0Candidate, vecBinsPt, cutValues)) {
          registry.fill(HIST("h_d0_preselection"), 1.0);
        } else {
          registry.fill(HIST("h_d0_preselection"), 2.0);;
        }
        // filling table
        distJetTable(axisDistance,
                     jet.pt(), jet.eta(), jet.phi(), jet.tracks_as<aod::JetTracks>().size(),
                     d0Candidate.pt(), d0Candidate.eta(), d0Candidate.phi(), d0Candidate.m(), d0Candidate.y(), d0Candidate.mlScores()[0], d0Candidate.mlScores()[1], d0Candidate.mlScores()[2]);

        break; // get out of candidates' loop after first HF particle is found in jet
      } // end of D0 candidates loop

    } // end of jets loop

  } // end of process function
  PROCESS_SWITCH(HfFragmentationFunction, processDataChargedSubstructure, "charged HF jet substructure", false);

  void processMcEfficiency(aod::JetMcCollisions const& mccollisions,
                           aod::JetCollisionsMCD const& collisions,
                           JetMCDTable const& mcdjets,
                           JetMCPTable const& mcpjets,
                           aod::CandidatesD0MCD const&,
                           aod::CandidatesD0MCP const&,
                           aod::JetTracks const&,
                           aod::JetParticles const&)
  {
    for (const auto& mccollision : mccollisions) {

      registry.fill(HIST("h_collision_counter"), 0.0);
      // skip collisions outside of |z| < vertexZCut
      if (std::abs(mccollision.posZ()) > vertexZCut) {
        continue;
      }
      registry.fill(HIST("h_collision_counter"), 1.0);

      // reconstructed collisions associated to same mccollision
      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {

        registry.fill(HIST("h_collision_counter"), 2.0);
        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
          continue;
        }
        registry.fill(HIST("h_collision_counter"), 3.0);

        // d0 detector level jets associated to the current same collision
        const auto d0mcdJetsPerCollision = mcdjets.sliceBy(d0MCDJetsPerCollisionPreslice, collision.globalIndex());
        for (const auto& mcdjet : d0mcdJetsPerCollision) {

          registry.fill(HIST("h_jet_counter"), 0.5);

          // obtain leading HF candidate in jet
          auto mcdd0cand = mcdjet.candidates_first_as<aod::CandidatesD0MCD>();

          if (mcdjet.has_matchedJetCand()) {
            registry.fill(HIST("h_jet_counter"), 1.5);
          }

          // reflection information for storage: D0 = +1, D0bar = -1, neither = 0
          int matchedFrom = 0;
          int decayChannel = 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;
          int selectedAs = 0;

          if (mcdd0cand.flagMcMatchRec() == decayChannel) { // matched to D0 on truth level
            matchedFrom = 1;
          } else if (mcdd0cand.flagMcMatchRec() == -decayChannel) { // matched to D0bar on truth level
            matchedFrom = -1;
          }
          // bitwise AND operation: Checks whether BIT(i) is set, regardless of other bits
          if (mcdd0cand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as D0
            selectedAs = 1;
          } else if (mcdd0cand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as D0bar
            selectedAs = -1;
          }

          // store data in MC detector level table
          mcddistJetTable(jetutilities::deltaR(mcdjet, mcdd0cand),
                          mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), mcdjet.tracks_as<aod::JetTracks>().size(),                                                         // detector level jet
                          mcdd0cand.pt(), mcdd0cand.eta(), mcdd0cand.phi(), mcdd0cand.m(), mcdd0cand.y(), (mcdd0cand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level D0 candidate
                          mcdjet.has_matchedJetCand(), mcdd0cand.mlScores()[0], mcdd0cand.mlScores()[1], mcdd0cand.mlScores()[2],                                     // // Machine Learning PID scores: background, prompt, non-prompt
                          matchedFrom, selectedAs);                                                                                                                   // D0 = +1, D0bar = -1, neither = 0
        }
      }

      // d0 particle level jets associated to same mccollision
      const auto d0mcpJetsPerMCCollision = mcpjets.sliceBy(d0MCPJetsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& mcpjet : d0mcpJetsPerMCCollision) {

        registry.fill(HIST("h_jet_counter"), 0.0);

        // obtain leading HF particle in jet
        auto mcpd0cand = mcpjet.candidates_first_as<aod::CandidatesD0MCP>();

        if (mcpjet.has_matchedJetCand()) {
          registry.fill(HIST("h_jet_counter"), 1.0);
        }

        // store data in MC detector level table (calculate angular distance in eta-phi plane on the fly)
        mcpdistJetTable(jetutilities::deltaR(mcpjet, mcpd0cand),
                        mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.tracks_as<aod::JetParticles>().size(),                                       // particle level jet
                        mcpd0cand.pt(), mcpd0cand.eta(), mcpd0cand.phi(), mcpd0cand.y(), (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt), // particle level D0
                        mcpjet.has_matchedJetCand());
      }
    }
  }
  PROCESS_SWITCH(HfFragmentationFunction, processMcEfficiency, "non-matched and matched MC HF and jets", false);

  void processMcChargedMatched(aod::JetMcCollisions const& mccollisions,
                               aod::JetCollisionsMCD const& collisions,
                               JetMCDTable const& mcdjets,
                               JetMCPTable const&,
                               aod::CandidatesD0MCD const&,
                               aod::CandidatesD0MCP const&,
                               aod::JetTracks const&,
                               aod::JetParticles const&)
  {
    for (const auto& mccollision : mccollisions) {

      registry.fill(HIST("h_collision_counter"), 0.0);

      // skip collisions outside of |z| < vertexZCut
      if (std::abs(mccollision.posZ()) > vertexZCut) {
        continue;
      }
      registry.fill(HIST("h_collision_counter"), 1.0);

      // reconstructed collisions associated to same mccollision
      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {

        registry.fill(HIST("h_collision_counter"), 2.0);
        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
          continue;
        }
        registry.fill(HIST("h_collision_counter"), 3.0);
        // d0 detector level jets associated to the current same collision
        const auto d0mcdJetsPerCollision = mcdjets.sliceBy(d0MCDJetsPerCollisionPreslice, collision.globalIndex());
        for (const auto& mcdjet : d0mcdJetsPerCollision) {

          registry.fill(HIST("h_jet_counter"), 0.5);

          // comparison with fill on bin on 2.5 for sanity check
          if (mcdjet.has_matchedJetCand()) {
            registry.fill(HIST("h_jet_counter"), 1.5);
          }

          // obtain leading HF candidate in jet
          auto mcdd0cand = mcdjet.candidates_first_as<aod::CandidatesD0MCD>();

          // reflection information for storage: D0 = +1, D0bar = -1, neither = 0
          int matchedFrom = 0;
          int decayChannel = 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;
          int selectedAs = 0;

          if (mcdd0cand.flagMcMatchRec() == decayChannel) { // matched to D0 on truth level
            matchedFrom = 1;
          } else if (mcdd0cand.flagMcMatchRec() == -decayChannel) { // matched to D0bar on truth level
            matchedFrom = -1;
          }
          // bitwise AND operation: Checks whether BIT(i) is set, regardless of other bits
          if (mcdd0cand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as D0
            selectedAs = 1;
          } else if (mcdd0cand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as D0bar
            selectedAs = -1;
          }

          // loop through detector level matched to current particle level
          for (const auto& mcpjet : mcdjet.matchedJetCand_as<JetMCPTable>()) {

            registry.fill(HIST("h_jet_counter"), 2.5);

            // obtain leading HF candidate in jet
            auto mcpd0cand = mcpjet.candidates_first_as<aod::CandidatesD0MCP>();

            // store matched particle and detector level data in one single table (calculate angular distance in eta-phi plane on the fly)
            matchJetTable(jetutilities::deltaR(mcpjet, mcpd0cand), mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.tracks_as<aod::JetParticles>().size(),             // particle level jet
                          mcpd0cand.pt(), mcpd0cand.eta(), mcpd0cand.phi(), mcpd0cand.y(), (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt),                // particle level D0
                          jetutilities::deltaR(mcdjet, mcdd0cand), mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), mcdjet.tracks_as<aod::JetTracks>().size(),                // detector level jet
                          mcdd0cand.pt(), mcdd0cand.eta(), mcdd0cand.phi(), mcdd0cand.m(), mcdd0cand.y(), (mcdd0cand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level D0
                          mcdd0cand.mlScores()[0], mcdd0cand.mlScores()[1], mcdd0cand.mlScores()[2],                                                                  // Machine Learning PID scores: background, prompt, non-prompt
                          matchedFrom, selectedAs);                                                                                                                   // D0 = +1, D0bar = -1, neither = 0
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfFragmentationFunction, processMcChargedMatched, "matched MC HF and jets", false);

  /*void processMcFullChargedMatched(aod::JetMcCollision const& mccollision,
                               aod::JetCollisionsMCD const&,
                               JetMCDTable const&,
                               JetMCPTable const& mcpjets,
                               aod::CandidatesD0MCD const&,
                               aod::CandidatesD0MCP const&,
                               aod::JetTracks const&,
                               aod::JetParticles const&)
  {
    // 1st option
    // loop over all particle level jets
    // apply selections to both particle and detector level collisions
    // check if jet has matched detector level jet
    // if yes, store data in one table
    // if no, store data in the same table with -1 entries for the detector level corresponding space

    // 2nd option
    // loop over all mccollisions
    // apply selections to particle level collision
    // loop over all particle level jets
    // check if jet has matched detector level jet
    // apply selections to detector level collision
    // if yes, store data in one table
    // if no, store data in the same table with -1 entries for the detector level corresponding space

    // 3rd option
    // under the corresponding mccollision subscription, loop over all particle level jets
    // apply selections to particle level collision
    // loop over all particle level jets
    // check if jet has matched detector level jet
    // apply selections to detector level collision
    // if yes, store data in one table
    // if no, store data in the same table with -1 entries for the detector level corresponding space

    registry.fill(HIST("h_collision_counter"), 0.0);

    // skip collisions outside of |z| < vertexZCut
    if (std::abs(mccollision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("h_collision_counter"), 1.0);
    
    for (const auto& mcpjet : mcpjets) {

      registry.fill(HIST("h_jet_counter"), 0.0);

      // obtain leading HF candidate in jet
      auto mcpd0cand = mcpjet.candidates_first_as<aod::CandidatesD0MCP>();

      if (mcpjet.has_matchedJetCand()) {
        registry.fill(HIST("h_jet_counter"), 1.0);

        for (const auto& mcdjet : mcpjet.matchedJetCand_as<JetMCDTable>()) {

          registry.fill(HIST("h_jet_counter"), 2.0);

          // apply collision selection
          auto const& collision = mcdjet.collision();
          std::cout << "Type of collision: " << typeid(collision).name() << std::endl;
          std::cout << "collision from mcdjet posZ: " << collision.posZ() << std::endl;
          //auto collision = mcdjet.tracks_as<aod::JetTracks>()[0];
          registry.fill(HIST("h_collision_counter"), 2.0);
          //if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
          //  continue;
          //}
          registry.fill(HIST("h_collision_counter"), 3.0);

          // obtain leading HF candidate in jet
          auto mcdd0cand = mcdjet.candidates_first_as<aod::CandidatesD0MCD>();

          // reflection information for storage: D0 = +1, D0bar = -1, neither = 0
          int matchedFrom = 0;
          int decayChannel = 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;
          int selectedAs = 0;

          if (mcdd0cand.flagMcMatchRec() == decayChannel) { // matched to D0 on truth level
            matchedFrom = 1;
          } else if (mcdd0cand.flagMcMatchRec() == -decayChannel) { // matched to D0bar on truth level
            matchedFrom = -1;
          }
          // bitwise AND operation: Checks whether BIT(i) is set, regardless of other bits
          if (mcdd0cand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as D0
            selectedAs = 1;
          } else if (mcdd0cand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as D0bar
            selectedAs = -1;
          }

          // store matched particle and detector level data in one single table (calculate angular distance in eta-phi plane on the fly)
          matchJetTable(jetutilities::deltaR(mcpjet, mcpd0cand), mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.tracks_as<aod::JetParticles>().size(),             // particle level jet
                        mcpd0cand.pt(), mcpd0cand.eta(), mcpd0cand.phi(), mcpd0cand.y(), (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt),                // particle level D0
                        jetutilities::deltaR(mcdjet, mcdd0cand), mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), mcdjet.tracks_as<aod::JetTracks>().size(),                // detector level jet
                        mcdd0cand.pt(), mcdd0cand.eta(), mcdd0cand.phi(), mcdd0cand.m(), mcdd0cand.y(), (mcdd0cand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level D0
                        mcdd0cand.mlScores()[0], mcdd0cand.mlScores()[1], mcdd0cand.mlScores()[2],                                                                  // ML PID scores: background, prompt, non-prompt
                        matchedFrom, selectedAs);                                                                                                                   // D0 = +1, D0bar = -1, neither = 0
          
        }
        
      } else {
        
        // store matched particle and detector level data in one single table (-2 for detector level entries since there is no matched)
        matchJetTable(jetutilities::deltaR(mcpjet, mcpd0cand), mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.tracks_as<aod::JetParticles>().size(),             // particle level jet
                      mcpd0cand.pt(), mcpd0cand.eta(), mcpd0cand.phi(), mcpd0cand.y(), (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt),                // particle level D0
                      -2, -2, -2, -2, -2,                                                                                                                         // detector level jet
                      -2, -2, -2, -2, -2, -2,                                                                                                                     // detector level D0
                      -2, -2, -2,                                                                                                                                 // ML PID scores: background, prompt, non-prompt
                      -2, -2);                                                                                                                                    // D0 = +1, D0bar = -1, neither = 0

      }

    }
    

    //_______________________________________________________________________________________________________________________________________________________
  }
  PROCESS_SWITCH(HfFragmentationFunction, processMcFullChargedMatched, "all particle level and only detector level matched MC HF and jets", false);

  void processMcFullChargedMatched2(aod::JetMcCollisions const& mccollisions,
                               aod::JetCollisionsMCD const& collisions,
                               JetMCDTable const& mcdjets,
                               JetMCPTable const& mcpjets,
                               aod::CandidatesD0MCD const&,
                               aod::CandidatesD0MCP const&,
                               aod::JetTracks const&,
                               aod::JetParticles const&)
  {
    for (const auto& mccollision : mccollisions) {

      registry.fill(HIST("h_collision_counter"), 0.0);

      // skip collisions outside of |z| < vertexZCut
      if (std::abs(mccollision.posZ()) > vertexZCut) {
        continue;
      }
      registry.fill(HIST("h_collision_counter"), 1.0);

      // reconstructed collisions associated to same mccollision
      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {

        registry.fill(HIST("h_collision_counter"), 2.0);
        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
          continue;
        }
        registry.fill(HIST("h_collision_counter"), 3.0);
        // d0 detector level jets associated to the current same collision
        const auto d0mcdJetsPerCollision = mcdjets.sliceBy(d0MCDJetsPerCollisionPreslice, collision.globalIndex());
        for (const auto& mcdjet : d0mcdJetsPerCollision) {

          registry.fill(HIST("h_jet_counter"), 0.5);

          // comparison with fill on bin on 2.5 for sanity check
          if (mcdjet.has_matchedJetCand()) {
            registry.fill(HIST("h_jet_counter"), 1.5);
          }

          // obtain leading HF candidate in jet
          auto mcdd0cand = mcdjet.candidates_first_as<aod::CandidatesD0MCD>();

          // reflection information for storage: D0 = +1, D0bar = -1, neither = 0
          int matchedFrom = 0;
          int decayChannel = 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;
          int selectedAs = 0;

          if (mcdd0cand.flagMcMatchRec() == decayChannel) { // matched to D0 on truth level
            matchedFrom = 1;
          } else if (mcdd0cand.flagMcMatchRec() == -decayChannel) { // matched to D0bar on truth level
            matchedFrom = -1;
          }
          // bitwise AND operation: Checks whether BIT(i) is set, regardless of other bits
          if (mcdd0cand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as D0
            selectedAs = 1;
          } else if (mcdd0cand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as D0bar
            selectedAs = -1;
          }

          // loop through detector level matched to current particle level
          for (const auto& mcpjet : mcdjet.matchedJetCand_as<JetMCPTable>()) {

            registry.fill(HIST("h_jet_counter"), 2.5);

            // obtain leading HF candidate in jet
            auto mcpd0cand = mcpjet.candidates_first_as<aod::CandidatesD0MCP>();

            // store matched particle and detector level data in one single table (calculate angular distance in eta-phi plane on the fly)
            matchJetTable(jetutilities::deltaR(mcpjet, mcpd0cand), mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.tracks_as<aod::JetParticles>().size(),             // particle level jet
                          mcpd0cand.pt(), mcpd0cand.eta(), mcpd0cand.phi(), mcpd0cand.y(), (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt),                // particle level D0
                          jetutilities::deltaR(mcdjet, mcdd0cand), mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), mcdjet.tracks_as<aod::JetTracks>().size(),                // detector level jet
                          mcdd0cand.pt(), mcdd0cand.eta(), mcdd0cand.phi(), mcdd0cand.m(), mcdd0cand.y(), (mcdd0cand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level D0
                          mcdd0cand.mlScores()[0], mcdd0cand.mlScores()[1], mcdd0cand.mlScores()[2],                                                                  // Machine Learning PID scores: background, prompt, non-prompt
                          matchedFrom, selectedAs);                                                                                                                   // D0 = +1, D0bar = -1, neither = 0
          }
        }
      }

      // d0 particle level jets associated to same mccollision
      const auto d0mcpJetsPerMCCollision = mcpjets.sliceBy(d0MCPJetsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& mcpjet : d0mcpJetsPerMCCollision) {

        registry.fill(HIST("h_jet_counter"), 0.0);

        // obtain leading HF particle in jet
        auto mcpd0cand = mcpjet.candidates_first_as<aod::CandidatesD0MCP>();

        if (mcpjet.has_matchedJetCand()) {
          registry.fill(HIST("h_jet_counter"), 1.0);
        } else {
          // store matched particle and detector level data in one single table (-2 for detector level entries since there is no matched)
          matchJetTable(jetutilities::deltaR(mcpjet, mcpd0cand), mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.tracks_as<aod::JetParticles>().size(),             // particle level jet
                        mcpd0cand.pt(), mcpd0cand.eta(), mcpd0cand.phi(), mcpd0cand.y(), (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt),                // particle level D0
                        -2, -2, -2, -2, -2,                                                                                                                         // detector level jet
                        -2, -2, -2, -2, -2, -2,                                                                                                                     // detector level D0
                        -2, -2, -2,                                                                                                                                 // ML PID scores: background, prompt, non-prompt
                        -2, -2);                                                                                                                                    // D0 = +1, D0bar = -1, neither = 0
        }

      }





    }
  }
  PROCESS_SWITCH(HfFragmentationFunction, processMcFullChargedMatched2, "matched MC HF and jets", false);*/

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfFragmentationFunction>(cfgc)};
}
