/**
 * Lambda_c hadron analysis
 * @file Efficiency_run3_style_detector_level.C
 * @brief Calculates pT dependent efficiency of HF selections using methodology from run 3
 * Input: 
 *      - series of backSub_%.jetptMin_to_%.jetptMax_jetpt.root files -> contain histograms such that the background+reflections contribution was removed
 *      - run2_style_efficiency_%.jetptMin_to_%.jetptMax_jetpt.root
 * Outputs: one run3_style_efficiency_%.jetptMin_to_%.jetptMax_jetpt.root -> contain pT,HF dependent run 2 style efficiency, to be used by Efficiency_run3_Style.C
 * 
 * @author: Christian Reckziegel
 * Date: February 2026
 */
#include "../commonUtilities.h"

using namespace std;

struct EfficiencyData {
    // Yield 2D histograms: pT,jet vs pT,HF: first = prompt D0 data distribution, second = non-prompt D0 distribution
    std::pair<TH2D*, TH2D*> hYieldTruth; // all particle level entries (denominator): will be kinematically corrected and folded
    std::pair<TH2D*, TH2D*> hYieldMeasured; // detector level entries (numerator): went over smearing effects and passed the selection cuts

    // Response matrices
    std::pair<RooUnfoldResponse, RooUnfoldResponse> response; // first = prompt D^{0}, second = non-prompt D^{0}
    std::pair<std::vector<TH2D*>, std::vector<TH2D*>> responseProjections; // response projections matrix: first = prompt D^{0}, second = non-prompt D^{0}

    //
    // Kinematic efficiency histograms: first = prompt D^{0}, second = non-prompt D^{0}
    //
    // Response range
    std::pair<TH2D*, TH2D*> hKEffResponseParticle; // response range, particle level entries
    std::pair<TH2D*, TH2D*> hKEffResponseDetector; // response range, detector level entries
    // Response range / total particle range
    std::pair<TH2D*, TH2D*> hKEffTruthTotalParticle; // total truth range, particle level entries
    std::pair<TH2D*, TH2D*> hKEffResponseParticle_Over_TotalParticle; // first = prompt D^{0}, second = non-prompt D^{0}
    // Response range / total detector range
    std::pair<TH2D*, TH2D*> hKEffRecoTotalDetector; // total reco range, detector level entries
    std::pair<TH2D*, TH2D*> hKEffResponseDetector_Over_TotalDetector; // first = prompt D^{0}, second = non-prompt D^{0}

    // All particle level entries with corrections applied (kinematically corrected and folded)
    std::pair<TH2D*, TH2D*> hYieldTruthCorrected;

    // Final distributions to be treated and divided at the end: first = prompt D^{0}, second = non-prompt D^{0}
    std::pair<TH2D*, TH2D*> hNumerator; // detector level entries: went over smearing effects and passed the selection cuts
    std::pair<TH2D*, TH2D*> hDenominator; // all particle level entries: will be kinematically corrected and folded

    // Denominator original and corrected pT histograms: projection of the 2D histograms
    std::pair<TH1D*, TH1D*> hHfPtYieldTruth;
    std::pair<TH1D*, TH1D*> hHfPtYieldTruthCorrected;

    // Run 3 style prompt particle level selection efficiency (fetched from file)
    TH1D* hEfficiency_prompt_part_run3style;
    TH1D* hEfficiency_nonprompt_part_run3style;
    // Run 2 style prompt selection efficiency (fetched from file)
    TH1D* hEfficiency_prompt_run2style;
    TH1D* hEfficiency_nonprompt_run2style;

    //
    // Final efficiency histograms (numerator / denominator): first = prompt D^{0}, second = non-prompt D^{0}
    //
    std::pair<TH1D*, TH1D*> hSelectionEfficiency; // efficiency = numerator / denominator

    // Investigation histograms
    TH1D* hBDTBackgroundScore;

    // Efficiency corrected data (after background subtraction and now efficiency)
    std::pair<TH3D*, TH2D*> hEfficiencyCorrected;
};

// Module to create histograms including interest variable
EfficiencyData createHistograms(const BinningStruct& binning) {

    // Create struct to store data
    EfficiencyData dataContainer;

    // Create 2D histograms for prompt D^{0}: raw and folded
    dataContainer.hYieldTruth.first = new TH2D("hYieldTruthPrompt", "Particle level data prompt yield distribution;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hYieldTruth.first->Sumw2();
    dataContainer.hYieldMeasured.first = new TH2D("hYieldMeasuredPrompt", "Detector level data prompt yield distribution;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{det} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    // Create 2D histograms for non-prompt D^{0}: raw and folded
    dataContainer.hYieldTruth.second = new TH2D("hYieldTruthNonPrompt", "Particle level data non-prompt yield distribution;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hYieldTruth.second->Sumw2();
    dataContainer.hYieldMeasured.second = new TH2D("hYieldMeasuredNonPrompt", "Detector level data non-prompt yield distribution;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{det} (GeV/#it{c})",
                                                binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());

    // Create RooUnfoldResponse objects for prompt and non-prompt D^{0}
    dataContainer.response.first = RooUnfoldResponse(dataContainer.hYieldMeasured.first, dataContainer.hYieldTruth.first); // prompt D^{0}
    dataContainer.response.second = RooUnfoldResponse(dataContainer.hYieldMeasured.second, dataContainer.hYieldTruth.second); // non-prompt D^{0}
    
    // Create projections of response matrix object for prompt and non-prompt D^{0}
    dataContainer.responseProjections.first.push_back(new TH2D("responseProjectionsPtJetPrompt", "Prompt D^{0}'s reponse matrix #it{p}_{T, ch. jet} projection;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, ch. jet}^{gen} (GeV/#it{c})", binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data()));
    dataContainer.responseProjections.first.push_back(new TH2D("responseProjectionsPtHFPrompt", "Prompt D^{0}'s reponse matrix p_{T,D^{0}} projection;#it{p}_{T, D^{0}}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data(), binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data()));
    dataContainer.responseProjections.second.push_back(new TH2D("responseProjectionsPtJetNonPrompt", "Non-prompt D^{0}'s reponse matrix #it{p}_{T, ch. jet} projection;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, ch. jet}^{gen} (GeV/#it{c})", binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data()));
    dataContainer.responseProjections.second.push_back(new TH2D("responseProjectionsPtHFNonPrompt", "Non-prompt D^{0}'s reponse matrix p_{T,D^{0}} projection;#it{p}_{T, D^{0}}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data(), binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data()));
    
    // Creating investigation histogram
    dataContainer.hBDTBackgroundScore = new TH1D("hBDTBackgroundScore", "Entries that didn't pass the cuts;BDT background score;Counts", 100, 0, 1);

    // Kinematic efficiency histograms: truth entries
    dataContainer.hKEffResponseParticle.first = new TH2D("hKEffResponseParticlePrompt", "Truth prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hKEffResponseParticle.second = new TH2D("hKEffResponseParticleNonPrompt", "Truth non-prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hKEffTruthTotalParticle.first = new TH2D("hKEffTruthTotalParticlePrompt", "Truth prompt D^{0}'s within total particle range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})",
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hKEffTruthTotalParticle.second = new TH2D("hKEffTruthTotalParticleNonPrompt", "Truth non-prompt D^{0}'s within total particle range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})",
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    // Kinematic efficiency histograms: reco entries
    dataContainer.hKEffResponseDetector.first = new TH2D("hKEffResponseDetectorPrompt", "Reconstructed prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    dataContainer.hKEffResponseDetector.second = new TH2D("hKEffResponseDetectorNonPrompt", "Reconstructed non-prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    dataContainer.hKEffRecoTotalDetector.first = new TH2D("hKEffRecoTotalDetectorPrompt", "Reconstructed prompt D^{0}'s within total detector range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}",
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    dataContainer.hKEffRecoTotalDetector.second = new TH2D("hKEffRecoTotalDetectorNonPrompt", "Reconstructed non-prompt D^{0}'s within total detector range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}",
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());

    // Final histograms to be treated and divided at the end
    dataContainer.hNumerator.first = new TH2D("hNumeratorPrompt", "Reconstructed prompt D^{0}'s (after selection cuts);#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    dataContainer.hNumerator.second = new TH2D("hNumeratorNonPrompt", "Reconstructed non-prompt D^{0}'s (after selection cuts);#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    dataContainer.hDenominator.first = new TH2D("hDenominatorPrompt", "All truth prompt D^{0}'s;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
            binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
            binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hDenominator.second = new TH2D("hDenominatorNonPrompt", "All truth non-prompt D^{0}'s;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
            binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
            binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    
    std::cout << "Histograms created.\n";

    return dataContainer;
}
//_____________________________________________________________________________________________________________________________________________________________________________________


//_____________________________________________________________________________________________________________________________________________________________________________________
// Modules to fill histograms (of matched data for response matrix and non-matched data for numerator and denominator distributions) from TFile data
void fillHistograms(TFile* fSimulated, TFile* fEffRun2Style, EfficiencyData& dataContainer, const BinningStruct& binning) {
    
    // Get run 3 style particle level prompt selection efficiency
    dataContainer.hEfficiency_prompt_part_run3style = (TH1D*) fEffRun2Style->Get("efficiency_prompt_run3style_particleLevel")->Clone("efficiency_prompt_run3style_particleLevelClone");
    dataContainer.hEfficiency_nonprompt_part_run3style = (TH1D*) fEffRun2Style->Get("efficiency_nonprompt_run3style_particleLevel")->Clone("efficiency_nonprompt_run3style_particleLevelClone");

    // Get run 2 style prompt selection efficiency
    dataContainer.hEfficiency_prompt_run2style = (TH1D*) fEffRun2Style->Get("efficiency_run2style_prompt")->Clone("hEfficiency_prompt_run2styleClone");
    dataContainer.hEfficiency_nonprompt_run2style = (TH1D*) fEffRun2Style->Get("efficiency_run2style_nonprompt")->Clone("hEfficiency_nonprompt_run2styleClone");

    // 0.1 - Accessing simulated data
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = binning.deltaRBinEdges_particle[binning.deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double jetptMin = binning.ptjetBinEdges_particle[0]; // on both levels jet
    const double jetptMax = binning.ptjetBinEdges_particle[binning.ptjetBinEdges_particle.size() - 1]; // on both levels jet
    const double MCPHfPtMincut = binning.ptHFBinEdges_particle[0]; // on particle level D0
    const double MCDHfPtMincut = binning.ptHFBinEdges_detector[0]; // on detector level D0
    const double MCPHfPtMaxcut = binning.ptHFBinEdges_particle[binning.ptHFBinEdges_particle.size() - 1]; // on particle level D0
    const double MCDHfPtMaxcut = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1]; // on detector level D0

    //
    // MC generator level tree and histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)fSimulated->Get("DF_merged/O2matchtable");
    
    // Check for correct access
    if (!tree) {
        std::cout << "Error opening tree.\n";
    }
    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    int MCPjetNConst;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    int MCDjetNConst, MCDhfMatchedFrom, MCDhfSelectedAs;
    bool MCDhfprompt;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;

    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetNConst); // double
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    MCPhfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&MCPhfY);
    tree->SetBranchAddress("fMcHfPrompt",&MCPhfprompt);
    // detector level branches
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fJetNConst",&MCDjetNConst); // int
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);
    tree->SetBranchAddress("fHfMlScore0",&MCDhfMlScore0);
    tree->SetBranchAddress("fHfMlScore1",&MCDhfMlScore1);
    tree->SetBranchAddress("fHfMlScore2",&MCDhfMlScore2);
    tree->SetBranchAddress("fHfMatchedFrom",&MCDhfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&MCDhfSelectedAs);

    // Histogram for efficiency weighting of response matrix
    TH1D* hEffWeight;
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double MCPDeltaR = MCPaxisDistance; // or use MCPaxisDistance
        double MCDDeltaR = MCDaxisDistance; // or use MCDaxisDistance
        
        // Generator level selection cuts
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);//, remove upper bound for particle level
        bool genHfPtRange = ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)); // remove entirely for particle level
        bool genDeltaRRange = ((MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < MCPDeltaRcut));
        bool genAcceptance = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut);
        bool genLevelRange = genJetPtRange && genDeltaRRange && genAcceptance && genHfPtRange;

        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = (MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut);
        bool recoDeltaRRange = (MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < MCDDeltaRcut);
        bool recoAcceptance = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut);
        bool recoLevelRange = recoJetPtRange && recoHfPtRange && recoDeltaRRange && recoAcceptance;
        
        // Get the threshold for this pT range: TODO - do NOT erase this, BDT cuts will be included in next hyperloop wagon
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isRealD0 = isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs);
        
        // 1 --- Selection efficiencies
        if (genLevelRange) {

            // Matched entries must be of real D0s, not reflections or combinatorial background
            if (MCDhfmatch) {
                if (isRealD0) {
                    // Fill particle level entry
                    if (MCPhfprompt) {
                        dataContainer.hYieldTruth.first->Fill(MCPjetPt, MCPhfPt); // prompt particle level denominator
                    } else{
                        dataContainer.hYieldTruth.second->Fill(MCPjetPt, MCPhfPt); // non-prompt particle level denonimator
                    }
                    // Fill detector level entry
                    if (recoLevelRange && passBDTcut) {
                        if (MCDhfprompt) {
                            // prompt detector level numerator
                            dataContainer.hYieldMeasured.first->Fill(MCDjetPt, MCDhfPt);
                        } else{
                            // non-prompt detector level numerator
                            dataContainer.hYieldMeasured.second->Fill(MCDjetPt, MCDhfPt);
                        }
                    } else if (recoLevelRange && !passBDTcut) {
                        dataContainer.hBDTBackgroundScore->Fill(MCDhfMlScore0);
                    }
                    
                }
            } else {// fill particle level even if not matched
                if (MCPhfprompt) {
                    dataContainer.hYieldTruth.first->Fill(MCPjetPt, MCPhfPt); // prompt particle level denominator
                } else{
                    dataContainer.hYieldTruth.second->Fill(MCPjetPt, MCPhfPt); // non-prompt particle level denonimator
                }
            }   
        }

        // 2 --- Response matrix: fill histograms considering jet pT and detector acceptance (response range)
        if (MCDhfmatch && isRealD0) {
            if (genLevelRange && recoLevelRange) {

                // fill 2D yields histograms
                if (MCPhfprompt) {

                    // Get efficiency estimate to weight the response matrix
                    hEffWeight = (TH1D*) fEffRun2Style->Get("efficiency_run2style_prompt");
                    // Get efficiency estimate to weight the response matrix: jet pT shape is influenced by D0 pT efficiency
                    int bin = hEffWeight->FindBin(MCDhfPt);
                    double estimatedEfficiency = hEffWeight->GetBinContent(bin);
                    if (estimatedEfficiency == 0) {
                        std::cout << "Warning: Prompt efficiency is zero for pT,HF = " << MCDhfPt << " GeV/c with bin " << bin << " of efficiency_prompt run 2 histogram. Setting it to 1. How to properly deal with these entries?" << std::endl;
                        estimatedEfficiency = 1; // Avoid division by zero
                    }

                    // Fill 4D RooUnfoldResponse object
                    dataContainer.response.first.Fill(MCDjetPt, MCDhfPt, MCPjetPt, MCPhfPt, 1 / estimatedEfficiency); // jet pT shape is influenced by D0 pT efficiency
                    dataContainer.responseProjections.first[0]->Fill(MCDjetPt, MCPjetPt,1 / estimatedEfficiency); // pT,jet projection, prompt D^{0}
                    dataContainer.responseProjections.first[1]->Fill(MCDhfPt, MCPhfPt, 1 / estimatedEfficiency); // pT,HF projection, prompt D^{0}
                } else{

                    // Get efficiency estimate to weight the response matrix
                    hEffWeight = (TH1D*) fEffRun2Style->Get("efficiency_run2style_nonprompt");
                    // Get efficiency estimate to weight the response matrix: jet pT shape is influenced by D0 pT efficiency
                    int bin = hEffWeight->FindBin(MCDhfPt);
                    double estimatedEfficiency = hEffWeight->GetBinContent(bin);
                    if (estimatedEfficiency == 0) {
                        //std::cout << "Warning: Prompt efficiency is zero for pT = " << MCDhfPt << " with bin " << bin << ". Setting it to 1." << std::endl;
                        estimatedEfficiency = 1; // Avoid division by zero
                    }

                    // Fill 4D RooUnfoldResponse object
                    dataContainer.response.second.Fill(MCDjetPt, MCDhfPt, MCPjetPt, MCPhfPt, 1 / estimatedEfficiency); // jet pT shape is influenced by D0 pT efficiency
                    dataContainer.responseProjections.second[0]->Fill(MCDjetPt, MCPjetPt,1 / estimatedEfficiency); // pT,jet projection, non-prompt D^{0}
                    dataContainer.responseProjections.second[1]->Fill(MCDhfPt, MCPhfPt, 1 / estimatedEfficiency); // pT,HF projection, non-prompt D^{0}
                }
            }
        }
        
        // 3 --- Kinematic efficiencies
        if (MCDhfmatch && isRealD0) {
            // Particle level kinematic efficiency
            if (genLevelRange) {
                // prompt D0s
                if (MCPhfprompt) {
                    // fill prompt 2D yield: total particle range (denominator)
                    dataContainer.hKEffTruthTotalParticle.first->Fill(MCPjetPt, MCPhfPt);
                    if (recoLevelRange) {
                        // fill prompt 2D yield: response range
                        dataContainer.hKEffResponseParticle.first->Fill(MCPjetPt, MCPhfPt);
                    }
                } else {
                    // non-prompt D0s
                    dataContainer.hKEffTruthTotalParticle.second->Fill(MCPjetPt, MCPhfPt);
                    if (recoLevelRange) {
                        // fill prompt 2D yield: particle level matched
                        dataContainer.hKEffResponseParticle.second->Fill(MCPjetPt, MCPhfPt);
                    }
                } 
            }

            // Detector level kinematic efficiency
            if (recoLevelRange) {
                if (MCPhfprompt) {
                    // fill prompt 2D yield: total detector range
                    dataContainer.hKEffRecoTotalDetector.first->Fill(MCDjetPt, MCDhfPt);
                    if (genLevelRange) {
                        // fill prompt 2D yield: response range
                        dataContainer.hKEffResponseDetector.first->Fill(MCDjetPt, MCDhfPt);
                    }
                } else {
                    dataContainer.hKEffRecoTotalDetector.second->Fill(MCDjetPt, MCDhfPt);
                    if (genLevelRange) {
                        // fill prompt 2D yield: detector level matched
                        dataContainer.hKEffResponseDetector.second->Fill(MCDjetPt, MCDhfPt);
                    }
                }
            }
        }
    }

    
    std::cout << "pT,HF histograms (particle and detector level, prompt and non-prompt), response matrix and kinematic correction histograms filled.\n";

}
//_____________________________________________________________________________________________________________________________________________________________________________________


void performEfficiencyCorrection(TFile* fBackSub, EfficiencyData& dataContainer, double jetptMin, double jetptMax) {
    
    //
    // Calculating pT,HF dependent efficiency distribution
    //
    
    // 0.1 - Calculate truth kinematic efficiency
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first = (TH2D*)dataContainer.hKEffResponseParticle.first->Clone("hKEffTruthPrompt");
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second = (TH2D*)dataContainer.hKEffResponseParticle.second->Clone("hKEffTruthNonPrompt");
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->Sumw2(); // necessary for correct error propagation
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->Sumw2();
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->Divide(dataContainer.hKEffTruthTotalParticle.first);
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->Divide(dataContainer.hKEffTruthTotalParticle.second);
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->SetTitle("#varepsilon_{kinematic} of particle level prompt D^{0};#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->SetTitle("#varepsilon_{kinematic} of particle level non-prompt D^{0};#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");
    
    // 0.2 - Calculate reco kinematic efficiency
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first = (TH2D*)dataContainer.hKEffResponseDetector.first->Clone("hKEffRecoPrompt");
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second = (TH2D*)dataContainer.hKEffResponseDetector.second->Clone("hKEffRecoNonPrompt");
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->Sumw2();
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->Sumw2();
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->Divide(dataContainer.hKEffRecoTotalDetector.first);
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->Divide(dataContainer.hKEffRecoTotalDetector.second);
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->SetTitle("#varepsilon_{kinematic} of detector level prompt D^{0};#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->SetTitle("#varepsilon_{kinematic} of detector level non-prompt D^{0};#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");

    // 1- Correct the particle level input histogram: remove entries from input distribution whose matched entry is outside of detector level range encoding
    dataContainer.hYieldTruthCorrected.first = (TH2D*)dataContainer.hYieldTruth.first->Clone("hYieldTruthCorrectedPrompt");
    dataContainer.hYieldTruthCorrected.second = (TH2D*)dataContainer.hYieldTruth.second->Clone("hYieldTruthCorrectedNonPrompt");
    dataContainer.hYieldTruthCorrected.first->Sumw2();
    dataContainer.hYieldTruthCorrected.second->Sumw2();
    dataContainer.hYieldTruthCorrected.first->Multiply(dataContainer.hKEffResponseParticle_Over_TotalParticle.first);
    dataContainer.hYieldTruthCorrected.second->Multiply(dataContainer.hKEffResponseParticle_Over_TotalParticle.second);

    // 2 - Fold the corrected distribution into the response matrix
    dataContainer.hYieldTruthCorrected.first = manualFolding(dataContainer.response.first, dataContainer.hYieldTruthCorrected.first, dataContainer.hYieldMeasured.first);
    dataContainer.hYieldTruthCorrected.second = manualFolding(dataContainer.response.second, dataContainer.hYieldTruthCorrected.second, dataContainer.hYieldMeasured.second);

    // 3 - Correct the detector level output histogram: add entries to the output distribution which would be present at outside ranges of the response matrix at particle level
    dataContainer.hYieldTruthCorrected.first->Divide(dataContainer.hKEffResponseDetector_Over_TotalDetector.first);
    dataContainer.hYieldTruthCorrected.second->Divide(dataContainer.hKEffResponseDetector_Over_TotalDetector.second);

    // 3.1 - Obtain pT projection histograms
    int minBin, maxBin;
    minBin = dataContainer.hYieldTruth.first->GetXaxis()->FindBin(jetptMin);
    maxBin = dataContainer.hYieldTruth.first->GetXaxis()->FindBin(jetptMax - 1e-6); // // Tiny epsilon to stay within range: This includes the nearest bin center for jetptMax, but if jetptMax lies between two bins, it might give slightly unintuitive results
    dataContainer.hHfPtYieldTruth.first = dataContainer.hYieldTruth.first->ProjectionY("hHfPtYieldTruthPrompt", minBin, maxBin);
    dataContainer.hHfPtYieldTruth.first->SetTitle("Denominator prompt D^{0} p_{T} distribution before correction; #it{p}_{T, D^{0}}^{gen}; Counts");
    minBin = dataContainer.hYieldTruth.second->GetXaxis()->FindBin(jetptMin);
    maxBin = dataContainer.hYieldTruth.second->GetXaxis()->FindBin(jetptMax - 1e-6);
    dataContainer.hHfPtYieldTruth.second = dataContainer.hYieldTruth.second->ProjectionY("hHfPtYieldTruthNonPrompt", minBin, maxBin);
    dataContainer.hHfPtYieldTruth.second->SetTitle("Denominator non-prompt D^{0} p_{T} distribution before correction; #it{p}_{T, D^{0}}^{gen}; Counts");
    minBin = dataContainer.hYieldTruthCorrected.first->GetXaxis()->FindBin(jetptMin);
    maxBin = dataContainer.hYieldTruthCorrected.first->GetXaxis()->FindBin(jetptMax - 1e-6);
    dataContainer.hHfPtYieldTruthCorrected.first = dataContainer.hYieldTruthCorrected.first->ProjectionY("hHfPtYieldTruthCorrectedPrompt", minBin, maxBin);
    dataContainer.hHfPtYieldTruthCorrected.first->SetTitle("Denominator prompt D^{0} p_{T} distribution after correction; #it{p}_{T, D^{0}}^{gen}; Counts");
    minBin = dataContainer.hYieldTruthCorrected.second->GetXaxis()->FindBin(jetptMin);
    maxBin = dataContainer.hYieldTruthCorrected.second->GetXaxis()->FindBin(jetptMax - 1e-6);
    dataContainer.hHfPtYieldTruthCorrected.second = dataContainer.hYieldTruthCorrected.second->ProjectionY("hHfPtYieldTruthCorrectedNonPrompt", minBin, maxBin);
    dataContainer.hHfPtYieldTruthCorrected.second->SetTitle("Denominator non-prompt D^{0} p_{T} distribution after correction; #it{p}_{T, D^{0}}^{gen}; Counts");

    // 4 - Calculate the efficiency distributions with the pT projections
    dataContainer.hSelectionEfficiency.first = dataContainer.hYieldMeasured.first->ProjectionY("hSelectionEfficiencyPrompt", minBin, maxBin);
    dataContainer.hSelectionEfficiency.second = dataContainer.hYieldMeasured.second->ProjectionY("hSelectionEfficiencyNonPrompt", minBin, maxBin);
    dataContainer.hSelectionEfficiency.first->SetTitle("Efficiency prompt D^{0} p_{T} distribution; #it{p}_{T, D^{0}}^{reco}; Efficiency#times Acceptance");
    dataContainer.hSelectionEfficiency.second->SetTitle("Efficiency non-prompt D^{0} p_{T} distribution; #it{p}_{T, D^{0}}^{reco}; Efficiency#times Acceptance");
    dataContainer.hSelectionEfficiency.first->Sumw2();
    dataContainer.hSelectionEfficiency.second->Sumw2();
    dataContainer.hSelectionEfficiency.first->Divide(dataContainer.hHfPtYieldTruthCorrected.first);
    dataContainer.hSelectionEfficiency.second->Divide(dataContainer.hHfPtYieldTruthCorrected.second);
    
    //
    // Correcting 3D background subtracted distribution with efficiency for each pT,HF bin
    //
    // Clone for later comparison
    TH3D* h3DBackgroundSubtractedCloneComparison = (TH3D*)fBackSub->Get("h3DBackgroundSubtracted")->Clone("h3DBackgroundSubtractedCloneComparison");
    cleanNaNs(h3DBackgroundSubtractedCloneComparison);
    TH2D* h2DBackgroundSubtractedCloneComparison = (TH2D*)h3DBackgroundSubtractedCloneComparison->Project3D("yx");
    TH1D* h1DBackgroundSubtractedCloneComparison = (TH1D*)h2DBackgroundSubtractedCloneComparison->ProjectionY("h1DBackgroundSubtractedCloneComparison",1,h2DBackgroundSubtractedCloneComparison->GetXaxis()->GetNbins());
    // 1 - Clone background subtracted distribution
    dataContainer.hEfficiencyCorrected.first = (TH3D*)fBackSub->Get("h3DBackgroundSubtracted")->Clone("h3DEfficiencyCorrected");
    cleanNaNs(dataContainer.hEfficiencyCorrected.first);
    dataContainer.hEfficiencyCorrected.first->SetTitle("Background subtracted, efficiency corrected");
    int xBins = dataContainer.hEfficiencyCorrected.first->GetXaxis()->GetNbins();
    int yBins = dataContainer.hEfficiencyCorrected.first->GetYaxis()->GetNbins();
    int zBins = dataContainer.hEfficiencyCorrected.first->GetZaxis()->GetNbins();
    for (int xBin = 1; xBin <= xBins; xBin++) {
        for (int yBin = 1; yBin <= yBins; yBin++) {
            for (int zBin = 1; zBin <= zBins; zBin++) { // pT,HF bins will be corrected by 1/prompt_efficiency
                
                // Get current content and error
                double content = dataContainer.hEfficiencyCorrected.first->GetBinContent(xBin, yBin, zBin);
                double error = dataContainer.hEfficiencyCorrected.first->GetBinError(xBin, yBin, zBin);
                
                // Get the pT,HF bin center from the Z axis of the 3D histogram
                double ptDcenter = dataContainer.hEfficiencyCorrected.first->GetZaxis()->GetBinCenter(zBin);

                // if (std::isnan(content)) {
                //     std::cout << "Warning: nan value encountered on content of pT,HF center with " << ptDcenter << " GeV/c\t (xBin, yBin, zBin) = (" << xBin << "," << yBin << "," << zBin << ")" << std::endl;
                //     std::cout << "Setting it to zero." << std::endl;
                //     content = 0;
                //     error = 0;
                // }
                
                // Find corresponding bin in efficiency histogram
                int effBin = dataContainer.hSelectionEfficiency.first->GetXaxis()->FindBin(ptDcenter);
                double eff = dataContainer.hSelectionEfficiency.first->GetBinContent(effBin);

                // Avoid divide by zero or nonsense values
                if (eff > 0) {
                    double correction = 1. / eff;
                    dataContainer.hEfficiencyCorrected.first->SetBinContent(xBin, yBin, zBin, content * correction);
                    dataContainer.hEfficiencyCorrected.first->SetBinError(xBin, yBin, zBin, error * correction);
                } else {
                    // Optionally warn about zero efficiency
                    std::cerr << "Zero or invalid efficiency at pT,HF = " << ptDcenter << " (eff bin = " << effBin << ")" << std::endl;
                    if (eff < 0) {
                        std::cout << "Warning: negative efficiencies!\tWhere do they come from?\tHow to deal with the corresponding entries?" << std::endl;
                    } else if (eff == 0) {
                        dataContainer.hEfficiencyCorrected.first->SetBinContent(xBin, yBin, zBin, 0.);
                        dataContainer.hEfficiencyCorrected.first->SetBinError(xBin, yBin, zBin, 0.);
                    }
                    
                    
                }
            }
        }
    }
    
    // 2 - Project the z axis (pT,HF) in order to obtain a 2D distribution of pT,jet (x axis) vs DeltaR (y axis)
    TH2D* h2D = (TH2D*)dataContainer.hEfficiencyCorrected.first->Project3D("yx");
    dataContainer.hEfficiencyCorrected.second = (TH2D*)h2D->Clone("h2DEfficiencyCorrected");
    dataContainer.hEfficiencyCorrected.second->SetTitle("Background subtracted, efficiency corrected");
    TH1D* h1D = (TH1D*)h2D->ProjectionY("h1D",1,h2D->GetXaxis()->GetNbins());
    TCanvas* cBefAftEfficiencyCorrection = new TCanvas("cBefAftEfficiencyCorrection","Before and after efficiency correction",1920,1080);
    cBefAftEfficiencyCorrection->cd();
    TLegend* legBefAft = new TLegend(0.6,0.6,0.8,0.75);
    h1DBackgroundSubtractedCloneComparison->SetLineColor(kRed+1);
    h1D->SetLineColor(kBlue+1);
    h1D->Draw();
    h1DBackgroundSubtractedCloneComparison->Draw("same");
    legBefAft->AddEntry(h1DBackgroundSubtractedCloneComparison,"Before correction","lp");
    legBefAft->AddEntry(h1D,"After correction","lp");
    legBefAft->Draw();

    std::cout << "Efficiency correction performed.\n";
}

void PlotHistograms(const EfficiencyData& dataContainer, double jetptMin, double jetptMax, const BinningStruct& binning) {

    //
    // 2D yield histograms
    //
    TCanvas* cYield = new TCanvas("cYield","2D yield histograms");
    cYield->Divide(2,2);
    cYield->SetCanvasSize(1800,1000);
    cYield->cd(1);
    //dataContainer.hYieldTruth.first->SetStats(0);
    dataContainer.hYieldTruth.first->Draw("COLZ");
    cYield->cd(2);
    //dataContainer.hYieldTruth.second->SetStats(0);
    dataContainer.hYieldTruth.second->Draw("COLZ");
    cYield->cd(3);
    //dataContainer.hYieldMeasured.first->SetStats(0);
    dataContainer.hYieldMeasured.first->Draw("COLZ");
    cYield->cd(4);
    //dataContainer.hYieldMeasured.second->SetStats(0);
    dataContainer.hYieldMeasured.second->Draw("COLZ");
    cYield->Update();

    //
    // Response matrix representation in 2D histogram
    //
    TCanvas* cResponse = new TCanvas("cResponse","Response matrices for all pT,HF bins");
    cResponse->SetCanvasSize(1800,1000);
    cResponse->Divide(2,2);
    std::pair<const TH2*, const TH2*> hresponse2D;
    hresponse2D.first = dataContainer.response.first.Hresponse();
    hresponse2D.second = dataContainer.response.second.Hresponse();
    std::pair<TH2D*, TH2D*> hresponse2DClone;
    hresponse2DClone.first = static_cast<TH2D*>(hresponse2D.first->Clone("hResponse2DPrompt"));
    hresponse2DClone.second = static_cast<TH2D*>(hresponse2D.second->Clone("hResponse2DNonPrompt"));
    hresponse2DClone.first->SetTitle("2D response matrix from 4D RooUnfoldResponse - prompt D^{0}'s;2D Reconstructed;2D Truth");
    hresponse2DClone.second->SetTitle("2D response matrix from 4D RooUnfoldResponse - non-prompt D^{0}'s;2D Reconstructed;2D Truth");
    cResponse->cd(1);
    hresponse2DClone.first->Draw("colz");
    cResponse->cd(2);
    hresponse2DClone.second->Draw("colz");

    // Separate response matrices in two canvases
    TCanvas* cResponsePrompt = new TCanvas("cResponsePrompt","Response matrix for prompt D^{0}'s");
    cResponsePrompt->SetCanvasSize(1800,1000);
    cResponsePrompt->cd();
    hresponse2DClone.first->Draw("colz");
    TCanvas* cResponseNonPrompt = new TCanvas("cResponseNonPrompt","Response matrix for non-prompt D^{0}'s");
    cResponseNonPrompt->SetCanvasSize(1800,1000);
    cResponseNonPrompt->cd();
    hresponse2DClone.second->Draw("colz");

    TCanvas* cResponseProjections = new TCanvas("cResponseProjections","Response matrix projections");
    cResponseProjections->SetCanvasSize(1800,1000);
    cResponseProjections->Divide(2,2);
    cResponseProjections->cd(1);
    dataContainer.responseProjections.first[0]->Draw("colz");
    cResponseProjections->cd(2);
    dataContainer.responseProjections.first[1]->Draw("colz");
    cResponseProjections->cd(3);
    dataContainer.responseProjections.second[0]->Draw("colz");
    cResponseProjections->cd(4);
    dataContainer.responseProjections.second[1]->Draw("colz");

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    
    // Kinematic efficiency histograms: truth entries
    TCanvas* cKEffTruthPrompt = new TCanvas("cKEffTruthPrompt","Kinematic efficiency histograms for truth prompt entries");
    cKEffTruthPrompt->SetCanvasSize(1800,1000);
    cKEffTruthPrompt->Divide(2,2);
    cKEffTruthPrompt->cd(1);
    dataContainer.hKEffTruthTotalParticle.first->Draw("colz");
    cKEffTruthPrompt->cd(2);
    dataContainer.hKEffResponseParticle.first->Draw("colz");
    cKEffTruthPrompt->cd(3);
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->Draw("text");
    TCanvas* cKEffTruthNonPrompt = new TCanvas("cKEffTruthNonPrompt","Kinematic efficiency histograms for truth non-prompt entries");
    cKEffTruthNonPrompt->SetCanvasSize(1800,1000);
    cKEffTruthNonPrompt->Divide(2,2);
    cKEffTruthNonPrompt->cd(1);
    dataContainer.hKEffTruthTotalParticle.second->Draw("colz");
    cKEffTruthNonPrompt->cd(2);
    dataContainer.hKEffResponseParticle.second->Draw("colz");
    cKEffTruthNonPrompt->cd(3);
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->Draw("text");

    // Kinematic efficiency histograms: reco entries
    TCanvas* cKEffRecoPrompt = new TCanvas("cKEffRecoPrompt","Kinematic efficiency histograms for reco prompt entries");
    cKEffRecoPrompt->SetCanvasSize(1800,1000);
    cKEffRecoPrompt->Divide(2,2);
    cKEffRecoPrompt->cd(1);
    dataContainer.hKEffRecoTotalDetector.first->Draw("colz");
    cKEffRecoPrompt->cd(2);
    dataContainer.hKEffResponseDetector.first->Draw("colz");
    cKEffRecoPrompt->cd(3);
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->Draw("text");
    TCanvas* cKEffRecoNonPrompt = new TCanvas("cKEffRecoNonPrompt","Kinematic efficiency histograms for reco non-prompt entries");
    cKEffRecoNonPrompt->SetCanvasSize(1800,1000);
    cKEffRecoNonPrompt->Divide(2,2);
    cKEffRecoNonPrompt->cd(1);
    dataContainer.hKEffRecoTotalDetector.second->Draw("colz");
    cKEffRecoNonPrompt->cd(2);
    dataContainer.hKEffResponseDetector.second->Draw("colz");
    cKEffRecoNonPrompt->cd(3);
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->Draw("text");

    // All kinematic efficiency histograms
    TCanvas* cKEffAll = new TCanvas("cKEffAll","All kinematic efficiency histograms for all entries");
    cKEffAll->SetCanvasSize(1800,1000);
    cKEffAll->Divide(2,2);
    cKEffAll->cd(1);
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->SetMarkerSize(1.); // default = 1.5
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->Draw("text");
    cKEffAll->cd(2);
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->SetMarkerSize(1.);
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->Draw("text");
    cKEffAll->cd(3);
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->SetMarkerSize(1.);
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->Draw("text");
    cKEffAll->cd(4);
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->SetMarkerSize(1.);
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->Draw("text");
    cKEffAll->Update();

    TCanvas* cYieldTruthCorrected = new TCanvas("cYieldTruthCorrected","Corrected denominator histograms and original ones");
    cYieldTruthCorrected->SetCanvasSize(1800,1000);
    cYieldTruthCorrected->Divide(2,2);
    cYieldTruthCorrected->cd(1);
    dataContainer.hYieldTruthCorrected.first->SetTitle("Corrected particle level data prompt yield distribution");
    dataContainer.hYieldTruthCorrected.first->Draw("colz");
    cYieldTruthCorrected->cd(2);
    dataContainer.hYieldTruthCorrected.second->SetTitle("Corrected particle level data non-prompt yield distribution");
    dataContainer.hYieldTruthCorrected.second->Draw("colz");
    cYieldTruthCorrected->cd(3);
    dataContainer.hYieldTruth.first->Draw("colz");
    cYieldTruthCorrected->cd(4);
    dataContainer.hYieldTruth.second->Draw("colz");

    TCanvas* cPtCorrectionComparison = new TCanvas("cPtCorrectionComparison","pT correction comparison");
    cPtCorrectionComparison->SetCanvasSize(1800,1000);
    cPtCorrectionComparison->Divide(2,2);
    cPtCorrectionComparison->cd(1);
    dataContainer.hHfPtYieldTruth.first->Draw("colz");
    cPtCorrectionComparison->cd(2);
    dataContainer.hHfPtYieldTruth.second->Draw("colz");
    cPtCorrectionComparison->cd(3);
    dataContainer.hHfPtYieldTruthCorrected.first->Draw("colz");
    cPtCorrectionComparison->cd(4);
    dataContainer.hHfPtYieldTruthCorrected.second->Draw("colz");

    TCanvas* cSelectionEfficiency = new TCanvas("cSelectionEfficiency","Efficiency histograms",1800,1000);
    cSelectionEfficiency->cd();
    dataContainer.hSelectionEfficiency.first->SetMarkerStyle(kFullCircle);
    dataContainer.hSelectionEfficiency.first->SetMarkerColor(30+1*10);
    dataContainer.hSelectionEfficiency.first->SetLineColor(30+1*10);
    dataContainer.hSelectionEfficiency.second->SetMarkerStyle(kFullCircle);
    dataContainer.hSelectionEfficiency.second->SetMarkerColor(30+2*10);
    dataContainer.hSelectionEfficiency.second->SetLineColor(30+2*10);
    dataContainer.hSelectionEfficiency.first->Draw();
    dataContainer.hSelectionEfficiency.second->Draw("same");
    TLegend* legendEff = new TLegend(0.65,0.25,0.8,0.38);
    legendEff->AddEntry(dataContainer.hSelectionEfficiency.first,"Prompt D^{0}", "lpe");
    legendEff->AddEntry(dataContainer.hSelectionEfficiency.second,"Non-prompt D^{0}", "lpe");
    legendEff->Draw();

    // Separating in two canvases
    TCanvas* cSelectionEfficiencyPrompt = new TCanvas("cSelectionEfficiencyPrompt","Efficiency histograms for prompt D^{0}'s");
    cSelectionEfficiencyPrompt->SetCanvasSize(1800,1000);
    cSelectionEfficiencyPrompt->cd();
    dataContainer.hSelectionEfficiency.first->Draw();
    TCanvas* cSelectionEfficiencyNonPrompt = new TCanvas("cSelectionEfficiencyNonPrompt","Efficiency histograms for non-prompt D^{0}'s");
    cSelectionEfficiencyNonPrompt->SetCanvasSize(1800,1000);
    cSelectionEfficiencyNonPrompt->cd();
    dataContainer.hSelectionEfficiency.second->Draw();

    TCanvas* cCorrectedData = new TCanvas("cCorrectedData","Background subtracted, efficiency corrected data");
    cCorrectedData->Divide(2,2);
    cCorrectedData->cd(1);
    dataContainer.hEfficiencyCorrected.first->Draw("colz");
    cCorrectedData->cd(2);
    dataContainer.hEfficiencyCorrected.second->Draw("colz");
    cCorrectedData->cd(3);
    dataContainer.hEfficiencyCorrected.second->ProjectionY("hEfficiencyCorrectedProjectionDeltaRPrompt",1,dataContainer.hEfficiencyCorrected.second->GetXaxis()->GetNbins())->Draw();

    TCanvas* cCorrectedData2DOnly = new TCanvas("cCorrectedData2DOnly","Background subtracted, efficiency corrected 2D data");
    cCorrectedData2DOnly->cd();
    dataContainer.hEfficiencyCorrected.second->Draw("colz");
    TCanvas* cCorrectedData1DOnly = new TCanvas("cCorrectedData1DOnly","Background subtracted, efficiency corrected 1D data");
    cCorrectedData1DOnly->cd();
    //dataContainer.hEfficiencyCorrected.second->SetMarkerStyle(kCircle);
    dataContainer.hEfficiencyCorrected.second->ProjectionY("hEfficiencyCorrectedProjectionDeltaRPrompt",1,dataContainer.hEfficiencyCorrected.second->GetXaxis()->GetNbins())->Draw();

    // 2D efficiency
    TCanvas* cSelectionEfficiencyPrompt2D = new TCanvas("cSelectionEfficiencyPrompt2D","2D prompt D0s efficiency", 1800, 1000);
    cSelectionEfficiencyPrompt2D->Divide(2,2);
    cSelectionEfficiencyPrompt2D->cd(1);
    dataContainer.hYieldTruthCorrected.first->Draw("colz");
    cSelectionEfficiencyPrompt2D->cd(2);
    dataContainer.hYieldMeasured.first->Draw("colz");
    cSelectionEfficiencyPrompt2D->cd(3);
    TH2D* hEfficiency2D = (TH2D*) dataContainer.hYieldMeasured.first->Clone("hEfficiency2D");
    hEfficiency2D->Divide(dataContainer.hYieldTruthCorrected.first);
    hEfficiency2D->SetTitle("Prompt D^{0}s 2D efficiency");
    hEfficiency2D->Draw("colz");
    cSelectionEfficiencyPrompt2D->cd(4);
    hEfficiency2D->Draw("text");

    TCanvas* cEfficiencyPerJetPt = new TCanvas("cEfficiencyPerJetPt"," Efficiency per jet pT", 1800, 1000);
    cEfficiencyPerJetPt->cd();
    TLegend* leg_per_jetpt = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t iPtJetBin = 0; iPtJetBin < binning.ptjetBinEdges_particle.size() - 1; iPtJetBin++) {
        // Get numerator
        TH1D* hNumerator = dataContainer.hYieldMeasured.first->ProjectionY(Form("hNumerator_%zu", iPtJetBin+1), iPtJetBin+1, iPtJetBin+1); // (x;y) = (pT,jet;pT,D0)
        hNumerator->Sumw2();
        // Get denominator
        TH1D* hDenominator = dataContainer.hYieldTruthCorrected.first->ProjectionY(Form("hDenominator_%zu", iPtJetBin), iPtJetBin+1, iPtJetBin+1); // (x;y) = (pT,jet;pT,D0)
        hDenominator->Sumw2();
        // Divide them
        TH1D* hEffJetPtInterval = (TH1D*) hNumerator->Clone(Form("hEffJetPtInterval_%zu",iPtJetBin+1));
        hEffJetPtInterval->Divide(hDenominator);

        // Set drawing style
        hEffJetPtInterval->SetTitle("Run 3 style prompt D^{0}s efficiencies");
        hEffJetPtInterval->SetLineColor(kBlack+iPtJetBin);
        hEffJetPtInterval->SetMarkerColor(kBlack+iPtJetBin);
        hEffJetPtInterval->SetMarkerStyle(kOpenCircle);
        hEffJetPtInterval->SetLineWidth(1);
        hEffJetPtInterval->GetYaxis()->SetRangeUser(0.,0.25);
        hEffJetPtInterval->Draw(iPtJetBin == 0 ? "" : "same");

        // Add to legend
        double iJetptMin = binning.ptjetBinEdges_particle[iPtJetBin];
        double iJetptMax = binning.ptjetBinEdges_particle[iPtJetBin+1];
        leg_per_jetpt->AddEntry(hEffJetPtInterval,Form("p_{T,jet}#in[%.0f-%.0f] GeV/c",iJetptMin,iJetptMax),"ple");
    }
    leg_per_jetpt->Draw();

    TLegend* legendEffTypes = new TLegend(0.65,0.45,0.9,0.68);
    TCanvas* cEfficiencyTypes = new TCanvas("cEfficiencyTypes","Efficiency types comparison",1800,1000);
    cEfficiencyTypes->cd();
    // Run 3 style particle level
    double minRangeY = std::min({dataContainer.hEfficiency_prompt_part_run3style->GetMinimum(), dataContainer.hEfficiency_prompt_run2style->GetMinimum(), dataContainer.hSelectionEfficiency.first->GetMinimum()});
    double maxRangeY = std::max({dataContainer.hEfficiency_prompt_part_run3style->GetMaximum(), dataContainer.hEfficiency_prompt_run2style->GetMaximum(), dataContainer.hSelectionEfficiency.first->GetMaximum()});
    dataContainer.hEfficiency_prompt_part_run3style->GetYaxis()->SetRangeUser(0.9 * minRangeY, 1.1 * maxRangeY);
    dataContainer.hEfficiency_prompt_part_run3style->SetMarkerStyle(kOpenCircle);
    dataContainer.hEfficiency_prompt_part_run3style->SetMarkerColor(kRed+2);
    dataContainer.hEfficiency_prompt_part_run3style->SetLineColor(kRed+2);
    dataContainer.hEfficiency_prompt_part_run3style->Draw();
    legendEffTypes->AddEntry(dataContainer.hEfficiency_prompt_part_run3style,"Run 3 style prompt particle level efficiency", "lpe");
    // Run 2 style
    dataContainer.hEfficiency_prompt_run2style->SetMarkerStyle(kOpenCircle);
    dataContainer.hEfficiency_prompt_run2style->SetMarkerColor(kBlue+2);
    dataContainer.hEfficiency_prompt_run2style->SetLineColor(kBlue+2);
    dataContainer.hEfficiency_prompt_run2style->Draw("same");
    legendEffTypes->AddEntry(dataContainer.hEfficiency_prompt_run2style,"Run 2 style prompt efficiency", "lpe");
    // Run 3 style detector level
    dataContainer.hSelectionEfficiency.first->SetMarkerStyle(kOpenCircle);
    dataContainer.hSelectionEfficiency.first->SetMarkerColor(kGreen+2);
    dataContainer.hSelectionEfficiency.first->SetLineColor(kGreen+2);
    dataContainer.hSelectionEfficiency.first->Draw("same");
    legendEffTypes->AddEntry(dataContainer.hSelectionEfficiency.first,"Run 3 style prompt detector level efficiency", "lpe");
    legendEffTypes->Draw();

    //
    // Storing in a single pdf file
    //
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../Images/2-Efficiency/Run3Style/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    
    cYield->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cResponse->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponseProjections->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cKEffAll->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cYieldTruthCorrected->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cPtCorrectionComparison->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cSelectionEfficiency->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cCorrectedData->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cCorrectedData2DOnly->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cSelectionEfficiencyPrompt2D->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cEfficiencyPerJetPt->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cEfficiencyTypes->Print(imagePath + Form("run3_style_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
    // in separate .pdf files
    cKEffAll->Print(imagePath + "Efficiency_run3style_kinematic_efficiency_" + sEmmaBins + ".pdf");

    std::cout << "Histograms plotted.\n";

}

void SaveData(const EfficiencyData& dataContainer, const BinningStruct& binning, double& jetptMin, double& jetptMax) {
    //
    // Open output file
    TFile* fOutput = new TFile(Form("selection_efficiency_run3style_%d_to_%d_jetpt_" + binning.dataPeriod + ".root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    
    // Response matrices
    dataContainer.response.first.Write("hResponsePrompt");
    dataContainer.response.second.Write("hResponseNonPrompt");

    // Truth kinematic efficiency histograms
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->Write();
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->Write();

    // Reco kinematic efficiency histograms
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->Write();
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->Write();

    // Run 3 style detector level selection efficiency histograms
    dataContainer.hSelectionEfficiency.first->Write();
    dataContainer.hSelectionEfficiency.second->Write();

    // Run 3 style particle level selection efficiency histograms
    dataContainer.hEfficiency_prompt_part_run3style->Write("selEffRun3Style_particleLevelPrompt");
    dataContainer.hEfficiency_nonprompt_part_run3style->Write("selEffRun3Style_particleLevelNonprompt");

    // Efficiency corrected data
    dataContainer.hEfficiencyCorrected.second->Write();

    // Also store the axes used for the histograms
    storeBinningInFile(fOutput, binning);

    fOutput->Close();
    delete fOutput;

    std::cout << "Data stored.\n";
}



void Efficiency_run3_style_detector_level(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Luminosity (for now arbitrary)
    double luminosity = 0;
    double luminosity_powheg = 0;
    double BR = 0.0623; // D0 -> KPi decay channel branching ratio = (6.23 +- 0.33) %

    // Load binning information
    TFile* fBinning = new TFile(Form("../1-SignalTreatment/Reflections/binningInfo.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    // Load binning from background subtracted file
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

    // Opening files
    TFile* fEffRun2Style = new TFile(Form("run2_style_efficiency_%d_to_%d_jetpt_%s.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax), binning.dataPeriod.Data()),"read");
    if (!fEffRun2Style || fEffRun2Style->IsZombie()) {
        std::cerr << "Error: Unable to open run 2 style efficiency data ROOT file." << std::endl;
    }
    TFile* fSimulated = new TFile("../" + binning.inputMC.second + "/AO2D_mergedDFs.root","read","read");
    TFile* fBackSub = new TFile(Form("../1-SignalTreatment/SideBand/full_merged_ranges_back_sub_%s.root",binning.dataPeriod.Data()),"read");
    if (!fSimulated || fSimulated->IsZombie()) {
        std::cerr << "Error: Unable to open full (not matched) simulated data ROOT file." << std::endl;
    }
    if (!fBackSub || fBackSub->IsZombie()) {
        std::cerr << "Error: Unable to open background subtracted data ROOT file." << std::endl;
    }
    

    EfficiencyData dataContainer = createHistograms(binning); // pT histograms
    
    // Fill matched histograms for corrections
    fillHistograms(fSimulated, fEffRun2Style, dataContainer, binning);

    // Calculate efficiency distribution
    performEfficiencyCorrection(fBackSub, dataContainer, jetptMin, jetptMax);

    // Plot the efficiency histogram and further corrected histograms
    PlotHistograms(dataContainer, jetptMin, jetptMax, binning);

    // Save corrected distributions to file
    SaveData(dataContainer, binning, jetptMin, jetptMax);

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    std::time_t now = std::time(nullptr);
    std::cout << "Finished at: " << std::ctime(&now);
}

int main(){
    Efficiency_run3_style_detector_level();
    return 0;
}
