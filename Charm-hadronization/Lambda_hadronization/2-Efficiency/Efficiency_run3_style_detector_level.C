/**
 * Lambda_c hadron analysis
 * @file Efficiency_run3_style_detector_level.C
 * @brief Calculates pT dependent efficiency of HF selections using methodology from run 3
 * Input: 
 *      - series of backSub_%.jetptMin_to_%.jetptMax_jetpt_with_reflections.root files -> contain histograms such that the background+reflections contribution was removed
 *      - run2_style_efficiency_%.jetptMin_to_%.jetptMax_jetpt.root
 * Outputs: one run3_style_efficiency_%.jetptMin_to_%.jetptMax_jetpt.root -> contain pT,HF dependent run 2 style efficiency, to be used by Efficiency_run3_Style.C
 * 
 * @author: Christian Reckziegel
 * Date: February 2026
 */
#include "../commonUtilities.h"

using namespace std;

struct EfficiencyData {
    // Yield 2D histograms: pT,jet vs pT,D: first = prompt D0 data distribution, second = non-prompt D0 distribution
    std::pair<TH2D*, TH2D*> hYieldTruth; // all particle level entries (denominator): will be kinematically corrected and folded
    std::pair<TH2D*, TH2D*> hYieldMeasured; // detector level entries (numerator): went over smearing effects and passed the selection cuts

    // Response matrices
    std::pair<RooUnfoldResponse, RooUnfoldResponse> response; // first = prompt #Lambda_{c}^{+}, second = non-prompt #Lambda_{c}^{+}
    std::pair<std::vector<TH2D*>, std::vector<TH2D*>> responseProjections; // response projections matrix: first = prompt #Lambda_{c}^{+}, second = non-prompt #Lambda_{c}^{+}

    //
    // Kinematic efficiency histograms: first = prompt #Lambda_{c}^{+}, second = non-prompt #Lambda_{c}^{+}
    //
    // Response range
    std::pair<TH2D*, TH2D*> hKEffResponseParticle; // response range, particle level entries
    std::pair<TH2D*, TH2D*> hKEffResponseDetector; // response range, detector level entries
    // Response range / total particle range
    std::pair<TH2D*, TH2D*> hKEffTruthTotalParticle; // total truth range, particle level entries
    std::pair<TH2D*, TH2D*> hKEffResponseParticle_Over_TotalParticle; // first = prompt #Lambda_{c}^{+}, second = non-prompt #Lambda_{c}^{+}
    // Response range / total detector range
    std::pair<TH2D*, TH2D*> hKEffRecoTotalDetector; // total reco range, detector level entries
    std::pair<TH2D*, TH2D*> hKEffResponseDetector_Over_TotalDetector; // first = prompt #Lambda_{c}^{+}, second = non-prompt #Lambda_{c}^{+}

    // All particle level entries with corrections applied (kinematically corrected and folded)
    std::pair<TH2D*, TH2D*> hYieldTruthCorrected;

    // Final distributions to be treated and divided at the end: first = prompt #Lambda_{c}^{+}, second = non-prompt #Lambda_{c}^{+}
    std::pair<TH2D*, TH2D*> hNumerator; // detector level entries: went over smearing effects and passed the selection cuts
    std::pair<TH2D*, TH2D*> hDenominator; // all particle level entries: will be kinematically corrected and folded

    // Denominator original and corrected pT histograms: projection of the 2D histograms
    std::pair<TH1D*, TH1D*> hHfPtYieldTruth;
    std::pair<TH1D*, TH1D*> hHfPtYieldTruthCorrected;

    //
    // Final efficiency histograms (numerator / denominator): first = prompt #Lambda_{c}^{+}, second = non-prompt #Lambda_{c}^{+}
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

    // Create 2D histograms for prompt #Lambda_{c}^{+}: raw and folded
    dataContainer.hYieldTruth.first = new TH2D("hYieldTruthPrompt", "Particle level data prompt yield distribution;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hYieldTruth.first->Sumw2();
    dataContainer.hYieldMeasured.first = new TH2D("hYieldMeasuredPrompt", "Detector level data prompt yield distribution;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{det} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    // Create 2D histograms for non-prompt #Lambda_{c}^{+}: raw and folded
    dataContainer.hYieldTruth.second = new TH2D("hYieldTruthNonPrompt", "Particle level data non-prompt yield distribution;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hYieldTruth.second->Sumw2();
    dataContainer.hYieldMeasured.second = new TH2D("hYieldMeasuredNonPrompt", "Detector level data non-prompt yield distribution;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{det} (GeV/#it{c})",
                                                binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());

    // Create RooUnfoldResponse objects for prompt and non-prompt #Lambda_{c}^{+}
    dataContainer.response.first = RooUnfoldResponse(dataContainer.hYieldMeasured.first, dataContainer.hYieldTruth.first); // prompt #Lambda_{c}^{+}
    dataContainer.response.second = RooUnfoldResponse(dataContainer.hYieldMeasured.second, dataContainer.hYieldTruth.second); // non-prompt #Lambda_{c}^{+}
    
    // Create projections of response matrix object for prompt and non-prompt #Lambda_{c}^{+}
    dataContainer.responseProjections.first.push_back(new TH2D("responseProjectionsPtJetPrompt", "Prompt #Lambda_{c}^{+}'s reponse matrix #it{p}_{T, ch. jet} projection;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, ch. jet}^{gen} (GeV/#it{c})", binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data()));
    dataContainer.responseProjections.first.push_back(new TH2D("responseProjectionsPtDPrompt", "Prompt #Lambda_{c}^{+}'s reponse matrix p_{T,#Lambda_{c}^{+}} projection;#it{p}_{T, #Lambda_{c}^{+}}^{det} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{gen} (GeV/#it{c})", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data(), binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data()));
    dataContainer.responseProjections.second.push_back(new TH2D("responseProjectionsPtJetNonPrompt", "Non-prompt #Lambda_{c}^{+}'s reponse matrix #it{p}_{T, ch. jet} projection;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, ch. jet}^{gen} (GeV/#it{c})", binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data()));
    dataContainer.responseProjections.second.push_back(new TH2D("responseProjectionsPtDNonPrompt", "Non-prompt #Lambda_{c}^{+}'s reponse matrix p_{T,#Lambda_{c}^{+}} projection;#it{p}_{T, #Lambda_{c}^{+}}^{det} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{gen} (GeV/#it{c})", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data(), binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data()));
    
    // Creating investigation histogram
    dataContainer.hBDTBackgroundScore = new TH1D("hBDTBackgroundScore", "Entries that didn't pass the cuts;BDT background score;Counts", 100, 0, 1);

    // Kinematic efficiency histograms: truth entries
    dataContainer.hKEffResponseParticle.first = new TH2D("hKEffResponseParticlePrompt", "Truth prompt #Lambda_{c}^{+}'s within response range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hKEffResponseParticle.second = new TH2D("hKEffResponseParticleNonPrompt", "Truth non-prompt #Lambda_{c}^{+}'s within response range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hKEffTruthTotalParticle.first = new TH2D("hKEffTruthTotalParticlePrompt", "Truth prompt #Lambda_{c}^{+}'s within total particle range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{gen} (GeV/#it{c})",
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hKEffTruthTotalParticle.second = new TH2D("hKEffTruthTotalParticleNonPrompt", "Truth non-prompt #Lambda_{c}^{+}'s within total particle range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{gen} (GeV/#it{c})",
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    // Kinematic efficiency histograms: reco entries
    dataContainer.hKEffResponseDetector.first = new TH2D("hKEffResponseDetectorPrompt", "Reconstructed prompt #Lambda_{c}^{+}'s within response range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, #Lambda_{c}^{+}}^{reco}", 
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    dataContainer.hKEffResponseDetector.second = new TH2D("hKEffResponseDetectorNonPrompt", "Reconstructed non-prompt #Lambda_{c}^{+}'s within response range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, #Lambda_{c}^{+}}^{reco}", 
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    dataContainer.hKEffRecoTotalDetector.first = new TH2D("hKEffRecoTotalDetectorPrompt", "Reconstructed prompt #Lambda_{c}^{+}'s within total detector range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, #Lambda_{c}^{+}}^{reco}",
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    dataContainer.hKEffRecoTotalDetector.second = new TH2D("hKEffRecoTotalDetectorNonPrompt", "Reconstructed non-prompt #Lambda_{c}^{+}'s within total detector range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, #Lambda_{c}^{+}}^{reco}",
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());

    // Final histograms to be treated and divided at the end
    dataContainer.hNumerator.first = new TH2D("hNumeratorPrompt", "Reconstructed prompt #Lambda_{c}^{+}'s (after selection cuts);#it{p}_{T, ch. jet}^{reco};#it{p}_{T, #Lambda_{c}^{+}}^{reco}", 
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    dataContainer.hNumerator.second = new TH2D("hNumeratorNonPrompt", "Reconstructed non-prompt #Lambda_{c}^{+}'s (after selection cuts);#it{p}_{T, ch. jet}^{reco};#it{p}_{T, #Lambda_{c}^{+}}^{reco}", 
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    dataContainer.hDenominator.first = new TH2D("hDenominatorPrompt", "All truth prompt #Lambda_{c}^{+}'s;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{gen} (GeV/#it{c})", 
            binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
            binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    dataContainer.hDenominator.second = new TH2D("hDenominatorNonPrompt", "All truth non-prompt #Lambda_{c}^{+}'s;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, #Lambda_{c}^{+}}^{gen} (GeV/#it{c})", 
            binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
            binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    
    std::cout << "Histograms created.\n";

    return dataContainer;
}
//_____________________________________________________________________________________________________________________________________________________________________________________


//_____________________________________________________________________________________________________________________________________________________________________________________
// Modules to fill histograms (of matched data for response matrix and non-matched data for numerator and denominator distributions) from TFile data
void fillHistograms(TFile* fSimulated, TFile* fEffRun2Style, EfficiencyData& dataContainer, const BinningStruct& binning) {
    
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
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY, MCPjetNConst;
    bool MCPhfprompt;
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
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2)); // or use MCPaxisDistance
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2)); // or use MCDaxisDistance

        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);
        bool genHfPtRange = (MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut);
        bool genDeltaRRange = (MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < MCPDeltaRcut);
        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < MCPDeltaRcut)); // currently used
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = (MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut);
        bool recoDeltaRRange = (MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < MCDDeltaRcut);
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < MCDDeltaRcut)); // currently used

        // Get the threshold for this pT range: TODO - do NOT erase this, BDT cuts will be included in next hyperloop wagon
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;

        // 1 --- Selection efficiencies
        // Fill histograms considering jet pT and ALICE's acceptance
        if (genLevelRange) {
            
            // Fill inclusive histogram
            //dataContainer.hMcpPt[0]->Fill(MCPhfPt);
            // fill prompt efficiency histogram
            if (MCPhfprompt) {
                dataContainer.hYieldTruth.first->Fill(MCPjetPt, MCPhfPt);
            } else{
                // fill non-prompt efficiency histogram
                dataContainer.hYieldTruth.second->Fill(MCPjetPt, MCPhfPt);
            }
        }
        // only compute matched detector level candidates, but compute all particle level ones
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isReflection = (MCDhfMatchedFrom != MCDhfSelectedAs) ? true : false;
        if (!MCDhfmatch || isReflection) { // !MCDhfmatch || isReflection
            continue;
        }
        // Fill histograms considering jet pT and detector acceptance
        if (recoLevelRange) {
            
            // Get the threshold for this pT range
            double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);

            // Fill histogram only if the BDT cut is passed
            if (MCDhfMlScore0 < maxBkgProb) {
                // Fill inclusive histogram
                //dataContainer.hMcdPt[0]->Fill(MCDhfPt);
                // fill prompt efficiency histogram
                if (MCDhfprompt) {
                    dataContainer.hYieldMeasured.first->Fill(MCPjetPt, MCPhfPt);
                } else{
                    // fill non-prompt efficiency histogram
                    dataContainer.hYieldMeasured.second->Fill(MCPjetPt, MCPhfPt);
                }
            } else {
                dataContainer.hBDTBackgroundScore->Fill(MCDhfMlScore0);
            } 
        }

        // 2 --- Response matrix: fill histograms considering jet pT and detector acceptance (response range)
        if (genLevelRange && recoLevelRange) {

            // fill 2D yields histograms
            if (MCPhfprompt) {

                // Get efficiency estimate to weight the response matrix
                hEffWeight = (TH1D*) fEffRun2Style->Get("efficiency_prompt");
                // Get efficiency estimate to weight the response matrix: jet pT shape is influenced by D0 pT efficiency
                int bin = hEffWeight->FindBin(MCDhfPt);
                double estimatedEfficiency = hEffWeight->GetBinContent(bin);
                if (estimatedEfficiency == 0) {
                    std::cout << "Warning: Prompt efficiency is zero for pT,HF = " << MCDhfPt << " GeV/c with bin " << bin << " of efficiency_prompt run 2 histogram. Setting it to 1." << std::endl;
                    estimatedEfficiency = 1; // Avoid division by zero
                }

                // Fill 4D RooUnfoldResponse object
                dataContainer.response.first.Fill(MCDjetPt, MCDhfPt, MCPjetPt, MCPhfPt, 1 / estimatedEfficiency); // jet pT shape is influenced by D0 pT efficiency
                dataContainer.responseProjections.first[0]->Fill(MCDjetPt, MCPjetPt,1 / estimatedEfficiency); // pT,jet projection, prompt #Lambda_{c}^{+}
                dataContainer.responseProjections.first[1]->Fill(MCDhfPt, MCPhfPt, 1 / estimatedEfficiency); // pT,D projection, prompt #Lambda_{c}^{+}
            } else{

                // Get efficiency estimate to weight the response matrix
                hEffWeight = (TH1D*) fEffRun2Style->Get("efficiency_nonprompt");
                // Get efficiency estimate to weight the response matrix: jet pT shape is influenced by D0 pT efficiency
                int bin = hEffWeight->FindBin(MCDhfPt);
                double estimatedEfficiency = hEffWeight->GetBinContent(bin);
                if (estimatedEfficiency == 0) {
                    //std::cout << "Warning: Prompt efficiency is zero for pT = " << MCDhfPt << " with bin " << bin << ". Setting it to 1." << std::endl;
                    estimatedEfficiency = 1; // Avoid division by zero
                }

                // Fill 4D RooUnfoldResponse object
                dataContainer.response.second.Fill(MCDjetPt, MCDhfPt, MCPjetPt, MCPhfPt, 1 / estimatedEfficiency); // jet pT shape is influenced by D0 pT efficiency
                dataContainer.responseProjections.second[0]->Fill(MCDjetPt, MCPjetPt,1 / estimatedEfficiency); // pT,jet projection, non-prompt #Lambda_{c}^{+}
                dataContainer.responseProjections.second[1]->Fill(MCDhfPt, MCPhfPt, 1 / estimatedEfficiency); // pT,D projection, non-prompt #Lambda_{c}^{+}
            }
        }
        // 3 --- Kinematic efficiencies
        if (genLevelRange && MCPhfprompt) {
            // fill prompt 2D yield: total particle range
            dataContainer.hKEffTruthTotalParticle.first->Fill(MCPjetPt, MCPhfPt);
            if (recoLevelRange) {
                // fill prompt 2D yield: response range
                dataContainer.hKEffResponseParticle.first->Fill(MCPjetPt, MCPhfPt);
                
            }
        // non-prompt #Lambda_{c}^{+}    
        } else if (genLevelRange && !MCPhfprompt) {
            dataContainer.hKEffTruthTotalParticle.second->Fill(MCPjetPt, MCPhfPt);
            if (recoLevelRange) {
                // fill prompt 2D yield: particle level matched
                dataContainer.hKEffResponseParticle.second->Fill(MCPjetPt, MCPhfPt);
                
            }
        }
        
        // Kinematic efficiency of reconstruction entries
        if (recoLevelRange && MCPhfprompt) {
            // fill prompt 2D yield: total detector range
            dataContainer.hKEffRecoTotalDetector.first->Fill(MCDjetPt, MCDhfPt);
            if (genLevelRange) {
                // fill prompt 2D yield: response range
                dataContainer.hKEffResponseDetector.first->Fill(MCDjetPt, MCDhfPt);
                
            }
        // non-prompt #Lambda_{c}^{+}    
        } else if (recoLevelRange && !MCPhfprompt) {
            dataContainer.hKEffRecoTotalDetector.second->Fill(MCDjetPt, MCDhfPt);
            if (genLevelRange) {
                // fill prompt 2D yield: detector level matched
                dataContainer.hKEffResponseDetector.second->Fill(MCDjetPt, MCDhfPt);
                
            }
        }
    }

    
    std::cout << "pT,HF histograms (particle and detector level, prompt and non-prompt), response matrix and kinematic correction histograms filled.\n";

}
//_____________________________________________________________________________________________________________________________________________________________________________________


void performEfficiencyCorrection(TFile* fBackSub, EfficiencyData& dataContainer, double jetptMin, double jetptMax) {
    
    //
    // Calculating pT,D dependent efficiency distribution
    //
    
    // 0.1 - Calculate truth kinematic efficiency
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first = (TH2D*)dataContainer.hKEffResponseParticle.first->Clone("hKEffTruthPrompt");
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second = (TH2D*)dataContainer.hKEffResponseParticle.second->Clone("hKEffTruthNonPrompt");
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->Sumw2(); // necessary for correct error propagation
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->Sumw2();
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->Divide(dataContainer.hKEffTruthTotalParticle.first);
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->Divide(dataContainer.hKEffTruthTotalParticle.second);
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->SetTitle("#varepsilon_{kinematic} of particle level prompt #Lambda_{c}^{+};#it{p}_{T, ch. jet};#it{p}_{T, #Lambda_{c}^{+}}");
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->SetTitle("#varepsilon_{kinematic} of particle level non-prompt #Lambda_{c}^{+};#it{p}_{T, ch. jet};#it{p}_{T, #Lambda_{c}^{+}}");
    
    // 0.2 - Calculate reco kinematic efficiency
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first = (TH2D*)dataContainer.hKEffResponseDetector.first->Clone("hKEffRecoPrompt");
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second = (TH2D*)dataContainer.hKEffResponseDetector.second->Clone("hKEffRecoNonPrompt");
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->Sumw2();
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->Sumw2();
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->Divide(dataContainer.hKEffRecoTotalDetector.first);
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->Divide(dataContainer.hKEffRecoTotalDetector.second);
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->SetTitle("#varepsilon_{kinematic} of detector level prompt #Lambda_{c}^{+};#it{p}_{T, ch. jet};#it{p}_{T, #Lambda_{c}^{+}}");
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->SetTitle("#varepsilon_{kinematic} of detector level non-prompt #Lambda_{c}^{+};#it{p}_{T, ch. jet};#it{p}_{T, #Lambda_{c}^{+}}");

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
    dataContainer.hHfPtYieldTruth.first->SetTitle("Denominator prompt #Lambda_{c}^{+} p_{T} distribution before correction; #it{p}_{T, #Lambda_{c}^{+}}^{gen}; Counts");
    minBin = dataContainer.hYieldTruth.second->GetXaxis()->FindBin(jetptMin);
    maxBin = dataContainer.hYieldTruth.second->GetXaxis()->FindBin(jetptMax - 1e-6);
    dataContainer.hHfPtYieldTruth.second = dataContainer.hYieldTruth.second->ProjectionY("hHfPtYieldTruthNonPrompt", minBin, maxBin);
    dataContainer.hHfPtYieldTruth.second->SetTitle("Denominator non-prompt #Lambda_{c}^{+} p_{T} distribution before correction; #it{p}_{T, #Lambda_{c}^{+}}^{gen}; Counts");
    minBin = dataContainer.hYieldTruthCorrected.first->GetXaxis()->FindBin(jetptMin);
    maxBin = dataContainer.hYieldTruthCorrected.first->GetXaxis()->FindBin(jetptMax - 1e-6);
    dataContainer.hHfPtYieldTruthCorrected.first = dataContainer.hYieldTruthCorrected.first->ProjectionY("hHfPtYieldTruthCorrectedPrompt", minBin, maxBin);
    dataContainer.hHfPtYieldTruthCorrected.first->SetTitle("Denominator prompt #Lambda_{c}^{+} p_{T} distribution after correction; #it{p}_{T, #Lambda_{c}^{+}}^{gen}; Counts");
    minBin = dataContainer.hYieldTruthCorrected.second->GetXaxis()->FindBin(jetptMin);
    maxBin = dataContainer.hYieldTruthCorrected.second->GetXaxis()->FindBin(jetptMax - 1e-6);
    dataContainer.hHfPtYieldTruthCorrected.second = dataContainer.hYieldTruthCorrected.second->ProjectionY("hHfPtYieldTruthCorrectedNonPrompt", minBin, maxBin);
    dataContainer.hHfPtYieldTruthCorrected.second->SetTitle("Denominator non-prompt #Lambda_{c}^{+} p_{T} distribution after correction; #it{p}_{T, #Lambda_{c}^{+}}^{gen}; Counts");

    // 4 - Calculate the efficiency distributions with the pT projections
    dataContainer.hSelectionEfficiency.first = dataContainer.hYieldMeasured.first->ProjectionY("hSelectionEfficiencyPrompt", minBin, maxBin);
    dataContainer.hSelectionEfficiency.second = dataContainer.hYieldMeasured.second->ProjectionY("hSelectionEfficiencyNonPrompt", minBin, maxBin);
    dataContainer.hSelectionEfficiency.first->SetTitle("Efficiency prompt #Lambda_{c}^{+} p_{T} distribution; #it{p}_{T, #Lambda_{c}^{+}}^{reco}; Efficiency#times Acceptance");
    dataContainer.hSelectionEfficiency.second->SetTitle("Efficiency non-prompt #Lambda_{c}^{+} p_{T} distribution; #it{p}_{T, #Lambda_{c}^{+}}^{reco}; Efficiency#times Acceptance");
    dataContainer.hSelectionEfficiency.first->Sumw2();
    dataContainer.hSelectionEfficiency.second->Sumw2();
    dataContainer.hSelectionEfficiency.first->Divide(dataContainer.hHfPtYieldTruthCorrected.first);
    dataContainer.hSelectionEfficiency.second->Divide(dataContainer.hHfPtYieldTruthCorrected.second);
    
    //
    // Correcting 3D background subtracted distribution with efficiency for each pT,D bin
    //
    // 1 - Clone background subtracted distribution
    dataContainer.hEfficiencyCorrected.first = (TH3D*)fBackSub->Get("h3DBackgroundSubtracted")->Clone("h3DEfficiencyCorrected");
    dataContainer.hEfficiencyCorrected.first->SetTitle("Background subtracted, efficiency corrected");
    int xBins = dataContainer.hEfficiencyCorrected.first->GetXaxis()->GetNbins();
    int yBins = dataContainer.hEfficiencyCorrected.first->GetYaxis()->GetNbins();
    int zBins = dataContainer.hEfficiencyCorrected.first->GetZaxis()->GetNbins();
    for (int xBin = 1; xBin <= xBins; xBin++) {
        for (int yBin = 1; yBin <= yBins; yBin++) {
            for (int zBin = 1; zBin <= zBins; zBin++) { // pT,D bins will be corrected by 1/prompt_efficiency
                
                // Get current content and error
                double content = dataContainer.hEfficiencyCorrected.first->GetBinContent(xBin, yBin, zBin);
                double error = dataContainer.hEfficiencyCorrected.first->GetBinError(xBin, yBin, zBin);

                // Get the pT,D bin center from the Z axis of the 3D histogram
                double ptDcenter = dataContainer.hEfficiencyCorrected.first->GetZaxis()->GetBinCenter(zBin);
                
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
                    std::cerr << "Zero or invalid efficiency at ptD = " << ptDcenter
                            << " (eff bin = " << effBin << ")" << std::endl;
                }
            }
        }
    }
    
    // 2 - Project the z axis (pT,D) in order to obtain a 2D distribution of pT,jet (x axis) vs DeltaR (y axis)
    TH2D* h2D = (TH2D*)dataContainer.hEfficiencyCorrected.first->Project3D("yx");
    dataContainer.hEfficiencyCorrected.second = (TH2D*)h2D->Clone("h2DEfficiencyCorrected");
    delete h2D;
    dataContainer.hEfficiencyCorrected.second->SetTitle("Background subtracted, efficiency corrected");

    std::cout << "Efficiency correction performed.\n";
}

void PlotHistograms(const EfficiencyData& dataContainer, double jetptMin, double jetptMax) {

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
    TCanvas* cResponse = new TCanvas("cResponse","Response matrices for all pT,D bins");
    cResponse->SetCanvasSize(1800,1000);
    cResponse->Divide(2,2);
    std::pair<const TH2*, const TH2*> hresponse2D;
    hresponse2D.first = dataContainer.response.first.Hresponse();
    hresponse2D.second = dataContainer.response.second.Hresponse();
    std::pair<TH2D*, TH2D*> hresponse2DClone;
    hresponse2DClone.first = static_cast<TH2D*>(hresponse2D.first->Clone("hResponse2DPrompt"));
    hresponse2DClone.second = static_cast<TH2D*>(hresponse2D.second->Clone("hResponse2DNonPrompt"));
    hresponse2DClone.first->SetTitle("2D response matrix from 4D RooUnfoldResponse - prompt #Lambda_{c}^{+}'s;2D Reconstructed;2D Truth");
    hresponse2DClone.second->SetTitle("2D response matrix from 4D RooUnfoldResponse - non-prompt #Lambda_{c}^{+}'s;2D Reconstructed;2D Truth");
    cResponse->cd(1);
    hresponse2DClone.first->Draw("colz");
    cResponse->cd(2);
    hresponse2DClone.second->Draw("colz");

    // Separate response matrices in two canvases
    TCanvas* cResponsePrompt = new TCanvas("cResponsePrompt","Response matrix for prompt #Lambda_{c}^{+}'s");
    cResponsePrompt->SetCanvasSize(1800,1000);
    cResponsePrompt->cd();
    hresponse2DClone.first->Draw("colz");
    TCanvas* cResponseNonPrompt = new TCanvas("cResponseNonPrompt","Response matrix for non-prompt #Lambda_{c}^{+}'s");
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

    TCanvas* cSelectionEfficiency = new TCanvas("cSelectionEfficiency","Efficiency histograms");
    cSelectionEfficiency->SetCanvasSize(1800,1000);
    cSelectionEfficiency->Divide(2,2);
    cSelectionEfficiency->cd(1);
    dataContainer.hSelectionEfficiency.first->Draw();
    cSelectionEfficiency->cd(2);
    dataContainer.hSelectionEfficiency.second->Draw();
    // Separating in two canvases
    TCanvas* cSelectionEfficiencyPrompt = new TCanvas("cSelectionEfficiencyPrompt","Efficiency histograms for prompt #Lambda_{c}^{+}'s");
    cSelectionEfficiencyPrompt->SetCanvasSize(1800,1000);
    cSelectionEfficiencyPrompt->cd();
    dataContainer.hSelectionEfficiency.first->Draw();
    TCanvas* cSelectionEfficiencyNonPrompt = new TCanvas("cSelectionEfficiencyNonPrompt","Efficiency histograms for non-prompt #Lambda_{c}^{+}'s");
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
    dataContainer.hEfficiencyCorrected.second->ProjectionY("hEfficiencyCorrectedProjectionDeltaRPrompt")->Draw();

    TCanvas* cCorrectedData2DOnly = new TCanvas("cCorrectedData2DOnly","Background subtracted, efficiency corrected 2D data");
    cCorrectedData2DOnly->cd();
    dataContainer.hEfficiencyCorrected.second->Draw("colz");
    TCanvas* cCorrectedData1DOnly = new TCanvas("cCorrectedData1DOnly","Background subtracted, efficiency corrected 1D data");
    cCorrectedData1DOnly->cd();
    //dataContainer.hEfficiencyCorrected.second->SetMarkerStyle(kCircle);
    dataContainer.hEfficiencyCorrected.second->ProjectionY("hEfficiencyCorrectedProjectionDeltaRPrompt")->Draw();

    //
    // Storing in a single pdf file
    //
    TString imagePath = "../Images/2-Efficiency/Run3Style/";
    cYield->Print(imagePath + Form("run3_style_efficiency_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cResponse->Print(imagePath + Form("run3_style_efficiency_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponseProjections->Print(imagePath + Form("run3_style_efficiency_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cKEffAll->Print(imagePath + Form("run3_style_efficiency_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cYieldTruthCorrected->Print(imagePath + Form("run3_style_efficiency_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cPtCorrectionComparison->Print(imagePath + Form("run3_style_efficiency_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cSelectionEfficiency->Print(imagePath + Form("run3_style_efficiency_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cCorrectedData->Print(imagePath + Form("run3_style_efficiency_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
    // in separate .pdf files
    cKEffAll->Print(imagePath + "Efficiency_run3style_kinematic_efficiency.pdf");

    std::cout << "Histograms plotted.\n";

}

void SaveData(const EfficiencyData& dataContainer, const BinningStruct& binning, double& jetptMin, double& jetptMax) {
    //
    // Open output file
    TFile* fOutput = new TFile(Form("selection_efficiency_run3style_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    
    // Response matrices
    dataContainer.response.first.Write("hResponsePrompt");
    dataContainer.response.second.Write("hResponseNonPrompt");

    // Truth kinematic efficiency histograms
    dataContainer.hKEffResponseParticle_Over_TotalParticle.first->Write();
    dataContainer.hKEffResponseParticle_Over_TotalParticle.second->Write();

    // Reco kinematic efficiency histograms
    dataContainer.hKEffResponseDetector_Over_TotalDetector.first->Write();
    dataContainer.hKEffResponseDetector_Over_TotalDetector.second->Write();

    // Selection efficiency histograms
    dataContainer.hSelectionEfficiency.first->Write();
    dataContainer.hSelectionEfficiency.second->Write();

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
    TFile* fBinning = new TFile(Form("../1-SignalTreatment/Sideband/full_merged_ranges_back_sub.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    // Load binning from background subtracted file
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

    // Opening files
    TFile* fEffRun2Style = new TFile(Form("run2_style_efficiency_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fEffRun2Style || fEffRun2Style->IsZombie()) {
        std::cerr << "Error: Unable to open run 2 style efficiency data ROOT file." << std::endl;
    }
    TFile* fSimulated = new TFile("../Data/MonteCarlo/Train_607573/AO2D_mergedDFs.root","read","read");
    TFile* fBackSub = new TFile(Form("../1-SignalTreatment/Sideband/full_merged_ranges_back_sub.root"),"read");
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
    PlotHistograms(dataContainer, jetptMin, jetptMax);

    // Save corrected distributions to file
    SaveData(dataContainer, binning, jetptMin, jetptMax);

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    Efficiency_run3_style_detector_level();
    return 0;
}
