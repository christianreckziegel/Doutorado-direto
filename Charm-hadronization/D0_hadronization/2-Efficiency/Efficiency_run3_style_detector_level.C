/*
 *
 *
 * Macro for plotting HF pT dependent efficiency 
 * and applying correction to BackgroundSubtraction.C
 * and SignalExtraction.C resulting distributions.
 * 
 * 
 * 
 * 
**/


using namespace std;

// calculate number of objects inside file
int HistogramCounter(TFile* file) {
    // Get the list of keys (i.e., objects) in the file
    TList* keys = file->GetListOfKeys();

    int numHistograms = 0;

    // Loop over all keys in the file
    for (int i = 0; i < keys->GetSize(); ++i) {
        TKey* key = (TKey*)keys->At(i);
        TObject* obj = file->Get(key->GetName());

        // Check if the object is a histogram
        if (obj->IsA()->InheritsFrom(TH1::Class())) {
            numHistograms++;
        }
    }

    std::cout << numHistograms << " histograms found in file.\n";
    return numHistograms;
}

// calculate delta phi such that 0 < delta phi < 2*pi
double DeltaPhi(double phi1, double phi2) {
    // Compute the absolute difference between phi1 and phi2
    double dphi = std::abs(phi1 - phi2); 
    if (dphi > M_PI) {
        // subtract 2pi if the difference if bigger than pi
        dphi = dphi - 2*M_PI;
    }

    return dphi;
}

struct EfficiencyData {
    // Yield 2D histograms: pT,jet vs pT,D: first = prompt D0 data distribution, second = non-prompt D0 distribution
    std::pair<TH2D*, TH2D*> hYieldTruth; // all particle level entries (denominator): will be kinematically corrected and folded
    std::pair<TH2D*, TH2D*> hYieldMeasured; // detector level entries (numerator): went over smearing effects and passed the selection cuts

    // Response matrices
    std::pair<RooUnfoldResponse, RooUnfoldResponse> response; // first = prompt D0s, second = non-prompt D0s
    std::pair<std::vector<TH2D*>, std::vector<TH2D*>> responseProjections; // response projections matrix: first = prompt D0s, second = non-prompt D0s

    //
    // Kinematic efficiency histograms: first = prompt D0s, second = non-prompt D0s
    //
    // Response range
    std::pair<TH2D*, TH2D*> hKEffResponseParticle; // response range, particle level entries
    std::pair<TH2D*, TH2D*> hKEffResponseDetector; // response range, detector level entries
    // Response range / total particle range
    std::pair<TH2D*, TH2D*> hKEffTruthTotalParticle; // total truth range, particle level entries
    std::pair<TH2D*, TH2D*> hKEffResponseParticle_Over_TotalParticle; // first = prompt D0s, second = non-prompt D0s
    // Response range / total detector range
    std::pair<TH2D*, TH2D*> hKEffRecoTotalDetector; // total reco range, detector level entries
    std::pair<TH2D*, TH2D*> hKEffResponseDetector_Over_TotalDetector; // first = prompt D0s, second = non-prompt D0s

    // All particle level entries with corrections applied (kinematically corrected and folded)
    std::pair<TH2D*, TH2D*> hYieldTruthCorrected;

    // Final distributions to be treated and divided at the end: first = prompt D0s, second = non-prompt D0s
    std::pair<TH2D*, TH2D*> hNumerator; // detector level entries: went over smearing effects and passed the selection cuts
    std::pair<TH2D*, TH2D*> hDenominator; // all particle level entries: will be kinematically corrected and folded

    // Denominator original and corrected pT histograms: projection of the 2D histograms
    std::pair<TH1D*, TH1D*> hHfPtYieldTruth;
    std::pair<TH1D*, TH1D*> hHfPtYieldTruthCorrected;

    //
    // Final efficiency histograms (numerator / denominator): first = prompt D0s, second = non-prompt D0s
    //
    std::pair<TH1D*, TH1D*> hSelectionEfficiency; // efficiency = numerator / denominator

    // Investigation histograms
    TH1D* hBDTBackgroundScore;

    // Efficiency corrected data (after background subtraction and now efficiency)
    std::pair<TH3D*, TH2D*> hEfficiencyCorrected;
};

// Module to create histograms including interest variable
EfficiencyData createHistograms(const std::vector<double>& ptjetBinEdges_particle, const std::vector<double>& ptDBinEdges_particle,
                                const std::vector<double>& ptjetBinEdges_detector, const std::vector<double>& ptDBinEdges_detector) {

    // Create struct to store data
    EfficiencyData histStruct;

    // Create 2D histograms for prompt D0s: raw and folded
    histStruct.hYieldTruth.first = new TH2D("hYieldTruthPrompt", "Particle level data prompt yield distribution;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data());
    histStruct.hYieldTruth.first->Sumw2();
    histStruct.hYieldMeasured.first = new TH2D("hYieldMeasuredPrompt", "Detector level data prompt yield distribution;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{det} (GeV/#it{c})", 
                                                ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data());
    // Create 2D histograms for non-prompt D0s: raw and folded
    histStruct.hYieldTruth.second = new TH2D("hYieldTruthNonPrompt", "Particle level data non-prompt yield distribution;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data());
    histStruct.hYieldTruth.second->Sumw2();
    histStruct.hYieldMeasured.second = new TH2D("hYieldMeasuredNonPrompt", "Detector level data non-prompt yield distribution;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{det} (GeV/#it{c})",
                                                ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data());

    // Create RooUnfoldResponse objects for prompt and non-prompt D0s
    histStruct.response.first = RooUnfoldResponse(histStruct.hYieldMeasured.first, histStruct.hYieldTruth.first); // prompt D0s
    histStruct.response.second = RooUnfoldResponse(histStruct.hYieldMeasured.second, histStruct.hYieldTruth.second); // non-prompt D0s

    // Create projections of response matrix object for prompt and non-prompt D0s
    histStruct.responseProjections.first.push_back(new TH2D("responseProjectionsPtJetPrompt", "Prompt D^{0}'s reponse matrix #it{p}_{T, ch. jet} projection;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, ch. jet}^{gen} (GeV/#it{c})", ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data()));
    histStruct.responseProjections.first.push_back(new TH2D("responseProjectionsPtDPrompt", "Prompt D^{0}'s reponse matrix p_{T,D^{0}} projection;#it{p}_{T, D^{0}}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data(), ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data()));
    histStruct.responseProjections.second.push_back(new TH2D("responseProjectionsPtJetNonPrompt", "Non-prompt D^{0}'s reponse matrix #it{p}_{T, ch. jet} projection;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, ch. jet}^{gen} (GeV/#it{c})", ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data()));
    histStruct.responseProjections.second.push_back(new TH2D("responseProjectionsPtDNonPrompt", "Non-prompt D^{0}'s reponse matrix p_{T,D^{0}} projection;#it{p}_{T, D^{0}}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data(), ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data()));
    
    // Creating investigation histogram
    histStruct.hBDTBackgroundScore = new TH1D("hBDTBackgroundScore", "Entries that didn't pass the cuts;BDT background score;Counts", 100, 0, 1);
    
    // Kinematic efficiency histograms: truth entries
    histStruct.hKEffResponseParticle.first = new TH2D("hKEffResponseParticlePrompt", "Truth prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data());
    histStruct.hKEffResponseParticle.second = new TH2D("hKEffResponseParticleNonPrompt", "Truth non-prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data());
    histStruct.hKEffTruthTotalParticle.first = new TH2D("hKEffTruthTotalParticlePrompt", "Truth prompt D^{0}'s within total particle range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})",
                                                ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data());
    histStruct.hKEffTruthTotalParticle.second = new TH2D("hKEffTruthTotalParticleNonPrompt", "Truth non-prompt D^{0}'s within total particle range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})",
                                                ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data());
    // Kinematic efficiency histograms: reco entries
    histStruct.hKEffResponseDetector.first = new TH2D("hKEffResponseDetectorPrompt", "Reconstructed prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
            ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data());
    histStruct.hKEffResponseDetector.second = new TH2D("hKEffResponseDetectorNonPrompt", "Reconstructed non-prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
            ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data());
    histStruct.hKEffRecoTotalDetector.first = new TH2D("hKEffRecoTotalDetectorPrompt", "Reconstructed prompt D^{0}'s within total detector range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}",
            ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
            ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data());
    histStruct.hKEffRecoTotalDetector.second = new TH2D("hKEffRecoTotalDetectorNonPrompt", "Reconstructed non-prompt D^{0}'s within total detector range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}",
            ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
            ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data());

    // Final histograms to be treated and divided at the end
    histStruct.hNumerator.first = new TH2D("hNumeratorPrompt", "Reconstructed prompt D^{0}'s (after selection cuts);#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
            ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data());
    histStruct.hNumerator.second = new TH2D("hNumeratorNonPrompt", "Reconstructed non-prompt D^{0}'s (after selection cuts);#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
            ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data());
    histStruct.hDenominator.first = new TH2D("hDenominatorPrompt", "All truth prompt D^{0}'s;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
            ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
            ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data());
    histStruct.hDenominator.second = new TH2D("hDenominatorNonPrompt", "All truth non-prompt D^{0}'s;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
            ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
            ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data());
    
    std::cout << "Histograms created.\n";

    return histStruct;
}
//_____________________________________________________________________________________________________________________________________________________________________________________
// Find the appropriate cut value based on pT
double GetBkgProbabilityCut(double pT, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    for (size_t i = 0; i < bdtPtCuts.size() - 1; ++i) {
        if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) {
            return bdtPtCuts[i].second;
        }
    }
    return 1.0; // Default: accept all if out of range
}

//_____________________________________________________________________________________________________________________________________________________________________________________
// Modules to fill histograms (of matched data for response matrix and non-matched data for numerator and denominator distributions) from TFile data
void fillMatchedHistograms(TFile* fSimulatedMCMatched, TFile* fEffRun2Style, EfficiencyData& histStruct, 
                           std::vector<double>& ptjetBinEdges_detector, std::vector<double>& ptDBinEdges_detector, std::vector<double>& deltaRBinEdges_detector,
                           std::vector<double>& ptjetBinEdges_particle, std::vector<double>& ptDBinEdges_particle, std::vector<double>& deltaRBinEdges_particle,
                           const std::vector<std::pair<double, double>>& bdtPtCuts) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = deltaRBinEdges_particle[deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = deltaRBinEdges_detector[deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double jetptMin = ptjetBinEdges_particle[0]; // on both levels jet
    const double jetptMax = ptjetBinEdges_particle[ptjetBinEdges_particle.size() - 1]; // on both levels jet
    const double MCPHfPtMincut = ptDBinEdges_particle[0]; // on particle level D0
    const double MCDHfPtMincut = ptDBinEdges_detector[0]; // on detector level D0
    const double MCPHfPtMaxcut = ptDBinEdges_particle[ptDBinEdges_particle.size() - 1]; // on particle level D0
    const double MCDHfPtMaxcut = ptDBinEdges_detector[ptDBinEdges_detector.size() - 1]; // on detector level D0

    //
    // MC generator level tree and histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)fSimulatedMCMatched->Get("DF_merged/O2matchtable");
    
    // Check for correct access
    if (!tree) {
        cout << "Error opening tree.\n";
    }

    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    int MCPjetnconst;//, MCPhfmatch;
    bool MCPhfprompt;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    int MCDjetnconst, MCDhfMatchedFrom, MCDhfSelectedAs;
    bool MCDhfprompt;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;

    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    MCPhfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&MCPhfY);
    tree->SetBranchAddress("fMcHfPrompt",&MCPhfprompt);
    //tree->SetBranchAddress("fMCHfMatch",&MCPhfmatch);
    // detector level branches
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    //tree->SetBranchAddress("fMCJetNConst",&jetnconst_float);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);
    //tree->SetBranchAddress("fHfMatch",&MCDhfmatch);
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

        // only compute matched detector level candidates, but compute all particle level ones
        bool isReflection = (MCDhfMatchedFrom != MCDhfSelectedAs) ? true : false;
        if (isReflection) {
            continue;
        }

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);
        bool genHfPtRange = (MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut);
        bool genDeltaRRange = (MCPDeltaR >= deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut);
        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPDeltaR >= 0.) && (MCPDeltaR < MCPDeltaRcut)); // currently used
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = (MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut);
        bool recoDeltaRRange = (MCDDeltaR >= deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut);
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDDeltaR >= 0.) && (MCDDeltaR < MCDDeltaRcut)); // currently used

        // Get the threshold for this pT range: TODO - do NOT erase this, BDT cuts will be included in next hyperloop wagon
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;

        // Fill histograms considering jet pT and detector acceptance (response range)
        //if (genLevelRange && recoLevelRange && passBDTcut) {
        if (genLevelRange && recoLevelRange) {

            // fill 2D yields histograms
            if (MCPhfprompt) {

                // Get efficiency estimate to weight the response matrix
                hEffWeight = (TH1D*) fEffRun2Style->Get("efficiency_prompt");
                // Get efficiency estimate to weight the response matrix: jet pT shape is influenced by D0 pT efficiency
                int bin = hEffWeight->FindBin(MCDhfPt);
                double estimatedEfficiency = hEffWeight->GetBinContent(bin);
                if (estimatedEfficiency == 0) {
                    //std::cout << "Warning: Prompt efficiency is zero for pT = " << MCDhfPt << " with bin " << bin << ". Setting it to 1." << std::endl;
                    estimatedEfficiency = 1; // Avoid division by zero
                }

                // Fill 4D RooUnfoldResponse object
                histStruct.response.first.Fill(MCDjetPt, MCDhfPt, MCPjetPt, MCPhfPt, 1 / estimatedEfficiency); // jet pT shape is influenced by D0 pT efficiency
                histStruct.responseProjections.first[0]->Fill(MCDjetPt, MCPjetPt,1 / estimatedEfficiency); // pT,jet projection, prompt D0s
                histStruct.responseProjections.first[1]->Fill(MCDhfPt, MCPhfPt, 1 / estimatedEfficiency); // pT,D projection, prompt D0s
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
                histStruct.response.second.Fill(MCDjetPt, MCDhfPt, MCPjetPt, MCPhfPt, 1 / estimatedEfficiency); // jet pT shape is influenced by D0 pT efficiency
                histStruct.responseProjections.second[0]->Fill(MCDjetPt, MCPjetPt,1 / estimatedEfficiency); // pT,jet projection, non-prompt D0s
                histStruct.responseProjections.second[1]->Fill(MCDhfPt, MCPhfPt, 1 / estimatedEfficiency); // pT,D projection, non-prompt D0s
            }

            
            
        }

        // Kinematic efficiency of truth entries
        if (genLevelRange && MCPhfprompt) {
            // fill prompt 2D yield: total particle range
            histStruct.hKEffTruthTotalParticle.first->Fill(MCPjetPt, MCPhfPt);
            if (recoLevelRange) {
                // fill prompt 2D yield: response range
                histStruct.hKEffResponseParticle.first->Fill(MCPjetPt, MCPhfPt);
                
            }
        // non-prompt D0s    
        } else if (genLevelRange && !MCPhfprompt) {
            histStruct.hKEffTruthTotalParticle.second->Fill(MCPjetPt, MCPhfPt);
            if (recoLevelRange) {
                // fill prompt 2D yield: particle level matched
                histStruct.hKEffResponseParticle.second->Fill(MCPjetPt, MCPhfPt);
                
            }
        }
        
        // Kinematic efficiency of reconstruction entries
        if (recoLevelRange && MCPhfprompt) {
            // fill prompt 2D yield: total detector range
            histStruct.hKEffRecoTotalDetector.first->Fill(MCDjetPt, MCDhfPt);
            if (genLevelRange) {
                // fill prompt 2D yield: response range
                histStruct.hKEffResponseDetector.first->Fill(MCDjetPt, MCDhfPt);
                
            }
        // non-prompt D0s    
        } else if (recoLevelRange && !MCPhfprompt) {
            histStruct.hKEffRecoTotalDetector.second->Fill(MCDjetPt, MCDhfPt);
            if (genLevelRange) {
                // fill prompt 2D yield: detector level matched
                histStruct.hKEffResponseDetector.second->Fill(MCDjetPt, MCDhfPt);
                
            }
        }
        
        
    }
    cout << "Response matrix and kinematic correction histograms filled.\n";

}

void fillNonMatchedHistograms(TFile* fSimulatedMCNonMatched, EfficiencyData& histStruct, double& jetptMin, double& jetptMax, double& hfptMin, double& hfptMax, 
                              std::vector<double>& deltaRBinEdges_particle, std::vector<double>& deltaRBinEdges_detector, const std::vector<std::pair<double, double>>& bdtPtCuts) {

    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = deltaRBinEdges_particle[deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = deltaRBinEdges_detector[deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double MCPHfPtMincut = hfptMin; // on particle level D0
    const double MCDHfPtMincut = hfptMin; // on detector level D0
    const double MCPHfPtMaxcut = hfptMax; // on particle level D0
    const double MCDHfPtMaxcut = hfptMax; // on detector level D0

    //
    // MC generator level tree and histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)fSimulatedMCNonMatched->Get("DF_merged/O2mcpjetdisttable");
    
    // Check for correct access
    if (!tree) {
        cout << "Error opening mcpjetdisttable tree.\n";
    }
    // defining variables for accessing the TTree
    float axisDistance, jetPt, jetEta, jetPhi;
    float hfPt, hfEta, hfPhi, hfMass, hfY;
    float jetnconst_float;
    bool hfprompt, hfmatch;

    tree->SetBranchAddress("fMcJetHfDist",&axisDistance);
    tree->SetBranchAddress("fMcJetPt",&jetPt);
    tree->SetBranchAddress("fMcJetEta",&jetEta);
    tree->SetBranchAddress("fMcJetPhi",&jetPhi);
    tree->SetBranchAddress("fMcJetNConst",&jetnconst_float);
    tree->SetBranchAddress("fMcHfPt",&hfPt);
    tree->SetBranchAddress("fMcHfEta",&hfEta);
    tree->SetBranchAddress("fMcHfPhi",&hfPhi);
    hfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&hfY);
    tree->SetBranchAddress("fMcHfPrompt",&hfprompt);
    tree->SetBranchAddress("fMcHfMatch",&hfmatch);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));

        bool genLevelRange = (abs(jetEta) < MCPetaCut) && (abs(hfY) < MCPyCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= deltaRBinEdges_particle[0]) && (deltaR < MCPDeltaRcut)) && ((hfPt >= MCPHfPtMincut) && (hfPt < MCPHfPtMaxcut));
        
        // Fill histograms considering jet pT and detector acceptance
        if (genLevelRange) {
            
            // fill prompt efficiency histogram
            if (hfprompt) {
                // fill prompt 2D yield: particle level matched
                histStruct.hYieldTruth.first->Fill(jetPt, hfPt);
            } else{
                // fill non-prompt 2D yield: particle level matched
                histStruct.hYieldTruth.second->Fill(jetPt, hfPt);
            }
            
        }
    }
    std::cout << "Non-matched generator level histograms filled.\n";

    //
    // MC detector level tree and histograms
    //
    // Accessing TTree
    tree = (TTree*)fSimulatedMCNonMatched->Get("DF_merged/O2mcdjetdisttable");
    // Check for correct access
    if (!tree) {
        cout << "Error opening mcdjetdisttable tree.\n";
    }
    // defining ML score variables for accessing the TTree
    float hfmlscore0, hfmlscore1, hfmlscore2;
    int jetnconst_int, hfMatchedFrom, hfSelectedAs;
    tree->SetBranchAddress("fJetHfDist",&axisDistance);
    tree->SetBranchAddress("fJetPt",&jetPt);
    tree->SetBranchAddress("fJetEta",&jetEta);
    tree->SetBranchAddress("fJetPhi",&jetPhi);
    tree->SetBranchAddress("fJetNConst",&jetnconst_int);
    tree->SetBranchAddress("fHfPt",&hfPt);
    tree->SetBranchAddress("fHfEta",&hfEta);
    tree->SetBranchAddress("fHfPhi",&hfPhi);
    tree->SetBranchAddress("fHfMass",&hfMass);
    tree->SetBranchAddress("fHfY",&hfY);
    tree->SetBranchAddress("fHfPrompt",&hfprompt);
    tree->SetBranchAddress("fHfMatch",&hfmatch);
    tree->SetBranchAddress("fHfMlScore0",&hfmlscore0);
    tree->SetBranchAddress("fHfMlScore1",&hfmlscore1);
    tree->SetBranchAddress("fHfMlScore2",&hfmlscore2);
    tree->SetBranchAddress("fHfMatchedFrom",&hfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&hfSelectedAs);

    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // only compute matched detector level candidates, but compute all particle level ones
        bool isReflection = (hfMatchedFrom != hfSelectedAs) ? true : false;
        if (!hfmatch || isReflection) {
            continue;
        }
        
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));
        bool recoLevelRange = (abs(jetEta) < MCDetaCut) && (abs(hfY) < MCDyCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= deltaRBinEdges_detector[0]) && (deltaR < MCDDeltaRcut)) && ((hfPt >= MCDHfPtMincut) && (hfPt < MCDHfPtMaxcut));
        // Get the threshold for this pT range
        double maxBkgProb = GetBkgProbabilityCut(hfPt, bdtPtCuts);
        bool passBDTcut = (hfmlscore0 < maxBkgProb) ? true : false;

        // Fill histograms considering jet pT and detector acceptance
        if (recoLevelRange && passBDTcut) {
            
            // fill prompt efficiency histogram
            if (hfprompt) {
                // fill prompt 2D yield: detector level matched
                histStruct.hYieldMeasured.first->Fill(jetPt, hfPt);
            } else{
                // fill non-prompt 2D yield: detector level matched
                histStruct.hYieldMeasured.second->Fill(jetPt, hfPt);
            }
            
        }
        
        
    }
    std::cout << "Non-matched detector level histograms filled.\n";

}
//_____________________________________________________________________________________________________________________________________________________________________________________

// Multiply hTruth and response matrix, The result will the the folded histogram
TH2D* manualFolding(RooUnfoldResponse response, TH2D* hTruth, TH2D* hMeasured) {
    
    // Create empty histogram for the folded data
    TH2D* hFolded = (TH2D*)hMeasured->Clone("hFolded");
    hFolded->Reset();
    hFolded->Sumw2();

    // Get the number of bins in measured and truth histograms
    int nBinsXMeasured = hFolded->GetNbinsX();
    int nBinsYMeasured = hFolded->GetNbinsY();
    int nBinsXTruth = hTruth->GetNbinsX();
    int nBinsYTruth = hTruth->GetNbinsY();
    
    // Debug: print a few values from the response matrix
    bool debug = false;
    if (debug) {
        std::cout << "Measured histogram: " << nBinsXMeasured << " x " << nBinsYMeasured << std::endl;
        std::cout << "Response matrix dimensions: " 
                << response.GetNbinsMeasured() << " x " << response.GetNbinsTruth() << std::endl;
    }
    

    // loop through detector level bins
    for (int iMeasured = 0; iMeasured < nBinsXMeasured; iMeasured++) {
        for (int jMeasured = 0; jMeasured < nBinsYMeasured; jMeasured++) {
            double foldedValue = 0;
            double foldedError2 = 0;

            // obtaining flattened 1D index through row-major ordering
            int index_x_measured = iMeasured + nBinsXMeasured*jMeasured;

            // calculating element iMeasured,jMeasured of folded 2D matrix
            for (int iTruth = 0; iTruth < nBinsXTruth; iTruth++) {
                for (int jTruth = 0; jTruth < nBinsYTruth; jTruth++) {
                    // obtaining flattened 1D index through row-major ordering
                    int index_x_truth = iTruth + nBinsXTruth*jTruth;
                    // calculating matrix element product
                    double truthValue = hTruth->GetBinContent(iTruth + 1,jTruth + 1);
                    double responseValue = response(index_x_measured, index_x_truth);
                    foldedValue += truthValue*responseValue;
                    foldedError2 += std::pow(hTruth->GetBinError(iTruth + 1,jTruth + 1),2) * std::pow(response(index_x_measured, index_x_truth),2);
                }
            }
            hFolded->SetBinContent(iMeasured + 1, jMeasured + 1, foldedValue);
            hFolded->SetBinError(iMeasured + 1, jMeasured + 1, std::sqrt(foldedError2));
        }
    }
    
    std::cout << "Folded manually with bin index flattening." << std::endl;

    return hFolded;
}

void performEfficiencyCorrection(TFile* fBackSub, EfficiencyData& histStruct, double jetptMin, double jetptMax) {
    
    //
    // Calculating pT,D dependent efficiency distribution
    //
    
    // 0.1 - Calculate truth kinematic efficiency
    histStruct.hKEffResponseParticle_Over_TotalParticle.first = (TH2D*)histStruct.hKEffResponseParticle.first->Clone("hKEffTruthPrompt");
    histStruct.hKEffResponseParticle_Over_TotalParticle.second = (TH2D*)histStruct.hKEffResponseParticle.second->Clone("hKEffTruthNonPrompt");
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->Sumw2(); // necessary for correct error propagation
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->Sumw2();
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->Divide(histStruct.hKEffTruthTotalParticle.first);
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->Divide(histStruct.hKEffTruthTotalParticle.second);
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->SetTitle("#varepsilon_{kinematic} of particle level prompt D0s;#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->SetTitle("#varepsilon_{kinematic} of particle level non-prompt D0s;#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");
    
    // 0.2 - Calculate reco kinematic efficiency
    histStruct.hKEffResponseDetector_Over_TotalDetector.first = (TH2D*)histStruct.hKEffResponseDetector.first->Clone("hKEffRecoPrompt");
    histStruct.hKEffResponseDetector_Over_TotalDetector.second = (TH2D*)histStruct.hKEffResponseDetector.second->Clone("hKEffRecoNonPrompt");
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->Sumw2();
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->Sumw2();
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->Divide(histStruct.hKEffRecoTotalDetector.first);
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->Divide(histStruct.hKEffRecoTotalDetector.second);
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->SetTitle("#varepsilon_{kinematic} of detector level prompt D0s;#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->SetTitle("#varepsilon_{kinematic} of detector level non-prompt D0s;#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");

    // 1- Correct the particle level input histogram: remove entries from input distribution whose matched entry is outside of detector level range encoding
    histStruct.hYieldTruthCorrected.first = (TH2D*)histStruct.hYieldTruth.first->Clone("hYieldTruthCorrectedPrompt");
    histStruct.hYieldTruthCorrected.second = (TH2D*)histStruct.hYieldTruth.second->Clone("hYieldTruthCorrectedNonPrompt");
    histStruct.hYieldTruthCorrected.first->Sumw2();
    histStruct.hYieldTruthCorrected.second->Sumw2();
    histStruct.hYieldTruthCorrected.first->Multiply(histStruct.hKEffResponseParticle_Over_TotalParticle.first);
    histStruct.hYieldTruthCorrected.second->Multiply(histStruct.hKEffResponseParticle_Over_TotalParticle.second);

    // 2 - Fold the corrected distribution into the response matrix
    histStruct.hYieldTruthCorrected.first = manualFolding(histStruct.response.first, histStruct.hYieldTruthCorrected.first, histStruct.hYieldMeasured.first);
    histStruct.hYieldTruthCorrected.second = manualFolding(histStruct.response.second, histStruct.hYieldTruthCorrected.second, histStruct.hYieldMeasured.second);

    // 3 - Correct the detector level output histogram: add entries to the output distribution which would be present at outside ranges of the response matrix at particle level
    histStruct.hYieldTruthCorrected.first->Divide(histStruct.hKEffResponseDetector_Over_TotalDetector.first);
    histStruct.hYieldTruthCorrected.second->Divide(histStruct.hKEffResponseDetector_Over_TotalDetector.second);

    // 3.1 - Obtain pT projection histograms
    int minBin, maxBin;
    minBin = histStruct.hYieldTruth.first->GetXaxis()->FindBin(jetptMin);
    maxBin = histStruct.hYieldTruth.first->GetXaxis()->FindBin(jetptMax - 1e-6); // // Tiny epsilon to stay within range: This includes the nearest bin center for jetptMax, but if jetptMax lies between two bins, it might give slightly unintuitive results
    histStruct.hHfPtYieldTruth.first = histStruct.hYieldTruth.first->ProjectionY("hHfPtYieldTruthPrompt", minBin, maxBin);
    histStruct.hHfPtYieldTruth.first->SetTitle("Denominator prompt D^{0} p_{T} distribution before correction; #it{p}_{T, D^{0}}^{gen}; Counts");
    minBin = histStruct.hYieldTruth.second->GetXaxis()->FindBin(jetptMin);
    maxBin = histStruct.hYieldTruth.second->GetXaxis()->FindBin(jetptMax - 1e-6);
    histStruct.hHfPtYieldTruth.second = histStruct.hYieldTruth.second->ProjectionY("hHfPtYieldTruthNonPrompt", minBin, maxBin);
    histStruct.hHfPtYieldTruth.second->SetTitle("Denominator non-prompt D^{0} p_{T} distribution before correction; #it{p}_{T, D^{0}}^{gen}; Counts");
    minBin = histStruct.hYieldTruthCorrected.first->GetXaxis()->FindBin(jetptMin);
    maxBin = histStruct.hYieldTruthCorrected.first->GetXaxis()->FindBin(jetptMax - 1e-6);
    histStruct.hHfPtYieldTruthCorrected.first = histStruct.hYieldTruthCorrected.first->ProjectionY("hHfPtYieldTruthCorrectedPrompt", minBin, maxBin);
    histStruct.hHfPtYieldTruthCorrected.first->SetTitle("Denominator prompt D^{0} p_{T} distribution after correction; #it{p}_{T, D^{0}}^{gen}; Counts");
    minBin = histStruct.hYieldTruthCorrected.second->GetXaxis()->FindBin(jetptMin);
    maxBin = histStruct.hYieldTruthCorrected.second->GetXaxis()->FindBin(jetptMax - 1e-6);
    histStruct.hHfPtYieldTruthCorrected.second = histStruct.hYieldTruthCorrected.second->ProjectionY("hHfPtYieldTruthCorrectedNonPrompt", minBin, maxBin);
    histStruct.hHfPtYieldTruthCorrected.second->SetTitle("Denominator non-prompt D^{0} p_{T} distribution after correction; #it{p}_{T, D^{0}}^{gen}; Counts");

    // 4 - Calculate the efficiency distributions with the pT projections
    histStruct.hSelectionEfficiency.first = histStruct.hYieldMeasured.first->ProjectionY("hSelectionEfficiencyPrompt", minBin, maxBin);
    histStruct.hSelectionEfficiency.second = histStruct.hYieldMeasured.second->ProjectionY("hSelectionEfficiencyNonPrompt", minBin, maxBin);
    histStruct.hSelectionEfficiency.first->SetTitle("Efficiency prompt D^{0} p_{T} distribution; #it{p}_{T, D^{0}}^{reco}; Efficiency#times Acceptance");
    histStruct.hSelectionEfficiency.second->SetTitle("Efficiency non-prompt D^{0} p_{T} distribution; #it{p}_{T, D^{0}}^{reco}; Efficiency#times Acceptance");
    histStruct.hSelectionEfficiency.first->Sumw2();
    histStruct.hSelectionEfficiency.second->Sumw2();
    histStruct.hSelectionEfficiency.first->Divide(histStruct.hHfPtYieldTruthCorrected.first);
    histStruct.hSelectionEfficiency.second->Divide(histStruct.hHfPtYieldTruthCorrected.second);
    
    //
    // Correcting 3D background subtracted distribution with efficiency for each pT,D bin
    //
    // 1 - Clone background subtracted distribution
    histStruct.hEfficiencyCorrected.first = (TH3D*)fBackSub->Get("h3DBackgroundSubtracted")->Clone("h3DEfficiencyCorrected");
    histStruct.hEfficiencyCorrected.first->SetTitle("Background subtracted, efficiency corrected");
    int xBins = histStruct.hEfficiencyCorrected.first->GetXaxis()->GetNbins();
    int yBins = histStruct.hEfficiencyCorrected.first->GetYaxis()->GetNbins();
    int zBins = histStruct.hEfficiencyCorrected.first->GetZaxis()->GetNbins();
    for (int xBin = 1; xBin <= xBins; xBin++) {
        for (int yBin = 1; yBin <= yBins; yBin++) {
            for (int zBin = 1; zBin <= zBins; zBin++) { // pT,D bins will be corrected by 1/prompt_efficiency
                
                // Get current content and error
                double content = histStruct.hEfficiencyCorrected.first->GetBinContent(xBin, yBin, zBin);
                double error = histStruct.hEfficiencyCorrected.first->GetBinError(xBin, yBin, zBin);

                // Get the pT,D bin center from the Z axis of the 3D histogram
                double ptDcenter = histStruct.hEfficiencyCorrected.first->GetZaxis()->GetBinCenter(zBin);
                
                // Find corresponding bin in efficiency histogram
                int effBin = histStruct.hSelectionEfficiency.first->GetXaxis()->FindBin(ptDcenter);
                double eff = histStruct.hSelectionEfficiency.first->GetBinContent(effBin);

                // Avoid divide by zero or nonsense values
                if (eff > 0) {
                    double correction = 1. / eff;
                    histStruct.hEfficiencyCorrected.first->SetBinContent(xBin, yBin, zBin, content * correction);
                    histStruct.hEfficiencyCorrected.first->SetBinError(xBin, yBin, zBin, error * correction);
                } else {
                    // Optionally warn about zero efficiency
                    std::cerr << "Zero or invalid efficiency at ptD = " << ptDcenter
                            << " (eff bin = " << effBin << ")" << std::endl;
                }
            }
        }
    }
    
    // 2 - Project the z axis (pT,D) in order to obtain a 2D distribution of pT,jet (x axis) vs DeltaR (y axis)
    TH2D* h2D = (TH2D*)histStruct.hEfficiencyCorrected.first->Project3D("yx");
    histStruct.hEfficiencyCorrected.second = (TH2D*)h2D->Clone("h2DEfficiencyCorrected");
    delete h2D;
    histStruct.hEfficiencyCorrected.second->SetTitle("Background subtracted, efficiency corrected");

    std::cout << "Efficiency correction performed.\n";
}

TCanvas* plotMCnet_template_response(TH2D* originalHistogram, double jetptMin, double jetptMax) {
    
    TCanvas* cMCnetFormat = new TCanvas("cMCnetFormat", "MCnet_format_response", 800, 600);
    cMCnetFormat->cd();
    //cMCnetFormat->SetRightMargin(0.2); // default ~0.1; increase to give space to color bar
    //cMCnetFormat->SetLeftMargin(0.13);  // decrease if there's too much space on the left
    double statBoxPos = gPad->GetUxmax(); // Height of the stat box
    gStyle->SetOptStat(0); // Turn off the default stats box
    gStyle->SetTextFont(42); // Helvetica-like, cleaner than default
    gPad->SetLogz();

    TH2D* histogram = (TH2D*)originalHistogram->Clone("hMCnetFormat");
    histogram->SetMarkerStyle(kFullCircle); //kFullDotMedium
    histogram->SetMarkerColor(kBlack);
    histogram->SetLineColor(kBlack);
    TString title = histogram->GetTitle();
    histogram->SetTitle("");
    //histogram->GetXaxis()->SetTitle("p_{T, ch jet}^{det} (GeV/#it{c})");
    //histogram->GetYaxis()->SetTitle("p_{T, ch jet}^{gen} (GeV/#it{c})");
    histogram->GetXaxis()->SetTitleSize(0.04);
    histogram->GetYaxis()->SetTitleSize(0.04);
    //histogram->GetYaxis()->SetTitleOffset(1.1);
    histogram->GetZaxis()->SetTitle("Counts");
    histogram->GetZaxis()->SetTitleOffset(1.01);
    histogram->SetMinimum(0);
    histogram->Draw("colz");
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    float topY = 0.83; // Top Y position for the text
    float bottomY = 0.15; // Bottom Y position for the text
    float deltaY = 0.06; // Vertical spacing between lines
    float middleX = 0.52; // Middle X position for the text
    float leftX = 0.13; // Left X position for the text
    float rightX = 0.87; // or 0.90, depending on your layout
    float textSize = 0.04; // Text size for the labels

    // Allign text to the left
    TLatex* latexLeft = new TLatex();
    latexLeft->SetTextColor(kWhite); // ROOT's predefined white color
    latexLeft->SetTextSizePixels(24);
    latexLeft->SetNDC(); // Enables normalized coordinates (0 to 1)
    latexLeft->SetTextSize(textSize); // default = 0.05
    latexLeft->DrawLatex(leftX, topY, "ALICE Performance");
    latexLeft->DrawLatex(leftX, topY-deltaY, "pp, #sqrt{#it{s}} = 13.6 TeV");
    latexLeft->DrawLatex(leftX, topY-2*deltaY, "Prompt D^{0}");

    // Allign text to the right
    TLatex* latexRight = new TLatex();
    latexRight->SetTextColor(kWhite); // ROOT's predefined white color
    latexRight->SetTextSizePixels(24);
    latexRight->SetNDC(); // Enables normalized coordinates (0 to 1)
    latexRight->SetTextSize(textSize); // default = 0.05
    latexRight->SetTextAlign(31); // Right-align text to that X position
    latexRight->DrawLatex(rightX, topY, "D^{0}#rightarrow K^{#minus}+#pi^{+} and ch. conj.");
    latexRight->DrawLatex(rightX, topY - deltaY, "in ch-particle jets, anti-#it{k}_{T}, #it{R} = 0.4");
    //latexRight->DrawLatex(rightX, topY - 2*deltaY, "1 < #it{p}_{T, D^{0}} (GeV/#it{c}) < 30, |#it{y}_{D^{0}}| #leq 0.8"); // middle GeV/c - Vit
    latexRight->DrawLatex(rightX, topY - 2*deltaY, "1 < #it{p}_{T, D^{0}} < 30 GeV/#it{c}, |#it{y}_{D^{0}}| #leq 0.8"); // right GeV/c - Raymond
    latexRight->DrawLatex(rightX, topY - 3*deltaY, "|#it{#eta}_{jet}| #leq 0.5");
    //latexRight->DrawLatex(rightX, topY - 3*deltaY, Form("%.0f < #it{p}_{T, ch. jet} < %.0f GeV/#it{c}, |#it{#eta}_{jet}| #leq 0.5", jetptMin, jetptMax));
    
    cMCnetFormat->Update();

    // Now get the palette and change its width
    TPaletteAxis* palette = (TPaletteAxis*)histogram->GetListOfFunctions()->FindObject("palette");
    if (palette) {
        //palette->SetX1NDC(0.88); // Left edge of color bar
        palette->SetX2NDC(0.93); // Right edge of color bar
        // Example: default might be 0.85–0.90, so this makes it a bit narrower
    }

    cMCnetFormat->Modified();
    cMCnetFormat->Update();

    return cMCnetFormat;
}

void PlotHistograms(const EfficiencyData& histStruct, double jetptMin, double jetptMax) {

    //
    // 2D yield histograms
    //
    TCanvas* cYield = new TCanvas("cYield","2D yield histograms");
    cYield->Divide(2,2);
    cYield->SetCanvasSize(1800,1000);
    cYield->cd(1);
    //histStruct.hYieldTruth.first->SetStats(0);
    histStruct.hYieldTruth.first->Draw("COLZ");
    cYield->cd(2);
    //histStruct.hYieldTruth.second->SetStats(0);
    histStruct.hYieldTruth.second->Draw("COLZ");
    cYield->cd(3);
    //histStruct.hYieldMeasured.first->SetStats(0);
    histStruct.hYieldMeasured.first->Draw("COLZ");
    cYield->cd(4);
    //histStruct.hYieldMeasured.second->SetStats(0);
    histStruct.hYieldMeasured.second->Draw("COLZ");
    cYield->Update();

    //
    // Response matrix representation in 2D histogram
    //
    TCanvas* cResponse = new TCanvas("cResponse","Response matrices for all pT,D bins");
    cResponse->SetCanvasSize(1800,1000);
    cResponse->Divide(2,2);
    std::pair<const TH2*, const TH2*> hresponse2D;
    hresponse2D.first = histStruct.response.first.Hresponse();
    hresponse2D.second = histStruct.response.second.Hresponse();
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
    histStruct.responseProjections.first[0]->Draw("colz");
    cResponseProjections->cd(2);
    histStruct.responseProjections.first[1]->Draw("colz");
    cResponseProjections->cd(3);
    histStruct.responseProjections.second[0]->Draw("colz");
    cResponseProjections->cd(4);
    histStruct.responseProjections.second[1]->Draw("colz");

    TCanvas* cMCnet1 = new TCanvas("cMCnet1","MCnet1");
    cMCnet1->SetCanvasSize(1800,1000);
    cMCnet1->cd();
    histStruct.responseProjections.first[0]->Draw("colz");
    TCanvas* cMCnet2 = new TCanvas("cMCnet2","MCnet2");
    cMCnet2->SetCanvasSize(1800,1000);
    cMCnet2->cd();
    histStruct.responseProjections.first[1]->Draw("colz");

    TCanvas* cMCnet_template_response = plotMCnet_template_response(histStruct.responseProjections.first[0], jetptMin, jetptMax);
    cMCnet_template_response->Draw();

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    
    // Kinematic efficiency histograms: truth entries
    TCanvas* cKEffTruthPrompt = new TCanvas("cKEffTruthPrompt","Kinematic efficiency histograms for truth prompt entries");
    cKEffTruthPrompt->SetCanvasSize(1800,1000);
    cKEffTruthPrompt->Divide(2,2);
    cKEffTruthPrompt->cd(1);
    histStruct.hKEffTruthTotalParticle.first->Draw("colz");
    cKEffTruthPrompt->cd(2);
    histStruct.hKEffResponseParticle.first->Draw("colz");
    cKEffTruthPrompt->cd(3);
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->Draw("text");
    TCanvas* cKEffTruthNonPrompt = new TCanvas("cKEffTruthNonPrompt","Kinematic efficiency histograms for truth non-prompt entries");
    cKEffTruthNonPrompt->SetCanvasSize(1800,1000);
    cKEffTruthNonPrompt->Divide(2,2);
    cKEffTruthNonPrompt->cd(1);
    histStruct.hKEffTruthTotalParticle.second->Draw("colz");
    cKEffTruthNonPrompt->cd(2);
    histStruct.hKEffResponseParticle.second->Draw("colz");
    cKEffTruthNonPrompt->cd(3);
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->Draw("text");

    // Kinematic efficiency histograms: reco entries
    TCanvas* cKEffRecoPrompt = new TCanvas("cKEffRecoPrompt","Kinematic efficiency histograms for reco prompt entries");
    cKEffRecoPrompt->SetCanvasSize(1800,1000);
    cKEffRecoPrompt->Divide(2,2);
    cKEffRecoPrompt->cd(1);
    histStruct.hKEffRecoTotalDetector.first->Draw("colz");
    cKEffRecoPrompt->cd(2);
    histStruct.hKEffResponseDetector.first->Draw("colz");
    cKEffRecoPrompt->cd(3);
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->Draw("text");
    TCanvas* cKEffRecoNonPrompt = new TCanvas("cKEffRecoNonPrompt","Kinematic efficiency histograms for reco non-prompt entries");
    cKEffRecoNonPrompt->SetCanvasSize(1800,1000);
    cKEffRecoNonPrompt->Divide(2,2);
    cKEffRecoNonPrompt->cd(1);
    histStruct.hKEffRecoTotalDetector.second->Draw("colz");
    cKEffRecoNonPrompt->cd(2);
    histStruct.hKEffResponseDetector.second->Draw("colz");
    cKEffRecoNonPrompt->cd(3);
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->Draw("text");

    // All kinematic efficiency histograms
    TCanvas* cKEffAll = new TCanvas("cKEffAll","All kinematic efficiency histograms for all entries");
    cKEffAll->SetCanvasSize(1800,1000);
    cKEffAll->Divide(2,2);
    cKEffAll->cd(1);
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->SetMarkerSize(1.); // default = 1.5
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->Draw("text");
    cKEffAll->cd(2);
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->SetMarkerSize(1.);
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->Draw("text");
    cKEffAll->cd(3);
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->SetMarkerSize(1.);
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->Draw("text");
    cKEffAll->cd(4);
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->SetMarkerSize(1.);
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->Draw("text");
    cKEffAll->Update();

    TCanvas* cYieldTruthCorrected = new TCanvas("cYieldTruthCorrected","Corrected denominator histograms and original ones");
    cYieldTruthCorrected->SetCanvasSize(1800,1000);
    cYieldTruthCorrected->Divide(2,2);
    cYieldTruthCorrected->cd(1);
    histStruct.hYieldTruthCorrected.first->SetTitle("Corrected particle level data prompt yield distribution");
    histStruct.hYieldTruthCorrected.first->Draw("colz");
    cYieldTruthCorrected->cd(2);
    histStruct.hYieldTruthCorrected.second->SetTitle("Corrected particle level data non-prompt yield distribution");
    histStruct.hYieldTruthCorrected.second->Draw("colz");
    cYieldTruthCorrected->cd(3);
    histStruct.hYieldTruth.first->Draw("colz");
    cYieldTruthCorrected->cd(4);
    histStruct.hYieldTruth.second->Draw("colz");

    TCanvas* cPtCorrectionComparison = new TCanvas("cPtCorrectionComparison","pT correction comparison");
    cPtCorrectionComparison->SetCanvasSize(1800,1000);
    cPtCorrectionComparison->Divide(2,2);
    cPtCorrectionComparison->cd(1);
    histStruct.hHfPtYieldTruth.first->Draw("colz");
    cPtCorrectionComparison->cd(2);
    histStruct.hHfPtYieldTruth.second->Draw("colz");
    cPtCorrectionComparison->cd(3);
    histStruct.hHfPtYieldTruthCorrected.first->Draw("colz");
    cPtCorrectionComparison->cd(4);
    histStruct.hHfPtYieldTruthCorrected.second->Draw("colz");

    TCanvas* cSelectionEfficiency = new TCanvas("cSelectionEfficiency","Efficiency histograms");
    cSelectionEfficiency->SetCanvasSize(1800,1000);
    cSelectionEfficiency->Divide(2,2);
    cSelectionEfficiency->cd(1);
    histStruct.hSelectionEfficiency.first->Draw();
    cSelectionEfficiency->cd(2);
    histStruct.hSelectionEfficiency.second->Draw();
    // Separating in two canvases
    TCanvas* cSelectionEfficiencyPrompt = new TCanvas("cSelectionEfficiencyPrompt","Efficiency histograms for prompt D^{0}'s");
    cSelectionEfficiencyPrompt->SetCanvasSize(1800,1000);
    cSelectionEfficiencyPrompt->cd();
    histStruct.hSelectionEfficiency.first->Draw();
    TCanvas* cSelectionEfficiencyNonPrompt = new TCanvas("cSelectionEfficiencyNonPrompt","Efficiency histograms for non-prompt D^{0}'s");
    cSelectionEfficiencyNonPrompt->SetCanvasSize(1800,1000);
    cSelectionEfficiencyNonPrompt->cd();
    histStruct.hSelectionEfficiency.second->Draw();

    TCanvas* cCorrectedData = new TCanvas("cCorrectedData","Background subtracted, efficiency corrected data");
    cCorrectedData->Divide(2,2);
    cCorrectedData->cd(1);
    histStruct.hEfficiencyCorrected.first->Draw("colz");
    cCorrectedData->cd(2);
    histStruct.hEfficiencyCorrected.second->Draw("colz");
    cCorrectedData->cd(3);
    histStruct.hEfficiencyCorrected.second->ProjectionY("hEfficiencyCorrectedProjectionDeltaRPrompt")->Draw();

    TCanvas* cCorrectedData2DOnly = new TCanvas("cCorrectedData2DOnly","Background subtracted, efficiency corrected 2D data");
    cCorrectedData2DOnly->cd();
    histStruct.hEfficiencyCorrected.second->Draw("colz");
    TCanvas* cCorrectedData1DOnly = new TCanvas("cCorrectedData1DOnly","Background subtracted, efficiency corrected 1D data");
    cCorrectedData1DOnly->cd();
    //histStruct.hEfficiencyCorrected.second->SetMarkerStyle(kCircle);
    histStruct.hEfficiencyCorrected.second->ProjectionY("hEfficiencyCorrectedProjectionDeltaRPrompt")->Draw();


    //
    // Storing images
    //
    TString imagePath = "../Images/2-Efficiency/Run3Style/";
    cYield->Update();
    cYield->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cYield->SaveAs(imagePath + "Efficiency_run3style_yields.png");
    cResponse->Update();
    cResponse->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cResponse->SaveAs(imagePath + "Efficiency_run3style_response.png");
    cResponsePrompt->Update();
    cResponsePrompt->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cResponsePrompt->SaveAs(imagePath + "Efficiency_run3style_response_prompt.png");
    cResponseNonPrompt->Update();
    cResponseNonPrompt->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cResponseNonPrompt->SaveAs(imagePath + "Efficiency_run3style_response_nonprompt.png");
    cResponseProjections->Update();
    cResponseProjections->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cResponseProjections->SaveAs(imagePath + "Efficiency_run3style_response_projections.png");
    cMCnet1->Update();
    cMCnet1->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cMCnet1->SaveAs(imagePath + "Efficiency_run3style_response_projections_MCnet1.png");
    cMCnet2->Update();
    cMCnet2->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cMCnet2->SaveAs(imagePath + "Efficiency_run3style_response_projections_MCnet2.png");
    cMCnet_template_response->Update();
    cMCnet_template_response->SaveAs(imagePath + "Efficiency_run3style_response_projections_MCnet_template_response.eps");
    cKEffTruthPrompt->Update();
    cKEffAll->Update();
    cKEffAll->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cKEffAll->SaveAs(imagePath + "Efficiency_run3style_kinematic_efficiency.png");
    cYieldTruthCorrected->Update();
    cYieldTruthCorrected->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cYieldTruthCorrected->SaveAs(imagePath + "Efficiency_run3style_yields_corrected.png");
    cPtCorrectionComparison->Update();
    cPtCorrectionComparison->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cPtCorrectionComparison->SaveAs(imagePath + "Efficiency_run3style_yields_pt_correction_comparison.png");
    cSelectionEfficiency->Update();
    cSelectionEfficiency->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cSelectionEfficiency->SaveAs(imagePath + "Efficiency_run3style_selection_efficiency.png");
    cSelectionEfficiencyPrompt->Update();
    cSelectionEfficiencyPrompt->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cSelectionEfficiencyPrompt->SaveAs(imagePath + "Efficiency_run3style_selection_efficiency_prompt.png");
    cSelectionEfficiencyNonPrompt->Update();
    cSelectionEfficiencyNonPrompt->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cSelectionEfficiencyNonPrompt->SaveAs(imagePath + "Efficiency_run3style_selection_efficiency_nonprompt.png");
    cCorrectedData->Update();
    cCorrectedData->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cCorrectedData->SaveAs(imagePath + "Efficiency_run3style_corrected_data.png");
    cCorrectedData2DOnly->Update();
    cCorrectedData2DOnly->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cCorrectedData2DOnly->SaveAs(imagePath + "Efficiency_run3style_corrected_data_2d_only.png");
    cCorrectedData1DOnly->Update();
    cCorrectedData1DOnly->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cCorrectedData1DOnly->SaveAs(imagePath + "Efficiency_run3style_corrected_data_1d_only.png");

    //
    // Storing in a single pdf file
    //
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

void SaveData(const EfficiencyData& histStruct, double jetptMin, double jetptMax, 
              const std::vector<double>& ptjetBinEdges_detector, const std::vector<double>& ptjetBinEdges_particle,
              const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& deltaRBinEdges_particle,
              const std::vector<double>& ptDBinEdges_detector, const std::vector<double>& ptDBinEdges_particle) {
    // Open output file
    TFile* outFile = new TFile(Form("selection_efficiency_run3style_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    
    // Response matrices
    histStruct.response.first.Write("hResponsePrompt");
    histStruct.response.second.Write("hResponseNonPrompt");

    // Truth kinematic efficiency histograms
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->Write();
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->Write();

    // Reco kinematic efficiency histograms
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->Write();
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->Write();

    // Selection efficiency histograms
    histStruct.hSelectionEfficiency.first->Write();
    histStruct.hSelectionEfficiency.second->Write();

    // Efficiency corrected data
    histStruct.hEfficiencyCorrected.second->Write();
    
    // Also store the axes used for the histograms
    // Create a directory for axes
    outFile->mkdir("axes");
    outFile->cd("axes");
    // Create TVectorD with same content
    TVectorD vecPtJet_detector(ptjetBinEdges_detector.size());
    for (size_t i = 0; i < ptjetBinEdges_detector.size(); ++i) {
        vecPtJet_detector[i] = ptjetBinEdges_detector[i];
    }
    vecPtJet_detector.Write("ptjetBinEdges_detector");
    TVectorD vecDeltaR_detector(deltaRBinEdges_detector.size());
    for (size_t i = 0; i < deltaRBinEdges_detector.size(); ++i) {
        vecDeltaR_detector[i] = deltaRBinEdges_detector[i];
    }
    vecDeltaR_detector.Write("deltaRBinEdges_detector");
    TVectorD vecPtD_detector(ptDBinEdges_detector.size());
    for (size_t i = 0; i < ptDBinEdges_detector.size(); ++i) {
        vecPtD_detector[i] = ptDBinEdges_detector[i];
    }
    vecPtD_detector.Write("ptDBinEdges_detector");
    TVectorD vecPtJet_particle(ptjetBinEdges_particle.size());
    for (size_t i = 0; i < ptjetBinEdges_particle.size(); ++i) {
        vecPtJet_particle[i] = ptjetBinEdges_particle[i];
    }
    vecPtJet_particle.Write("ptjetBinEdges_particle");
    TVectorD vecDeltaR_particle(deltaRBinEdges_particle.size());
    for (size_t i = 0; i < deltaRBinEdges_particle.size(); ++i) {
        vecDeltaR_particle[i] = deltaRBinEdges_particle[i];
    }
    vecDeltaR_particle.Write("deltaRBinEdges_particle");
    TVectorD vecPtD_particle(ptDBinEdges_particle.size());
    for (size_t i = 0; i < ptDBinEdges_particle.size(); ++i) {
        vecPtD_particle[i] = ptDBinEdges_particle[i];
    }
    vecPtD_particle.Write("ptDBinEdges_particle");
    // Return to root directory (optional)
    outFile->cd();

    outFile->Close();
    delete outFile;
    
    std::cout << "Data stored.\n";
}

std::vector<double> LoadBinning(TFile* fInput, const char* pathInFile) {
    auto* vec = (TVectorD*)fInput->Get(pathInFile);
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", pathInFile));
    }
    return std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
}

void Efficiency_run3_style_detector_level(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Luminosity (for now arbitrary)
    double luminosity = 0;
    double luminosity_powheg = 0;
    double BR = 0.0393; // D0 -> KPi decay channel branching ratio = (3.93 +- 0.04) %

    

    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.03 (from first of pair)
    //std::vector<std::pair<double, double>> bdtPtCuts = {
    //    {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}
    //};
    // Dataset: JE_HF_LHC24g5_All_D0
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.4}, {8, 0.6}, {10, 0.8}, {12, 0.8}, {16, 1.0}, {50, 1.0}
    };
    
    // opening files
    //TFile* fSimulatedMCNonMatched = new TFile("../SimulatedData/Hyperloop_output/McEfficiency/New_with_reflections/Merged_AO2D_HF_LHC24d3a_All.root","read"); // previous dataset used
    //TFile* fSimulatedMCMatched = new TFile("../SimulatedData/Hyperloop_output/McChargedMatched/mergedMatched.root","read"); // previous dataset used
    TFile* fSimulatedMCNonMatched = new TFile("../SimulatedData/Hyperloop_output/Train_runs/410602_Eff/AO2D_mergedDFs.root","read");
    TFile* fSimulatedMCMatched = new TFile("../SimulatedData/Hyperloop_output/Train_runs/410603_Match/AO2D_mergedDFs.root","read");
    TFile* fBackSub = new TFile(Form("../1-SignalTreatment/SideBand/full_merged_ranges_back_sub.root"),"read");
    if (!fSimulatedMCNonMatched || fSimulatedMCNonMatched->IsZombie()) {
        std::cerr << "Error: Unable to open full (not matched) simulated data ROOT file." << std::endl;
    }
    if (!fSimulatedMCMatched || fSimulatedMCMatched->IsZombie()) {
        std::cerr << "Error: Unable to open matched simulated data ROOT file." << std::endl;
    }

    // Load pTjet bin edges
    std::vector<double> ptjetBinEdges_detector = LoadBinning(fBackSub, "axes/ptjetBinEdges_detector");
    std::vector<double> ptjetBinEdges_particle = ptjetBinEdges_detector;
    double jetptMin = ptjetBinEdges_detector[0]; // GeV
    double jetptMax = ptjetBinEdges_detector[ptjetBinEdges_detector.size() - 1]; // GeV
    // Load ΔR bin edges
    std::vector<double> deltaRBinEdges_detector = LoadBinning(fBackSub, "axes/deltaRBinEdges_detector");
    std::vector<double> deltaRBinEdges_particle = deltaRBinEdges_detector;
    double minDeltaR = deltaRBinEdges_detector[0];
    double maxDeltaR = deltaRBinEdges_detector[deltaRBinEdges_detector.size() - 1];
    // Load pTD bin edges
    std::vector<double> ptDBinEdges_detector = LoadBinning(fBackSub, "axes/ptDBinEdges_detector");
    std::vector<double> ptDBinEdges_particle = ptDBinEdges_detector;
    double hfptMin = ptDBinEdges_detector[0]; //ptDBinEdges[0] - should start from 0 or from the lowest pT,D value?
    double hfptMax = ptDBinEdges_detector[ptDBinEdges_detector.size() - 1];

    TFile* fEffRun2Style = new TFile(Form("run2_style_efficiency_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fEffRun2Style || fEffRun2Style->IsZombie()) {
        std::cerr << "Error: Unable to open run 2 style efficiency data ROOT file." << std::endl;
    }

    EfficiencyData histStruct = createHistograms(ptjetBinEdges_particle, ptDBinEdges_particle, ptjetBinEdges_detector, ptDBinEdges_detector); // pT histograms
    
    // Fill matched histograms for corrections
    fillMatchedHistograms(fSimulatedMCMatched, fEffRun2Style, histStruct, ptjetBinEdges_detector, ptDBinEdges_detector, deltaRBinEdges_detector, ptjetBinEdges_particle, ptDBinEdges_particle, deltaRBinEdges_particle, bdtPtCuts);
    // Fill non-matched histograms for final numerator and denominator division
    fillNonMatchedHistograms(fSimulatedMCNonMatched, histStruct, jetptMin, jetptMax, hfptMin, hfptMax, deltaRBinEdges_particle, deltaRBinEdges_detector, bdtPtCuts);

    // Calculate efficiency distribution
    performEfficiencyCorrection(fBackSub, histStruct, jetptMin, jetptMax);

    // Plot the efficiency histogram and further corrected histograms
    PlotHistograms(histStruct, jetptMin, jetptMax);

    // Save corrected distributions to file
    SaveData(histStruct, jetptMin, jetptMax, 
             ptjetBinEdges_detector, ptjetBinEdges_particle,
             deltaRBinEdges_detector, deltaRBinEdges_particle,
             ptDBinEdges_detector, ptDBinEdges_particle);

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
