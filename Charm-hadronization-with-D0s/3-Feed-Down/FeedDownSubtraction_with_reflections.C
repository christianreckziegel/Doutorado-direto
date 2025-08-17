/*
 *
 *
 * Macro for performing B feed-down subtraction 
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

    return numHistograms;
}
// calculate number of background subtracted histograms inside file
int HistogramCounter(TFile* file, bool byName = true) {
    TList* keys = file->GetListOfKeys();
    std::set<TString> uniqueHistograms;
    int numHistograms = 0;

    for (int i = 0; i < keys->GetSize(); ++i) {
        TKey* key = (TKey*)keys->At(i);
        TString name = key->GetName();

        if (name.BeginsWith("h_back_subtracted_")) {
            // Only consider the latest cycle (ignore duplicates)
            if (uniqueHistograms.find(name) == uniqueHistograms.end()) {
                TObject* obj = file->Get(name);
                if (obj && obj->IsA()->InheritsFrom(TH1::Class())) {
                    uniqueHistograms.insert(name);
                    numHistograms++;
                }
            }
        }
    }

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

struct FeedDownData {
    // Step 1: POWHEG + PYTHIA particle level data of non-prompt D0s: deltaR vs pT,jet vs pT,D
    TH3D* hPowheg = nullptr;
    
    // Step 2: apply efficiency (pT,D dependent) correction
    std::pair<TH1D*, TH1D*> hSelectionEfficiency; // first = prompt D0s, second = non-prompt D0s
    TH3D* hPowhegEffCorrected = nullptr; // pT,D-efficiency-corrected 3D histogram

    // Step 3: projection to 2D (jet pT vs DeltaR)
    TH2D* hProjected2D = nullptr;

    // Step 4: folding
    std::vector<TH2D*> hKineEffParticle = {nullptr, nullptr, nullptr}; // particle level (addition): [0] = numerator, [1] = denominator, [2] = efficiency
    std::vector<TH2D*> hKineEffDetector = {nullptr, nullptr, nullptr}; // detector level (removal): [0] = numerator, [1] = denominator, [2] = efficiency
    RooUnfoldResponse response; // response matrix for folding
    TH2D* hProjected2DKinCorrected = nullptr; // 3D histogram projection with particle level kinematic efficiency correction
    TH2D* hMeasuredTemplate = nullptr; // template for measured data (jet pT vs DeltaR), used for hFolded binning
    TH2D* hFolded2D = nullptr; // folded 2D histogram (jet pT vs DeltaR)
    TH2D* hFolded2DKinCorrected = nullptr; // folded 2D histogram with detector level kinematic correction

    // Step 5: final scaled histogram (jet pT vs DeltaR)
    TH2D* hFinalScaled = nullptr;

    // Metadata
    double lumiData = 1; // luminosity of the data in 1/pb
    double lumiMC = 1; // luminosity of the MC in 1/pb
    double branchingRatio = 0.0393; // D0 -> KPi decay channel branching ratio = (3.93 +- 0.04) %

    // Step 6: B subtracted data distribution
    TH2D* hBFedDownData;
};

// Obtain the bin edges of a histogram (useful for asymmetrical bin sizes)
std::vector<double> getBinEdges(const TAxis* axis) {
    
    if (!axis) {
        std::cerr << "Error: null axis pointer.\n";
    }
    
    int nBins = axis->GetNbins();

    // vector for storing bin edges
    std::vector<double> binEdges(nBins + 1);

    for (int iBin = 0; iBin <= nBins; iBin++) {
        binEdges[iBin] = axis->GetBinLowEdge(iBin + 1);
    }

    return binEdges;
}

// Module to create TH2D histograms including interest variable: ptjetBinEdges_particle, deltaRBinEdges_particle, ptDBinEdges_particle
FeedDownData createHistograms(const std::vector<double>& ptjetBinEdges_particle, const std::vector<double>& deltaRBinEdges_particle, const std::vector<double>& ptDBinEdges_particle,
                              const std::vector<double>& ptjetBinEdges_detector, const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& ptDBinEdges_detector) {
                              //const double& jetptMin, const double& jetptMax) {
    // Create struct to store data
    FeedDownData dataContainer;

    // Create 3D histogram
    dataContainer.hPowheg = new TH3D("hPowheg", "POWHEG + PYTHIA;p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR;p_{T,D}^{gen}", ptjetBinEdges_particle.size()-1, ptjetBinEdges_particle.data(), 
                                                                                                                          deltaRBinEdges_particle.size()-1, deltaRBinEdges_particle.data(), 
                                                                                                                          ptDBinEdges_particle.size()-1, ptDBinEdges_particle.data());
    dataContainer.hPowheg->SetMarkerColor(30);
    dataContainer.hPowheg->SetLineColor(30); // 30 = pastel green
    dataContainer.hPowheg->SetMarkerStyle(kCircle);
    dataContainer.hPowheg->Sumw2();
    dataContainer.hPowheg->SetStats(0);

    
    //
    // Matching histograms for folding process
    //
    dataContainer.hKineEffParticle[0] = new TH2D("hKineEffParticleNumerator", "Particle level kinematic efficiency numerator (non-prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                                                  deltaRBinEdges_particle.size() - 1, deltaRBinEdges_particle.data());
    dataContainer.hKineEffParticle[1] = new TH2D("hKineEffParticleDenominator", "Particle level kinematic efficiency denominator (non-prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                                                  deltaRBinEdges_particle.size() - 1, deltaRBinEdges_particle.data());
    dataContainer.hKineEffDetector[0] = new TH2D("hKineEffDetectorNumerator", "Detector level kinematic efficiency numerator (non-prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                                                  deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data());
    dataContainer.hKineEffDetector[1] = new TH2D("hKineEffDetectorDenominator", "Detector level kinematic efficiency denominator (non-prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                                                  deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data());
    
    dataContainer.hMeasuredTemplate = new TH2D("hMeasuredTemplate", "Folded 2D histogram;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}", 
                                                                                  ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                                                  deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data());
    std::cout << "Histograms created.\n";

    return dataContainer;
}

// Module to fill POWHEG+PYTHIA 3D data histogram (all data in file is non-prompt and on particle level)
void fillNonMatchedHistograms(TFile* fPowheg, TFile* fEfficiency, FeedDownData& dataContainer, double& jetptMin, double& jetptMax, double& hfptMin, double& hfptMax,
                              const std::vector<double>& deltaRBinEdges_particle, const std::vector<double>& deltaRBinEdges_detector) {
    
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
    TTree* tree = (TTree*)fPowheg->Get("tree_D0");

    // Check for correct access
    if (!tree) {
        cout << "Error opening POWHEG tree.\n";
    }
    
    // defining variables for accessing data on TTree
    double hfPt, hfEta, hfPhi, hfY;
    double jetPt, jetEta, jetPhi, axisDistance;

    tree->SetBranchAddress("pt_cand",&hfPt);
    tree->SetBranchAddress("eta_cand",&hfEta);
    tree->SetBranchAddress("phi_cand",&hfPhi);
    tree->SetBranchAddress("y_cand",&hfY);
    tree->SetBranchAddress("pt_jet",&jetPt);
    tree->SetBranchAddress("eta_jet",&jetEta);
    tree->SetBranchAddress("phi_jet",&jetPhi);
    tree->SetBranchAddress("delta_r_jet",&axisDistance);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));
        bool genLevelRange = (abs(jetEta) < MCPetaCut) && (abs(hfY) < MCPyCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= deltaRBinEdges_particle[0]) && (deltaR < MCPDeltaRcut)) && ((hfPt >= MCPHfPtMincut) && (hfPt < MCPHfPtMaxcut));
        
        // Fill 2D histogram considering jet pT and detector acceptance
        if (genLevelRange) {
            
            dataContainer.hPowheg->Fill(jetPt, axisDistance, hfPt);
            
        }
        
    }
    std::cout << "Generator level (POWHEG+PYTHIA) histograms filled.\n";

    // Access prompt efficiency
    dataContainer.hSelectionEfficiency.first = (TH1D*)fEfficiency->Get("hSelectionEfficiencyPrompt");
    dataContainer.hSelectionEfficiency.second = (TH1D*)fEfficiency->Get("hSelectionEfficiencyNonPrompt");

}

// Get the optimal BDT score cut for the corresponding pT,D of the entry
double GetBkgProbabilityCut(double pT, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    for (size_t i = 0; i < bdtPtCuts.size() - 1; ++i) {
        if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) {
            return bdtPtCuts[i].second;
        }
    }
    return 1.0; // Default: accept all if out of range
}
// Module to fill 2D response matrix and kinematic efficiency histograms
void fillMatchedHistograms(TFile* fSimulatedMCMatched, FeedDownData& dataContainer, double& jetptMin, double& jetptMax, double& hfptMin, double& hfptMax,
                           const std::vector<double>& deltaRBinEdges_particle, const std::vector<double>& deltaRBinEdges_detector, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    
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
    
    // Create response matrix for folding
    std::cout << "Creating response matrix" << std::endl;
    if (dataContainer.hKineEffParticle[0] == nullptr || dataContainer.hKineEffDetector[0] == nullptr) {
        std::cerr << "Error: Kinematic efficiency histograms are not initialized.\n";
        return;
    }
    dataContainer.response = RooUnfoldResponse(dataContainer.hKineEffDetector[0], dataContainer.hKineEffParticle[0]);
    std::cout << "Response matrix created.\n";

    //___________________________________________
    //
    // Fill response matrix and kinematic efficiency histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)fSimulatedMCMatched->Get("DF_merged/O2matchtable");
    // Check for correct access
    if (!tree) {
        cout << "Error opening O2 matching tree.\n";
    }

    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
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
    // detector level branches
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);
    tree->SetBranchAddress("fHfMlScore0",&MCDhfMlScore0);
    tree->SetBranchAddress("fHfMlScore1",&MCDhfMlScore1);
    tree->SetBranchAddress("fHfMlScore2",&MCDhfMlScore2);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // Apply non-prompt selection (i.e., only B → D0)
        if (MCDhfprompt) {
            continue;
        }

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        // Get the threshold for this pT range
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;

        // Fill response matrix and kinematic efficiency histograms
        if (genLevelRange && recoLevelRange && passBDTcut) {
            
            // Find the bin corresponding to the given pT,D value
            int bin = dataContainer.hSelectionEfficiency.first->FindBin(MCDhfPt);
            // Get the efficiency value from the bin content
            double efficiency_prompt = dataContainer.hSelectionEfficiency.first->GetBinContent(bin);
            // Fill 4D RooUnfoldResponse object (jet pT shape is influenced by D0 pT efficiency)
            dataContainer.response.Fill(MCDjetPt, MCDDeltaR, MCPjetPt, MCPDeltaR, 1./efficiency_prompt);
            
            // Fill kinematic efficiency numerator histograms
            dataContainer.hKineEffParticle[0]->Fill(MCPjetPt, MCPDeltaR);
            dataContainer.hKineEffDetector[0]->Fill(MCDjetPt, MCDDeltaR);
        }
        if (genLevelRange) {
            // Fill kinematic efficiency denominator histogram for full particle level range
            dataContainer.hKineEffParticle[1]->Fill(MCPjetPt, MCPDeltaR);
        }
        if (recoLevelRange && passBDTcut) {
            // Fill kinematic efficiency denominator histogram for full detector level range
            dataContainer.hKineEffDetector[1]->Fill(MCDjetPt, MCDDeltaR);
        }
        
        
    }

    // Correct underflow and overflow bins of kinematic efficiency numerator and denominator histograms (wrong concept, no need)
    
    
    std::cout << "Response matrix filled.\n";
}

// Manual folding of TH2D data using a RooUnfoldResponse object.
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

/**
 * @brief Remove entries in region not treated by response matrix from the folding input data.
 *
 * This helper function
 *
 * @param dataContainer Container with total and intersection range data and the POWHEG data to be corrected.
 *
 * @return Corrected POWHEG truth input 2D distribution
 *
 * @note 
 *
 * @see smearGeneratorData() [Correct truth POWHEG data before folding.]
 */
void removeOutsideData(FeedDownData& dataContainer) {

    // Calculate particle level kinematic efficiency histograms
    dataContainer.hKineEffParticle[2] = (TH2D*)dataContainer.hKineEffParticle[0]->Clone("hKineEffParticleEfficiency");
    dataContainer.hKineEffParticle[2]->Divide(dataContainer.hKineEffParticle[1]); // A = A / B:  A = A->Divide(B)

    // Copy POWHEG efficiency corrected 2D histogram hProjected2D
    dataContainer.hProjected2DKinCorrected = (TH2D*)dataContainer.hProjected2D->Clone("hProjected2DKinCorrected");
    dataContainer.hProjected2DKinCorrected->SetTitle("Folded 2D histogram with #varepsilon_{kin}^{particle} correction;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR^{gen}");

    // Apply kinematic efficiency correction
    dataContainer.hProjected2DKinCorrected->Multiply(dataContainer.hKineEffParticle[2]);

    //return dataContainer.hProjected2DKinCorrected;
}

/**
 * @brief Add entries of the folding output data from region not treated by response matrix .
 *
 * This helper function
 *
 * @param dataContainer Container with total and intersection range data and the folded data to be corrected.
 *
 * @return Corrected Folded detector level output 2D distribution
 *
 * @note 
 *
 * @see smearGeneratorData() [Correct detector level folded data.]
 */
void addOutsideData(FeedDownData& dataContainer) {

    // Calculate detector level kinematic efficiency histograms
    dataContainer.hKineEffDetector[2] = (TH2D*)dataContainer.hKineEffDetector[0]->Clone("hKineEffDetectorEfficiency");
    dataContainer.hKineEffDetector[2]->Divide(dataContainer.hKineEffDetector[1]); // A = A / B:  A = A->Divide(B)

    // Copy POWHEG folded 2D histogram hProjected2D
    dataContainer.hFolded2DKinCorrected = (TH2D*)dataContainer.hFolded2D->Clone("hFolded2DKinCorrected");
    dataContainer.hFolded2DKinCorrected->SetTitle("Folded 2D histogram with #varepsilon_{kin}^{detector} correction;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}");

    // Apply kinematic efficiency correction
    dataContainer.hFolded2DKinCorrected->Divide(dataContainer.hKineEffDetector[2]);

    //return dataContainer.hFolded2DKinCorrected;

}

// Get POWHEG and data luminosities
void getLuminosities(TFile* fPowheg, TFile* fData, FeedDownData& dataContainer, double luminosity_powheg, double luminosity) {
    // Accessing total cross section value stored in first bin (in mb)
    TH1D* xSection_powheg = dynamic_cast<TH1D*>(fPowheg->Get("fHistXsection"));
    double crossSecPowheg = xSection_powheg->GetBinContent(1);

    // Accessing number of events
    TTree* tree_D0 = dynamic_cast<TTree*>(fPowheg->Get("tree_D0"));
    double numOfEventsPowheg = tree_D0->GetEntries();

    // integrated POWHEG luminosity 
    luminosity_powheg = numOfEventsPowheg/crossSecPowheg;
    dataContainer.lumiMC = luminosity_powheg; // Store in dataContainer for later use
    std::cout << "POWHEG luminosity = " << luminosity_powheg << " mb^-1" << std::endl;

    // Calculating measured luminosity
    // folder: jet-luminosity-calculator
    // Histogram name: counter
    // bin labels: "BC+TVX", "Coll+TVX", "Coll+TVX+VtxZ+Sel8"
    TH1D* hDataLumi = (TH1D*)fData->Get("jet-luminosity-calculator/counter");
    int bin = hDataLumi->GetXaxis()->FindBin("BC+TVX");
    double bcTVX = hDataLumi->GetBinContent(bin); // BC+TVX
    std::cout << "On bin number " << bin <<", BC+TVX = " << bcTVX << std::endl;
    bin = hDataLumi->GetXaxis()->FindBin("Coll+TVX+VtxZ+Sel8");
    double selection = hDataLumi->GetBinContent(bin); // Coll+TVX+VtxZ+Sel8
    std::cout << "On bin number " << bin <<", Coll+TVX+VtxZ+Sel8 = " << selection << std::endl;
    bin = hDataLumi->GetXaxis()->FindBin("Coll+TVX");
    double collTVX = hDataLumi->GetBinContent(bin); // Coll+TVX
    std::cout << "On bin number " << bin <<", Coll+TVX = " << collTVX << std::endl;
    double triggered = bcTVX * selection / collTVX; // number of TVX triggered BC that correspond to your selections and your train efficiencies
    double runLuminosity = 1.0/0.0594e6; // luminosity value for the runs (// in mb⁻¹?)
    double dataLuminosity = triggered * runLuminosity;
    std::cout << "triggered = " << triggered << "; runLuminosity = " << runLuminosity << std::endl;
    std::cout << "BC+TVX = " << bcTVX << ", Coll+TVX = " << collTVX << ", Selection = " << selection << ", runLuminosity = " << runLuminosity << std::endl;
    dataContainer.lumiData = dataLuminosity; // Store in dataContainer for later use
    std::cout << "Measured luminosity = " << dataLuminosity << " mb^-1" << std::endl;
}

// Module for folding particle level data from POWHEG simulation
void smearGeneratorData(FeedDownData& dataContainer, double& luminosity_powheg, TFile* fEfficiency, double& luminosity, double& BR) {
    //
    // 0th step: clone 3D histogram in order to save the smeared at the end (while keeping the original)
    //
    dataContainer.hPowhegEffCorrected = (TH3D*)dataContainer.hPowheg->Clone("hPowhegEffCorrected");
    dataContainer.hPowhegEffCorrected->SetTitle("1st step: POWHEG distribution efficiency corrected;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR;p_{T,D}^{gen}");

    //
    // 1st step: apply efficiency (pT,D dependent) correction -> scale by efficiency ratio
    //
    for (int xBin = 1; xBin <= dataContainer.hPowhegEffCorrected->GetNbinsX(); xBin++) {
        for (int yBin = 1; yBin <= dataContainer.hPowhegEffCorrected->GetNbinsY(); yBin++) {
            for (int zBin = 1; zBin <= dataContainer.hPowhegEffCorrected->GetNbinsZ(); zBin++) { //Inner loop over z-axis bins: this is scaled axis (pT,D bins)

                // Bin center of the Z axis (i.e., pT,D)
                double ptDcenter = dataContainer.hPowhegEffCorrected->GetZaxis()->GetBinCenter(zBin);

                // Find corresponding bin in the efficiency histograms
                int effBinPrompt = dataContainer.hSelectionEfficiency.first->FindBin(ptDcenter);
                int effBinNonPrompt = dataContainer.hSelectionEfficiency.second->FindBin(ptDcenter);

                // Get efficiency values
                double effPrompt = dataContainer.hSelectionEfficiency.first->GetBinContent(effBinPrompt);
                double effNonPrompt = dataContainer.hSelectionEfficiency.second->GetBinContent(effBinNonPrompt);

                // Get 3D histogram content
                double binContent = dataContainer.hPowhegEffCorrected->GetBinContent(xBin, yBin, zBin);
                double binError = dataContainer.hPowhegEffCorrected->GetBinError(xBin, yBin, zBin);

                // Rescale the bin content by the ratio of efficiencies
                if (effPrompt > 0) {
                    double efficiencyRatio = effNonPrompt / effPrompt;
                    dataContainer.hPowhegEffCorrected->SetBinContent(xBin, yBin, zBin, binContent * efficiencyRatio);
                    // Rescale the bin error accordingly
                    dataContainer.hPowhegEffCorrected->SetBinError(xBin, yBin, zBin, binError * efficiencyRatio);
                } else {
                    double ptDlow = dataContainer.hPowhegEffCorrected->GetZaxis()->GetBinLowEdge(zBin);
                    double ptDhigh = dataContainer.hPowhegEffCorrected->GetZaxis()->GetBinUpEdge(zBin);

                    std::cout << "Warning: Efficiency for pT,D bin " << zBin 
                              << " (range: " << ptDlow << " - " << ptDhigh << " GeV/c, center: " << ptDcenter 
                              << ") is zero or negative(effPrompt = " << effPrompt << "). Skipping scaling." << std::endl;

                }

            }
            
        }
        
    }
    
    //
    // 2nd step: obtain the 2D Delta R vs pT,jet projection for all pT,D bins before folding
    // 
    dataContainer.hProjected2D = (TH2D*)dataContainer.hPowhegEffCorrected->Project3D("yxy"); // xy = pT,jet vs Delta R
    dataContainer.hProjected2D->SetName("hProjected2D");
    dataContainer.hProjected2D->SetTitle("2nd step: projected POWHEG efficiency corrected histogram;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR");

    //
    // 3.1th step: remove outside of response range data in POWHEG (apply particle levelkinematic efficiency)
    //
    removeOutsideData(dataContainer);

    //
    // 3.2th step: fold Delta R vs pT,jet distribution using detector response matrix of non-prompt D0 jets
    //
    dataContainer.hFolded2D = manualFolding(dataContainer.response, dataContainer.hProjected2DKinCorrected, dataContainer.hMeasuredTemplate);
    delete dataContainer.hMeasuredTemplate;

    //
    // 3.3th step: add outside of response range data in folded data (apply 1 / detector level kinematic efficiency)
    //
    addOutsideData(dataContainer);

    //
    // 4th step: scale by 1 / POWHEG integrated luminosity, the measured integrated luminosity and BR of D0 decay channel
    //
    dataContainer.hFinalScaled = (TH2D*)dataContainer.hFolded2DKinCorrected->Clone("hFinalScaled");
    dataContainer.hFinalScaled->Scale(dataContainer.branchingRatio * dataContainer.lumiData / dataContainer.lumiMC);

    std::cout << "Generator data smeared.\n";
}

// Module to subtract non-prompt D0 jets from prompt efficiency corrected distribution
void feedDown(TFile* fEfficiency, FeedDownData& dataContainer, const double jetptMin, const double jetptMax) {
    
    dataContainer.hBFedDownData = (TH2D*)fEfficiency->Get("h2DEfficiencyCorrected")->Clone("hBFedDownData");
    dataContainer.hBFedDownData->SetTitle("Background subtracted, efficiency corrected, fed down data distribution");

    dataContainer.hBFedDownData->Add(dataContainer.hFinalScaled,-1);
    
    std::cout << "Feed-down subtracted from experimental data.\n";
}

void plotHistograms(const FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax) {
    cout << "Plotting histograms...\n";

    // Changing color palette to a bigger one
    gStyle->SetPalette(kRainBow);

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    //
    // 3D histograms
    //
    TCanvas* cPowheg = new TCanvas("cPowheg","POWHEG data");
    cPowheg->SetCanvasSize(1800,1000);
    cPowheg->cd();
    dataContainer.hPowheg->Draw("colz");
    
    //
    // Response matrix
    //
    TCanvas* cResponse = new TCanvas("cResponse","Response matrix");
    cResponse->SetCanvasSize(1800,1000);
    cResponse->cd();
    const TH2* hresponse2D = dataContainer.response.Hresponse();
    TH2D* hresponse2DClone = static_cast<TH2D*>(hresponse2D->Clone("hResponse2D"));
    hresponse2DClone->SetTitle("2D response matrix from 4D RooUnfoldResponse - non-prompt D^{0}'s;2D Reconstructed;2D Truth");
    hresponse2DClone->Draw("colz");
    //
    // Kinematic efficiency histograms
    //
    TCanvas* cKinEfficiency = new TCanvas("cKineEfficiency","Kinematic efficiency");
    cKinEfficiency->SetCanvasSize(1800,1000);
    cKinEfficiency->Divide(3,2);
    cKinEfficiency->cd(1);
    dataContainer.hKineEffParticle[0]->Draw("colz"); // numerator
    cKinEfficiency->cd(2);
    dataContainer.hKineEffParticle[1]->Draw("colz"); // denominator
    cKinEfficiency->cd(3);
    dataContainer.hKineEffParticle[2]->Draw("text"); // efficiency
    cKinEfficiency->cd(4);
    dataContainer.hKineEffDetector[0]->Draw("colz"); // numerator
    cKinEfficiency->cd(5);
    dataContainer.hKineEffDetector[1]->Draw("colz"); // denominator
    cKinEfficiency->cd(6);
    dataContainer.hKineEffDetector[2]->Draw("text"); // efficiency
    
    //
    // Folded data
    //
    TCanvas* cFolded = new TCanvas("cFolded","Folded data");
    cFolded->SetCanvasSize(1800,1000);
    cFolded->Divide(2,2);
    cFolded->cd(1);
    dataContainer.hProjected2DKinCorrected->Draw("colz");
    cFolded->cd(2);
    dataContainer.hFolded2D->Draw("colz");
    cFolded->cd(3);
    dataContainer.hFolded2DKinCorrected->Draw("colz");

    TCanvas* cFedDownData = new TCanvas("cFedDownData","Background subtracted, efficiency corrected, fed down data distribution");
    cFedDownData->Divide(2,2);
    cFedDownData->cd(1);
    dataContainer.hBFedDownData->Draw("colz");
    cFedDownData->cd(2);
    dataContainer.hBFedDownData->ProjectionY("hBFedDownDataDeltaR")->Draw();

    //
    // Storing images
    //
    TString imagePath = "../Images/3-Feed-Down/";
    cPowheg->Update();
    cPowheg->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cPowheg->SaveAs(imagePath + "powheg_3d.png");
    cResponse->Update();
    cResponse->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cResponse->SaveAs(imagePath + "response_matrix.png");
    cKinEfficiency->Update();
    cKinEfficiency->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cKinEfficiency->SaveAs(imagePath + "kinematic_efficiencies.png");
    cFolded->Update();
    cFolded->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cFolded->SaveAs(imagePath + "folded_stages.png");
    cFedDownData->Update();
    cFedDownData->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cFedDownData->SaveAs(imagePath + "fed_down_data.png");

    //
    // Storing in a single pdf file
    //
    cPowheg->Print(imagePath + Form("feeddown_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cResponse->Print(imagePath + Form("feeddown_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cKinEfficiency->Print(imagePath + Form("feeddown_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFolded->Print(imagePath + Form("feeddown_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFedDownData->Print(imagePath + Form("feeddown_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

}

void saveData(const FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax, 
              const std::vector<double>& ptjetBinEdges_detector, const std::vector<double>& ptjetBinEdges_particle,
              const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& deltaRBinEdges_particle,
              const std::vector<double>& ptDBinEdges_detector, const std::vector<double>& ptDBinEdges_particle){
    // Open output file
    TFile* outFile = new TFile(Form("outputFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");

    // Store B -> D0 contribution histogram
    dataContainer.hFinalScaled->Write();

    // Store fed down data histogram
    dataContainer.hBFedDownData->Write();
    
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
    
    cout << "Data stored in file" << Form("outputFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)) << endl;
}

std::vector<double> LoadBinning(TFile* fInput, const char* pathInFile) {
    auto* vec = (TVectorD*)fInput->Get(pathInFile);
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", pathInFile));
    }
    return std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
}

void FeedDownSubtraction_with_reflections(){

    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Luminosity (for now arbitrary)
    double luminosity = 0;
    double luminosity_powheg = 0;
    double BR = 0.0393; // D0 -> KPi decay channel branching ratio = (3.93 +- 0.04) %

    TFile* fAxes = new TFile(Form("../1-SignalTreatment/SideBand/full_merged_ranges_back_sub.root"),"read");
    if (!fAxes || fAxes->IsZombie()) {
        std::cerr << "Error: Unable to open simulated data ROOT file." << std::endl;
    }
    // Load pTjet bin edges
    std::vector<double> ptjetBinEdges_detector = LoadBinning(fAxes, "axes/ptjetBinEdges_detector");
    double jetptMin = ptjetBinEdges_detector[0]; // GeV
    double jetptMax = ptjetBinEdges_detector[ptjetBinEdges_detector.size() - 1]; // GeV
    // Load ΔR bin edges
    std::vector<double> deltaRBinEdges_detector = LoadBinning(fAxes, "axes/deltaRBinEdges_detector");
    double minDeltaR = deltaRBinEdges_detector[0];
    double maxDeltaR = deltaRBinEdges_detector[deltaRBinEdges_detector.size() - 1];
    // Load pTD bin edges
    std::vector<double> ptDBinEdges_detector = LoadBinning(fAxes, "axes/ptDBinEdges_detector");
    double hfptMin = ptDBinEdges_detector[0]; //ptDBinEdges[0] - should start from 0 or from the lowest pT,D value?
    double hfptMax = ptDBinEdges_detector[ptDBinEdges_detector.size() - 1];
    fAxes->Close();

    TFile* fEfficiency = new TFile(Form("../2-Efficiency/selection_efficiency_run3style_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fEfficiency || fEfficiency->IsZombie()) {
        std::cerr << "Error: Unable to open estimated selection efficiency ROOT file." << std::endl;
    }
    std::vector<double> ptjetBinEdges_particle = LoadBinning(fEfficiency, "axes/ptjetBinEdges_particle");
    std::vector<double> deltaRBinEdges_particle = LoadBinning(fEfficiency, "axes/deltaRBinEdges_particle");
    std::vector<double> ptDBinEdges_particle = LoadBinning(fEfficiency, "axes/ptDBinEdges_particle");

    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.03 (from first of pair)
    //std::vector<std::pair<double, double>> bdtPtCuts = {
    //    {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}
    //};
    // Dataset: JE_HF_LHC24g5_All_D0
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.4}, {8, 0.6}, {10, 0.8}, {12, 0.8}, {16, 1.0}, {50, 1.0}
    };

    // Opening files
    TFile* fPowheg = new TFile("../SimulatedData/POWHEG/trees_powheg_fd_central.root","read");
    TFile* fSimulatedMCMatched = new TFile("../SimulatedData/Hyperloop_output/Train_runs/410603_Match/AO2D_mergedDFs.root","read");
    TFile* fData = new TFile("../ExperimentalData/Hyperloop_output/HF_LHC23_pass4_Thin_small_2P3PDstar_DATA_newMLModel/AnalysisResults.root","read");
    if (!fPowheg || fPowheg->IsZombie()) {
        std::cerr << "Error: Unable to open POWHEG simulation data ROOT file." << std::endl;
    }
    if (!fSimulatedMCMatched || fSimulatedMCMatched->IsZombie()) {
        std::cerr << "Error: Unable to open O2 MC matched ROOT file." << std::endl;
    }
    if (!fData || fData->IsZombie()) {
        std::cerr << "Error: Unable to open AnalysisResults.root data ROOT file." << std::endl;
    }
    FeedDownData dataContainer = createHistograms(ptjetBinEdges_particle, deltaRBinEdges_particle, ptDBinEdges_particle,
                                                  ptjetBinEdges_detector, deltaRBinEdges_detector, ptDBinEdges_detector);

    // Fill histograms with POWHEG simulation data and estimated efficiency
    fillNonMatchedHistograms(fPowheg, fEfficiency, dataContainer, jetptMin, jetptMax, hfptMin, hfptMax, deltaRBinEdges_particle, deltaRBinEdges_detector);
    // Fill histograms with matched MC simulation data and response matrix
    fillMatchedHistograms(fSimulatedMCMatched, dataContainer, jetptMin, jetptMax, hfptMin, hfptMax, deltaRBinEdges_particle, deltaRBinEdges_detector, bdtPtCuts);

    // Get POWHEG and data luminosities
    getLuminosities(fPowheg, fData, dataContainer, luminosity_powheg, luminosity);

    // Fold data using two methods
    smearGeneratorData(dataContainer, luminosity_powheg, fEfficiency, luminosity, BR);

    // Subtract non-prompt distribution from prompt efficiency corrected ones
    feedDown(fEfficiency, dataContainer, jetptMin, jetptMax);

    // Plot the efficiency histogram and further corrected histograms
    plotHistograms(dataContainer, jetptMin, jetptMax);

    // Save corrected distributions to file
    saveData(dataContainer, jetptMin, jetptMax, 
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
    FeedDownSubtraction_with_reflections();
    return 0;
}
