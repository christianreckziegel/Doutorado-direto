/*
 *
 *
 * Macro for performing unfolding on prompt D0s measured 
 * data and applying correction to FeedDownSubtraction.C
 * resulting distributions.
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

struct UnfoldData {
    // Step 2: apply efficiency (pT,D dependent) correction
    std::pair<TH1D*, TH1D*> hSelectionEfficiency;                           // first = prompt D0s, second = non-prompt D0s

    // Step: background subtracted, efficiency corrected, B fed-down 2D distribution (pT,jet vs DeltaR)
    TH2D* hBFedDownData;

    // Step 4: bayesian unfolding
    std::vector<TH2D*> hKineEffParticle = {nullptr, nullptr, nullptr};      // particle level (addition): [0] = numerator, [1] = denominator, [2] = efficiency
    std::vector<TH2D*> hKineEffDetector = {nullptr, nullptr, nullptr};      // detector level (removal): [0] = numerator, [1] = denominator, [2] = efficiency
    TH2D* hBFedDownDataKinCorrected = nullptr;                            // 2D histogram with detector level kinematic efficiency correction
    RooUnfoldResponse* response;                                             // response matrix for folding
    std::vector<TH2D*> hResponse2D = {nullptr, nullptr};                     // response projections matrix: first = DeltaR, second = pT,jet
    std::vector<RooUnfoldBayes*> unfold;                                     // unfolding objects, there are iterationNumber unfolding objects
    TH2D* hMeasuredTemplate = nullptr;                                      // template for measured data (jet pT vs DeltaR), used for hUnfolded binning
    TH2D* hTruthTemplate = nullptr;
    std::vector<TH2D*> hUnfolded;                                         // unfolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects
    std::vector<TH2D*> hUnfoldedKinCorrected;                             // unfolded 2D histogram with detector level kinematic correction, there are iterationNumber unfolding objects

    // Refolding test
    std::vector<TH2D*> hRefolded;                                           // refolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects
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

// Module to create TH2D histograms including interest variable
UnfoldData createHistograms(const std::vector<double>& ptjetBinEdges_particle, const std::vector<double>& deltaRBinEdges_particle, const std::vector<double>& ptDBinEdges_particle,
                              const std::vector<double>& ptjetBinEdges_detector, const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& ptDBinEdges_detector,
                            int& iterationNumber) {
                              //const double& jetptMin, const double& jetptMax) {
    // Create struct to store data
    UnfoldData dataContainer;
    
    // Create kinematic efficiency histograms
    dataContainer.hKineEffParticle[0] = new TH2D("hKineEffParticleNumerator", "Particle level kinematic efficiency numerator (prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                                                  deltaRBinEdges_particle.size() - 1, deltaRBinEdges_particle.data());
    dataContainer.hKineEffParticle[1] = new TH2D("hKineEffParticleDenominator", "Particle level kinematic efficiency denominator (prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                                                  deltaRBinEdges_particle.size() - 1, deltaRBinEdges_particle.data());
    dataContainer.hKineEffDetector[0] = new TH2D("hKineEffDetectorNumerator", "Detector level kinematic efficiency numerator (prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                                                  deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data());
    dataContainer.hKineEffDetector[1] = new TH2D("hKineEffDetectorDenominator", "Detector level kinematic efficiency denominator (prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                                                  deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data());
    
    // Create template histograms used for response matrix creation
    dataContainer.hMeasuredTemplate = new TH2D("hMeasuredTemplate", "Measured template;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}", 
                                                                                  ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                                                  deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data());
    dataContainer.hTruthTemplate = new TH2D("hTruthTemplate", "Truth template;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR^{gen}", 
                                                                                  ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                                                  deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data());

    // Create response matrices
    dataContainer.response = new RooUnfoldResponse(dataContainer.hKineEffDetector[0], dataContainer.hKineEffParticle[0]);

    dataContainer.hResponse2D[0] = new TH2D("hResponseDeltaR", "Response matrix projection on #DeltaR;#DeltaR^{reco};#DeltaR^{gen}", 
                                                                                  deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data(), 
                                                                                  deltaRBinEdges_particle.size() - 1, deltaRBinEdges_particle.data());
    dataContainer.hResponse2D[1] = new TH2D("hResponsePtJet", "Response matrix projection on p_{T,jet};p_{T,jet}^{reco} (GeV/#it{c});p_{T,jet}^{gen} (GeV/#it{c})", 
                                                                                  ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                                                  ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data());

    // Reserve space for the unfolding objects with the number of iterations
    dataContainer.unfold.resize(iterationNumber);

    // Reserve space for the unfolded TH2D* histograms with the number of iterations (and point to nullptr)
    dataContainer.hUnfolded.resize(iterationNumber, nullptr);

    // And also for the kinematic efficiency corrected version of the previous mentioned (and point to nullptr)
    dataContainer.hUnfoldedKinCorrected.resize(iterationNumber, nullptr);

    std::cout << "Template histograms and response object created." << std::endl;
    return dataContainer;
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
/**
 * @brief Module to fill histograms from O2 matching task simulation.
 * 
 * This helper function builds two kinds of objects:
 * - the response matrix content
 * - obtain the selection efficiency histograms
 *
 * @param fFeedDown The ROOT file with output from Feed-down subtracted step.
 * @param fEfficiency The ROOT file with run 3 style detector level selection efficiency.
 * @param dataContainer Container with already instanciated histograms.
 * @param jetptMin Min jet pT cut.
 * @param jetptMax Max jet pT cut.
 *
 *
 * @see createHistograms() [Instanciate histograms.]
 */
void fillHistograms(TFile* fFeedDown, TFile* fEfficiency, TFile* fSimulatedMCMatched, UnfoldData& dataContainer, double& jetptMin, double& jetptMax, double& hfptMin, double& hfptMax, 
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

    // Access the background subtracted, efficiency corrected, fed-down distribution
    dataContainer.hBFedDownData = (TH2D*)fFeedDown->Get("hBFedDownData");
    if (dataContainer.hBFedDownData) {
        dataContainer.hBFedDownData = (TH2D*)dataContainer.hBFedDownData->Clone("hBFedDownData");
        std::cout << "Background subtracted, efficiency corrected, fed-down distribution acquired." << std::endl;
    } else {
        std::cout << "Error: failed to acquire background subtracted, efficiency corrected, fed-down distribution acquired." << std::endl;
    }

    // Access prompt efficiency
    dataContainer.hSelectionEfficiency.first = (TH1D*)fEfficiency->Get("hSelectionEfficiencyPrompt");
    dataContainer.hSelectionEfficiency.second = (TH1D*)fEfficiency->Get("hSelectionEfficiencyNonPrompt");
    if (dataContainer.hSelectionEfficiency.first || dataContainer.hSelectionEfficiency.second) {
        std::cout << "Efficiency histograms (run 3 style detector level) acquired." << std::endl;
    } else {
        std::cout << "Error: failed to acquire efficiency histograms (run 3 style detector level)." << std::endl;
    }

    //
    // Response matrix and kinematic efficiency histograms
    //
    // Accessing TTree
    TTree* tree;
    int nEntries;

    //
    // O2 matching task measured and truth histograms
    //
    // Accessing TTree
    tree = (TTree*)fSimulatedMCMatched->Get("DF_merged/O2matchtable");
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

    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // Apply prompt selection (i.e., only c → D0, not B → D0)
        if (!MCDhfprompt) {
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
            dataContainer.response->Fill(MCDjetPt, MCDDeltaR, MCPjetPt, MCPDeltaR, 1./efficiency_prompt);
            // Fill response matrix projections
            dataContainer.hResponse2D[0]->Fill(MCDDeltaR, MCPDeltaR, 1./efficiency_prompt);
            dataContainer.hResponse2D[1]->Fill(MCDjetPt, MCPjetPt, 1./efficiency_prompt);
            
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
    std::cout << "Response object and kinematic efficiency histograms filled." << std::endl;
}

// Detector level kinematic efficiency: removal of data
void removeOutsideData(UnfoldData& dataContainer) {
    // Calculate detector level kinematic efficiency histograms
    dataContainer.hKineEffDetector[2] = (TH2D*)dataContainer.hKineEffDetector[0]->Clone("hKineEffDetectorEfficiency");
    dataContainer.hKineEffDetector[2]->Sumw2();
    dataContainer.hKineEffDetector[2]->Divide(dataContainer.hKineEffDetector[1]); // A = A / B:  A = A->Divide(B)

    // Copy background subtracted, efficiency corrected, fed-down subtracted 2D histogram
    dataContainer.hBFedDownDataKinCorrected = (TH2D*)dataContainer.hBFedDownData->Clone("hBFedDownDataKinCorrected");
    dataContainer.hBFedDownDataKinCorrected->SetTitle("Fed-down with #varepsilon_{kin}^{detector} correction;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}");

    // Apply kinematic efficiency correction
    dataContainer.hBFedDownDataKinCorrected->Sumw2();
    dataContainer.hBFedDownDataKinCorrected->Multiply(dataContainer.hKineEffDetector[2]);

    //return dataContainer.hBFedDownDataKinCorrected;
    std::cout << "Detector level kinematic efficiency applied (removal)." << std::endl;
}

// Particle level kinematic efficiency: addition of data
TH2D* addOutsideData(TH2D* hKineEffParticle, TH2D* hUnfolded, int& iterationNumber) {

    // Copy unfolded 2D histogram
    TH2D* hUnfoldedKinCorrected = (TH2D*)hUnfolded->Clone(Form("hUnfoldedKinCor_iter%d", iterationNumber));
    hUnfoldedKinCorrected->SetTitle(Form("Unfolded with #varepsilon_{kin}^{particle} correction with %d iterations", iterationNumber));
    // Apply kinematic efficiency correction
    hUnfoldedKinCorrected->Sumw2();
    hKineEffParticle->Sumw2();
    hUnfoldedKinCorrected->Divide(hKineEffParticle); // A = A / B:  A = A->Divide(B)

    std::cout << "Particle level kinematic efficiency applied (addition)." << std::endl;
    return hUnfoldedKinCorrected;

}

void performUnfolding(UnfoldData& dataContainer, int& iterationNumber) {

    // Correct distribution with detector level kinematic efficiency (remove entries)
    removeOutsideData(dataContainer);

    // Calculate particle level kinematic efficiency histograms
    dataContainer.hKineEffParticle[2] = (TH2D*)dataContainer.hKineEffParticle[0]->Clone("hKineEffParticleEfficiency");
    dataContainer.hKineEffParticle[2]->Sumw2();
    dataContainer.hKineEffParticle[2]->Divide(dataContainer.hKineEffParticle[1]); // A = A / B:  A = A->Divide(B)

    // Unfold in multiple iterations
    for (int iIter = 1; iIter <= iterationNumber; iIter++) {
        dataContainer.unfold[iIter - 1] = new RooUnfoldBayes(dataContainer.response, dataContainer.hBFedDownDataKinCorrected, iIter);
        dataContainer.hUnfolded[iIter - 1] = (TH2D*) dataContainer.unfold[iIter - 1]->Hreco(RooUnfold::kCovariance); // kCovariance = 2 = Errors from the square root of of the covariance matrix given by the unfolding
        dataContainer.hUnfolded[iIter - 1]->SetTitle(Form("Immediately unfolded 2D histogram with %d iterations;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR^{gen}", iIter));
        dataContainer.hUnfolded[iIter - 1]->SetName(Form("hUnfolded_iter%d", iIter));
        dataContainer.hUnfolded[iIter - 1]->Sumw2();

        // Correct distribution with particle level kinematic efficiency (add entries)
        dataContainer.hUnfoldedKinCorrected[iIter - 1] = addOutsideData(dataContainer.hKineEffParticle[2],dataContainer.hUnfolded[iIter - 1], iIter);
        
    }
    std::cout << "Unfolding procedure performed." << std::endl;
}

void convergenceTest(UnfoldData& dataContainer) {
    //
    std::cout << "Convergence test performed." << std::endl;
}

// Manual folding of TH2D data using a RooUnfoldResponse object.
TH2D* manualFolding(RooUnfoldResponse* response, TH2D* hTruth, TH2D* hMeasured) {
    
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
                << response->GetNbinsMeasured() << " x " << response->GetNbinsTruth() << std::endl;
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
                    double responseValue = (*response)(index_x_measured, index_x_truth);
                    foldedValue += truthValue * responseValue;
                    foldedError2 += std::pow(hTruth->GetBinError(iTruth + 1,jTruth + 1),2) * std::pow((*response)(index_x_measured, index_x_truth),2);
                }
            }
            hFolded->SetBinContent(iMeasured + 1, jMeasured + 1, foldedValue);
            hFolded->SetBinError(iMeasured + 1, jMeasured + 1, std::sqrt(foldedError2));
        }
    }
    
    std::cout << "Folded manually with bin index flattening." << std::endl;

    return hFolded;
}

void refoldingTest(UnfoldData& dataContainer) {
    
    // Ensure the vector has the right size and initialized with nullptrs
    size_t nIter = dataContainer.hUnfoldedKinCorrected.size();
    dataContainer.hRefolded.resize(nIter, nullptr);

    // Loop through unfolded histograms of different iterations number
    for (size_t iHisto = 0; iHisto < nIter; iHisto++) {
        // Copy unfolded histogram (with detector level kinematic efficiency correction)
        dataContainer.hRefolded[iHisto] = (TH2D*)dataContainer.hUnfoldedKinCorrected[iHisto]->Clone(Form("hRefolded_iter%zu",iHisto+1));
        dataContainer.hRefolded[iHisto]->SetTitle("Refolded with particle and detector level #varepsilon_{kin} corrections;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}");
        dataContainer.hRefolded[iHisto]->Sumw2(); // Ensure errors are propagated correctly
        dataContainer.hKineEffParticle[2]->Sumw2();

        // Apply detector level kinematic efficiency correction (multiply)
        dataContainer.hRefolded[iHisto]->Multiply(dataContainer.hKineEffParticle[2]);

        // Fold it with the response matrix
        dataContainer.hRefolded[iHisto] = manualFolding(dataContainer.response, dataContainer.hRefolded[iHisto], dataContainer.hMeasuredTemplate);

        // Apply particle level kinematic efficiency correction (divide)
        //dataContainer.hRefolded[iHisto]->Divide(dataContainer.hKineEffDetector[2]);

        dataContainer.hRefolded[iHisto]->SetTitle(Form("Refolded with %zu iterations", iHisto+1));
        // Clean up intermediate clone (optional if not reused)
        //delete hRefoldInput;

    }
    
    

    std::cout << "Refolding test performed." << std::endl;
}

void closureTest(UnfoldData& dataContainer) {
    // Build MC input sample (20%): matched data with detector level to test and particle level to compare with

    // Build MC correction sample (80%): build all MC level correction objects (efficiencies and response matrices)

    // 1 - Unfolding closure test

    // 2 - Unfolding closure test + sideband subtraction correction steps

    // 3 - Unfolding closure test + sideband subtraction correction steps + feed-down subtraction steps

    std::cout << "Closure test performed." << std::endl;
}

void plotHistograms(const UnfoldData& dataContainer, const double& jetptMin, const double& jetptMax) {
    cout << "Plotting histograms...\n";

    gStyle->SetPalette(kRainbow);

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    TCanvas* cKinEff = new TCanvas("cKinEff","Kinematic efficiencies");
    cKinEff->Divide(2,2);
    cKinEff->cd(1);
    dataContainer.hKineEffDetector[2]->Draw("text");
    cKinEff->cd(2);
    dataContainer.hKineEffParticle[2]->Draw("text");

    TCanvas* cResponse = new TCanvas("cResponse","Response matrix 2D representation");
    cResponse->Divide(2,2);
    cResponse->cd(1);
    const TH2* hresponse2D = dataContainer.response->Hresponse();
    TH2D* hresponse2DClone = static_cast<TH2D*>(hresponse2D->Clone("hResponse2D"));
    hresponse2DClone->SetTitle("2D response matrix from 4D RooUnfoldResponse - prompt D^{0}'s;2D Reconstructed;2D Truth");
    hresponse2DClone->Draw("colz");
    cResponse->cd(2);
    dataContainer.hResponse2D[0]->Draw("colz");
    cResponse->cd(3);
    dataContainer.hResponse2D[1]->Draw("colz");

    TCanvas* cUnfoldedIter = new TCanvas("cUnfoldedIter","Unfolded histograms for each iteration");
    cUnfoldedIter->cd();
    TLegend* lUnfoldedIter = new TLegend(0.6,0.57,0.7,0.77);
    std::vector<TH1D*> hUnfoldedKinCorrectedProj(dataContainer.hUnfoldedKinCorrected.size());
    int secondBin = 2;
    int lastButOneBin = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetNbins() - 1; // penultimate bin
    for (size_t iIter = 0; iIter < dataContainer.hUnfoldedKinCorrected.size(); iIter++) {
        // Project Y axis (DeltaR) using X axis (pTjet) range [secondBin, oneBeforeLastBin] (excluding the padding bins)
        hUnfoldedKinCorrectedProj[iIter] = dataContainer.hUnfoldedKinCorrected[iIter]->ProjectionY(Form("hProjIter_%zu", iIter),secondBin, lastButOneBin); // bins specified in X (pT,jet) dimension
        hUnfoldedKinCorrectedProj[iIter]->SetLineColor(kBlack + iIter);
        lUnfoldedIter->AddEntry(hUnfoldedKinCorrectedProj[iIter],Form("Iteration %zu", iIter+1), "le");
        if (iIter == 0) {
            hUnfoldedKinCorrectedProj[iIter]->Draw();
        } else {
            hUnfoldedKinCorrectedProj[iIter]->Draw("same");
        }
        
    }
    lUnfoldedIter->Draw();
    double reportingJetPtMin = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetBinLowEdge(secondBin);
    double reportingJetPtMax = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));

    TCanvas* cSteps = new TCanvas("cSteps","Unfolding steps");
    cSteps->Divide(4,2);
    cSteps->cd(1);
    dataContainer.hBFedDownData->Draw("colz");
    cSteps->cd(2);
    dataContainer.hBFedDownData->ProjectionY("hBFedDownDataProjY")->Draw();
    cSteps->cd(3);
    dataContainer.hBFedDownDataKinCorrected->Draw("colz");
    cSteps->cd(4);
    dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjY")->Draw();
    cSteps->cd(5);
    dataContainer.hUnfolded[7]->Draw("colz");
    cSteps->cd(6);
    dataContainer.hUnfolded[7]->ProjectionY("hUnfoldedProjY")->Draw();
    cSteps->cd(7);
    dataContainer.hUnfoldedKinCorrected[7]->Draw("colz");
    cSteps->cd(8);
    dataContainer.hUnfoldedKinCorrected[7]->ProjectionY("hUnfoldedKinCorrectedProjY")->Draw();
    
    TCanvas* cRefoldedIter = new TCanvas("cRefoldedIter","Refolded histograms for each iteration");
    cRefoldedIter->cd();
    TLegend* lRefoldedIter = new TLegend(0.5,0.57,0.85,0.87);
    for (size_t iHisto = 0; iHisto < dataContainer.hRefolded.size(); iHisto++)
    {
        dataContainer.hRefolded[iHisto]->SetLineColor(kBlack + iHisto);
        lRefoldedIter->AddEntry(dataContainer.hRefolded[iHisto],Form("Refolded Iteration %zu with #epsilon_{kin}^{part} correction", iHisto+1), "le");
        if (iHisto == 0) {
            dataContainer.hRefolded[iHisto]->ProjectionY(Form("hRefolded_iter%zu_Proj",iHisto+1))->Sumw2();
            dataContainer.hRefolded[iHisto]->ProjectionY(Form("hRefolded_iter%zu_Proj",iHisto+1))->Draw();
        } else {
            dataContainer.hRefolded[iHisto]->ProjectionY(Form("hRefolded_iter%zu_Proj",iHisto+1))->Sumw2();
            dataContainer.hRefolded[iHisto]->ProjectionY(Form("hRefolded_iter%zu_Proj",iHisto+1))->Draw("same");
        }
    }
    dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjY")->SetLineWidth(2);
    lRefoldedIter->AddEntry(dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjY"),"Original data with #epsilon_{kin}^{det} correction", "le");
    dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjY")->Draw("same");
    lRefoldedIter->Draw();

    TCanvas* cRefoldedPtRanges = new TCanvas("cRefoldedPtRanges", "Refolded histograms for last iteration in different pT,jet ranges");
    cRefoldedPtRanges->cd();
    TLegend* lRefoldedPtRanges = new TLegend(0.5,0.57,0.85,0.87);
    // Get last iteration refolded histogram
    auto* hLastRefolded = dataContainer.hRefolded.back();
    // Define bin ranges
    int bin7  = hLastRefolded->GetXaxis()->FindBin(7.0);
    int bin15 = hLastRefolded->GetXaxis()->FindBin(15.0);
    int bin30 = hLastRefolded->GetXaxis()->FindBin(30.0);
    // Create projections for refolded
    auto* hProj_7_15 = hLastRefolded->ProjectionY("hRefolded_iter_8_Proj_7_15GeV", bin7, bin15);
    hProj_7_15->SetLineColor(kBlue);
    hProj_7_15->Draw();
    lRefoldedPtRanges->AddEntry(hProj_7_15, "Refolded: 7 < p_{T,jet} < 15", "le");
    auto* hProj_15_30 = hLastRefolded->ProjectionY("hRefolded_iter_8_Proj_15_30GeV", bin15+1, bin30);
    hProj_15_30->SetLineColor(kRed);
    hProj_15_30->Draw("same");
    lRefoldedPtRanges->AddEntry(hProj_15_30, "Refolded: 15 < p_{T,jet} < 30", "le");
    // Create projections for original (BFedDown) data
    bin7  = dataContainer.hBFedDownDataKinCorrected->GetXaxis()->FindBin(7.0);
    bin15 = dataContainer.hBFedDownDataKinCorrected->GetXaxis()->FindBin(15.0);
    bin30 = dataContainer.hBFedDownDataKinCorrected->GetXaxis()->FindBin(30.0);
    auto* hOrig_7_15 = dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjY_7_15GeV", bin7, bin15);
    hOrig_7_15->SetLineColor(kBlue);
    hOrig_7_15->SetLineStyle(2); // dashed line for distinction
    hOrig_7_15->Draw("same");
    lRefoldedPtRanges->AddEntry(hOrig_7_15, "Original: 7 < p_{T,jet} < 15", "le");
    auto* hOrig_15_30 = dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjY_15_30GeV", bin15+1, bin30);
    hOrig_15_30->SetLineColor(kRed);
    hOrig_15_30->SetLineStyle(2);
    hOrig_15_30->Draw("same");
    lRefoldedPtRanges->AddEntry(hOrig_15_30, "Original: 15 < p_{T,jet} < 30", "le");
    lRefoldedPtRanges->Draw();

    //
    // Storing images
    //
    TString imagePath = "../Images/4-Unfolding/";
    cKinEff->Update();
    cKinEff->SaveAs(imagePath + "Unfolding_kin_efficiencies.png");
    cResponse->Update();
    cResponse->SaveAs(imagePath + "Unfolding_response.png");
    cUnfoldedIter->Update();
    cUnfoldedIter->SaveAs(imagePath + "Unfolding_unfolded_iterations.png");
    cRefoldedIter->Update();
    cRefoldedIter->SaveAs(imagePath + "Unfolding_refolded_comparison.png");
    cRefoldedPtRanges->Update();
    cRefoldedPtRanges->SaveAs(imagePath + "Unfolding_refolded_pt_ranges.png");
    
    //
    // Storing in a single pdf file
    //
    cKinEff->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cResponse->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cUnfoldedIter->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cRefoldedIter->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cRefoldedPtRanges->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
    //cKinEfficiency->Print(imagePath + Form("feeddown_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    //cFolded->Print(imagePath + Form("feeddown_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    //cFedDownData->Print(imagePath + Form("feeddown_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

}

void saveData(const UnfoldData& dataContainer, const double& jetptMin, const double& jetptMax) {
    // Open output file
    TFile* outFile = new TFile(Form("unfolding_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");

    // Store each histogram in file

    // Selection efficiency
    dataContainer.hSelectionEfficiency.first->Write();

    // Kinematic efficiency
    dataContainer.hKineEffParticle[2]->Write();
    dataContainer.hKineEffDetector[2]->Write();

    // 2D response matrices projections
    dataContainer.hResponse2D[0]->Write();
    dataContainer.hResponse2D[1]->Write();

    // Create subdirectories
    TDirectory* dirUnfolded = outFile->mkdir("Unfolded");
    TDirectory* dirRefolded = outFile->mkdir("Refolded");

    for (size_t iHisto = 0; iHisto < dataContainer.hUnfoldedKinCorrected.size(); iHisto++) {
        // Unfolded histograms in "Unfolded" folder
        dirUnfolded->cd();
        dataContainer.hUnfoldedKinCorrected[iHisto]->Write();
        // Refolded histograms in "Refolded" folder
        dirRefolded->cd();
        dataContainer.hRefolded[iHisto]->Write();
    }
    
    
    outFile->Close();
    delete outFile;
    
    cout << "Data stored in file" << Form("unfolding_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)) << endl;
}

void Unfolding(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // Number of unfolding procedure iterations
    int iterationNumber = 8;

    // jet pT cuts
    std::vector<double> ptjetBinEdges_particle = {5., 7., 15., 30., 50.};
    std::vector<double> ptjetBinEdges_detector = {5., 7., 15., 30., 50.};
    double jetptMin = ptjetBinEdges_particle[0]; // GeV
    double jetptMax = ptjetBinEdges_particle[ptjetBinEdges_particle.size() - 1]; // GeV

    // deltaR histogram
    std::vector<double> deltaRBinEdges_particle = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15,0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5}; // default = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4} chosen by Nima
    std::vector<double> deltaRBinEdges_detector = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15,0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5}; // default = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4} chosen by Nima
    double minDeltaR = deltaRBinEdges_particle[0];
    double maxDeltaR = deltaRBinEdges_particle[deltaRBinEdges_particle.size() - 1];
    
    // pT,D histograms
    std::vector<double> ptDBinEdges_particle = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.}; // default pT,D = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.}
    std::vector<double> ptDBinEdges_detector = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.}; // default pT,D = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.}
    double hfptMin = ptDBinEdges_particle[0];
    double hfptMax = ptDBinEdges_particle[ptDBinEdges_particle.size() - 1];

    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.03 (from first of pair)
    //std::vector<std::pair<double, double>> bdtPtCuts = {
    //    {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}
    //};
    // Dataset: JE_HF_LHC24g5_All_D0
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.4}, {8, 0.6}, {10, 0.8}, {12, 0.8}, {16, 1.0}, {50, 1.0}
    };

    // Opening files
    TFile* fSimulatedMCMatched = new TFile("../SimulatedData/Hyperloop_output/Train_runs/410603_Match/AO2D_mergedDFs.root","read");
    TFile* fEfficiency = new TFile(Form("../2-Efficiency/selection_efficiency_run3style_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fData = new TFile("../ExperimentalData/Hyperloop_output/HF_LHC23_pass4_Thin_small_2P3PDstar_DATA_newMLModel/AnalysisResults.root","read");
    TFile* fFeedDown = new TFile(Form("../3-Feed-Down/outputFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fSimulatedMCMatched || fSimulatedMCMatched->IsZombie()) {
        std::cerr << "Error: Unable to open O2 MC matched ROOT file." << std::endl;
    }
    if (!fEfficiency || fEfficiency->IsZombie()) {
        std::cerr << "Error: Unable to open estimated selection efficiency ROOT file." << std::endl;
    }
    if (!fData || fData->IsZombie()) {
        std::cerr << "Error: Unable to open AnalysisResults.root data ROOT file." << std::endl;
    }
    if (!fFeedDown || fFeedDown->IsZombie()) {
        std::cerr << "Error: Unable to open AnalysisResults.root data ROOT file." << std::endl;
    }
    
    //UnfoldData dataContainer = createHistograms(deltaRBinEdges, ptDBinEdges, jetptMin, jetptMax);
    UnfoldData dataContainer = createHistograms(ptjetBinEdges_particle, deltaRBinEdges_particle, ptDBinEdges_particle,
                                                ptjetBinEdges_detector, deltaRBinEdges_detector, ptDBinEdges_detector,
                                                iterationNumber);

    // Fill histograms with POWHEG simulation data
    fillHistograms(fFeedDown, fEfficiency, fSimulatedMCMatched, dataContainer, jetptMin, jetptMax, hfptMin, hfptMax,
                   deltaRBinEdges_particle, deltaRBinEdges_detector, bdtPtCuts);

    performUnfolding(dataContainer, iterationNumber);

    refoldingTest(dataContainer);

    // Plot the efficiency histogram and further corrected histograms
    plotHistograms(dataContainer, jetptMin, jetptMax);

    // Save corrected distributions to file
    saveData(dataContainer, jetptMin, jetptMax);


    

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    Unfolding();
    return 0;
}
