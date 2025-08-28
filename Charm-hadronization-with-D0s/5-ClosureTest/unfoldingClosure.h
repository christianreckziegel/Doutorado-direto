/*
 * Macro for performing unfolding in second closure test
 * 
 * 
 * 
 * Author: Christian Reckziegel
**/

using namespace std;

struct UnfoldData {
    // Step 2: apply efficiency (pT,D dependent) correction
    std::pair<TH1D*, TH1D*> hSelectionEfficiency;                           // first = prompt D0s, second = non-prompt D0s

    // Step 3: background subtracted, efficiency corrected, B fed-down 2D distribution (pT,jet vs DeltaR)
    TH2D* hBFedDownData;

    TH2D* hEfficiencyCorrectedInput;

    // Step 4: bayesian unfolding
    std::vector<TH2D*> hKineEffParticle = {nullptr, nullptr, nullptr};      // particle level (addition): [0] = numerator, [1] = denominator, [2] = efficiency
    std::vector<TH2D*> hKineEffDetector = {nullptr, nullptr, nullptr};      // detector level (removal): [0] = numerator, [1] = denominator, [2] = efficiency
    TH2D* hBFedDownDataKinCorrected = nullptr;                            // 2D histogram with detector level kinematic efficiency correction
    RooUnfoldResponse* response;                                             // response matrix for folding
    std::vector<RooUnfoldBayes*> unfold;                                     // unfolding objects, there are iterationNumber unfolding objects
    TH2D* hMeasuredTemplate = nullptr;                                      // template for measured data (jet pT vs DeltaR), used for hUnfolded binning
    TH2D* hTruthTemplate = nullptr;
    std::vector<TH2D*> hUnfolded;                                         // unfolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects
    std::vector<TH2D*> hUnfoldedKinCorrected;                             // unfolded 2D histogram with detector level kinematic correction, there are iterationNumber unfolding objects

    // Refolding test
    std::vector<TH2D*> hRefolded;                                           // refolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects
};

UnfoldData calculateKinematics(TFile* fClosureInputMatched, std::vector<TH1D*>& hSelEff_run3style, const BinningStruct& binningStruct, const std::vector<std::pair<double, double>>& bdtPtCuts, int& iterationNumber) {
    UnfoldData dataContainer;
    
    // 1 ----- Create histograms
    // Create kinematic efficiency histograms
    dataContainer.hKineEffParticle[0] = new TH2D("hKineEffParticleNumerator", "Particle level kinematic efficiency numerator (prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  binningStruct.ptjetBinEdges_particle.size() - 1, binningStruct.ptjetBinEdges_particle.data(), 
                                                                                  binningStruct.deltaRBinEdges_particle.size() - 1, binningStruct.deltaRBinEdges_particle.data());
    dataContainer.hKineEffParticle[1] = new TH2D("hKineEffParticleDenominator", "Particle level kinematic efficiency denominator (prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  binningStruct.ptjetBinEdges_particle.size() - 1, binningStruct.ptjetBinEdges_particle.data(), 
                                                                                  binningStruct.deltaRBinEdges_particle.size() - 1, binningStruct.deltaRBinEdges_particle.data());
    dataContainer.hKineEffDetector[0] = new TH2D("hKineEffDetectorNumerator", "Detector level kinematic efficiency numerator (prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  binningStruct.ptjetBinEdges_detector.size() - 1, binningStruct.ptjetBinEdges_detector.data(), 
                                                                                  binningStruct.deltaRBinEdges_detector.size() - 1, binningStruct.deltaRBinEdges_detector.data());
    dataContainer.hKineEffDetector[1] = new TH2D("hKineEffDetectorDenominator", "Detector level kinematic efficiency denominator (prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  binningStruct.ptjetBinEdges_detector.size() - 1, binningStruct.ptjetBinEdges_detector.data(), 
                                                                                  binningStruct.deltaRBinEdges_detector.size() - 1, binningStruct.deltaRBinEdges_detector.data());
    // Create template histograms used for response matrix creation
    dataContainer.hMeasuredTemplate = new TH2D("hMeasuredTemplate", "Measured template;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}", 
                                                                                  binningStruct.ptjetBinEdges_detector.size() - 1, binningStruct.ptjetBinEdges_detector.data(), 
                                                                                  binningStruct.deltaRBinEdges_detector.size() - 1, binningStruct.deltaRBinEdges_detector.data());
    dataContainer.hTruthTemplate = new TH2D("hTruthTemplate", "Truth template;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR^{gen}", 
                                                                                  binningStruct.ptjetBinEdges_detector.size() - 1, binningStruct.ptjetBinEdges_detector.data(), 
                                                                                  binningStruct.deltaRBinEdges_detector.size() - 1, binningStruct.deltaRBinEdges_detector.data());
    // Create response matrices
    dataContainer.response = new RooUnfoldResponse(dataContainer.hKineEffDetector[0], dataContainer.hKineEffParticle[0]);
    // Reserve space for the unfolding objects with the number of iterations
    dataContainer.unfold.resize(iterationNumber);
    // Reserve space for the unfolded TH2D* histograms with the number of iterations (and point to nullptr)
    dataContainer.hUnfolded.resize(iterationNumber, nullptr);
    // And also for the kinematic efficiency corrected version of the previous mentioned (and point to nullptr)
    dataContainer.hUnfoldedKinCorrected.resize(iterationNumber, nullptr);

    // 2 ----- Fill histograms
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = binningStruct.deltaRBinEdges_particle[binningStruct.deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = binningStruct.deltaRBinEdges_detector[binningStruct.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double jetptMin = binningStruct.ptjetBinEdges_particle[0]; // on both levels jet
    const double jetptMax = binningStruct.ptjetBinEdges_particle[binningStruct.ptjetBinEdges_particle.size() - 1]; // on both levels jet
    const double MCPHfPtMincut = binningStruct.ptDBinEdges_particle[0]; // on particle level D0
    const double MCDHfPtMincut = binningStruct.ptDBinEdges_detector[0]; // on detector level D0
    const double MCPHfPtMaxcut = binningStruct.ptDBinEdges_particle[binningStruct.ptDBinEdges_particle.size() - 1]; // on particle level D0
    const double MCDHfPtMaxcut = binningStruct.ptDBinEdges_detector[binningStruct.ptDBinEdges_detector.size() - 1]; // on detector level D0
    // Access prompt and non-prompt efficiency
    dataContainer.hSelectionEfficiency.first = hSelEff_run3style[0];
    dataContainer.hSelectionEfficiency.second = hSelEff_run3style[1];
    TTree* tree = (TTree*)fClosureInputMatched->Get("CorrectionTree");
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
        
        // Apply prompt selection (i.e., only c → D0, not B → D0)
        if (!MCDhfprompt) {
            continue;
        }

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= binningStruct.deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= binningStruct.deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
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

    return dataContainer;
}
// Detector level kinematic efficiency: removal of data
void removeOutsideData(UnfoldData& dataContainer) {
    // Calculate detector level kinematic efficiency histograms
    dataContainer.hKineEffDetector[2] = (TH2D*)dataContainer.hKineEffDetector[0]->Clone("hKineEffDetectorEfficiency");
    dataContainer.hKineEffDetector[2]->Sumw2();
    dataContainer.hKineEffDetector[2]->Divide(dataContainer.hKineEffDetector[1]); // A = A / B:  A = A->Divide(B)

    // Copy background subtracted, efficiency corrected, fed-down subtracted 2D histogram
    //dataContainer.hEfficiencyCorrected = (TH2D*)dataContainer.hEfficiencyCorrected->Clone("hEfficiencyCorrectedUnfoldingInput");
    dataContainer.hEfficiencyCorrectedInput->SetTitle("Fed-down with #varepsilon_{kin}^{detector} correction;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}");

    // Apply kinematic efficiency correction
    dataContainer.hEfficiencyCorrectedInput->Sumw2();
    dataContainer.hEfficiencyCorrectedInput->Multiply(dataContainer.hKineEffDetector[2]);

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

std::vector<TH2D*> performUnfolding(UnfoldData& dataContainer, TH2D* hEfficiencyCorrected, int& iterationNumber) {
    
    // Obtain input distribution
    dataContainer.hEfficiencyCorrectedInput = (TH2D*)hEfficiencyCorrected->Clone("hEfficiencyCorrectedUnfoldingInput");
    // Correct distribution with detector level kinematic efficiency (remove entries)
    removeOutsideData(dataContainer);

    // Calculate particle level kinematic efficiency histograms
    dataContainer.hKineEffParticle[2] = (TH2D*)dataContainer.hKineEffParticle[0]->Clone("hKineEffParticleEfficiency");
    dataContainer.hKineEffParticle[2]->Sumw2();
    dataContainer.hKineEffParticle[2]->Divide(dataContainer.hKineEffParticle[1]); // A = A / B:  A = A->Divide(B)

    // Unfold in multiple iterations
    for (int iIter = 1; iIter <= iterationNumber; iIter++) {
        dataContainer.unfold[iIter - 1] = new RooUnfoldBayes(dataContainer.response, dataContainer.hEfficiencyCorrectedInput, iIter);
        dataContainer.hUnfolded[iIter - 1] = (TH2D*) dataContainer.unfold[iIter - 1]->Hreco(RooUnfold::kCovariance); // kCovariance = 2 = Errors from the square root of of the covariance matrix given by the unfolding
        dataContainer.hUnfolded[iIter - 1]->SetTitle(Form("Immediately unfolded 2D histogram with %d iterations;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR^{gen}", iIter));
        dataContainer.hUnfolded[iIter - 1]->SetName(Form("hUnfolded_iter%d", iIter));
        dataContainer.hUnfolded[iIter - 1]->Sumw2();

        // Correct distribution with particle level kinematic efficiency (add entries)
        dataContainer.hUnfoldedKinCorrected[iIter - 1] = addOutsideData(dataContainer.hKineEffParticle[2],dataContainer.hUnfolded[iIter - 1], iIter);
        
    }
    std::cout << "Unfolding procedure performed." << std::endl;

    return dataContainer.hUnfoldedKinCorrected;
}

std::vector<TH2D*> UnfoldingClosure(TFile* fClosureInputMatched, std::vector<TH1D*>& hSelEff_run3style, TH2D* hEfficiencyCorrected, const BinningStruct& binningStruct, const std::vector<std::pair<double, double>>& bdtPtCuts) {

    int iterationNumber = 8;

    // Calculate response matrix and kinematic efficiencies
    UnfoldData dataContainer = calculateKinematics(fClosureInputMatched, hSelEff_run3style, binningStruct, bdtPtCuts, iterationNumber);

    // Unfold with "iterationNumber" iterations
    std::vector<TH2D*> hUnfKinCorrected = performUnfolding(dataContainer, hEfficiencyCorrected, iterationNumber);

    return hUnfKinCorrected;
}
