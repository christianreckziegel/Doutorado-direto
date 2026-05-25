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
    std::pair<TH1D*, TH1D*> hSelectionEfficiency_run3_particleLevel;        // first = prompt D0s, second = non-prompt D0s

    // Step 3: input (background subtracted, efficiency corrected, B fed-down) 2D distribution (pT,jet vs DeltaR) with corrections applied to be unfolded
    TH2D* hBFedDownData;

    TH2D* hCorrectedInput;

    // Step 4: bayesian unfolding
    std::vector<TH2D*> hKineEffParticle = {nullptr, nullptr, nullptr};      // particle level (addition): [0] = numerator, [1] = denominator, [2] = efficiency
    std::vector<TH2D*> hKineEffDetector = {nullptr, nullptr, nullptr};      // detector level (removal): [0] = numerator, [1] = denominator, [2] = efficiency
    RooUnfoldResponse* response;                                            // response matrix for folding
    std::vector<TH2D*> hResponse2D = {nullptr, nullptr};                    // response projections matrix: first = DeltaR, second = pT,jet
    std::vector<RooUnfoldBayes*> unfold;                                    // unfolding objects, there are iterationNumber unfolding objects
    TH2D* hMeasuredTemplate = nullptr;                                      // template for measured data (jet pT vs DeltaR), used for hUnfolded binning
    TH2D* hTruthTemplate = nullptr;
    std::vector<TH2D*> hUnfolded;                                           // unfolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects
    std::vector<TH2D*> hUnfoldedKinCorrected;                               // unfolded 2D histogram with detector level kinematic correction, there are iterationNumber unfolding objects

    // Refolding test
    std::vector<TH2D*> hRefolded;                                           // refolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects
};

UnfoldData calculateKinematics(TFile* fClosureInput, const std::vector<TH1D*>& hSelEff_run3style, const std::vector<TH1D*>& hSelEff_run3_particleLevel, const BinningStruct& binning, int& iterationNumber) {
    UnfoldData dataContainer;
    
    // 1 ----- Create histograms
    // Create kinematic efficiency histograms
    dataContainer.hKineEffParticle[0] = new TH2D("hKineEffParticleNumerator", "Particle level kinematic efficiency numerator (prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hKineEffParticle[1] = new TH2D("hKineEffParticleDenominator", "Particle level kinematic efficiency denominator (prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hKineEffDetector[0] = new TH2D("hKineEffDetectorNumerator", "Detector level kinematic efficiency numerator (prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());
    dataContainer.hKineEffDetector[1] = new TH2D("hKineEffDetectorDenominator", "Detector level kinematic efficiency denominator (prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());
    // Create template histograms used for response matrix creation
    dataContainer.hMeasuredTemplate = new TH2D("hMeasuredTemplate", "Measured template;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());
    dataContainer.hTruthTemplate = new TH2D("hTruthTemplate", "Truth template;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR^{gen}", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());
    // Create response matrices
    dataContainer.response = new RooUnfoldResponse(dataContainer.hKineEffDetector[0], dataContainer.hKineEffParticle[0]);
    dataContainer.hResponse2D[0] = new TH2D("hResponseDeltaR", "Response matrix projection on #DeltaR;#DeltaR^{reco};#DeltaR^{gen}", 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hResponse2D[1] = new TH2D("hResponsePtJet", "Response matrix projection on p_{T,jet};p_{T,jet}^{reco} (GeV/#it{c});p_{T,jet}^{gen} (GeV/#it{c})", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data());
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
    const double MCPDeltaRcut = binning.deltaRBinEdges_particle[binning.deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double jetptMin = binning.ptjetBinEdges_particle[0]; // on both levels jet
    const double jetptMax = binning.ptjetBinEdges_particle[binning.ptjetBinEdges_particle.size() - 1]; // on both levels jet
    const double MCPHfPtMincut = binning.ptHFBinEdges_particle[0]; // on particle level D0
    const double MCDHfPtMincut = binning.ptHFBinEdges_detector[0]; // on detector level D0
    const double MCPHfPtMaxcut = binning.ptHFBinEdges_particle[binning.ptHFBinEdges_particle.size() - 1]; // on particle level D0
    const double MCDHfPtMaxcut = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1]; // on detector level D0
    
    // Access prompt and non-prompt efficiency
    dataContainer.hSelectionEfficiency.first = hSelEff_run3style[0];
    dataContainer.hSelectionEfficiency.second = hSelEff_run3style[1];
    dataContainer.hSelectionEfficiency_run3_particleLevel.first = hSelEff_run3_particleLevel[0];
    dataContainer.hSelectionEfficiency_run3_particleLevel.second = hSelEff_run3_particleLevel[1];

    // Accessing TTree
    TTree* tree = (TTree*)fClosureInput->Get("CorrectionTree");
    // Check for correct access
    if (!tree) {
        cout << "Error opening O2 matching tree.\n";
    }
    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    int MCPjetNConst;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    int MCDhfMatchedFrom, MCDhfSelectedAs;
    int MCDjetNConst;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;

    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetNConst); // float
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

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // Only prompt, real (not reflections) and matched candidates
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isPrompt = MCDhfprompt || MCPhfprompt;                                 // both need to be non-prompt in order to reject the candidate
        if (!MCPhfprompt || !MCDhfmatch || !isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs)) {
            continue;
        }
  
        // Generator level selection cuts
        double MCPDeltaR = MCPaxisDistance;
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);
        bool genHfPtRange = ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        bool genDeltaRRange = ((MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < MCPDeltaRcut));
        bool genAcceptance = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut);
        bool genLevelRange = genAcceptance && genJetPtRange && genDeltaRRange && genHfPtRange; // --> this is new!
        // Reconstruction level selection cuts
        double MCDDeltaR = MCDaxisDistance;
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb);
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        bool recoDeltaRRange = ((MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < MCDDeltaRcut));
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && recoJetPtRange && recoHfPtRange && recoDeltaRRange;
        
        // Find the bin corresponding to the given pT,D value
        int bin = dataContainer.hSelectionEfficiency.first->FindBin(MCDhfPt);
        // Get the efficiency value from the bin content
        double efficiency_prompt = dataContainer.hSelectionEfficiency.first->GetBinContent(bin);

        // Fill response matrix and kinematic efficiency histograms
        if (genLevelRange && recoLevelRange && passBDTcut) { // use matched selected sample to describe kinematic smearing requires correction with the inverse of selection efficiency

            // Fill 4D RooUnfoldResponse object (jet pT shape is influenced by D0 pT efficiency)
            dataContainer.response->Fill(MCDjetPt, MCDDeltaR, MCPjetPt, MCPDeltaR, 1./efficiency_prompt);
            
            // Fill response matrix projections
            dataContainer.hResponse2D[0]->Fill(MCDDeltaR, MCPDeltaR, 1./efficiency_prompt);
            dataContainer.hResponse2D[1]->Fill(MCDjetPt, MCPjetPt, 1./efficiency_prompt);

            
        }
        // No need for passBDTCut condition in the following, since we want to fill the kinematic efficiency histograms for all the entries that are inside the kinematic range,
        // independently of whether they pass the BDT cut or not, since the kinematic efficiency is meant to correct for the geometrical and kinematic acceptance of the detector,
        // not for the BDT selection efficiency (which is already corrected for in the fed-down subtraction step).
        int particleBin = dataContainer.hSelectionEfficiency_run3_particleLevel.first->FindBin(MCPhfPt);  // ← use MCPhfPt
        double eff_particle = dataContainer.hSelectionEfficiency_run3_particleLevel.first->GetBinContent(particleBin);
        if (genLevelRange && recoLevelRange) {
            // Fill kinematic efficiency numerator histograms
            dataContainer.hKineEffParticle[0]->Fill(MCPjetPt, MCPDeltaR);
            dataContainer.hKineEffDetector[0]->Fill(MCDjetPt, MCDDeltaR);
        }
        if (genLevelRange) {
            
            // Fill kinematic efficiency denominator histogram for full particle level range
            dataContainer.hKineEffParticle[1]->Fill(MCPjetPt, MCPDeltaR);
        }
        if (recoLevelRange) {
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
    //dataContainer.hEfficiencyCorrectedData = (TH2D*)dataContainer.hEfficiencyCorrectedData->Clone("hEfficiencyCorrectedDataUnfoldingInput");
    dataContainer.hCorrectedInput->SetTitle("Fed-down with #varepsilon_{kin}^{detector} correction;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}");

    // Apply kinematic efficiency correction
    dataContainer.hCorrectedInput->Sumw2();
    dataContainer.hCorrectedInput->Multiply(dataContainer.hKineEffDetector[2]); // test with perfect kinematic efficiency

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
    hUnfoldedKinCorrected->Divide(hKineEffParticle); // A = A / B:  A = A->Divide(B), // test with perfect kinematic efficiency

    std::cout << "Particle level kinematic efficiency applied (addition)." << std::endl;
    return hUnfoldedKinCorrected;

}

std::vector<TH2D*> performUnfolding(UnfoldData& dataContainer, TH2D* hEfficiencyCorrectedData, int& iterationNumber) {
    
    // Obtain input distribution
    dataContainer.hCorrectedInput = (TH2D*)hEfficiencyCorrectedData->Clone("hUnfoldingInput");
    // Correct distribution with detector level kinematic efficiency (remove entries)
    removeOutsideData(dataContainer);

    // Calculate particle level kinematic efficiency histograms
    dataContainer.hKineEffParticle[2] = (TH2D*)dataContainer.hKineEffParticle[0]->Clone("hKineEffParticleEfficiency");
    dataContainer.hKineEffParticle[2]->Sumw2();
    dataContainer.hKineEffParticle[2]->Divide(dataContainer.hKineEffParticle[1]); // A = A / B:  A = A->Divide(B)

    // Unfold in multiple iterations
    for (int iIter = 1; iIter <= iterationNumber; iIter++) {
        dataContainer.unfold[iIter - 1] = new RooUnfoldBayes(dataContainer.response, dataContainer.hCorrectedInput, iIter);
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

std::vector<TH1D*> convergenceTest(const std::vector<TH2D*>& hUnfKinCorrected) {
    // To be implemented if needed
    std::vector<TH1D*> hConvergenceTest;

    // Compute ratios between successive iterations
    for (size_t iHisto = 1; iHisto < hUnfKinCorrected.size(); iHisto++) {

        // Project 2D histograms to 1D (along X axis)
        TH1D* hCurrent = hUnfKinCorrected[iHisto]->ProjectionY(Form("projCurrent_%zu", iHisto));
        TH1D* hPrevious = hUnfKinCorrected[iHisto-1]->ProjectionY(Form("projPrevious_%zu", iHisto-1));

        // Clone the current histogram to hold the ratio
        TH1D* hRatio = (TH1D*)hCurrent->Clone(Form("ratioIter_%zu", iHisto));
        hRatio->SetTitle(Form("Iteration %zu / %zu ratio", iHisto, iHisto-1));
        
        // Divide by the previous iteration
        hRatio->Divide(hPrevious);
        
        // Optional: style
        hRatio->SetLineColor(kBlack + iHisto);
        hRatio->SetMarkerStyle(20);
        
        hConvergenceTest.push_back(hRatio);
    }

    // Quantitative check using RMS deviation
    std::vector<double> rmsValues;
    for (size_t iHisto = 0; iHisto < hConvergenceTest.size(); iHisto++) {
        TH1D* hRatio = hConvergenceTest[iHisto];
        double rms = 0;
        int nBins = hRatio->GetNbinsX();
        for (int iBin = 1; iBin <= nBins; iBin++) {
            double val = hRatio->GetBinContent(iBin) - 1.0; // deviation from 1
            rms += val * val;
        }
        rms = sqrt(rms / nBins);
        rmsValues.push_back(rms);
        std::cout << "Iteration " << iHisto+1 << "/" << iHisto << " RMS deviation from 1: " << rms << std::endl;
    }

    // Plot
    TCanvas* cRatio = new TCanvas("cRatioUnfoldingConvergence", "Unfolding Convergence Ratios",1920,1080);
    cRatio->cd();

    for (size_t iHisto = 0; iHisto < hConvergenceTest.size(); iHisto++) {
        if (iHisto == 0) {
            hConvergenceTest[iHisto]->Draw("E");      // draw the first ratio
        } else {
            hConvergenceTest[iHisto]->Draw("E SAME");        // draw others on the same canvas
        }
    }

    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t iHisto = 0; iHisto < hConvergenceTest.size(); iHisto++) {
        legend->AddEntry(hConvergenceTest[iHisto], Form("Iteration %zu / %zu's RMS: %.3f", iHisto+1, iHisto, rmsValues[iHisto]), "l");
    }
    legend->Draw();

    return hConvergenceTest;

}


void plotHistograms(const UnfoldData& dataContainer, const BinningStruct& binning) {
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    
    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);
    
    // Kinematic efficiencies
    TCanvas* cKinematicEffPart = new TCanvas("cKinematicEffPart","Particle level kinematic efficiency",1920,1080);
    cKinematicEffPart->cd();
    dataContainer.hKineEffParticle[2]->SetTitle("Particle level kinematic efficiency (prompt D^{0}'s)");
    dataContainer.hKineEffParticle[2]->Draw("text");
    TCanvas* cKinematicEffDet = new TCanvas("cKinematicEffDet","Detector level kinematic efficiency",1920,1080);
    cKinematicEffDet->cd();
    dataContainer.hKineEffDetector[2]->SetTitle("Detector level kinematic efficiency (prompt D^{0}'s)");
    dataContainer.hKineEffDetector[2]->Draw("text");

    TCanvas* cReponseProjDeltaR = new TCanvas("cReponseProjDeltaR","DeltaR response matrix projection",1920,1080);
    cReponseProjDeltaR->cd();
    dataContainer.hResponse2D[0]->Draw("colz");
    TCanvas* cReponseProjPtjet = new TCanvas("cReponseProjPtjet","Jet transverse momentum response matrix projection",1920,1080);
    cReponseProjPtjet->cd();
    dataContainer.hResponse2D[1]->Draw("colz");

    TCanvas* cFullyCorrected2D = new TCanvas("cFullyCorrected2D","Fully corrected deltaR 2D distribution",1920,1080);
    dataContainer.hUnfoldedKinCorrected[3]->Draw("colz"); // 4th iteration is considered good enough through convergence test result?

    TCanvas* cFullyCorrected1D = new TCanvas("cFullyCorrected1D","Fully corrected deltaR 1D distribution iterations",1920,1080);
    TLegend* lUnfoldedIter = new TLegend(0.6,0.57,0.7,0.77);
    std::vector<TH1D*> hUnfoldedKinCorrectedProj(dataContainer.hUnfoldedKinCorrected.size());
    int secondBin = 2;
    int lastButOneBin = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetNbins() - 1; // penultimate bin
    for (size_t iIter = 0; iIter < dataContainer.hUnfoldedKinCorrected.size(); iIter++) {
        cFullyCorrected1D->cd();
        hUnfoldedKinCorrectedProj[iIter] = dataContainer.hUnfoldedKinCorrected[iIter]->ProjectionY(Form("hUnfoldedKinCorrected_iter%zu",iIter), secondBin, lastButOneBin);
        hUnfoldedKinCorrectedProj[iIter]->SetLineColor(kBlack + iIter);
        lUnfoldedIter->AddEntry(hUnfoldedKinCorrectedProj[iIter],Form("Iteration %zu", iIter+1), "le");
        if (iIter == 0) {
            hUnfoldedKinCorrectedProj[iIter]->GetYaxis()->SetRangeUser(0.,hUnfoldedKinCorrectedProj[iIter]->GetMaximum()*1.2);
            hUnfoldedKinCorrectedProj[iIter]->SetTitle("Fully corrected distribution projection");
            hUnfoldedKinCorrectedProj[iIter]->GetYaxis()->SetTitle("dN");
            hUnfoldedKinCorrectedProj[iIter]->Draw();
        } else {
            hUnfoldedKinCorrectedProj[iIter]->Draw("same");
        }
    }
    double reportingJetPtMin = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetBinLowEdge(secondBin);
    double reportingJetPtMax = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));

    
    //
    // Storing images
    //
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../Images/5-ClosureTest/Second/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    cKinematicEffPart->Print(imagePath + Form("closureTest2_unfolding_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cKinematicEffDet->Print(imagePath + Form("closureTest2_unfolding_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cReponseProjDeltaR->Print(imagePath + Form("closureTest2_unfolding_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cReponseProjPtjet->Print(imagePath + Form("closureTest2_unfolding_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFullyCorrected2D->Print(imagePath + Form("closureTest2_unfolding_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFullyCorrected1D->Print(imagePath + Form("closureTest2_unfolding_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
}

UnfoldData UnfoldingClosure(TFile* fClosureInput, TH2D* hEfficiencyCorrectedData, const std::vector<TH1D*>& hSelEff_run3style, const std::vector<TH1D*>& hSelEff_run3_particleLevel, const BinningStruct& binning) {

    int iterationNumber = 8;

    // Calculate response matrix and kinematic efficiencies
    UnfoldData dataContainer = calculateKinematics(fClosureInput, hSelEff_run3style, hSelEff_run3_particleLevel, binning, iterationNumber);

    // Unfold with "iterationNumber" iterations
    std::vector<TH2D*> hUnfKinCorrected = performUnfolding(dataContainer, hEfficiencyCorrectedData, iterationNumber);

    // Convergence test
    std::vector<TH1D*> hConvergenceTest = convergenceTest(hUnfKinCorrected);

    // Plot response matrix projectios and histograms
    plotHistograms(dataContainer, binning);

    return dataContainer;
}
