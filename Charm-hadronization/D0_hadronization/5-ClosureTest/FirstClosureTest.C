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

#include "../commonUtilities.h"

using namespace std;

struct ClosureTestData1 {

    // Input objects
    TH2D* hInputParticle = nullptr;
    TH2D* hInputDetector = nullptr;
    TH2D* hInputDetectorRaw = nullptr;

    // Correction objects
    RooUnfoldResponse* response;                                        // response matrix
    std::vector<TH2D*> hResponseProj = {nullptr, nullptr};              // response matrix projections -> 0 = deltaR, 1 = pT,jet
    std::vector<TH2D*> hResponseProjDeltaR;                             // DeltaR response matrix projections for each pT,jet bin
    std::vector<TH2D*> hKineEffParticle = {nullptr, nullptr, nullptr};  // particle level kinematic efficiency -> 0 = numerator, 1 = denominator, 2 = efficiency
    std::vector<TH2D*> hKineEffDetector = {nullptr, nullptr, nullptr};  // detector level kinematic efficiency -> 0 = numerator, 1 = denominator, 2 = efficiency
    
    // Unfolding objects
    std::vector<RooUnfoldBayes*> unfold;                                // unfolding objects, there are iterationNumber unfolding objects
    std::vector<TH2D*> hUnfolded;                                       // unfolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects
};

ClosureTestData1 createAnalysisObjects(TFile* fClosureInput, double& jetptMin, double& jetptMax, double& hfptMin, double& hfptMax, const BinningStruct& binning) {
    // Create empty struct to store data
    ClosureTestData1 dataContainer;
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = binning.deltaRBinEdges_particle[binning.deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double MCPHfPtMincut = hfptMin; // on particle level D0
    const double MCDHfPtMincut = hfptMin; // on detector level D0
    const double MCPHfPtMaxcut = hfptMax; // on particle level D0
    const double MCDHfPtMaxcut = hfptMax; // on detector level D0
    
    // Create 2D (detector level, prompt D0's, matched to particle level) input distribution pT,jet vs DeltaR
    dataContainer.hInputDetector = new TH2D("hInputDetector", "Detector level prompt D^{0} jets distribution; pT,jet (GeV); #DeltaR", 
                        binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                        binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());
    dataContainer.hInputDetector->Sumw2();
    // Create 2D (matched particle level, prompt D0's, matched to the previous detector level distribution) input distribution pT,jet vs DeltaR
    dataContainer.hInputParticle = new TH2D("hInputParticle", "Particle level prompt D^{0} jets distribution; pT,jet (GeV); #DeltaR", 
                        binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                        binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hInputParticle->Sumw2();

    // Create kinematic efficiency histograms
    dataContainer.hKineEffParticle[0] = new TH2D("hKineEffParticleNumerator", "Particle level kinematic efficiency numerator (non-prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hKineEffParticle[1] = new TH2D("hKineEffParticleDenominator", "Particle level kinematic efficiency denominator (non-prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hKineEffDetector[0] = new TH2D("hKineEffDetectorNumerator", "Detector level kinematic efficiency numerator (non-prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());
    dataContainer.hKineEffDetector[1] = new TH2D("hKineEffDetectorDenominator", "Detector level kinematic efficiency denominator (non-prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());

    // Create response matrix for unfolding
    dataContainer.response = new RooUnfoldResponse(dataContainer.hKineEffDetector[0], dataContainer.hKineEffParticle[0]);
    dataContainer.hResponseProj[0] = new TH2D("hResponseProj_pTjet","Response matrix projection on p_{T,jet};p_{T,jet}^{reco ch} (GeV/#it{c});p_{T,jet}^{truth ch} (GeV/#it{c})",
                                             binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(),
                                             binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data());
    dataContainer.hResponseProj[1] = new TH2D("hResponseProj_deltaR","Response matrix projection on #DeltaR;#DeltaR^{reco};#DeltaR^{gen}",
                                             binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data(),
                                             binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hResponseProjDeltaR.resize(binning.ptjetBinEdges_detector.size() - 1);

    //_______________________________________________Correction data_______________________________________________
    //
    // Accessing TTree for correction data
    //
    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    int MCPjetNConst;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;
    int MCDjetNConst;
    // Accessing TTree
    TTree* tree = (TTree*)fClosureInput->Get("CorrectionTree");
    // Check for correct access
    if (!tree) {
        cout << "Error opening correction data tree.\n";
    }
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
    int MCDhfMatchedFrom, MCDhfSelectedAs;
    tree->SetBranchAddress("fHfMatchedFrom",&MCDhfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&MCDhfSelectedAs);

    //
    // Filling objects
    //
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // Apply prompt selection (i.e., only c → D0)
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        if (!MCDhfprompt || !isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs) || !MCDhfmatch) {
            continue;
        }

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= binning.deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= binning.deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        // Get the threshold for this pT range
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;

        // Fill response matrix and kinematic efficiency histograms
        if (genLevelRange && recoLevelRange && passBDTcut) {
            
            // Fill 4D RooUnfoldResponse object (jet pT shape is influenced by D0 pT efficiency)
            dataContainer.response->Fill(MCDjetPt, MCDDeltaR, MCPjetPt, MCPDeltaR);
            dataContainer.hResponseProj[0]->Fill(MCDjetPt, MCPjetPt);
            dataContainer.hResponseProj[1]->Fill(MCDDeltaR, MCPDeltaR);
            
            // Fill kinematic efficiency numerator histograms
            dataContainer.hKineEffParticle[0]->Fill(MCPjetPt, MCPDeltaR);
            dataContainer.hKineEffDetector[0]->Fill(MCDjetPt, MCDDeltaR);

            // if (pT in the iJetPtBin) {
            //     dataContainer.hResponseProjDeltaR[iJetPt]->Fill(MCDDeltaR, MCPDeltaR);
            // }
            
        }
        if (genLevelRange && passBDTcut) {
            // Fill kinematic efficiency denominator histogram for full particle level range
            dataContainer.hKineEffParticle[1]->Fill(MCPjetPt, MCPDeltaR);
            
        }
        if (recoLevelRange && passBDTcut) {
            // Fill kinematic efficiency denominator histogram for full detector level range
            dataContainer.hKineEffDetector[1]->Fill(MCDjetPt, MCDDeltaR);
        }
    }

    // Calculate kinematic efficiency
    dataContainer.hKineEffParticle[0]->Sumw2(); // numerator
    dataContainer.hKineEffParticle[1]->Sumw2(); // denominator
    dataContainer.hKineEffParticle[2] = (TH2D*)dataContainer.hKineEffParticle[0]->Clone("hKineEffParticleEfficiency");
    dataContainer.hKineEffParticle[2]->Divide(dataContainer.hKineEffParticle[1]); // A = A / B:  A = A->Divide(B)
    dataContainer.hKineEffDetector[0]->Sumw2(); // numerator
    dataContainer.hKineEffDetector[1]->Sumw2(); // denominator
    dataContainer.hKineEffDetector[2] = (TH2D*)dataContainer.hKineEffDetector[0]->Clone("hKineEffDetectorEfficiency");
    dataContainer.hKineEffDetector[2]->Divide(dataContainer.hKineEffDetector[1]); // A = A / B:  A = A->Divide(B)
    //_______________________________________________Input data_______________________________________________
    //
    // Accessing TTree for correction data
    //
    tree = (TTree*)fClosureInput->Get("InputTree");
    // Check for correct access
    if (!tree) {
        cout << "Error opening correction data tree.\n";
    }
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
    //int MCDhfMatchedFrom, MCDhfSelectedAs;
    tree->SetBranchAddress("fHfMatchedFrom",&MCDhfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&MCDhfSelectedAs);

    //
    // Filling objects
    //
    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // Apply prompt (excluding reflections) selection (i.e., only c → D0)
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        if (!MCDhfprompt || !isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs) || !MCDhfmatch) {
            continue;
        }

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= binning.deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= binning.deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        // Get the threshold for this pT range
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;

        if (genLevelRange && passBDTcut) {
            // Fill input distribution for particle level
            dataContainer.hInputParticle->Fill(MCPjetPt, MCPDeltaR);
            //std::cout << "Filling particle level input: pT,jet = " << MCPjetPt << ", DeltaR = " << MCPDeltaR << std::endl;
        }
        if (recoLevelRange && passBDTcut) {
            // Fill input distribution for detector level
            dataContainer.hInputDetector->Fill(MCDjetPt, MCDDeltaR);
            //std::cout << "Filling detector level input: pT,jet = " << MCDjetPt << ", DeltaR = " << MCDDeltaR << std::endl;
        }
    }

    dataContainer.hInputDetectorRaw = (TH2D*) dataContainer.hInputDetector->Clone("hInputDetectorRaw");
    return dataContainer;

}

void unfoldInputDetector(ClosureTestData1& dataContainer, int& iterationNumber) {
    // Correct input with detector level kinematic efficiency
    dataContainer.hInputDetector->Multiply(dataContainer.hKineEffDetector[2]); // A = A * B:  A = A->Multiply(B)
    //std::cout << "Input detector level histogram has " << dataContainer.hInputDetector->GetEntries() << " entries after kinematic efficiency correction." << std::endl;

    // Unfold with multiple iterations
    // Unfold in multiple iterations
    for (int iIter = 1; iIter <= iterationNumber; iIter++) {
        dataContainer.unfold.push_back(new RooUnfoldBayes(dataContainer.response, dataContainer.hInputDetector, iIter));
        //dataContainer.unfold[iIter - 1] = new RooUnfoldBayes(dataContainer.response, dataContainer.hInputDetector, iIter);
        dataContainer.hUnfolded.push_back((TH2D*) dataContainer.unfold[iIter - 1]->Hreco(RooUnfold::kCovariance)); // kCovariance = 2 = Errors from the square root of of the covariance matrix given by the unfolding
        dataContainer.hUnfolded[iIter - 1]->SetTitle(Form("Unfolded 2D histogram with %d iterations;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR^{gen}", iIter));
        dataContainer.hUnfolded[iIter - 1]->SetName(Form("hUnfolded_iter%d", iIter));
        dataContainer.hUnfolded[iIter - 1]->Sumw2();

        // Correct distribution with particle level kinematic efficiency (add entries)
        dataContainer.hUnfolded[iIter - 1]->Divide(dataContainer.hKineEffParticle[2]); // A = A / B:  A = A->Divide(B)
        
    }
}

void plotHistograms(const ClosureTestData1& dataContainer, const double& jetptMin, const double& jetptMax, const BinningStruct& binning) {
    cout << "Plotting histograms...\n";

    gStyle->SetPalette(kRainbow);

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    TH1D* hUnfoldedProjNotNorm = nullptr; // to store the last iteration unfolded distribution before normalization for plotting

    // --- Raw counts canvas with ratio pad ---
    TCanvas* cFirstClosureTest = new TCanvas("cFirstClosureTest", "First Closure Test: Unfolding", 1800, 1000);
    TPad* padTop    = new TPad("padTop",    "top pad",    0, 0.3, 1, 1.0);
    TPad* padBottom = new TPad("padBottom", "bottom pad", 0, 0.0, 1, 0.3);
    padTop->SetBottomMargin(0.0);
    padBottom->SetTopMargin(0.0);
    padBottom->SetBottomMargin(0.35);
    padTop->SetLeftMargin(0.12);
    padBottom->SetLeftMargin(0.12);
    padTop->SetTickx();    padTop->SetTicky();
    padBottom->SetTickx(); padBottom->SetTicky();
    padTop->Draw();
    padBottom->Draw();

    // --- Normalized canvas with ratio pad ---
    TCanvas* cFirstClosureTestNorm = new TCanvas("cFirstClosureTestNorm", "First Closure Test: Unfolding, self-normalized", 1800, 1000);
    TPad* padTopNorm    = new TPad("padTopNorm",    "top pad norm",    0, 0.3, 1, 1.0);
    TPad* padBottomNorm = new TPad("padBottomNorm", "bottom pad norm", 0, 0.0, 1, 0.3);
    padTopNorm->SetBottomMargin(0.0);
    padBottomNorm->SetTopMargin(0.0);
    padBottomNorm->SetBottomMargin(0.35);
    padTopNorm->SetLeftMargin(0.12);
    padBottomNorm->SetLeftMargin(0.12);
    padTopNorm->SetTickx();    padTopNorm->SetTicky();
    padBottomNorm->SetTickx(); padBottomNorm->SetTicky();
    padTopNorm->Draw();
    padBottomNorm->Draw();

    TLegend* lUnfoldedIter     = new TLegend(0.6, 0.57, 0.88, 0.87);
    TLegend* lUnfoldedIterNorm = new TLegend(0.6, 0.57, 0.88, 0.87);
    std::vector<TH1D*> hUnfoldedProj(dataContainer.hUnfolded.size());
    std::vector<TH1D*> hUnfoldedProjNorm(dataContainer.hUnfolded.size());

    int secondBin = 2;
    int lastButOneBin = dataContainer.hUnfolded[0]->GetXaxis()->GetNbins() - 1;

    for (size_t iIter = 0; iIter < dataContainer.hUnfolded.size(); iIter++) {

        hUnfoldedProj[iIter] = dataContainer.hUnfolded[iIter]->ProjectionY(Form("hProjIter_%zu", iIter), secondBin, lastButOneBin);
        hUnfoldedProj[iIter]->SetLineColor(kBlack + iIter);
        lUnfoldedIter->AddEntry(hUnfoldedProj[iIter], Form("Iteration %zu", iIter+1), "le");

        hUnfoldedProjNorm[iIter] = (TH1D*) hUnfoldedProj[iIter]->Clone(
                                        Form("hProjIterNorm_%zu", iIter));
        lUnfoldedIterNorm->AddEntry(hUnfoldedProjNorm[iIter], Form("Iteration %zu", iIter+1), "le");

        double integral = hUnfoldedProjNorm[iIter]->Integral();
        if (integral != 0) hUnfoldedProjNorm[iIter]->Scale(1.0 / integral, "width");

        // suppress X axis labels on top pads
        hUnfoldedProj[iIter]->GetXaxis()->SetLabelSize(0);
        hUnfoldedProj[iIter]->GetXaxis()->SetTitleSize(0);
        hUnfoldedProjNorm[iIter]->GetXaxis()->SetLabelSize(0);
        hUnfoldedProjNorm[iIter]->GetXaxis()->SetTitleSize(0);

        // save last iteration for use elsewhere if needed
        if (iIter == dataContainer.hUnfolded.size() - 1) {
            hUnfoldedProjNotNorm = (TH1D*) hUnfoldedProj[iIter]->Clone("hUnfoldedProjNotNorm");
        }

        padTop->cd();
        if (iIter == 0) {
            hUnfoldedProj[iIter]->SetTitle("Closure test 1: unfolding");
            hUnfoldedProj[iIter]->GetYaxis()->SetTitle("dN");
            hUnfoldedProj[iIter]->Draw();
        } else {
            hUnfoldedProj[iIter]->Draw("same");
        }

        padTopNorm->cd();
        if (iIter == 0) {
            hUnfoldedProjNorm[iIter]->SetTitle("Closure test 1: unfolding, all histograms self-normalized");
            hUnfoldedProjNorm[iIter]->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
            hUnfoldedProjNorm[iIter]->Draw();
        } else {
            hUnfoldedProjNorm[iIter]->Draw("same");
        }
    }

    // --- Truth reference ---
    TH1D* hInputParticleProj = dataContainer.hInputParticle->ProjectionY(
                                "hInputParticleProj", secondBin, lastButOneBin);
    hInputParticleProj->SetLineColor(kRed);
    hInputParticleProj->SetLineStyle(2);
    hInputParticleProj->SetLineWidth(2);
    hInputParticleProj->GetXaxis()->SetLabelSize(0);
    hInputParticleProj->GetXaxis()->SetTitleSize(0);

    TH1D* hInputParticleProjNorm = (TH1D*) hInputParticleProj->Clone("hInputParticleProjNorm");
    double integralTruth = hInputParticleProjNorm->Integral();
    if (integralTruth != 0) hInputParticleProjNorm->Scale(1.0 / integralTruth, "width");
    hInputParticleProjNorm->GetXaxis()->SetLabelSize(0);
    hInputParticleProjNorm->GetXaxis()->SetTitleSize(0);

    double reportingJetPtMin = dataContainer.hUnfolded[0]->GetXaxis()->GetBinLowEdge(secondBin);
    double reportingJetPtMax = dataContainer.hUnfolded[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);

    padTop->cd();
    hInputParticleProj->Draw("same");
    lUnfoldedIter->AddEntry(hInputParticleProj, "Input particle level", "le");
    lUnfoldedIter->Draw();
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c",
                                    reportingJetPtMin, reportingJetPtMax));

    padTopNorm->cd();
    hInputParticleProjNorm->Draw("same");
    lUnfoldedIterNorm->AddEntry(hInputParticleProjNorm, "Input particle level", "le");
    lUnfoldedIterNorm->Draw();
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));

    // --- Ratio pads ---
    auto drawRatioPad = [&](TPad* pad, const std::vector<TH1D*>& hVec, TH1D* hRef) {
        pad->cd();
        for (size_t iIter = 0; iIter < hVec.size(); iIter++) {
            TH1D* hRatio = (TH1D*) hVec[iIter]->Clone(Form("hRatioCT1_iter%zu_%s", iIter+1, pad->GetName()));
            hRatio->Divide(hRef);
            hRatio->SetLineColor(kBlack + iIter);
            hRatio->SetTitle("");
            hRatio->GetYaxis()->SetTitle("#frac{Unfolded_{i}}{Truth}");
            hRatio->GetYaxis()->SetTitleSize(0.10);
            hRatio->GetYaxis()->SetLabelSize(0.09);
            hRatio->GetXaxis()->SetTitleSize(0.12);
            hRatio->GetXaxis()->SetLabelSize(0.10);
            hRatio->GetYaxis()->SetTitleOffset(0.45);
            hRatio->GetXaxis()->SetTitle("#DeltaR^{gen}");
            if (iIter == 0) {
                if (binning.useEmmaYeatsBins) {
                    hRatio->GetYaxis()->SetRangeUser(0.9, 1.1);
                } else {
                    hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
                }
                //hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
                hRatio->Draw("E");
            } else {
                hRatio->Draw("E same");
            }
        }
        TLine* unity = new TLine(hRef->GetXaxis()->GetXmin(), 1.0,hRef->GetXaxis()->GetXmax(), 1.0);
        unity->SetLineStyle(2);
        unity->SetLineColor(kBlack);
        unity->Draw("same");
    };

    drawRatioPad(padBottom,     hUnfoldedProj,     hInputParticleProj);
    drawRatioPad(padBottomNorm, hUnfoldedProjNorm, hInputParticleProjNorm);

    TCanvas* cRatiosTruthToUnfolded = new TCanvas("cRatiosTruthToUnfolded", "Ratios of truth and unfolded distributions", 1800, 1000);
    TLegend* lRatios = new TLegend(0.2, 0.2, 0.3, 0.4);
    for (size_t iIter = 0; iIter < hUnfoldedProj.size(); iIter++) {
        // TH1D* hRatio = (TH1D*) hUnfoldedProj[iIter]->Clone(Form("hRatioTruthToUnfolded_iter%zu", iIter+1));
        TH1D* hRatio = (TH1D*) hInputParticleProj->Clone(Form("hRatioTruthToUnfolded_iter%zu", iIter+1));
        hRatio->Divide(hUnfoldedProj[iIter]);
        hRatio->SetLineStyle(kSolid);
        hRatio->SetLineWidth(1);
        hRatio->SetTitle(Form("Ratio of input particle level to unfolded distribution, for each iteration;#DeltaR;Input particle level / Unfolded"));
        hRatio->GetYaxis()->SetTitle("#frac{Input reference}{Fully corrected}");
        hRatio->SetLineColor(kBlack + iIter);
        lRatios->AddEntry(hRatio, Form("%zu iterations", iIter+1), "le");
        cRatiosTruthToUnfolded->cd();
        hRatio->Draw(iIter == 0 ? "" : "same");
    }
    TLine* lUnity = new TLine(hInputParticleProj->GetXaxis()->GetXmin(), 1.0,hInputParticleProj->GetXaxis()->GetXmax(), 1.0);
    lUnity->SetLineStyle(2);
    lUnity->SetLineColor(kBlack);
    lUnity->Draw("same");
    lRatios->Draw();

    
    TCanvas* cFullRange = new TCanvas("cFullRange", "Closure test 1: unfolding, projection in full pT,jet range");
    cFullRange->Divide(2,2);
    cFullRange->cd(1);
    int binMin = 1;
    int binMax = dataContainer.hUnfolded[0]->GetXaxis()->GetNbins();
    double fullReportingJetPtMin = dataContainer.hUnfolded[0]->GetXaxis()->GetBinLowEdge(binMin);
    double fullReportingJetPtMax = dataContainer.hUnfolded[0]->GetXaxis()->GetBinUpEdge(binMax);
    TH1D* hUnfoldedProjFullRange = dataContainer.hUnfolded.back()->ProjectionY("hUnfoldedProjFullRange", binMin, binMax);
    hUnfoldedProjFullRange->Draw();
    latex->DrawLatex(0.15, 0.85,Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c",fullReportingJetPtMin, fullReportingJetPtMax));
    cFullRange->cd(2);
    dataContainer.hUnfolded.back()->Draw("colz");
    latex->DrawLatex(0.15, 0.85,Form("Unfolded and kinematically corrected distribution, p_{T,jet} #in [%.1f, %.1f] GeV/c",jetptMin, jetptMax));
    cFullRange->cd(3);
    hInputParticleProj->SetTitle("Not normalized input particle level vs fully corrected unfolded distribution;#DeltaR;dN");
    hInputParticleProj->Draw();
    hUnfoldedProjNotNorm->Draw("same");
    latex->DrawLatex(0.15, 0.85,Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c",reportingJetPtMin, reportingJetPtMax));
    cFullRange->cd(4);
    TH1D* hInputParticleProjBinwidth = (TH1D*) hInputParticleProj->Clone("hInputParticleProjBinwidth");
    hInputParticleProjBinwidth->Scale(1.0, "width");
    hInputParticleProjBinwidth->SetTitle("Input particle level normalized to bin width;#DeltaR;dN/d#DeltaR");
    TH1D* hUnfoldedProjNotNormBinwidth = (TH1D*) hUnfoldedProjNotNorm->Clone("hUnfoldedProjNotNormBinwidth");
    hUnfoldedProjNotNormBinwidth->Scale(1.0, "width");
    hUnfoldedProjNotNormBinwidth->SetTitle("Fully corrected unfolded distribution normalized to bin width;#DeltaR;dN/d#DeltaR");
    hInputParticleProjBinwidth->Draw();
    hUnfoldedProjNotNormBinwidth->Draw("same");
    latex->DrawLatex(0.15, 0.85,Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c",reportingJetPtMin, reportingJetPtMax));
    

    TCanvas* cCorrectionObjects = new TCanvas("cCorrectionObjects", "Correction objects");
    cCorrectionObjects->Divide(2,2);
    dataContainer.hKineEffParticle[2]->SetTitle("Particle level kinematic efficiency;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR^{gen}");
    cCorrectionObjects->cd(1);
    dataContainer.hKineEffParticle[2]->Draw("text");
    dataContainer.hKineEffDetector[2]->SetTitle("Detector level kinematic efficiency;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}");
    cCorrectionObjects->cd(2);
    dataContainer.hKineEffDetector[2]->Draw("text");
    cCorrectionObjects->cd(3);
    dataContainer.hResponseProj[0]->Draw("colz");
    cCorrectionObjects->cd(4);
    dataContainer.hResponseProj[1]->Draw("colz");

    TCanvas* cResponseObjects = new TCanvas("cResponseObjects","Response matrices and projections");
    cResponseObjects->Divide(2,2);
    cResponseObjects->cd(1);
    dataContainer.hResponseProj[0]->Draw("colz");
    cResponseObjects->cd(2);
    dataContainer.hResponseProj[1]->Draw("colz");
    cResponseObjects->cd(3);
    TH1D* hResponseProjDeltaRReco = dataContainer.hResponseProj[1]->ProjectionX("hResponseProjDeltaRReco",2,dataContainer.hResponseProj[1]->GetYaxis()->GetNbins() - 1);
    TH1D* hResponseProjDeltaRGen = dataContainer.hResponseProj[1]->ProjectionY("hResponseProjDeltaRGen",2,dataContainer.hResponseProj[1]->GetXaxis()->GetNbins() - 1);
    hResponseProjDeltaRReco->SetLineColor(kGreen+1);
    hResponseProjDeltaRGen->SetLineColor(kRed+1);
    hResponseProjDeltaRReco->Draw();
    hResponseProjDeltaRGen->Draw("same");
    TLegend* legResponseProjDeltaR = new TLegend(0.6,0.57,0.7,0.77);
    legResponseProjDeltaR->AddEntry(hResponseProjDeltaRReco, "Reco", "lp");
    legResponseProjDeltaR->AddEntry(hResponseProjDeltaRGen,"Gen", "lp");
    legResponseProjDeltaR->Draw();

    TCanvas* cInputs = new TCanvas("cInputs","Input particle and detector level distributions with and w/o normalization");
    cInputs->Divide(2,2);
    TH1D* hInputDetectorRange = dataContainer.hInputDetectorRaw->ProjectionY("hInputDetectorRange",2,dataContainer.hInputDetectorRaw->GetXaxis()->GetNbins() - 1);
    hInputDetectorRange->GetYaxis()->SetTitle("dN");
    hInputDetectorRange->SetTitle("Input detector level");
    TH1D* hInputDetectorRangeNorm = (TH1D*) hInputDetectorRange->Clone("hInputDetectorRangeNorm");
    hInputDetectorRangeNorm->SetTitle("Input detector level normalized");
    hInputDetectorRangeNorm->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
    hInputDetectorRangeNorm->Scale(1 / hInputDetectorRangeNorm->Integral(), "width");
    TH1D* hInputParticleRange = dataContainer.hInputParticle->ProjectionY("hInputParticleRange",2,dataContainer.hInputParticle->GetXaxis()->GetNbins() - 1);
    hInputParticleRange->GetYaxis()->SetTitle("dN");
    hInputParticleRange->SetTitle("Input particle level flag");
    TH1D* hInputParticleRangeNorm = (TH1D*) hInputParticleRange->Clone("hInputParticleRangeNorm");
    hInputParticleRangeNorm->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
    hInputParticleRangeNorm->SetTitle("Input particle level normalized");
    hInputParticleRangeNorm->Scale(1 / hInputParticleRangeNorm->Integral(), "width");
    cInputs->cd(1);
    hInputDetectorRange->Draw();
    cInputs->cd(2);
    hInputDetectorRangeNorm->Draw();
    cInputs->cd(3);
    hInputParticleRange->Draw();
    cInputs->cd(4);
    hInputParticleRangeNorm->Draw();


    //
    // Storing images
    //
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../Images/5-ClosureTest/First/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    cFirstClosureTest->Update();
    cFirstClosureTest->SaveAs(imagePath + "ClosureTest1_unfolding.png");    
    
    //cKinEff->Print(imagePath + Form("unfolding" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cCorrectionObjects->Print(imagePath + Form("closureTest1_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cFirstClosureTest->Print(imagePath + Form("closureTest1_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFirstClosureTestNorm->Print(imagePath + Form("closureTest1_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFullRange->Print(imagePath + Form("closureTest1_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cInputs->Print(imagePath + Form("closureTest1_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cRatiosTruthToUnfolded->Print(imagePath + Form("closureTest1_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
    //cUnfoldedIter->Print(imagePath + Form("unfolding" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

}

void saveData(const ClosureTestData1& dataContainer, const double& jetptMin, const double& jetptMax, const BinningStruct& binning) {
    // Open output file
    TFile* outFile = new TFile(Form("FirstOutput/closure_test_1_results_%d_to_%d_jetpt_" + binning.dataPeriod + ".root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");

    // store each histogram in file
    //dataContainer.hSBUnfolded->Write();
    dataContainer.hInputParticle->Write();
    dataContainer.hInputDetector->Write(); // with kinematic efficiency correction
    dataContainer.hResponseProj[0]->Write();
    dataContainer.hResponseProj[1]->Write();
    
    // Return to root directory (optional)
    outFile->cd();

    // Also store the axes used for the histograms
    storeBinningInFile(outFile, binning);

    outFile->Close();
    delete outFile;
    
    std::cout << "Data stored in file " << Form("closure_test_1_results_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)) << endl;
}

// 1 - Unfolding closure test
void FirstClosureTest(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // Number of unfolding procedure iterations
    int iterationNumber = 8;

    // Load binning from reflections file
    TFile* fBinning = new TFile(Form("../1-SignalTreatment/Reflections/binningInfo.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    double minDeltaR = binning.deltaRBinEdges_detector[0];
    double maxDeltaR = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1];
    double hfptMin = binning.ptHFBinEdges_detector[0]; //ptHFBinEdges[0] - should start from 0 or from the lowest pT,D value?
    double hfptMax = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1];

    TFile* fEfficiency = new TFile(Form("../2-Efficiency/selection_efficiency_run3style_%d_to_%d_jetpt_" + binning.dataPeriod + ".root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fEfficiency || fEfficiency->IsZombie()) {
        std::cerr << "Error: Unable to open estimated selection efficiency ROOT file." << std::endl;
    }

    // Opening files
    TFile* fClosureInput = new TFile("mc_closure_input_data_" + binning.dataPeriod + ".root","read");
    if (!fClosureInput || fClosureInput->IsZombie()) {
        std::cerr << "Error: Unable to open 1st closure input ROOT file." << std::endl;
    }

    // Create response matrix with correction sample (without efficiency scaling)
    ClosureTestData1 dataContainer = createAnalysisObjects(fClosureInput, jetptMin, jetptMax, hfptMin, hfptMax, binning);
    
    // Unfold the detector level distribution (with particle and detector level kinematic efficiency corrections)
    unfoldInputDetector(dataContainer, iterationNumber);

    // Compare the unfolded distribution with the particle level distribution

    // Plot the efficiency histogram and further corrected histograms
    plotHistograms(dataContainer, jetptMin, jetptMax, binning);

    // Save corrected distributions to file
    saveData(dataContainer, jetptMin, jetptMax, binning);

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
    FirstClosureTest();
    return 0;
}
