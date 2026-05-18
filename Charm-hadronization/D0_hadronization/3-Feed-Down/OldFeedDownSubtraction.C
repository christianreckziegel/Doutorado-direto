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
#include "../commonUtilities.h"

using namespace std;

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
    TH2D* hResponseDeltaR;
    TH2D* hResponsePtJet;
    TH2D* hProjected2DKinCorrected = nullptr; // 3D histogram projection with particle level kinematic efficiency correction
    TH2D* hMeasuredTemplate = nullptr; // template for measured data (jet pT vs DeltaR), used for hFolded binning
    TH2D* hFolded2D = nullptr; // folded 2D histogram (jet pT vs DeltaR)
    TH2D* hFolded2DKinCorrected = nullptr; // folded 2D histogram with detector level kinematic correction

    // Step 5: final scaled histogram (jet pT vs DeltaR)
    TH2D* hFinalScaled = nullptr;

    // Metadata
    double lumiData = 1; // luminosity of the data in 1/pb
    double lumiMC = 1; // luminosity of the MC in 1/pb
    const double branchingRatio = 0.03947; // D0 -> KPi decay channel branching ratio = (3.947 +- 0.04) %

    // Step 6: B subtracted data distribution
    TH2D* hBFedDownData;
};

FeedDownData createHistograms(const BinningStruct& binning) {
    // Create struct to store data
    FeedDownData dataContainer;

    //
    // Create 3D POWHEG data histogram
    //
    dataContainer.hPowheg = new TH3D("hPowheg", "POWHEG + PYTHIA;p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR;p_{T,D}^{gen}", binning.ptjetBinEdges_particle.size()-1, binning.ptjetBinEdges_particle.data(), 
                                                                                                                          binning.deltaRBinEdges_particle.size()-1, binning.deltaRBinEdges_particle.data(), 
                                                                                                                          binning.ptHFBinEdges_particle.size()-1, binning.ptHFBinEdges_particle.data());
    dataContainer.hPowheg->SetMarkerColor(30);
    dataContainer.hPowheg->SetLineColor(30); // 30 = pastel green
    dataContainer.hPowheg->SetMarkerStyle(kCircle);
    dataContainer.hPowheg->Sumw2();
    dataContainer.hPowheg->SetStats(0);

    
    //
    // Matching histograms for folding process: kinematic efficiencies numerator and denominator for particle level and detector level
    //
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
    
    dataContainer.hMeasuredTemplate = new TH2D("hMeasuredTemplate", "Folded 2D histogram;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco}", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());
    //
    // Response matrix projections: DeltaR and pT,jet
    //
    dataContainer.hResponseDeltaR = new TH2D("hResponseDeltaR","Response matrix #DeltaR projection;#DeltaR^{reco};#DeltaR^{gen}",
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data(),
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hResponsePtJet = new TH2D("hResponsePtJet","Response matrix p_{T,jet} projection;p_{T,jet}^{reco ch};p_{T,jet}^{gen}",
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(),
                                                                                  binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data());
    std::cout << "Histograms created.\n";

    // Create response matrix for folding
    std::cout << "Creating response matrix" << std::endl;
    if (dataContainer.hKineEffParticle[0] == nullptr || dataContainer.hKineEffDetector[0] == nullptr) {
        std::cerr << "Error: Kinematic efficiency histograms are not initialized.\n";
        return dataContainer;
    }
    dataContainer.response = RooUnfoldResponse(dataContainer.hKineEffDetector[0], dataContainer.hKineEffParticle[0]);
    std::cout << "Response matrix created.\n";

    return dataContainer;
};

// Fill histograms with data from files: 
void fillHistograms(TFile* fEfficiency, TFile* fPowhegNonPrompt, TFile* fSimulatedMCMatched, TFile* fLumi, FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax, const BinningStruct& binning) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level HF
    const double MCDyCut = 0.8; // on detector level HF
    const double MCPDeltaRcut = binning.deltaRBinEdges_particle[binning.deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double MCPHfPtMincut = binning.ptHFBinEdges_particle[0]; // on particle level HF
    const double MCDHfPtMincut = binning.ptHFBinEdges_detector[0]; // on detector level HF
    const double MCPHfPtMaxcut = binning.ptHFBinEdges_particle[binning.ptHFBinEdges_particle.size() - 1]; // on particle level HF
    const double MCDHfPtMaxcut = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1]; // on detector level HF

    // Defining variables for accessing POWHEG tree
    double PowaxisDistance, PowjetPt, PowjetEta, PowjetPhi;
    double PowhfPt, PowhfEta, PowhfPhi, PowhfY;
    // Defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    int MCDhfMatchedFrom, MCDhfSelectedAs;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;
    int MCPjetNConst;
    int MCDjetNConst;

    // (i) Fetch efficiency distributions from previous step
    dataContainer.hSelectionEfficiency.first = (TH1D*)fEfficiency->Get("hSelectionEfficiencyPrompt")->Clone("hSelectionEfficiencyPrompt_clone");
    if (!dataContainer.hSelectionEfficiency.first) {
        std::cout << "Error fetching prompt D0 selection efficiency histogram from file.\n";
    }
    dataContainer.hSelectionEfficiency.second = (TH1D*)fEfficiency->Get("hSelectionEfficiencyNonPrompt")->Clone("hSelectionEfficiencyNonPrompt_clone");
    if (!dataContainer.hSelectionEfficiency.second) {
        std::cout << "Error fetching non-prompt D0 selection efficiency histogram from file.\n";
    }

    // (ii) Access POWHEG+PYTHIA non-prompt D0 data and fill 3D histogram
    TTree* tree = (TTree*)fPowhegNonPrompt->Get("tree_D0");
    if (!tree) {
        cout << "Error opening POWHEG tree.\n";
    }
    tree->SetBranchAddress("pt_cand",&PowhfPt);
    tree->SetBranchAddress("eta_cand",&PowhfEta);
    tree->SetBranchAddress("phi_cand",&PowhfPhi);
    tree->SetBranchAddress("y_cand",&PowhfY);
    tree->SetBranchAddress("pt_jet",&PowjetPt);
    tree->SetBranchAddress("eta_jet",&PowjetEta);
    tree->SetBranchAddress("phi_jet",&PowjetPhi);
    tree->SetBranchAddress("delta_r_jet",&PowaxisDistance);
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double deltaR = PowaxisDistance;
        bool genLevelRange = (abs(PowjetEta) < MCPetaCut) && (abs(PowhfY) < MCPyCut) && ((PowjetPt >= jetptMin) && (PowjetPt < jetptMax)) && ((deltaR >= binning.deltaRBinEdges_particle[0]) && (deltaR < MCPDeltaRcut)) && ((PowhfPt >= MCPHfPtMincut) && (PowhfPt < MCPHfPtMaxcut));
        
        // Fill 2D histogram considering jet pT and detector acceptance
        if (genLevelRange) {
            dataContainer.hPowheg->Fill(PowjetPt, PowaxisDistance, PowhfPt);
            // use Emma's run 2 cuts in case useEmmaYeatsBins is true
            if (!binning.useEmmaYeatsBins || passEmmaCut(PowjetPt, PowhfPt)) {
                dataContainer.hPowheg->Fill(PowjetPt, PowaxisDistance, PowhfPt);
            }
        }
        
    }
    std::cout << "Generator level (POWHEG+PYTHIA) histograms filled.\n";
    
    // (iii) Access simulated MC matched data and fill response matrix histograms and numerator and denominator histograms for kinematic efficiency calculation
    tree = (TTree*)fSimulatedMCMatched->Get("DF_merged/O2matchtable");
    if (!tree) {
        cout << "Error opening O2 matching tree.\n";
    }
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance); // particle level branches
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetNConst);
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    MCPhfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&MCPhfY);
    tree->SetBranchAddress("fMcHfPrompt",&MCPhfprompt);
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance); // detector level branches
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fJetNConst",&MCDjetNConst);
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
    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // Apply non-prompt selection (i.e., only B → D0)
        //bool isReflection = (MCDhfMatchedFrom != MCDhfSelectedAs) ? true : false;
        if (MCDhfprompt || !isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs) || (MCDjetNConst < 0)) {
            continue;
        }

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));
        // Define gen-level and reco-level selection criteria
        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= binning.deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= binning.deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        // Get the threshold for this pT range
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;

        // Fill response matrix and kinematic efficiency histograms
        if (genLevelRange && recoLevelRange && passBDTcut) {
            
            // Find the bin corresponding to the given pT,D value
            int bin = dataContainer.hSelectionEfficiency.first->FindBin(MCDhfPt);
            // Get the efficiency value from the bin content
            double efficiency_prompt = dataContainer.hSelectionEfficiency.first->GetBinContent(bin);
            // Fill 4D RooUnfoldResponse object (jet pT shape is influenced by D0 pT efficiency)
            dataContainer.response.Fill(MCDjetPt, MCDDeltaR, MCPjetPt, MCPDeltaR, 1./efficiency_prompt);
            dataContainer.hResponseDeltaR->Fill(MCDDeltaR,MCPDeltaR,1./efficiency_prompt);
            dataContainer.hResponsePtJet->Fill(MCDjetPt,MCPjetPt,1./efficiency_prompt);
            
            // Fill kinematic efficiency numerator histograms
            dataContainer.hKineEffParticle[0]->Fill(MCPjetPt, MCPDeltaR);
            dataContainer.hKineEffDetector[0]->Fill(MCDjetPt, MCDDeltaR);
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
    std::cout << "Response matrix and kinematic efficiencies histograms filled.\n";

}

std::vector<double> retrieveLuminosityFromFile(TFile* fPowhegNonPrompt, TFile* fLumi, FeedDownData& dataContainer) {
    
    // Define temporary histogram to access luminosity values from file
    TH1D* hLumi = nullptr;
    //
    // ----- POWHEG luminosity calculation
    //
    hLumi = (TH1D*) fPowhegNonPrompt->Get("fHistXsection");
    double crossSecPowheg = hLumi->GetBinContent(1);
    TTree* tree_D0 = dynamic_cast<TTree*>(fPowhegNonPrompt->Get("tree_D0")); // number of events in POWHEG file
    double numOfEventsPowheg = tree_D0->GetEntries();
    double powhegLuminosity = numOfEventsPowheg / crossSecPowheg;
    dataContainer.lumiMC = powhegLuminosity; // Store in dataContainer for later use
    //
    // ----- Data luminosity calculation
    //
    // Bunch crossings counting
    TTree* tLumiBunchCrossing = (TTree*) fLumi->Get("DF_merged/O2bccount"); // branch -> fCountsWithTVX
    Int_t bcTVX = 0, bcTVXSum = 0;
    tLumiBunchCrossing->SetBranchAddress("fCountsWithTVX", &bcTVX);
    std::cout << "Entries in BC+TVX tree: " << tLumiBunchCrossing->GetEntries() << std::endl;
    for (size_t iEntry = 0; iEntry < tLumiBunchCrossing->GetEntries(); iEntry++) {
        tLumiBunchCrossing->GetEntry(iEntry);
        bcTVXSum += bcTVX; // Sum over all entries to get total number of TVX triggered BC
    }
    // Collisions counting
    TTree* tLumiCollision = (TTree*) fLumi->Get("DF_merged/O2collcount");
    Int_t selection, collTVX = 0;
    Int_t selectionSum = 0, collTVXSum = 0;
    tLumiCollision->SetBranchAddress("fCountsWithTVXAndZVertexAndSel8", &selection);
    tLumiCollision->SetBranchAddress("fCountsWithTVX", &collTVX);
    for (size_t iEntry = 0; iEntry < tLumiCollision->GetEntries(); iEntry++) {
        tLumiCollision->GetEntry(iEntry);
        selectionSum += selection; // Sum over all entries to get total number of selected collisions
        collTVXSum += collTVX; // Sum over all entries to get total number of TVX triggered BC
    }

    // Number of TVX triggered BC that correspond to your selections and your train efficiencies
    double triggered = bcTVXSum * (selectionSum / collTVXSum);
    double runLuminosity = 1.0 / 0.0595e6; // luminosity value for the runs (// in mb⁻¹?)
    double dataLuminosity = triggered * runLuminosity;
    dataContainer.lumiData = dataLuminosity; // Store in dataContainer for later use

    std::cout << "Luminosities calculated." << std::endl;
    return {powhegLuminosity, dataLuminosity};
}

// Particle level kinematic efficiency
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

// Detector level kinematic efficiency
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

// Module for folding particle level data from POWHEG simulation
void smearGeneratorData(FeedDownData& dataContainer, TFile* fEfficiency) {
    std::cout << "Entries in prompt efficiency histogram: " << dataContainer.hSelectionEfficiency.first->GetEntries() << std::endl;
    // dataContainer.hSelectionEfficiency.first = (TH1D*)fEfficiency->Get("hSelectionEfficiencyPrompt");
    // dataContainer.hSelectionEfficiency.second = (TH1D*)fEfficiency->Get("hSelectionEfficiencyNonPrompt");

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
    dataContainer.hProjected2D = (TH2D*)dataContainer.hPowhegEffCorrected->Project3D("yx"); // xy = pT,jet vs Delta R
    dataContainer.hProjected2D->SetName("hProjected2D");
    dataContainer.hProjected2D->SetTitle("2nd step: projected POWHEG efficiency corrected histogram;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR");

    //
    // 3.1th step: remove outside of response range data in POWHEG (apply particle level kinematic efficiency)
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

void plotHistograms(const FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax, const BinningStruct& binning) {
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
    gPad->SetLogz();
    TCanvas* cResponseDeltaR = new TCanvas("cResponseDeltaR","cResponseDeltaR",1800,1000);
    cResponseDeltaR->cd();
    dataContainer.hResponseDeltaR->Draw("colz");
    gPad->SetLogz();
    TCanvas* cResponsePtJet = new TCanvas("cResponsePtJet","cResponsePtJet",1800,1000);
    cResponsePtJet->cd();
    dataContainer.hResponsePtJet->Draw("colz");
    gPad->SetLogz();
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
    TCanvas* cKinEffParticle = new TCanvas("cKinEffParticle","cKinEffParticle");
    cKinEffParticle->cd();
    dataContainer.hKineEffParticle[2]->SetTitle("Particle level kinematic efficiency (non-prompt D^{0}'s)");
    dataContainer.hKineEffParticle[2]->Draw("text");
    TCanvas* cKinEffDetector = new TCanvas("cKinEffDetector","cKinEffDetector");
    cKinEffDetector->cd();
    dataContainer.hKineEffDetector[2]->SetTitle("Detector level kinematic efficiency (non-prompt D^{0}'s)");
    dataContainer.hKineEffDetector[2]->Draw("text");
    
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
    gPad->SetLogz();
    cFedDownData->cd(2);
    TH1D* hBFedDownDataDeltaR =dataContainer.hBFedDownData->ProjectionY("hBFedDownDataDeltaR");
    hBFedDownDataDeltaR->Draw();
    cFedDownData->cd(3);
    double reportedJetptMin = dataContainer.hBFedDownData->GetXaxis()->GetBinLowEdge(2);
    double reportedJetptMax = dataContainer.hBFedDownData->GetXaxis()->GetBinUpEdge(dataContainer.hBFedDownData->GetXaxis()->GetLast()-1);
    TH1D* hBFedDownData1DRange = dataContainer.hBFedDownData->ProjectionY("hBFedDownData1DRange", 2, dataContainer.hBFedDownData->GetXaxis()->GetLast()-1);
    hBFedDownData1DRange->GetYaxis()->SetTitle("dN");
    hBFedDownData1DRange->Draw();
    latex->DrawLatex(0.2,0.2,Form("p_{T,jet} = [%.1f,%.1f] GeV/c",reportedJetptMin,reportedJetptMax));
    cFedDownData->cd(4);
    TH1D* hBFedDownData1DRangeNorm = (TH1D*) hBFedDownData1DRange->Clone("hBFedDownData1DRangeNorm");
    hBFedDownData1DRangeNorm->Scale(1 / hBFedDownData1DRangeNorm->Integral(), "width");
    hBFedDownData1DRangeNorm->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
    hBFedDownData1DRangeNorm->Draw();
    latex->DrawLatex(0.2,0.2,Form("p_{T,jet} = [%.1f,%.1f] GeV/c",reportedJetptMin,reportedJetptMax));
    TCanvas* cFedDownData2DOnly = new TCanvas("cFedDownData2DOnly","Background subtracted, efficiency corrected, fed down 2D data distribution");
    cFedDownData2DOnly->cd();
    dataContainer.hBFedDownData->Draw("colz");
    TCanvas* cFedDownData1DOnly = new TCanvas("cFedDownData1DOnly","Background subtracted, efficiency corrected, fed down 1D data distribution");
    cFedDownData1DOnly->cd();
    //dataContainer.hBFedDownData->ProjectionY("hBFedDownDataDeltaR")->Draw();
    hBFedDownDataDeltaR->Draw();

    //
    // Storing images
    //
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../Images/3-Feed-Down/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    
    cPowheg->Update();
    cPowheg->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cPowheg->SaveAs(imagePath + "powheg_3d.png");
    cResponse->Update();
    cResponse->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cResponse->SaveAs(imagePath + "response_matrix.png");
    cKinEfficiency->Update();
    cKinEfficiency->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cKinEfficiency->SaveAs(imagePath + "kinematic_efficiencies.png");
    cKinEffParticle->Update();
    cKinEffParticle->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cKinEffParticle->SaveAs(imagePath + "kinematic_efficiency_particle.png");
    cKinEffDetector->Update();
    cKinEffDetector->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cKinEffDetector->SaveAs(imagePath + "kinematic_efficiency_detector.png");
    cFolded->Update();
    cFolded->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cFolded->SaveAs(imagePath + "folded_stages.png");
    cFedDownData->Update();
    cFedDownData->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cFedDownData->SaveAs(imagePath + "fed_down_data.png");
    cFedDownData2DOnly->Update();
    cFedDownData2DOnly->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cFedDownData2DOnly->SaveAs(imagePath + "fed_down_2d_data.png");
    cFedDownData1DOnly->Update();
    cFedDownData1DOnly->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cFedDownData1DOnly->SaveAs(imagePath + "fed_down_1d_data.png");
    cResponseDeltaR->Update();
    cResponseDeltaR->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cResponseDeltaR->SaveAs(imagePath + "response_matrix_proj_deltar.png");
    cResponsePtJet->Update();
    cResponsePtJet->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cResponsePtJet->SaveAs(imagePath + "response_matrix_proj_ptjet.png");

    //
    // Storing in a single pdf file
    //
    cPowheg->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cResponse->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cKinEfficiency->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cKinEffParticle->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cKinEffDetector->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFolded->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFedDownData->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponseDeltaR->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponsePtJet->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFedDownData2DOnly->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
    // in separate .pdf file
    cKinEffParticle->Print(imagePath + "kinematic_efficiency_particle_" + sEmmaBins + ".pdf");
    cKinEffDetector->Print(imagePath + "kinematic_efficiency_detector_" + sEmmaBins + ".pdf");

}

void saveData(const FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax, 
              const BinningStruct& binning) {
    // Open output file
    TFile* outFile = new TFile(Form("outputFeedDown_%d_to_%d_jetpt_" + binning.dataPeriod + ".root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");

    // Store B -> D0 contribution histogram
    dataContainer.hFinalScaled->Write();

    // Store fed down data histogram
    dataContainer.hBFedDownData->Write();
    
    // Also store the axes used for the histograms
    storeBinningInFile(outFile, binning);

    // Return to root directory (optional)
    outFile->cd();

    outFile->Close();
    delete outFile;
    
    cout << "Data stored in file" << Form("outputFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)) << endl;
}

void FeedDownSubtraction(){

    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Load binning from reflections file
    TFile* fBinning = new TFile(Form("../1-SignalTreatment/Reflections/binningInfo.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

    // Open files
    TFile* fEfficiency = new TFile(Form("../2-Efficiency/selection_efficiency_run3style_%d_to_%d_jetpt_" + binning.dataPeriod + ".root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fEfficiency || fEfficiency->IsZombie()) {
        std::cerr << "Error: Unable to open estimated selection efficiency ROOT file." << std::endl;
    }
    TFile* fPowhegNonPrompt = new TFile("../Data/MonteCarlo/POWHEG/OldRun2/trees_powheg_fd_central.root","read");
    if (!fPowhegNonPrompt || fPowhegNonPrompt->IsZombie()) {
        std::cerr << "Error: Unable to open POWHEG non-prompt data ROOT file." << std::endl;
    }
    TFile* fSimulatedMCMatched = new TFile("../" + binning.inputMC.second + "/AO2D_mergedDFs.root","read");
    if (!fSimulatedMCMatched || fSimulatedMCMatched->IsZombie()) {
        std::cerr << "Error: Unable to open simulated MC matched data ROOT file." << std::endl;
    }
    TFile* fLumi = new TFile("../" + binning.inputDATA.second + "/AO2D_mergedDFs.root","read");
    if (!fLumi || fLumi->IsZombie()) {
        std::cerr << "Error: Unable to open luminosity data ROOT file." << std::endl;
    }

    // 1 - Create histograms
    FeedDownData dataContainer = createHistograms(binning);

    // 2 - Fill data
    fillHistograms(fEfficiency, fPowhegNonPrompt, fSimulatedMCMatched, fLumi, dataContainer, jetptMin, jetptMax, binning);

    // 3 - Calculate luminosities for scaling
    std::vector<double> lumiValues = retrieveLuminosityFromFile(fPowhegNonPrompt, fLumi, dataContainer);

    // 4 - Smearing and folding of POWHEG data to detector level
    smearGeneratorData(dataContainer, fEfficiency);

    // 5 - Subtract feed-down contribution from measured data distribution
    feedDown(fEfficiency, dataContainer, jetptMin, jetptMax);

    // 6 - Plot histograms
    plotHistograms(dataContainer, jetptMin, jetptMax, binning);

    // 7 - Save histograms to output file
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
    FeedDownSubtraction();
    return 0;
}
