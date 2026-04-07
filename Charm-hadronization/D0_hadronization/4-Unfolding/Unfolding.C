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

struct UnfoldData {
    // Step 1: apply efficiency (pT,D dependent) correction
    std::pair<TH1D*, TH1D*> hSelectionEfficiency;                           // first = prompt D0s, second = non-prompt D0s

    // Step 2: background subtracted, efficiency corrected, B fed-down 2D distribution (pT,jet vs DeltaR)
    TH2D* hBFedDownData;

    // Step 3: bayesian unfolding
    std::vector<TH2D*> hKineEffParticle = {nullptr, nullptr, nullptr};      // particle level (addition): [0] = numerator, [1] = denominator, [2] = efficiency
    std::vector<TH2D*> hKineEffDetector = {nullptr, nullptr, nullptr};      // detector level (removal): [0] = numerator, [1] = denominator, [2] = efficiency
    TH2D* hBFedDownDataKinCorrected = nullptr;                              // 2D histogram with detector level kinematic efficiency correction
    RooUnfoldResponse* response;                                            // response matrix for folding
    std::vector<TH2D*> hResponse2D = {nullptr, nullptr};                    // response projections matrix: first = DeltaR, second = pT,jet
    std::vector<TH2D*> hResponse2DJetptRange = {nullptr, nullptr, nullptr};                    // response projections matrix: first = DeltaR, second = pT,jet
    std::vector<RooUnfoldBayes*> unfold;                                    // unfolding objects, there are iterationNumber unfolding objects
    TH2D* hMeasuredTemplate = nullptr;                                      // template for measured data (jet pT vs DeltaR), used for hUnfolded binning
    TH2D* hTruthTemplate = nullptr;
    std::vector<TH2D*> hUnfolded;                                           // unfolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects
    std::vector<TH2D*> hUnfoldedKinCorrected;                               // unfolded 2D histogram with detector level kinematic correction, there are iterationNumber unfolding objects

    // Step 4: Refolding test
    std::vector<TH2D*> hRefolded;                                           // refolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects

    // Step 5: Convergence test
    std::vector<TH1D*> hConvergenceTest;                                    // between sucessive iterations
    std::vector<TH1D*> hConvergenceTest2;                                   // between unfolded and original fed-down subtracted data

    // Investigation of migration problem
    std::vector<TH1D*> hDeltaRMigration = {nullptr, nullptr};                // particle-level, detector-level
    std::vector<TH1D*> hJetptMigration = {nullptr, nullptr};                // particle-level, detector-level
    std::vector<TH1D*> hJetetaMigration = {nullptr, nullptr};                // particle-level, detector-level
    std::vector<TH1D*> hJetphiMigration = {nullptr, nullptr};                // particle-level, detector-level
    std::vector<TH1D*> hJetNConstMigration = {nullptr, nullptr};                // particle-level, detector-level
    std::vector<TH1D*> hHfptMigration = {nullptr, nullptr};                // particle-level, detector-level
    std::vector<TH1D*> hHfetaMigration = {nullptr, nullptr};                // particle-level, detector-level
    std::vector<TH1D*> hHfphiMigration = {nullptr, nullptr};                // particle-level, detector-level
};

// Module to create TH2D histograms including interest variable
UnfoldData createHistograms(const BinningStruct& binning, const int& iterationNumber) {
                              //const double& jetptMin, const double& jetptMax) {
    // Create struct to store data
    UnfoldData dataContainer;
    
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
    // for jetpt in [10;20] GeV/c
    dataContainer.hResponse2DJetptRange[0] = new TH2D("hResponseDeltaRJetptRange0", "Response matrix projection on #DeltaR in p_{T,jet}^{gen}#in[10,20] GeV/c;#DeltaR^{reco};#DeltaR^{gen}", 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hResponse2DJetptRange[1] = new TH2D("hResponsePtJetJetptRange1", "Response matrix projection on #DeltaR in p_{T,jet}^{reco}#in[10,20] GeV/c;#DeltaR^{reco};#DeltaR^{gen}", 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hResponse2DJetptRange[2] = new TH2D("hResponsePtJetJetptRange2", "Response matrix projection on #DeltaR in p_{T,jet}^{gen},p_{T,jet}^{reco}#in[10,20] GeV/c;#DeltaR^{reco};#DeltaR^{gen}", 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    
    // Investigation of migration problem
    dataContainer.hDeltaRMigration[0] = new TH1D("hDeltaRMigrationMCP", "Migration entries p_{T,jet}#in[10,20]GeV/c;#DeltaR^{gen};dN", 400, 0., 0.2);
    dataContainer.hDeltaRMigration[1] = new TH1D("hDeltaRMigrationMCD", "Migration entries p_{T,jet}#in[10,20]GeV/c;#DeltaR^{reco};dN", 400, 0., 0.2);

    dataContainer.hJetptMigration[0] = new TH1D("hJetptMigrationMCP", "Migration entries p_{T,jet}#in[10,20]GeV/c;p_{T,jet}^{gen};dN", 50, 10., 20.);
    dataContainer.hJetptMigration[1] = new TH1D("hJetptMigrationMCD", "Migration entries p_{T,jet}#in[10,20]GeV/c;p_{T,jet}^{reco};dN", 50, 10., 20.);

    dataContainer.hJetetaMigration[0] = new TH1D("hJetetaMigrationMCP", "Migration entries p_{T,jet}#in[10,20]GeV/c;#eta_{jet}^{gen};dN", 20, -0.4, 0.4);
    dataContainer.hJetetaMigration[1] = new TH1D("hJetetaMigrationMCD", "Migration entries p_{T,jet}#in[10,20]GeV/c;#eta_{jet}^{reco};dN", 20, -0.4, 0.4);

    dataContainer.hJetphiMigration[0] = new TH1D("hJetphiMigrationMCP", "Migration entries p_{T,jet}#in[10,20]GeV/c;#phi_{jet}^{gen};dN", 50, 0., 2 * TMath::Pi());
    dataContainer.hJetphiMigration[1] = new TH1D("hJetphiMigrationMCD", "Migration entries p_{T,jet}#in[10,20]GeV/c;#phi_{jet}^{reco};dN", 50, 0., 2 * TMath::Pi());

    dataContainer.hJetNConstMigration[0] = new TH1D("hJetNConstMigrationMCP", "Migration entries p_{T,jet}#in[10,20]GeV/c;N_{const}^{gen};dN", 25, 0., 25);
    dataContainer.hJetNConstMigration[1] = new TH1D("hJetNConstMigrationMCD", "Migration entries p_{T,jet}#in[10,20]GeV/c;N_{const}^{reco};dN", 25, 0., 25);

    dataContainer.hHfptMigration[0] = new TH1D("hHfptMigrationMCP", "Migration entries p_{T,jet}#in[10,20]GeV/c;p_{T,HF}^{gen};dN", 50, 0., 20.);
    dataContainer.hHfptMigration[1] = new TH1D("hHfptMigrationMCD", "Migration entries p_{T,jet}#in[10,20]GeV/c;p_{T,HF}^{reco};dN", 50, 0., 20.);

    dataContainer.hHfetaMigration[0] = new TH1D("hHfetaMigrationMCP", "Migration entries p_{T,jet}#in[10,20]GeV/c;#eta_{HF}^{gen};dN", 20, -0.4, 0.4);
    dataContainer.hHfetaMigration[1] = new TH1D("hHfetaMigrationMCD", "Migration entries p_{T,jet}#in[10,20]GeV/c;#eta_{HF}^{reco};dN", 20, -0.4, 0.4);
    
    dataContainer.hHfphiMigration[0] = new TH1D("hHfphiMigrationMCP", "Migration entries p_{T,jet}#in[10,20]GeV/c;#phi_{HF}^{gen};dN", 50, 0., 2 * TMath::Pi());
    dataContainer.hHfphiMigration[1] = new TH1D("hHfphiMigrationMCD", "Migration entries p_{T,jet}#in[10,20]GeV/c;#phi_{HF}^{reco};dN", 50, 0., 2 * TMath::Pi());

    // Reserve space for the unfolding objects with the number of iterations
    dataContainer.unfold.resize(iterationNumber);

    // Reserve space for the unfolded TH2D* histograms with the number of iterations (and point to nullptr)
    dataContainer.hUnfolded.resize(iterationNumber, nullptr);

    // And also for the kinematic efficiency corrected version of the previous mentioned (and point to nullptr)
    dataContainer.hUnfoldedKinCorrected.resize(iterationNumber, nullptr);

    std::cout << "Template histograms and response object created." << std::endl;
    return dataContainer;
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
void fillHistograms(TFile* fFeedDown, TFile* fEfficiency, TFile* fSimulatedMCMatched, UnfoldData& dataContainer, const BinningStruct& binning) {

    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = binning.deltaRBinEdges_particle[binning.deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double jetptMin = binning.ptjetBinEdges_detector[0];
    const double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    const double hfptMin = binning.ptHFBinEdges_detector[0];
    const double hfptMax = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1];
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

    // Access prompt and non-prompt efficiency
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
    float MCPjetNConst;
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

    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // Apply prompt selection (i.e., only c → D0, not B → D0)
        bool isReflection = (MCDhfMatchedFrom != MCDhfSelectedAs) ? true : false;
        if (!MCDhfprompt || isReflection || (MCDjetNConst < 0)) {
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
        
        bool recoInside = (MCDDeltaR >= 0.) && (MCDDeltaR < 0.01);
        bool genOutside = (MCPDeltaR >= 0.01);
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
            if ((MCPjetPt >= 10.) && (MCPjetPt < 20.)) {
                dataContainer.hResponse2DJetptRange[0]->Fill(MCDDeltaR, MCPDeltaR, 1./efficiency_prompt);
            }
            if ((MCDjetPt >= 10.) && (MCDjetPt < 20.)) {
                //std::cout << "MCDDeltaR = " << MCDDeltaR << "\tMCPDeltaR = " << MCPDeltaR << "\tprompt_effiency = " << efficiency_prompt << std::endl;
                dataContainer.hResponse2DJetptRange[1]->Fill(MCDDeltaR, MCPDeltaR, 1./efficiency_prompt);
            }
            if ((MCPjetPt >= 10.) && (MCPjetPt < 20.) && (MCDjetPt >= 10.) && (MCDjetPt < 20.) && recoInside && genOutside) { //  && recoInside && genOutside
                dataContainer.hResponse2DJetptRange[2]->Fill(MCDDeltaR, MCPDeltaR, 1./efficiency_prompt);
                dataContainer.hDeltaRMigration[0]->Fill(MCPDeltaR); // or ->Fill(MCPDeltaR, 1./efficiency_prompt);
                dataContainer.hDeltaRMigration[1]->Fill(MCDDeltaR);

                dataContainer.hJetptMigration[0]->Fill(MCPjetPt);
                dataContainer.hJetptMigration[1]->Fill(MCDjetPt);

                dataContainer.hJetetaMigration[0]->Fill(MCPjetEta);
                dataContainer.hJetetaMigration[1]->Fill(MCDjetEta);

                dataContainer.hJetphiMigration[0]->Fill(MCPjetPhi);
                dataContainer.hJetphiMigration[1]->Fill(MCDjetPhi);

                dataContainer.hJetNConstMigration[0]->Fill(MCPjetNConst);
                dataContainer.hJetNConstMigration[1]->Fill(MCDjetNConst);

                dataContainer.hHfptMigration[0]->Fill(MCPhfPt);
                dataContainer.hHfptMigration[1]->Fill(MCDhfPt);

                dataContainer.hHfetaMigration[0]->Fill(MCPhfEta);
                dataContainer.hHfetaMigration[1]->Fill(MCDhfEta);

                dataContainer.hHfphiMigration[0]->Fill(MCPhfPhi);
                dataContainer.hHfphiMigration[1]->Fill(MCDhfPhi);
            }
            
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
        dataContainer.hRefolded[iHisto] = manualFolding(*dataContainer.response, dataContainer.hRefolded[iHisto], dataContainer.hMeasuredTemplate);

        // Apply particle level kinematic efficiency correction (divide)
        //dataContainer.hRefolded[iHisto]->Divide(dataContainer.hKineEffDetector[2]);

        dataContainer.hRefolded[iHisto]->SetTitle(Form("Refolded with %zu iterations", iHisto+1));
        // Clean up intermediate clone (optional if not reused)
        //delete hRefoldInput;

    }
    
    

    std::cout << "Refolding test performed." << std::endl;
}

std::vector<TH1D*> convergenceTest(UnfoldData& dataContainer, const std::vector<TH2D*>& hUnfoldedKinCorrected) {
    // First convergence test: ratio between successive iterations
    //std::vector<TH1D*> hConvergenceTest;

    // Compute ratios between successive iterations
    for (size_t iHisto = 1; iHisto < hUnfoldedKinCorrected.size(); iHisto++) {

        // Project 2D histograms to 1D (along X axis)
        TH1D* hCurrent = hUnfoldedKinCorrected[iHisto]->ProjectionY(Form("projCurrent_%zu", iHisto));
        TH1D* hPrevious = hUnfoldedKinCorrected[iHisto-1]->ProjectionY(Form("projPrevious_%zu", iHisto-1));

        // Clone the current histogram to hold the ratio
        TH1D* hRatio = (TH1D*)hCurrent->Clone(Form("ratioIter_%zu", iHisto));
        hRatio->SetTitle(Form("Iteration %zu / %zu ratio", iHisto, iHisto-1));
        
        // Divide by the previous iteration
        hRatio->Divide(hPrevious);
        
        // Optional: style
        hRatio->SetLineColor(kBlack + iHisto);
        //hRatio->SetMarkerStyle(20);
        
        dataContainer.hConvergenceTest.push_back(hRatio);
    }

    // Quantitative check using RMS deviation
    std::vector<double> rmsValues;
    for (size_t iHisto = 0; iHisto < dataContainer.hConvergenceTest.size(); iHisto++) {
        TH1D* hRatio = dataContainer.hConvergenceTest[iHisto];
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

    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t iHisto = 0; iHisto < dataContainer.hConvergenceTest.size(); iHisto++) {
        legend->AddEntry(dataContainer.hConvergenceTest[iHisto], Form("Iteration %zu / %zu's RMS: %.5f", iHisto+1, iHisto, rmsValues[iHisto]), "l");
    }
    legend->Draw();

    std::cout << "Convergence test performed." << std::endl;

    // Second convergence test: ratio to original fed-down data
    TCanvas* cRatio2 = new TCanvas("cRatioUnfoldingConvergence2", "Unfolding Convergence Ratios to Original Fed-Down Data", 1200, 800);
    gStyle->SetOptStat(0);
    cRatio2->cd();
    TH1D* hFedDown1D = dataContainer.hBFedDownDataKinCorrected->ProjectionY("OriginalFromFeedDown");
    TH1D* hUnfoldedKinCorrectedFirstIter = hUnfoldedKinCorrected[0]->ProjectionY("FirstIterProjection");
    //std::vector<TH1D*> hConvergenceTest2;
    // Compute ratios between successive iterations
    for (size_t iHisto = 1; iHisto < hUnfoldedKinCorrected.size(); iHisto++) {

        // Project 2D histograms to 1D (along X axis)
        TH1D* hCurrent = hUnfoldedKinCorrected[iHisto]->ProjectionY(Form("projCurrent_%zu", iHisto));
        //TH1D* hPrevious = hUnfoldedKinCorrected[iHisto-1]->ProjectionY(Form("projPrevious_%zu", iHisto-1));

        // Clone the current histogram to hold the ratio
        TH1D* hRatio = (TH1D*)hCurrent->Clone(Form("ratio2Iter_%zu", iHisto));
        hRatio->SetTitle(Form("Convergence test: iterations divided by the first iteration"));
        
        // Divide by the previous iteration
        hRatio->Divide(hUnfoldedKinCorrectedFirstIter);
        
        // Optional: style
        hRatio->SetLineColor(kBlack + iHisto);
        //hRatio->SetMarkerStyle(20);
        
        dataContainer.hConvergenceTest2.push_back(hRatio);
        dataContainer.hConvergenceTest2[dataContainer.hConvergenceTest2.size()-1]->Draw((iHisto == 1) ? "" : "SAME");
    }
    TLegend* origlegend = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t iHisto = 0; iHisto < dataContainer.hConvergenceTest2.size(); iHisto++) {
        origlegend->AddEntry(dataContainer.hConvergenceTest2[iHisto], Form("Iteration %zu / original", iHisto+1), "l");
    }
    origlegend->Draw();

    return dataContainer.hConvergenceTest2;

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
    TCanvas* cKinEffParticle = new TCanvas("cKinEffParticle","Particle level kinematic efficiency");
    cKinEffParticle->cd();
    // gPad->SetLogx();
    dataContainer.hKineEffParticle[2]->SetTitle("Particle level kinematic efficiency (prompt D^{0}'s)");
    dataContainer.hKineEffParticle[2]->Draw("text");
    TCanvas* cKinEffDetector = new TCanvas("cKinEffDetector","Detector level kinematic efficiency");
    cKinEffDetector->cd();
    // gPad->SetLogx();
    dataContainer.hKineEffDetector[2]->SetTitle("Detector level kinematic efficiency (prompt D^{0}'s)");
    dataContainer.hKineEffDetector[2]->Draw("text");

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
    TCanvas* cResponseDeltaR = new TCanvas("cResponseDeltaR","Response matrix DeltaR projection");
    cResponseDeltaR->cd();
    dataContainer.hResponse2D[0]->Draw("colz");
    TCanvas* cResponsePtJet = new TCanvas("cResponsePtJet","Response matrix pT,jet projection");
    cResponsePtJet->cd();
    dataContainer.hResponse2D[1]->Draw("colz");

    TCanvas* cResponseJetptRanges = new TCanvas("cResponseJetptRanges","Response matrix 2D representation for jet pt in 10-20 GeV");
    cResponseJetptRanges->Divide(2,2);
    cResponseJetptRanges->cd(1);
    dataContainer.hResponse2DJetptRange[0]->Draw("colz");
    cResponseJetptRanges->cd(2);
    dataContainer.hResponse2DJetptRange[1]->Draw("colz");
    cResponseJetptRanges->cd(3);
    dataContainer.hResponse2DJetptRange[2]->Draw("colz");
    cResponseJetptRanges->cd(4);
    TH1D* hResponse2DJetptRangeReco = (TH1D*) dataContainer.hResponse2DJetptRange[2]->ProjectionX("hResponse2DJetptRangeReco");
    TH1D* hResponse2DJetptRangeGen = (TH1D*) dataContainer.hResponse2DJetptRange[2]->ProjectionY("hResponse2DJetptRangeGen");
    hResponse2DJetptRangeReco->SetLineWidth(2);
    hResponse2DJetptRangeReco->SetLineColor(kGreen+1);
    hResponse2DJetptRangeGen->SetLineWidth(2);
    hResponse2DJetptRangeGen->SetLineColor(kRed+1);
    if (hResponse2DJetptRangeReco->GetMaximum() > hResponse2DJetptRangeGen->GetMaximum()) {
        hResponse2DJetptRangeReco->Draw();
        hResponse2DJetptRangeGen->Draw("same");
    } else {
        hResponse2DJetptRangeGen->Draw();
        hResponse2DJetptRangeReco->Draw("same");
    }
    TLegend* lResponse2DJetptRange = new TLegend(0.6,0.57,0.7,0.77);
    lResponse2DJetptRange->AddEntry(hResponse2DJetptRangeReco,"Reco", "le");
    lResponse2DJetptRange->AddEntry(hResponse2DJetptRangeGen,"Gen", "le");
    lResponse2DJetptRange->Draw();

    std::vector<TLegend*> legMigration;
    TCanvas* cMigrationEntries = new TCanvas("cMigrationEntries","Migration entries");
    cMigrationEntries->Divide(3,3);
    cMigrationEntries->cd(1);
    dataContainer.hDeltaRMigration[0]->SetLineWidth(2);
    dataContainer.hDeltaRMigration[0]->SetLineColor(kRed+1);
    dataContainer.hDeltaRMigration[1]->SetLineWidth(2);
    dataContainer.hDeltaRMigration[1]->SetLineColor(kGreen+1);
    dataContainer.hDeltaRMigration[1]->Draw();
    dataContainer.hDeltaRMigration[0]->Draw("same");
    legMigration.push_back(new TLegend(0.6,0.57,0.7,0.77));
    legMigration[0]->AddEntry(dataContainer.hDeltaRMigration[0],"Gen", "le");
    legMigration[0]->AddEntry(dataContainer.hDeltaRMigration[1],"Reco", "le");
    legMigration[0]->Draw();
    cMigrationEntries->cd(2);
    dataContainer.hJetptMigration[0]->SetLineWidth(2);
    dataContainer.hJetptMigration[0]->SetLineColor(kRed+1);
    dataContainer.hJetptMigration[1]->SetLineWidth(2);
    dataContainer.hJetptMigration[1]->SetLineColor(kGreen+1);
    dataContainer.hJetptMigration[1]->Draw();
    dataContainer.hJetptMigration[0]->Draw("same");
    legMigration.push_back(new TLegend(0.6,0.57,0.7,0.77));
    legMigration[1]->AddEntry(dataContainer.hJetptMigration[0],"Gen", "le");
    legMigration[1]->AddEntry(dataContainer.hJetptMigration[1],"Reco", "le");
    legMigration[1]->Draw();
    cMigrationEntries->cd(3);
    dataContainer.hJetetaMigration[0]->SetLineWidth(2);
    dataContainer.hJetetaMigration[0]->SetLineColor(kRed+1);
    dataContainer.hJetetaMigration[1]->SetLineWidth(2);
    dataContainer.hJetetaMigration[1]->SetLineColor(kGreen+1);
    dataContainer.hJetetaMigration[1]->Draw();
    dataContainer.hJetetaMigration[0]->Draw("same");
    legMigration.push_back(new TLegend(0.6,0.57,0.7,0.77));
    legMigration[2]->AddEntry(dataContainer.hJetetaMigration[0],"Gen", "le");
    legMigration[2]->AddEntry(dataContainer.hJetetaMigration[1],"Reco", "le");
    legMigration[2]->Draw();
    cMigrationEntries->cd(4);
    dataContainer.hJetphiMigration[0]->SetLineWidth(2);
    dataContainer.hJetphiMigration[0]->SetLineColor(kRed+1);
    dataContainer.hJetphiMigration[1]->SetLineWidth(2);
    dataContainer.hJetphiMigration[1]->SetLineColor(kGreen+1);
    dataContainer.hJetphiMigration[1]->Draw();
    dataContainer.hJetphiMigration[0]->Draw("same");
    legMigration.push_back(new TLegend(0.6,0.57,0.7,0.77));
    legMigration[3]->AddEntry(dataContainer.hJetphiMigration[0],"Gen", "le");
    legMigration[3]->AddEntry(dataContainer.hJetphiMigration[1],"Reco", "le");
    legMigration[3]->Draw();
    cMigrationEntries->cd(5);
    dataContainer.hJetNConstMigration[0]->SetLineWidth(2);
    dataContainer.hJetNConstMigration[0]->SetLineColor(kRed+1);
    dataContainer.hJetNConstMigration[1]->SetLineWidth(2);
    dataContainer.hJetNConstMigration[1]->SetLineColor(kGreen+1);
    dataContainer.hJetNConstMigration[1]->Draw();
    dataContainer.hJetNConstMigration[0]->Draw("same");
    legMigration.push_back(new TLegend(0.6,0.57,0.7,0.77));
    legMigration[4]->AddEntry(dataContainer.hJetNConstMigration[0],"Gen", "le");
    legMigration[4]->AddEntry(dataContainer.hJetNConstMigration[1],"Reco", "le");
    legMigration[4]->Draw();
    cMigrationEntries->cd(6);
    dataContainer.hHfptMigration[0]->SetLineWidth(2);
    dataContainer.hHfptMigration[0]->SetLineColor(kRed+1);
    dataContainer.hHfptMigration[1]->SetLineWidth(2);
    dataContainer.hHfptMigration[1]->SetLineColor(kGreen+1);
    dataContainer.hHfptMigration[1]->Draw();
    dataContainer.hHfptMigration[0]->Draw("same");
    legMigration.push_back(new TLegend(0.6,0.57,0.7,0.77));
    legMigration[5]->AddEntry(dataContainer.hHfptMigration[0],"Gen", "le");
    legMigration[5]->AddEntry(dataContainer.hHfptMigration[1],"Reco", "le");
    legMigration[5]->Draw();
    cMigrationEntries->cd(7);
    dataContainer.hHfetaMigration[0]->SetLineWidth(2);
    dataContainer.hHfetaMigration[0]->SetLineColor(kRed+1);
    dataContainer.hHfetaMigration[0]->Draw();
    dataContainer.hHfetaMigration[1]->SetLineWidth(2);
    dataContainer.hHfetaMigration[1]->SetLineColor(kGreen+1);
    dataContainer.hHfetaMigration[1]->Draw();
    dataContainer.hHfetaMigration[0]->Draw("same");
    legMigration.push_back(new TLegend(0.6,0.57,0.7,0.77));
    legMigration[6]->AddEntry(dataContainer.hHfetaMigration[0],"Gen", "le");
    legMigration[6]->AddEntry(dataContainer.hHfetaMigration[1],"Reco", "le");
    legMigration[6]->Draw();
    cMigrationEntries->cd(8);
    dataContainer.hHfphiMigration[0]->SetLineWidth(2);
    dataContainer.hHfphiMigration[0]->SetLineColor(kRed+1);
    dataContainer.hHfphiMigration[0]->Draw();
    dataContainer.hHfphiMigration[1]->SetLineWidth(2);
    dataContainer.hHfphiMigration[1]->SetLineColor(kGreen+1);
    dataContainer.hHfphiMigration[1]->Draw();
    dataContainer.hHfphiMigration[0]->Draw("same");
    legMigration.push_back(new TLegend(0.6,0.57,0.7,0.77));
    legMigration[7]->AddEntry(dataContainer.hHfphiMigration[0],"Gen", "le");
    legMigration[7]->AddEntry(dataContainer.hHfphiMigration[1],"Reco", "le");
    legMigration[7]->Draw();

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
        hUnfoldedKinCorrectedProj[iIter]->GetYaxis()->SetTitle("dN");
        if (iIter == 0) {
            hUnfoldedKinCorrectedProj[iIter]->GetYaxis()->SetRangeUser(0.,hUnfoldedKinCorrectedProj[iIter]->GetMaximum()*1.2);
        }
        lUnfoldedIter->AddEntry(hUnfoldedKinCorrectedProj[iIter],Form("Iteration %zu", iIter+1), "le");
        if (iIter == 0) {
            hUnfoldedKinCorrectedProj[iIter]->SetTitle("Unfolded with #epsilon^{particle}_{kin} correction");
            hUnfoldedKinCorrectedProj[iIter]->Draw();
        } else {
            hUnfoldedKinCorrectedProj[iIter]->Draw("same");
        }
        
    }
    lUnfoldedIter->Draw();
    double reportingJetPtMin = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetBinLowEdge(secondBin);
    double reportingJetPtMax = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));

    TCanvas* cUnfoldedIterNorm = new TCanvas("cUnfoldedIterNorm","Unfolded histograms for each iteration self-normalized");
    cUnfoldedIterNorm->cd();
    TLegend* lUnfoldedIterNorm = new TLegend(0.6,0.57,0.7,0.77);
    std::vector<TH1D*> hUnfoldedKinCorrectedProjNorm(dataContainer.hUnfoldedKinCorrected.size());
    secondBin = 2;
    lastButOneBin = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetNbins() - 1; // penultimate bin
    for (size_t iIter = 0; iIter < dataContainer.hUnfoldedKinCorrected.size(); iIter++) {
        // Project Y axis (DeltaR) using X axis (pTjet) range [secondBin, oneBeforeLastBin] (excluding the padding bins)
        hUnfoldedKinCorrectedProjNorm[iIter] = dataContainer.hUnfoldedKinCorrected[iIter]->ProjectionY(Form("hProjIterNorm_%zu", iIter),secondBin, lastButOneBin); // bins specified in X (pT,jet) dimension
        hUnfoldedKinCorrectedProjNorm[iIter]->SetLineColor(kBlack + iIter);
        hUnfoldedKinCorrectedProjNorm[iIter]->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
        hUnfoldedKinCorrectedProjNorm[iIter]->Scale(1.0 / hUnfoldedKinCorrectedProjNorm[iIter]->Integral(),"width"); // self-normalization
        lUnfoldedIterNorm->AddEntry(hUnfoldedKinCorrectedProjNorm[iIter],Form("Iteration %zu", iIter+1), "le");
        if (iIter == 0) {
            hUnfoldedKinCorrectedProjNorm[iIter]->SetTitle("Unfolded with #epsilon^{particle}_{kin} correction, self-normalized");
            hUnfoldedKinCorrectedProjNorm[iIter]->GetYaxis()->SetRangeUser(0.,hUnfoldedKinCorrectedProjNorm[iIter]->GetMaximum()*1.2);
            hUnfoldedKinCorrectedProjNorm[iIter]->Draw();
        } else {
            hUnfoldedKinCorrectedProjNorm[iIter]->Draw("same");
        }
    }
    lUnfoldedIterNorm->Draw();
    reportingJetPtMin = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetBinLowEdge(secondBin);
    reportingJetPtMax = dataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
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

    // Plot input distribution fully fed-down, self normalized with bin width for later comparison
    TCanvas* cFeddownSelfNorm = new TCanvas("cFeddownSelfNorm","Input distribution fully fed-down, self normalized with bin width");
    cFeddownSelfNorm->Divide(2,2);
    std::vector<TH1D*> hFeddownSelfNorm;
    hFeddownSelfNorm.push_back(dataContainer.hBFedDownData->ProjectionY("hFeddownSelfNorm_0",2,dataContainer.hBFedDownData->GetXaxis()->GetNbins() - 1));
    hFeddownSelfNorm[0]->GetYaxis()->SetTitle("dN");
    hFeddownSelfNorm.push_back((TH1D*) hFeddownSelfNorm[0]->Clone("hFeddownSelfNorm_1"));
    hFeddownSelfNorm[1]->Scale(1.0, "width");
    hFeddownSelfNorm[1]->GetYaxis()->SetTitle("#frac{dN}{d#DeltaR}");
    hFeddownSelfNorm.push_back((TH1D*) hFeddownSelfNorm[1]->Clone("hFeddownSelfNorm_2"));
    hFeddownSelfNorm[2]->Scale(1 / hFeddownSelfNorm[2]->Integral());
    hFeddownSelfNorm[2]->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
    cFeddownSelfNorm->cd(1);
    hFeddownSelfNorm[0]->Draw();
    cFeddownSelfNorm->cd(2);
    hFeddownSelfNorm[1]->Draw();
    cFeddownSelfNorm->cd(3);
    hFeddownSelfNorm[2]->Draw();

    TCanvas* cStepsReportingRange = new TCanvas("cStepsReportingRange","Unfolding steps for reporting range");
    cStepsReportingRange->Divide(4,2);
    cStepsReportingRange->cd(1);
    dataContainer.hBFedDownData->Draw("colz");
    cStepsReportingRange->cd(2);
    TH1D* hBFedDownDataRange = dataContainer.hBFedDownData->ProjectionY("hBFedDownDataProjYRange",2,dataContainer.hBFedDownData->GetXaxis()->GetNbins() - 1);
    hBFedDownDataRange->Draw();
    cStepsReportingRange->cd(3);
    dataContainer.hBFedDownDataKinCorrected->Draw("colz");
    cStepsReportingRange->cd(4);
    TH1D* hBFedDownDataKinCorrectedRange = dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjYRange",2,dataContainer.hBFedDownDataKinCorrected->GetXaxis()->GetNbins() - 1);
    hBFedDownDataKinCorrectedRange->Draw();
    cStepsReportingRange->cd(5);
    dataContainer.hUnfolded[7]->Draw("colz");
    cStepsReportingRange->cd(6);
    TH1D* hUnfolded_7Range = dataContainer.hUnfolded[7]->ProjectionY("hUnfoldedProjYRange",2,dataContainer.hUnfolded[7]->GetXaxis()->GetNbins() - 1);
    hUnfolded_7Range->Draw();
    cStepsReportingRange->cd(7);
    dataContainer.hUnfoldedKinCorrected[7]->Draw("colz");
    cStepsReportingRange->cd(8);
    TH1D* hUnfoldedKinCorrected_7Range = dataContainer.hUnfoldedKinCorrected[7]->ProjectionY("hUnfoldedKinCorrectedProjYRange",2,dataContainer.hUnfoldedKinCorrected[7]->GetXaxis()->GetNbins() - 1);
    hUnfoldedKinCorrected_7Range->Draw();
    
    TCanvas* cRefoldedIter = new TCanvas("cRefoldedIter","Refolded histograms for each iteration");
    TPad* pad1 = new TPad("pad1","top pad",0,0.3,1,1);
    TPad* pad2 = new TPad("pad2","bottom pad",0,0,1,0.3);
    pad2->SetBottomMargin(0.3);
    pad1->SetBottomMargin(0.0);
    pad2->SetTopMargin(0.0);
    pad1->SetTickx();
    pad2->SetTickx();
    pad1->SetTicky();
    pad2->SetTicky();
    pad1->SetLeftMargin(0.12);
    pad2->SetLeftMargin(0.12);
    pad1->Draw();
    pad2->Draw();
    // cRefoldedIter->cd();
    TLegend* lRefoldedIter = new TLegend(0.5,0.57,0.85,0.87);
    std::vector<TH1D*> projRefolded; // to hold ratios with respect to the fed-down distribution for each unfolding iteration
    for (size_t iHisto = 0; iHisto < dataContainer.hRefolded.size(); iHisto++) {
        TH1D* hProj = dataContainer.hRefolded[iHisto]->ProjectionY(Form("hRefolded_iter%zu_Proj", iHisto+1));
        hProj->GetYaxis()->SetTitle("dN");
        hProj->Sumw2();
        hProj->SetFillColorAlpha(kBlack + iHisto, 0.3); // 0 = transparent, 1 = solid
        hProj->SetLineColor(kBlack + iHisto);
        projRefolded.push_back(hProj);
        lRefoldedIter->AddEntry(hProj,Form("Refolded Iteration %zu with #epsilon_{kin}^{part} correction", iHisto+1), "le");
        if (iHisto == 0) {
            hProj->SetTitle("Refolding test");
            hProj->GetYaxis()->SetRangeUser(1e-6, hProj->GetMaximum()*1.2);
            // hProj->GetYaxis()->SetNdivisions(510);   // top pad
            hProj->GetYaxis()->SetLabelOffset(0.02);
            hProj->GetXaxis()->SetLabelSize(0);
            hProj->GetXaxis()->SetTitleSize(0);
            hProj->GetYaxis()->ChangeLabel(1, -1, 0); 
            pad1->cd();
            hProj->Draw();
        } else {
            pad1->cd();
            hProj->Draw("same");
        }
        //lRefoldedIter->AddEntry(hProj, Form("Refolded Iteration %zu", iHisto+1), "le");
    }
    TH1D* hDataProj = dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjY");
    hDataProj->SetLineWidth(2);
    lRefoldedIter->AddEntry(hDataProj,"Original data with #epsilon_{kin}^{det} correction", "le");
    hDataProj->Draw("same");
    lRefoldedIter->Draw();
    pad2->cd(); // now draw the ratio
    for (size_t iHisto = 0; iHisto < projRefolded.size(); iHisto++) {
        TH1D* hRatio = (TH1D*) projRefolded[iHisto]->Clone(Form("hRatio_iter%zu", iHisto+1));
        hRatio->Divide(hDataProj);
        hRatio->SetFillColorAlpha(kBlack + iHisto, 0.3); // 0 = transparent, 1 = solid
        hRatio->SetLineColor(kBlack + iHisto);
        hRatio->SetTitle("");
        hRatio->GetYaxis()->SetTitle("#frac{Refolded_{i}}{Original}");
        int nDiv = hRatio->GetYaxis()->GetNdivisions();
        int nMajor = nDiv / 100; // ROOT encoding: nDiv = 100*N1 + 10*N2 + N3
        hRatio->GetYaxis()->ChangeLabel(nMajor, -1, 0);
        // hRatio->GetYaxis()->SetNdivisions(505);
        hRatio->GetYaxis()->SetTitleSize(0.08);
        hRatio->GetYaxis()->SetLabelSize(0.08);
        hRatio->GetXaxis()->SetTitleSize(0.1);
        hRatio->GetXaxis()->SetLabelSize(0.1);
        hRatio->GetYaxis()->SetTitleOffset(0.5);
        // hRatio->GetYaxis()->SetTitleSize(0.08);
        if (iHisto == 0) {
            hRatio->GetYaxis()->SetRangeUser(0.7, 2.0); // typical range
            hRatio->Draw();
        } else {
            hRatio->Draw("same");
        }
    }
    TLine* line = new TLine(hDataProj->GetXaxis()->GetXmin(),1.0,hDataProj->GetXaxis()->GetXmax(),1.0); // add unity line
    line->SetLineStyle(2);
    line->Draw();

    TCanvas* cRefoldedPtRanges = new TCanvas("cRefoldedPtRanges", "Refolded histograms for last iteration in different pT,jet ranges");
    cRefoldedPtRanges->cd();
    TLegend* lRefoldedPtRanges = new TLegend(0.5,0.57,0.85,0.87);
    // Get last iteration refolded histogram
    auto* hLastRefolded = dataContainer.hRefolded.back();
    // Define bin ranges
    int bin7  = hLastRefolded->GetXaxis()->FindBin(7.0);
    int bin10  = hLastRefolded->GetXaxis()->FindBin(10.0);
    int bin15 = hLastRefolded->GetXaxis()->FindBin(15.0);
    int bin30 = hLastRefolded->GetXaxis()->FindBin(30.0);
    // Create projections for refolded
    auto* hProj_7_10 = hLastRefolded->ProjectionY("hRefolded_iter_8_Proj_7_10GeV", bin7, bin10);
    hProj_7_10->SetLineColor(kBlue);
    hProj_7_10->Draw();
    lRefoldedPtRanges->AddEntry(hProj_7_10, "Refolded: 7 < p_{T,jet} < 10", "le");
    auto* hProj_10_15 = hLastRefolded->ProjectionY("hRefolded_iter_8_Proj_10_15GeV", bin10+1, bin15);
    hProj_10_15->SetLineColor(kBlack);
    hProj_10_15->Draw("same");
    lRefoldedPtRanges->AddEntry(hProj_10_15, "Refolded: 10 < p_{T,jet} < 15", "le");
    auto* hProj_15_30 = hLastRefolded->ProjectionY("hRefolded_iter_8_Proj_15_30GeV", bin15+1, bin30);
    hProj_15_30->SetLineColor(kRed);
    hProj_15_30->Draw("same");
    lRefoldedPtRanges->AddEntry(hProj_15_30, "Refolded: 15 < p_{T,jet} < 30", "le");
    // Create projections for original (BFedDown) data
    bin7  = dataContainer.hBFedDownDataKinCorrected->GetXaxis()->FindBin(7.0);
    bin10 = dataContainer.hBFedDownDataKinCorrected->GetXaxis()->FindBin(10.0);
    bin15 = dataContainer.hBFedDownDataKinCorrected->GetXaxis()->FindBin(15.0);
    bin30 = dataContainer.hBFedDownDataKinCorrected->GetXaxis()->FindBin(30.0);
    auto* hOrig_7_10 = dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjY_7_10GeV", bin7, bin10);
    hOrig_7_10->SetLineColor(kBlue);
    hOrig_7_10->SetLineStyle(2); // dashed line for distinction
    hOrig_7_10->Draw("same");
    lRefoldedPtRanges->AddEntry(hOrig_7_10, "Original: 7 < p_{T,jet} < 10", "le");
    auto* hOrig_10_15 = dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjY_10_15GeV", bin10+1, bin15);
    hOrig_10_15->SetLineColor(kBlack);
    hOrig_10_15->SetLineStyle(2); // dashed line for distinction
    hOrig_10_15->Draw("same");
    lRefoldedPtRanges->AddEntry(hOrig_10_15, "Original: 10 < p_{T,jet} < 15", "le");
    auto* hOrig_15_30 = dataContainer.hBFedDownDataKinCorrected->ProjectionY("hBFedDownDataKinCorrectedProjY_15_30GeV", bin15+1, bin30);
    hOrig_15_30->SetLineColor(kRed);
    hOrig_15_30->SetLineStyle(2);
    hOrig_15_30->Draw("same");
    lRefoldedPtRanges->AddEntry(hOrig_15_30, "Original: 15 < p_{T,jet} < 30", "le");
    lRefoldedPtRanges->Draw();

    // Self normalized final iteration distribution
    TCanvas* cSelfNorm = new TCanvas("cSelfNorm","Self-normalized unfolded distribution for last iteration");
    cSelfNorm->cd();
    TH2D* hLastUnfolded = dataContainer.hUnfoldedKinCorrected.back();
    bin7  = hLastUnfolded->GetXaxis()->FindBin(7.0);
    bin10 = hLastUnfolded->GetXaxis()->FindBin(10.0);
    bin15  = hLastUnfolded->GetXaxis()->FindBin(15.0);
    bin30 = hLastUnfolded->GetXaxis()->FindBin(30.0);
    TH1D* hLastUnfoldedProjY = hLastUnfolded->ProjectionY("hLastUnfoldedProjY");
    hLastUnfoldedProjY->SetTitle("Self-normalized unfolded distribution for last iteration;#DeltaR^{gen};#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
    hLastUnfoldedProjY->Sumw2();
    double integral = hLastUnfoldedProjY->Integral();
    if (integral != 0) {
        hLastUnfoldedProjY->Scale(1.0 / integral, "width");
    }
    hLastUnfoldedProjY->Draw();
    // Multiple pT,jet ranges
    TCanvas* cSelfNormRanges = new TCanvas("cSelfNormRanges","Self-normalized unfolded distribution for last iteration, multiple ranges");
    cSelfNormRanges->cd();
    TH1D* hLastUnfoldedProjY_7_10 = hLastUnfolded->ProjectionY("hLastUnfoldedProjY_7_10GeV", bin7, bin10);
    hLastUnfoldedProjY_7_10->SetTitle("Self-normalized unfolded distribution for last iteration;#DeltaR^{gen};#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
    hLastUnfoldedProjY_7_10->SetLineColor(kBlue);
    hLastUnfoldedProjY_7_10->Sumw2();
    double integral_7_10 = hLastUnfoldedProjY_7_10->Integral();
    if (integral_7_10 != 0) {
        hLastUnfoldedProjY_7_10->Scale(1.0 / integral_7_10, "width");
    }
    hLastUnfoldedProjY_7_10->Draw();
    TH1D* hLastUnfoldedProjY_10_15 = hLastUnfolded->ProjectionY("hLastUnfoldedProjY_10_15GeV", bin10+1, bin15);
    hLastUnfoldedProjY_10_15->SetLineColor(kBlack);
    hLastUnfoldedProjY_10_15->Sumw2();
    double integral_10_15 = hLastUnfoldedProjY_10_15->Integral();
    if (integral_10_15 != 0) {
        hLastUnfoldedProjY_10_15->Scale(1.0 / integral_10_15, "width");
    }
    hLastUnfoldedProjY_10_15->Draw("same");
    TH1D* hLastUnfoldedProjY_15_30 = hLastUnfolded->ProjectionY("hLastUnfoldedProjY_15_30GeV", bin15+1, bin30);
    hLastUnfoldedProjY_15_30->SetLineColor(kRed);
    hLastUnfoldedProjY_15_30->Sumw2();
    double integral_15_30 = hLastUnfoldedProjY_15_30->Integral();
    if (integral_15_30 != 0) {
        hLastUnfoldedProjY_15_30->Scale(1.0 / integral_15_30, "width");
    }
    hLastUnfoldedProjY_15_30->Draw("same");
    TLegend* lSelfNormRanges = new TLegend(0.5,0.57,0.85,0.87);
    lSelfNormRanges->AddEntry(hLastUnfoldedProjY_7_10, "7 < p_{T,jet} < 10", "le");
    lSelfNormRanges->AddEntry(hLastUnfoldedProjY_10_15, "10 < p_{T,jet} < 15", "le");
    lSelfNormRanges->AddEntry(hLastUnfoldedProjY_15_30, "15 < p_{T,jet} < 30", "le");
    lSelfNormRanges->Draw();

    // Convergence test plots
    // Second convergence test: ratio to original fed-down data
    TCanvas* cConvTest = new TCanvas("cConvTest", "Convergence Ratios to Original Fed-Down Data", 1200, 800);
    gStyle->SetOptStat(0);
    cConvTest->cd();
    TH1D* hFedDown1D = dataContainer.hBFedDownDataKinCorrected->ProjectionY("OriginalFromFeedDown_convtest");
    TLegend* origlegend = new TLegend(0.7, 0.7, 0.9, 0.9);
    dataContainer.hConvergenceTest2[0]->GetYaxis()->SetRangeUser(0.0, 2.0);
    dataContainer.hConvergenceTest2[0]->GetYaxis()->SetTitle("#frac{Unfolded_{i}}{Unfolded_{0}}");
    for (size_t iHisto = 0; iHisto < dataContainer.hConvergenceTest2.size(); iHisto++) {
        dataContainer.hConvergenceTest2[iHisto]->Draw((iHisto == 0) ? "" : "SAME");
        origlegend->AddEntry(dataContainer.hConvergenceTest2[iHisto], Form("Iteration %zu / Iteration 0", iHisto+1), "l");
    }
    origlegend->Draw();

    TCanvas* cConvTestRatio = new TCanvas("cRatioUnfoldingConvergence", "Unfolding Convergence Ratios", 1200, 800);
    TLegend* legConvTestRatio = new TLegend(0.3, 0.15, 0.5, 0.35);
    cConvTestRatio->cd();
    for (size_t iHisto = 0; iHisto < dataContainer.hConvergenceTest.size(); iHisto++) {
        if (iHisto == 0) {
            legConvTestRatio->AddEntry(dataContainer.hConvergenceTest[iHisto], Form("Iteration %zu / Iteration %zu", iHisto+2, iHisto+1), "l");
            dataContainer.hConvergenceTest[iHisto]->SetTitle("Consecutive ratios of unfolded distributions with different # of iterations;#DeltaR;#frac{Unfolded_{i+1}}{Unfolded_{i}}");
            dataContainer.hConvergenceTest[iHisto]->GetYaxis()->SetRangeUser(0.0, 2.0);
            dataContainer.hConvergenceTest[iHisto]->Draw("E");      // draw the first ratio
        } else {
            legConvTestRatio->AddEntry(dataContainer.hConvergenceTest[iHisto], Form("Iteration %zu / Iteration %zu", iHisto+2, iHisto+1), "l");
            dataContainer.hConvergenceTest[iHisto]->Draw("E SAME");        // draw others on the same canvas
        }
    }
    legConvTestRatio->Draw();

    //
    // Storing images
    //
    TString imagePath = "../Images/4-Unfolding/";
    cKinEff->Update();
    cKinEff->SaveAs(imagePath + "Unfolding_kin_efficiencies.png");
    cKinEffParticle->Update();
    cKinEffParticle->SaveAs(imagePath + "Unfolding_kin_efficiency_particle.png");
    cKinEffDetector->Update();
    cKinEffDetector->SaveAs(imagePath + "Unfolding_kin_efficiency_detector.png");
    cResponse->Update();
    cResponse->SaveAs(imagePath + "Unfolding_response.png");
    cResponseDeltaR->Update();
    cResponseDeltaR->SaveAs(imagePath + "Unfolding_response_deltar.png");
    cResponsePtJet->Update();
    cResponsePtJet->SaveAs(imagePath + "Unfolding_response_ptjet.png");
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
    cKinEffParticle->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cKinEffDetector->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponse->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponseDeltaR->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponsePtJet->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponseJetptRanges->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cStepsReportingRange->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cMigrationEntries->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cUnfoldedIter->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cUnfoldedIterNorm->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cRefoldedIter->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cRefoldedPtRanges->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cSelfNorm->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cSelfNormRanges->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cConvTest->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cConvTestRatio->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
    //cKinEfficiency->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    //cFolded->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    //cFedDownData->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.'pdf)",jetptMin,jetptMax));
    // in separate .pdf file
    cKinEffParticle->Print(imagePath + "Unfolding_kin_efficiency_particle.pdf");
    cKinEffDetector->Print(imagePath + "Unfolding_kin_efficiency_detector.pdf");
    cUnfoldedIter->Print(imagePath + "Unfolding_unfolded_iterations.pdf");
    cRefoldedIter->Print(imagePath + "Unfolding_refolded_comparison.pdf");
    

}

void saveData(const UnfoldData& dataContainer, const double& jetptMin, const double& jetptMax, const BinningStruct& binning) {
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
    // Return to root directory (optional)
    outFile->cd();

    // Also store the axes used for the histograms
    storeBinningInFile(outFile, binning);

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

    TFile* fEfficiency = new TFile(Form("../2-Efficiency/selection_efficiency_run3style_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fEfficiency || fEfficiency->IsZombie()) {
        std::cerr << "Error: Unable to open estimated selection efficiency ROOT file." << std::endl;
    }

    // Opening files
    TFile* fSimulatedMCMatched = new TFile("../Data/MonteCarlo/Train_645447/AO2D_mergedDFs.root","read");
    TFile* fFeedDown = new TFile(Form("../3-Feed-Down/outputFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fSimulatedMCMatched || fSimulatedMCMatched->IsZombie()) {
        std::cerr << "Error: Unable to open O2 MC matched ROOT file." << std::endl;
    }
    if (!fFeedDown || fFeedDown->IsZombie()) {
        std::cerr << "Error: Unable to open AnalysisResults.root data ROOT file." << std::endl;
    }
    
    //UnfoldData dataContainer = createHistograms(deltaRBinEdges, ptHFBinEdges, jetptMin, jetptMax);
    UnfoldData dataContainer = createHistograms(binning, iterationNumber);

    // Fill histograms with POWHEG simulation data
    fillHistograms(fFeedDown, fEfficiency, fSimulatedMCMatched, dataContainer, binning);

    // Perform unfolding with the chosen method and number of iterations
    performUnfolding(dataContainer, iterationNumber);

    // Fold the resulting distribution and compare with input detector level distribution
    refoldingTest(dataContainer);
    
    // Verify convergence and choose best iteration number
    std::vector<TH1D*> hConvergents = convergenceTest(dataContainer, dataContainer.hUnfoldedKinCorrected);

    // Plot the efficiency histogram and further corrected histograms
    plotHistograms(dataContainer, jetptMin, jetptMax);

    // Save corrected distributions to file
    saveData(dataContainer, jetptMin, jetptMax, binning);


    

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
