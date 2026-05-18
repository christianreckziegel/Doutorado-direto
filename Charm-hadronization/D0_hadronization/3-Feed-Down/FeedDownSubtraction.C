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
    
    // Step 2: apply  run 3 selection efficiency (pT,D dependent) correction
    std::pair<TH1D*, TH1D*> hSelEffRun3Style_detectorLevel;                 // first = prompt D0s, second = non-prompt D0s
    std::pair<TH1D*, TH1D*> hSelEffRun3Style_particleLevel;                 // first = prompt D0s, second = non-prompt D0s
    TH3D* hPowhegSelEffCorrected = nullptr;                                 // corrected by non-prompt D0 particle level selection efficiency

    // Step 3: particle level kinematic efficiency
    std::vector<TH3D*> hKineEffParticle = {nullptr, nullptr, nullptr};      // particle level (addition): [0] = numerator, [1] = denominator, [2] = efficiency
    TH3D* hPowhegSelKinEffCorrected = nullptr;                              // corrected by non-prompt D0 particle level selection efficiency AND by the appropriate particle level kinematic efficiency

    // Step 4: folding
    RooUnfoldResponse response;                                             // 6D response matrix for folding
    TH2D* hResponseDeltaR;                                                  // 2D response matrix projection on DeltaR
    TH2D* hResponsePtJet;                                                   // 2D response matrix projection on jet pT
    TH2D* hResponsePtHF;                                                    // 2D response matrix projection on HF pT
    TH3D* hMeasuredTemplate = nullptr;                                      // template of measured data (jet pT vs DeltaR vs HF pT), used for hFolded binning and response matrix creation
    TH3D* hFolded = nullptr;                                                // folded 3D histogram (jet pT vs DeltaR vs HF pT)

    // Step 5: detector level kinematic efficiency
    std::vector<TH3D*> hKineEffDetector = {nullptr, nullptr, nullptr};      // detector level (removal): [0] = numerator, [1] = denominator, [2] = efficiency
    TH3D* hFoldedKinCorrected = nullptr;                                    // folded 3D histogram with detector level kinematic correction

    // Step 6: non-prompt D0s corrected by prompt run 3 detector level selection efficiency
    TH3D* hFinalPromptSelEffCorrected = nullptr;

    // Step 7: scale to luminosity and branching ratio
    TH3D* hFinalScaled = nullptr;

    // Step 8: projection to 2D (jet pT vs DeltaR)
    TH2D* hFinalScaled2D = nullptr;

    // Metadata
    double lumiData = 1; // luminosity of the data in 1/pb
    double lumiMC = 1; // luminosity of the MC in 1/pb
    const double branchingRatio = 0.03947; // D0 -> KPi decay channel branching ratio = (3.947 +- 0.04) %

    // Total data distribution before feed-down subtraction, to be used for comparison with the final background subtracted distribution and for calculating the feed-down fraction
    TH2D* hTotalDataBeforeFeedDown = nullptr;

    // Step 9: remove non-prompt D0 contribution from data distribution
    TH2D* hBFedDownData;
};

FeedDownData createHistograms(const BinningStruct& binning) {
    // Create struct to store data
    FeedDownData dataContainer;

    //
    // Create 3D POWHEG data histogram
    //
    dataContainer.hPowheg = new TH3D("hPowheg", "POWHEG + PYTHIA;p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR^{gen};p_{T,D}^{gen}", binning.ptjetBinEdges_particle.size()-1, binning.ptjetBinEdges_particle.data(), 
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
    dataContainer.hKineEffParticle[0] = new TH3D("hKineEffParticleNumerator", "Particle level kinematic efficiency numerator (non-prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR^{gen};p_{T,D}^{gen}", 
                                                                                  binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data(),
                                                                                  binning.ptHFBinEdges_particle.size()-1, binning.ptHFBinEdges_particle.data());
    dataContainer.hKineEffParticle[1] = new TH3D("hKineEffParticleDenominator", "Particle level kinematic efficiency denominator (non-prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR^{gen};p_{T,D}^{gen}", 
                                                                                  binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data(),
                                                                                  binning.ptHFBinEdges_particle.size()-1, binning.ptHFBinEdges_particle.data());
    dataContainer.hKineEffDetector[0] = new TH3D("hKineEffDetectorNumerator", "Detector level kinematic efficiency numerator (non-prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR^{reco};p_{T,D}^{reco}", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data(),
                                                                                  binning.ptHFBinEdges_detector.size()-1, binning.ptHFBinEdges_detector.data());
    dataContainer.hKineEffDetector[1] = new TH3D("hKineEffDetectorDenominator", "Detector level kinematic efficiency denominator (non-prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR^{reco};p_{T,D}^{reco}", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data(),
                                                                                  binning.ptHFBinEdges_detector.size()-1, binning.ptHFBinEdges_detector.data());

    dataContainer.hMeasuredTemplate = new TH3D("hMeasuredTemplate", "Folded 3D histogram;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco};p_{T,D}^{reco}", 
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data(),
                                                                                  binning.ptHFBinEdges_detector.size()-1, binning.ptHFBinEdges_detector.data());
    //
    // Response matrix projections: DeltaR and pT,jet
    //
    dataContainer.hResponseDeltaR = new TH2D("hResponseDeltaR","Response matrix #DeltaR projection;#DeltaR^{reco};#DeltaR^{gen}",
                                                                                  binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data(),
                                                                                  binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    dataContainer.hResponsePtJet = new TH2D("hResponsePtJet","Response matrix p_{T,jet} projection;p_{T,jet}^{reco ch};p_{T,jet}^{gen}",
                                                                                  binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(),
                                                                                  binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data());
    dataContainer.hResponsePtHF = new TH2D("hResponsePtHF","Response matrix p_{T,HF} projection;p_{T,HF}^{reco};p_{T,HF}^{gen}",
                                                                                  binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data(),
                                                                                  binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());

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
void fillHistograms(TFile* fEfficiency, TFile* fPowhegNonPrompt, TFile* fSimulatedMCMatched, FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax, const BinningStruct& binning) {
    
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
    dataContainer.hSelEffRun3Style_detectorLevel.first = (TH1D*)fEfficiency->Get("hSelectionEfficiencyPrompt")->Clone("hSelEffRun3Style_detectorLevelPrompt");
    if (!dataContainer.hSelEffRun3Style_detectorLevel.first) {
        std::cout << "Error fetching prompt D0 selection efficiency histogram from file.\n";
    }
    dataContainer.hSelEffRun3Style_detectorLevel.second = (TH1D*)fEfficiency->Get("hSelectionEfficiencyNonPrompt")->Clone("hSelEffRun3Style_detectorLevelNonprompt");
    if (!dataContainer.hSelEffRun3Style_detectorLevel.second) {
        std::cout << "Error fetching non-prompt D0 selection efficiency histogram from file.\n";
    }
    dataContainer.hSelEffRun3Style_particleLevel.first = (TH1D*)fEfficiency->Get("selEffRun3Style_particleLevelPrompt")->Clone("hSelEffRun3Style_particleLevelPrompt");
    if (!dataContainer.hSelEffRun3Style_particleLevel.first) {
        std::cout << "Error fetching prompt D0 selection efficiency histogram from file.\n";
    }
    dataContainer.hSelEffRun3Style_particleLevel.second = (TH1D*)fEfficiency->Get("selEffRun3Style_particleLevelNonprompt")->Clone("hSelEffRun3Style_particleLevelNonprompt");
    if (!dataContainer.hSelEffRun3Style_particleLevel.second) {
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
            // if (!binning.useEmmaYeatsBins || passEmmaCut(PowjetPt, PowhfPt)) { // ToDo: remove Emma correspondent bins only at detector level, after all corrections, before luminosity scaling
            //     dataContainer.hPowheg->Fill(PowjetPt, PowaxisDistance, PowhfPt);
            // }
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
        bool recoAcceptance = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut);
        bool recoLevelRange = recoAcceptance && recoJetPtRange && recoHfPtRange && recoDeltaRRange;

        // Apply non-prompt selection (i.e., only B → D0)
        //bool isReflection = (MCDhfMatchedFrom != MCDhfSelectedAs) ? true : false;
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isRealD0 = MCDhfmatch && isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs);

        // Only real non-prompt D0s with existing match counterpart on particle level and detector level
        if (!MCDhfprompt && MCDhfmatch && isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs)) {
            // Response matrix and kinematic efficiencies numerators
            if (genLevelRange && recoLevelRange && passBDTcut) {
                // Fetch prompt D0 run 3 style detector level selection efficiency value corresponding to the given pT,D from the histogram
                int bin = dataContainer.hSelEffRun3Style_detectorLevel.first->FindBin(MCDhfPt);
                double efficiency_prompt = dataContainer.hSelEffRun3Style_detectorLevel.first->GetBinContent(bin);
                // Obs.: previously a weigth of 1./efficiency_prompt was used, is that correct?
                // Fill 4D RooUnfoldResponse object (jet pT shape is influenced by D0 pT efficiency)
                dataContainer.response.Fill(MCDjetPt, MCDDeltaR, MCDhfPt, MCPjetPt, MCPDeltaR, MCPhfPt, 1./ efficiency_prompt);
                dataContainer.hResponseDeltaR->Fill(MCDDeltaR,MCPDeltaR, 1./ efficiency_prompt);
                dataContainer.hResponsePtJet->Fill(MCDjetPt,MCPjetPt, 1./ efficiency_prompt);
                dataContainer.hResponsePtHF->Fill(MCDhfPt,MCPhfPt, 1./ efficiency_prompt);

                // Fill kinematic efficiency numerator histograms
                dataContainer.hKineEffParticle[0]->Fill(MCPjetPt, MCPDeltaR, MCPhfPt);
                dataContainer.hKineEffDetector[0]->Fill(MCDjetPt, MCDDeltaR, MCDhfPt);
            }
            // Particle level kinematic efficiency denominator (full particle level range)
            if (genLevelRange) {
                dataContainer.hKineEffParticle[1]->Fill(MCPjetPt, MCPDeltaR, MCPhfPt);
            }
            // Detector level kinematic efficiency denominator (full detector level range)
            if (recoLevelRange && passBDTcut) {
                dataContainer.hKineEffDetector[1]->Fill(MCDjetPt, MCDDeltaR, MCDhfPt);
            }
        }
    }
    std::cout << "Response matrix and kinematic efficiencies histograms filled.\n";

}

std::vector<double> retrieveLuminosityFromFile(TFile* fPowhegNonPrompt, TFile* fLumi, FeedDownData& dataContainer) {
    
    // Define temporary histogram to access luminosity values from file
    TH1D* hLumi = nullptr;

    // ---------------------------------------------------------
    // 1. POWHEG MC Luminosity Calculation
    // ---------------------------------------------------------
    
    // Properly fetch as a TProfile
    TProfile* pLumi = dynamic_cast<TProfile*>(fPowhegNonPrompt->Get("fHistXsection"));
    if (!pLumi) {
        std::cerr << "Error: fHistXsection is missing or not a TProfile!" << std::endl;
        return {0.0, 0.0};
    }

    // 1. Get the average cross section (Mean Y = 0.4263 mb)
    double sigmaPowheg_mb = pLumi->GetMean(2); // 2 specifies the Y-axis mean
    
    // 2. Get the total number of generated events (Entries = 2e7)
    double nEventsPowheg = pLumi->GetEntries(); // it's not the same number of entries on the D0 TTree

    // 3. Calculate MC Luminosity in mb^-1
    // L = 20,000,000 / 0.4263 mb = ~4.69e7 mb^-1
    double powhegLuminosity_mb = nEventsPowheg / sigmaPowheg_mb; 
    
    // 4. Convert to ub^-1 for consistency with data (1 mb^-1 = 1000 ub^-1)
    double powhegLuminosity_ub = powhegLuminosity_mb / 1000.0;

    // ---------------------------------------------------------
    // 2. Data Luminosity Calculation (Run 3 O2 Framework style)
    // ---------------------------------------------------------
    
    // Bunch crossings counting
    TTree* tLumiBunchCrossing = dynamic_cast<TTree*>(fLumi->Get("DF_merged/O2bccount"));
    Int_t bcTVX = 0;
    long long bcTVXSum = 0; // Use long long to avoid integer overflow loops
    tLumiBunchCrossing->SetBranchAddress("fCountsWithTVX", &bcTVX);
    
    for (Long64_t iEntry = 0; iEntry < tLumiBunchCrossing->GetEntries(); iEntry++) {
        tLumiBunchCrossing->GetEntry(iEntry);
        bcTVXSum += bcTVX;
    }

    // Collisions counting
    TTree* tLumiCollision = dynamic_cast<TTree*>(fLumi->Get("DF_merged/O2collcount"));
    Int_t selection = 0, collTVX = 0;
    long long selectionSum = 0, collTVXSum = 0;
    tLumiCollision->SetBranchAddress("fCountsWithTVXAndZVertexAndSel8", &selection);
    tLumiCollision->SetBranchAddress("fCountsWithTVX", &collTVX);
    
    for (Long64_t iEntry = 0; iEntry < tLumiCollision->GetEntries(); iEntry++) {
        tLumiCollision->GetEntry(iEntry);
        selectionSum += selection;
        collTVXSum += collTVX;
    }

    // Interaction cross section for TVX trigger in pp @ 13.6 TeV (in microbarns)
    const double csTVX_ub = 0.0594e6;
    
    // Train selection efficiency correction
    double efficiencyFactor = 0.0;
    if (collTVXSum > 0) {
        efficiencyFactor = static_cast<double>(selectionSum) / static_cast<double>(collTVXSum);
    }

    // L = (N_BC / sigma_trigger) * efficiency
    double dataLuminosity_ub = (static_cast<double>(bcTVXSum) / csTVX_ub) * efficiencyFactor;


    // ---------------------------------------------------------
    // 3. Save and Return
    // ---------------------------------------------------------
    dataContainer.lumiData = dataLuminosity_ub;
    dataContainer.lumiMC = powhegLuminosity_ub;

    std::cout << "========================================" << std::endl;
    std::cout << "MC Luminosity:   " << powhegLuminosity_ub << " ub^-1" << std::endl;
    std::cout << "Data Luminosity: " << dataLuminosity_ub << " ub^-1" << std::endl;
    std::cout << "========================================" << std::endl;

    return {powhegLuminosity_ub, dataLuminosity_ub};
}

// Particle level kinematic efficiency
void removeOutsideData(FeedDownData& dataContainer) {

    // Calculate particle level kinematic efficiency histograms
    dataContainer.hKineEffParticle[2] = (TH3D*)dataContainer.hKineEffParticle[0]->Clone("hKineEffParticleEfficiency");
    dataContainer.hKineEffParticle[2]->Divide(dataContainer.hKineEffParticle[1]); // A = A / B:  A = A->Divide(B)

    // Copy POWHEG efficiency corrected 3D histogram hPowhegEffCorrected
    dataContainer.hPowhegSelKinEffCorrected = (TH3D*)dataContainer.hPowhegSelEffCorrected->Clone("hPowhegSelKinEffCorrected");
    dataContainer.hPowhegSelKinEffCorrected->SetTitle("POWHEG distribution with #varepsilon_{non-prompt, sel}^{particle} and #varepsilon_{kin}^{particle} corrections");

    // Apply kinematic efficiency correction
    dataContainer.hPowhegSelKinEffCorrected->Multiply(dataContainer.hKineEffParticle[2]);
}

// Detector level kinematic efficiency
void addOutsideData(FeedDownData& dataContainer) {

    // Calculate detector level kinematic efficiency histograms
    dataContainer.hKineEffDetector[2] = (TH3D*)dataContainer.hKineEffDetector[0]->Clone("hKineEffDetectorEfficiency");
    dataContainer.hKineEffDetector[2]->Divide(dataContainer.hKineEffDetector[1]); // A = A / B:  A = A->Divide(B)

    // Copy POWHEG folded 2D histogram hFolded
    dataContainer.hFoldedKinCorrected = (TH3D*)dataContainer.hFolded->Clone("hFoldedKinCorrected");
    dataContainer.hFoldedKinCorrected->SetTitle("Folded 2D histogram with #varepsilon_{kin}^{detector} correction;p_{T,jet}^{reco} (GeV/#it{c});#DeltaR^{reco};p_{T,D^{0}}^{reco}");

    // Apply kinematic efficiency correction
    dataContainer.hFoldedKinCorrected->Divide(dataContainer.hKineEffDetector[2]);

}

// Module for folding particle level data from POWHEG simulation
void smearGeneratorData(FeedDownData& dataContainer, TFile* fEfficiency) {

    //
    // 1st step: correct by non-prompt particle level selection efficiency
    //
    dataContainer.hPowhegSelEffCorrected = (TH3D*)dataContainer.hPowheg->Clone("hPowhegSelEffCorrected");
    dataContainer.hPowhegSelEffCorrected->SetTitle("POWHEG distribution with #varepsilon_{non-prompt, sel}^{particle} correction");
    for (int xBin = 1; xBin <= dataContainer.hPowhegSelEffCorrected->GetNbinsX(); xBin++) {
        for (int yBin = 1; yBin <= dataContainer.hPowhegSelEffCorrected->GetNbinsY(); yBin++) {
            for (int zBin = 1; zBin <= dataContainer.hPowhegSelEffCorrected->GetNbinsZ(); zBin++) { //Inner loop over z-axis bins: this is scaled axis (pT,D bins)

                // Bin center of the Z axis (i.e., pT,D)
                double ptDcenter = dataContainer.hPowhegSelEffCorrected->GetZaxis()->GetBinCenter(zBin);

                // Find corresponding bin in the efficiency histograms
                int effBinPrompt = dataContainer.hSelEffRun3Style_particleLevel.first->FindBin(ptDcenter);
                int effBinNonPrompt = dataContainer.hSelEffRun3Style_particleLevel.second->FindBin(ptDcenter);

                // Get efficiency values
                double effPrompt = dataContainer.hSelEffRun3Style_particleLevel.first->GetBinContent(effBinPrompt);
                double effNonPrompt = dataContainer.hSelEffRun3Style_particleLevel.second->GetBinContent(effBinNonPrompt);

                // Get 3D histogram content
                double binContent = dataContainer.hPowhegSelEffCorrected->GetBinContent(xBin, yBin, zBin);
                double binError = dataContainer.hPowhegSelEffCorrected->GetBinError(xBin, yBin, zBin);

                // Rescale the bin content by the ratio of efficiencies
                if (effPrompt > 0) {
                    //double efficiencyRatio = effNonPrompt / effPrompt;
                    dataContainer.hPowhegSelEffCorrected->SetBinContent(xBin, yBin, zBin, binContent * effNonPrompt);
                    // Rescale the bin error accordingly
                    dataContainer.hPowhegSelEffCorrected->SetBinError(xBin, yBin, zBin, binError * effNonPrompt);
                } else {
                    double ptDlow = dataContainer.hPowhegSelEffCorrected->GetZaxis()->GetBinLowEdge(zBin);
                    double ptDhigh = dataContainer.hPowhegSelEffCorrected->GetZaxis()->GetBinUpEdge(zBin);

                    std::cout << "Warning: Efficiency for pT,D bin " << zBin 
                              << " (range: " << ptDlow << " - " << ptDhigh << " GeV/c, center: " << ptDcenter 
                              << ") is zero or negative(effNonPrompt = " << effNonPrompt << "). Skipping scaling." << std::endl;
                }
            }
        }
    }

    //
    // 2nd step: calculate particle level kinematic efficiency and apply to the distribution (remove outside of response range data in POWHEG)
    //
    removeOutsideData(dataContainer);

    //
    // 3rd step: fold pT,jet x DeltaR x pT,HF distribution using response matrix of selected non-prompt D0 jets
    //
    dataContainer.hFolded = manualFolding(dataContainer.response, dataContainer.hPowhegSelKinEffCorrected, dataContainer.hMeasuredTemplate);
    delete dataContainer.hMeasuredTemplate;

    //
    // 4th step: calculate detector level kinematic efficiency and apply to the distribution (add outside of response range data in folded distribution)
    //
    addOutsideData(dataContainer);
    
    //
    // 5th step: scale by 1 / efficiency_prompt (one over the detector level prompt d0 selection efficiency)
    // Obs.: also remove entries that doesn't pass Emma's cuts
    dataContainer.hFinalPromptSelEffCorrected = (TH3D*)dataContainer.hFoldedKinCorrected->Clone("hFinalPromptSelEffCorrected");
    for (int xBin = 1; xBin <= dataContainer.hFinalPromptSelEffCorrected->GetNbinsX(); xBin++) {
        for (int yBin = 1; yBin <= dataContainer.hFinalPromptSelEffCorrected->GetNbinsY(); yBin++) {
            for (int zBin = 1; zBin <= dataContainer.hFinalPromptSelEffCorrected->GetNbinsZ(); zBin++) { //Inner loop over z-axis bins: this is scaled axis (pT,D bins)

                // Bin center of the X axis (i.e., jet pT)
                double jetPtCenter = dataContainer.hFinalPromptSelEffCorrected->GetXaxis()->GetBinCenter(xBin);

                // Bin center of the Z axis (i.e., pT,D)
                double ptDcenter = dataContainer.hFinalPromptSelEffCorrected->GetZaxis()->GetBinCenter(zBin);

                // Find corresponding bin in the efficiency histograms
                int effBinPrompt = dataContainer.hSelEffRun3Style_detectorLevel.first->FindBin(ptDcenter);
                int effBinNonPrompt = dataContainer.hSelEffRun3Style_detectorLevel.second->FindBin(ptDcenter);

                // Get efficiency values
                double effPrompt = dataContainer.hSelEffRun3Style_detectorLevel.first->GetBinContent(effBinPrompt);
                double effNonPrompt = dataContainer.hSelEffRun3Style_detectorLevel.second->GetBinContent(effBinNonPrompt);

                // Get 3D histogram content
                double binContent = dataContainer.hFinalPromptSelEffCorrected->GetBinContent(xBin, yBin, zBin);
                double binError = dataContainer.hFinalPromptSelEffCorrected->GetBinError(xBin, yBin, zBin);

                double weightEmmaCut = passEmmaCut(jetPtCenter, ptDcenter) ? 1.0 : 0.0; // ToDo: remove Emma correspondent bins only at detector level, after all corrections, before luminosity scaling
                binContent *= weightEmmaCut;
                binError *= weightEmmaCut;

                // Rescale the bin content by the ratio of efficiencies
                if (effPrompt > 0) {
                    dataContainer.hFinalPromptSelEffCorrected->SetBinContent(xBin, yBin, zBin, binContent / effPrompt);
                    // Rescale the bin error accordingly
                    dataContainer.hFinalPromptSelEffCorrected->SetBinError(xBin, yBin, zBin, binError / effPrompt);
                } else {
                    double ptDlow = dataContainer.hFinalPromptSelEffCorrected->GetZaxis()->GetBinLowEdge(zBin);
                    double ptDhigh = dataContainer.hFinalPromptSelEffCorrected->GetZaxis()->GetBinUpEdge(zBin);

                    std::cout << "Warning: Efficiency for pT,D bin " << zBin 
                              << " (range: " << ptDlow << " - " << ptDhigh << " GeV/c, center: " << ptDcenter 
                              << ") is zero or negative(effNonPrompt = " << effNonPrompt << "). Skipping scaling." << std::endl;
                }
            }
        }
    }

    //
    // 6th step: scale by luminosity and branching ratio
    //
    dataContainer.hFinalScaled = (TH3D*)dataContainer.hFinalPromptSelEffCorrected->Clone("hFinalScaled");
    dataContainer.hFinalScaled->Scale(dataContainer.branchingRatio * dataContainer.lumiData / dataContainer.lumiMC);

    //
    // 7th step: project the final 3D distribution to obtain the final 2D distribution in jet pT and Delta R, after all corrections and folding
    //
    dataContainer.hFinalScaled2D = (TH2D*)dataContainer.hFinalScaled->Project3D("yx"); // (x;y) = (pT,jet;DeltaR)
    dataContainer.hFinalScaled2D->SetName("hFinalScaled2D");
    dataContainer.hFinalScaled2D->SetTitle("Pure non-prompt fully corrected distribution (as found in data)");

    std::cout << "Generator data smeared.\n";
}

// Module to subtract non-prompt D0 jets from prompt efficiency corrected distribution
void feedDown(TFile* fEfficiency, FeedDownData& dataContainer, const double jetptMin, const double jetptMax) {
    
    dataContainer.hTotalDataBeforeFeedDown = (TH2D*)fEfficiency->Get("h2DEfficiencyCorrected")->Clone("hTotalDataBeforeFeedDown");
    dataContainer.hBFedDownData = (TH2D*)fEfficiency->Get("h2DEfficiencyCorrected")->Clone("hBFedDownData");
    dataContainer.hBFedDownData->SetTitle("Background subtracted, efficiency corrected, fed down data distribution");

    dataContainer.hBFedDownData->Add(dataContainer.hFinalScaled2D,-1);
    
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
    TCanvas* cPowheg = new TCanvas("cPowheg","POWHEG data", 1800, 1000);
    cPowheg->Divide(2,2);
    cPowheg->cd(1);
    dataContainer.hPowheg->Draw("colz");
    cPowheg->cd(2);
    dataContainer.hPowheg->Project3D("yx")->Draw("colz");
    cPowheg->cd(3);
    dataContainer.hPowheg->Project3D("yx");
    TH2D* hPowheg2D = (TH2D*)dataContainer.hPowheg->Project3D("yx");
    hPowheg2D->ProjectionY("hPowhegProjY")->Draw();

    //
    // Response matrix
    //
    TCanvas* cResponse = new TCanvas("cResponse","Response matrix", 1800, 1000);
    cResponse->cd();
    const TH2* hresponse2D = dataContainer.response.Hresponse();
    TH2D* hresponse2DClone = static_cast<TH2D*>(hresponse2D->Clone("hResponse2D"));
    hresponse2DClone->SetTitle("2D response matrix from 6D RooUnfoldResponse - non-prompt D^{0}'s;3D Reconstructed;3D Truth");
    hresponse2DClone->SetStats(0);
    hresponse2DClone->Draw("colz");
    gPad->SetLogz();
    TCanvas* cResponseDeltaR = new TCanvas("cResponseDeltaR","cResponseDeltaR",1800,1000);
    cResponseDeltaR->cd();
    dataContainer.hResponseDeltaR->SetStats(0);
    dataContainer.hResponseDeltaR->Draw("colz");
    gPad->SetLogz();
    TCanvas* cResponsePtJet = new TCanvas("cResponsePtJet","cResponsePtJet",1800,1000);
    cResponsePtJet->cd();
    dataContainer.hResponsePtJet->SetStats(0);
    dataContainer.hResponsePtJet->Draw("colz");
    gPad->SetLogz();
    TCanvas* cResponsePtHF = new TCanvas("cResponsePtHF","cResponsePtHF",1800,1000);
    cResponsePtHF->cd();
    dataContainer.hResponsePtHF->SetStats(0);
    dataContainer.hResponsePtHF->Draw("colz");
    gPad->SetLogz();

    //
    // Kinematic efficiency histograms
    //
    TCanvas* cKinEffParticle = new TCanvas("cKinEffParticle","cKinEffParticle", 1800,1000);
    cKinEffParticle->Divide(2,2);
    cKinEffParticle->cd(1);
    TH2D* hNumPart_xy = (TH2D*)dataContainer.hKineEffParticle[0]->Project3D("xy");
    TH2D* hDenPart_xy = (TH2D*)dataContainer.hKineEffParticle[1]->Project3D("xy");
    hNumPart_xy->Divide(hDenPart_xy); // now values will be between 0 and 1
    hNumPart_xy->SetTitle("Particle level kinematic efficiency (non-prompt D^{0}'s) xy projection");
    hNumPart_xy->Draw("text");
    cKinEffParticle->cd(2);
    TH2D* hNumPart_yz = (TH2D*)dataContainer.hKineEffParticle[0]->Project3D("yz");
    TH2D* hDenPart_yz = (TH2D*)dataContainer.hKineEffParticle[1]->Project3D("yz");
    hNumPart_yz->Divide(hDenPart_yz); // now values will be between 0 and 1
    hNumPart_yz->SetTitle("Particle level kinematic efficiency (non-prompt D^{0}'s) yz projection");
    hNumPart_yz->Draw("text");
    cKinEffParticle->cd(3);
    TH2D* hNumPart_zx = (TH2D*)dataContainer.hKineEffParticle[0]->Project3D("zx");
    TH2D* hDenPart_zx = (TH2D*)dataContainer.hKineEffParticle[1]->Project3D("zx");
    hNumPart_zx->Divide(hDenPart_zx); // now values will be between 0 and 1
    hNumPart_zx->SetTitle("Particle level kinematic efficiency (non-prompt D^{0}'s) zx projection");
    hNumPart_zx->Draw("text");

    TCanvas* cKinEffDetector = new TCanvas("cKinEffDetector","cKinEffDetector", 1800,1000);
    cKinEffDetector->Divide(2,2);
    cKinEffDetector->cd(1);
    TH2D* hNumDet_xy = (TH2D*)dataContainer.hKineEffDetector[0]->Project3D("xy");
    TH2D* hDenDet_xy = (TH2D*)dataContainer.hKineEffDetector[1]->Project3D("xy");
    hNumDet_xy->Divide(hDenDet_xy); // now values will be between 0 and 1
    hNumDet_xy->SetTitle("Detector level kinematic efficiency (non-prompt D^{0}'s) xy projection");
    hNumDet_xy->Draw("text");
    cKinEffDetector->cd(2);
    TH2D* hNumDet_yz = (TH2D*)dataContainer.hKineEffDetector[0]->Project3D("yz");
    TH2D* hDenDet_yz = (TH2D*)dataContainer.hKineEffDetector[1]->Project3D("yz");
    hNumDet_yz->Divide(hDenDet_yz); // now values will be between 0 and 1
    hNumDet_yz->SetTitle("Detector level kinematic efficiency (non-prompt D^{0}'s) yz projection");
    hNumDet_yz->Draw("text");
    cKinEffDetector->cd(3);
    TH2D* hNumDet_zx = (TH2D*)dataContainer.hKineEffDetector[0]->Project3D("zx");
    TH2D* hDenDet_zx = (TH2D*)dataContainer.hKineEffDetector[1]->Project3D("zx");
    hNumDet_zx->Divide(hDenDet_zx); // now values will be between 0 and 1
    hNumDet_zx->SetTitle("Detector level kinematic efficiency (non-prompt D^{0}'s) zx projection");
    hNumDet_zx->Draw("text");

    //
    // Folded data
    //
    TCanvas* cFolded = new TCanvas("cFolded","Folded data", 1800, 1000);
    cFolded->Divide(2,2);
    cFolded->cd(1);
    dataContainer.hFolded->Project3D("yx")->SetStats(0);
    dataContainer.hFolded->Project3D("yx")->Draw("colz");
    cFolded->cd(2);
    dataContainer.hFoldedKinCorrected->Project3D("yx")->SetStats(0);
    dataContainer.hFoldedKinCorrected->Project3D("yx")->Draw("colz");
    cFolded->cd(3);
    dataContainer.hFinalPromptSelEffCorrected->Project3D("yx")->SetStats(0);
    dataContainer.hFinalPromptSelEffCorrected->Project3D("yx")->Draw("colz");
    cFolded->cd(4);
    dataContainer.hFinalScaled2D->SetStats(0);
    dataContainer.hFinalScaled2D->Draw("colz");
    gPad->SetLogz();

    TCanvas* cFedDownData = new TCanvas("cFedDownData","Background subtracted, efficiency corrected, fed down data distribution", 1800, 1000);
    cFedDownData->Divide(2,2);
    cFedDownData->cd(1);
    dataContainer.hBFedDownData->Draw("colz");
    gPad->SetLogz();
    cFedDownData->cd(2);
    dataContainer.hBFedDownData->ProjectionY("hBFedDownDataProjY")->Draw("colz");
    
    std::vector<TCanvas*> cNonpromptFraction(binning.ptjetBinEdges_detector.size() - 1);
    for (size_t iJetptBin = 0; iJetptBin < binning.ptjetBinEdges_detector.size() - 1; iJetptBin++) {
        // Get non-prompt projection
        TH1D* hNonpromptJetptBin = (TH1D*) dataContainer.hFinalScaled2D->ProjectionY(Form("hNonpromptJetptBin_%zu",iJetptBin), iJetptBin + 1, iJetptBin + 1);
        hNonpromptJetptBin->SetTitle(Form("Non-prompt fraction for jet pT bin %.0f - %.0f GeV/c", binning.ptjetBinEdges_detector[iJetptBin], binning.ptjetBinEdges_detector[iJetptBin + 1]));
        hNonpromptJetptBin->SetLineColor(kRed+2);
        hNonpromptJetptBin->SetMarkerStyle(kOpenCircle);
        hNonpromptJetptBin->SetMarkerColor(kRed+2);

        // Get total projection
        TH1D* hTotalJetptBin = (TH1D*) dataContainer.hTotalDataBeforeFeedDown->ProjectionY(Form("hTotalJetptBin_%zu",iJetptBin), iJetptBin + 1, iJetptBin + 1);
        hTotalJetptBin->SetTitle(Form("Total data for jet pT bin %.0f - %.0f GeV/c", binning.ptjetBinEdges_detector[iJetptBin], binning.ptjetBinEdges_detector[iJetptBin + 1]));
        hTotalJetptBin->SetLineColor(kBlue+2);
        hTotalJetptBin->SetMarkerStyle(kOpenCircle);
        hTotalJetptBin->SetMarkerColor(kBlue+2);
        hTotalJetptBin->GetYaxis()->SetRangeUser(std::min(hTotalJetptBin->GetMinimum(), hNonpromptJetptBin->GetMinimum()) * 0.8, std::max(hTotalJetptBin->GetMaximum(), hNonpromptJetptBin->GetMaximum()) * 1.2);

        TLegend* legend = new TLegend(0.6, 0.7, 0.8, 0.8);
        legend->AddEntry(hNonpromptJetptBin, "Non-prompt fraction", "lpe");
        legend->AddEntry(hTotalJetptBin, "Total data", "lpe");
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);

        cNonpromptFraction[iJetptBin] = new TCanvas(Form("cNonpromptFraction_%zu", iJetptBin), Form("Non-prompt fraction for jet pT bin %.0f - %.0f GeV/c", binning.ptjetBinEdges_detector[iJetptBin], binning.ptjetBinEdges_detector[iJetptBin + 1]), 1800, 1000);
        cNonpromptFraction[iJetptBin]->cd();
        gPad->SetLogy();
        hTotalJetptBin->Draw();
        hNonpromptJetptBin->Draw("same");
        latex->DrawLatex(0.2, 0.8, Form("Jet pT bin: %.0f - %.0f GeV/c", binning.ptjetBinEdges_detector[iJetptBin], binning.ptjetBinEdges_detector[iJetptBin + 1]));
        legend->Draw();
    }
    // Full jet pT range non-prompt fraction
    cNonpromptFraction.push_back(new TCanvas("cNonpromptFractionFull", "Non-prompt fraction for full jet pT range", 1800, 1000));
    cNonpromptFraction.back()->cd();
    TH1D* hNonpromptFull = (TH1D*) dataContainer.hFinalScaled2D->ProjectionY("hNonpromptFull");
    hNonpromptFull->SetTitle("Non-prompt fraction for full jet pT range"); 
    hNonpromptFull->SetLineColor(kRed+2);
    hNonpromptFull->SetMarkerStyle(kOpenCircle);
    hNonpromptFull->SetMarkerColor(kRed+2);
    TH1D* hTotalFull = (TH1D*) dataContainer.hTotalDataBeforeFeedDown->ProjectionY("hTotalFull");
    hTotalFull->SetTitle("Total data for full jet pT range");
    hTotalFull->SetLineColor(kBlue+2);
    hTotalFull->SetMarkerStyle(kOpenCircle);
    hTotalFull->SetMarkerColor(kBlue+2);
    hTotalFull->GetYaxis()->SetRangeUser(std::min(hTotalFull->GetMinimum(), hNonpromptFull->GetMinimum()) * 0.8, std::max(hTotalFull->GetMaximum(), hNonpromptFull->GetMaximum()) * 1.2);
    TLegend* legendFull = new TLegend(0.6, 0.7, 0.8, 0.8);
    legendFull->AddEntry(hNonpromptFull, "Non-prompt fraction", "lpe");
    legendFull->AddEntry(hTotalFull, "Total data", "lpe");
    legendFull->SetBorderSize(0);
    legendFull->SetFillStyle(0);
    gPad->SetLogy();
    hTotalFull->Draw();
    hNonpromptFull->Draw("same");
    latex->DrawLatex(0.2, 0.8, Form("Jet pT bin: %.0f - %.0f GeV/c", binning.ptjetBinEdges_detector[0], binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1]));
    legendFull->Draw();
    
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

    //
    // Storing in a single pdf file
    //
    cPowheg->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cResponse->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponseDeltaR->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponsePtJet->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponsePtHF->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cKinEffParticle->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cKinEffDetector->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFolded->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cFedDownData->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    for (size_t iJetptBin = 0; iJetptBin < cNonpromptFraction.size() - 1; iJetptBin++) {
        cNonpromptFraction[iJetptBin]->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    }
    cNonpromptFraction.back()->Print(imagePath + Form("feeddown_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

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
    TFile* fLumi = new TFile("../" + binning.inputDATA.second + "/AO2D_mergedDFs.root","read"); // previously: "AO2D_mergedDFs"
    if (!fLumi || fLumi->IsZombie()) {
        std::cerr << "Error: Unable to open luminosity data ROOT file." << std::endl;
    }

    // 1 - Create histograms
    FeedDownData dataContainer = createHistograms(binning);

    // 2 - Fill data
    fillHistograms(fEfficiency, fPowhegNonPrompt, fSimulatedMCMatched, dataContainer, jetptMin, jetptMax, binning);

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
