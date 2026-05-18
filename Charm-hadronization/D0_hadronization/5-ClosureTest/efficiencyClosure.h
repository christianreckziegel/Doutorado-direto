/*
 * Macro for performing selection efficiency estimation and correction to second closure test
 * 
 * 
 * 
 * Author: Christian Reckziegel
**/

using namespace std;

struct EfficiencyData {
    //
    // Run 2 style efficiency
    //
    std::vector<TH1D*> hMcpPt;
    std::vector<TH1D*> hMcdPt;
    std::vector<TH1D*> hSelEff_run2style;
    TH3D* hMcpPt_vs_ptJet_vs_deltaR;
    TH3D* hMcdPt_vs_ptJet_vs_deltaR;
    TH3D* hSelEff_run2style3d;
    std::vector<TH1D*> hSelEff_run2style_per_jetpt;
    std::vector<TH1D*> hSelEff_run2style_per_deltaR;
    
    //
    // Run 3 style efficiency
    //
    
    // Yield 2D histograms: pT,jet vs pT,D: first = prompt D0 data distribution, second = non-prompt D0 distribution
    std::pair<TH2D*, TH2D*> hYieldTruth; // all particle level entries (denominator): will be kinematically corrected and folded
    std::pair<TH2D*, TH2D*> hYieldMeasured; // detector level entries (numerator): went over smearing effects and passed the selection cuts

    // Response matrices
    std::pair<RooUnfoldResponse, RooUnfoldResponse> response; // first = prompt D0s, second = non-prompt D0s
    std::pair<std::vector<TH2D*>, std::vector<TH2D*>> responseProjections; // response projections matrix: first = prompt D^{0}, second = non-prompt D^{0}

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
    std::pair<TH2D*, TH2D*> hNumeratorRun3Style; // detector level entries: went over smearing effects and passed the selection cuts
    std::pair<TH2D*, TH2D*> hDenominatorRun3Style; // all particle level entries: will be kinematically corrected and folded

    // Denominator original and corrected pT histograms: projection of the 2D histograms
    std::pair<TH1D*, TH1D*> hHfPtYieldTruth;
    std::pair<TH1D*, TH1D*> hHfPtYieldTruthCorrected;

    //
    // Final efficiency histograms (numerator / denominator): first = prompt D0s, second = non-prompt D0s
    //
    std::pair<TH1D*, TH1D*> hSelectionEfficiency; // efficiency = numerator / denominator
    std::vector<TH1D*> hSelEff_run3style = {nullptr, nullptr};

    // Investigation histograms
    TH1D* hBDTBackgroundScore;

    // Efficiency corrected data (after background subtraction and now efficiency)
    std::pair<TH3D*, TH2D*> hEfficiencyDataBeforeCorrection;
    //std::pair<TH3D*, TH2D*> hEfficiencyDataBeforeCorrectionWithKinEff; // ToDo
    std::pair<TH3D*, TH2D*> hEfficiencyCorrectedData;
    //std::pair<TH3D*, TH2D*> hEfficiencyCorrectedDataWithKinEff; // ToDo
};

// Estimate the selection efficiency from detector level MC data (run 2 style)
EfficiencyData run2StyleEfficiency(TFile* fClosureInput, const BinningStruct& binning) { // modified testing version
    
    EfficiencyData histStruct;

    // The output object of efficiencies: inclusive, prompt only, non-prompt only
    std::vector<TH1D*> hSelEff_run2style;

    // 1 ----- Create pT,D histograms for efficiency calculation (3 cases: inclusive = 0, prompt only = 1, non-prompt only = 2)
    int histCaseNum = 3;
    for (size_t i = 0; i < histCaseNum; ++i) {
        histStruct.hMcpPt.push_back(new TH1D(Form("mcp_pt_%zu",i), ";#it{p}_{T,D}^{truth};dN^{truth}", binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data()));
        histStruct.hMcpPt[i]->SetMarkerColor(kBlack);
        histStruct.hMcpPt[i]->SetLineColor(kBlack);
        histStruct.hMcpPt[i]->SetMarkerStyle(kOpenCircle);
        histStruct.hMcpPt[i]->Sumw2();
        histStruct.hMcpPt[i]->SetStats(0);
        histStruct.hMcdPt.push_back(new TH1D(Form("mcd_pt_%zu",i), ";#it{p}_{T,D}^{reco};dN^{reco}", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data()));
        histStruct.hMcdPt[i]->SetMarkerColor(kBlue);
        histStruct.hMcdPt[i]->SetLineColor(kBlue);
        histStruct.hMcdPt[i]->SetMarkerStyle(kFullCircle);
        histStruct.hMcdPt[i]->Sumw2();
        histStruct.hMcdPt[i]->SetStats(0);
    }
    // 3D run 2 style efficiency pT,HF vs DeltaR (check deltaR dependence)
    histStruct.hMcpPt_vs_ptJet_vs_deltaR = new TH3D("hMcpPt_vs_ptJet_vs_deltaR", "Run 2 style prompt 3D efficiency;#it{p}_{T,D}^{truth};#it{p}_{T,jet}^{truth};#DeltaR^{truth}",binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data(), binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    histStruct.hMcpPt_vs_ptJet_vs_deltaR->Sumw2();
    histStruct.hMcpPt_vs_ptJet_vs_deltaR->SetStats(0);
    histStruct.hMcdPt_vs_ptJet_vs_deltaR = new TH3D("hMcdPt_vs_ptJet_vs_deltaR", "Run 2 style prompt 3D efficiency;#it{p}_{T,D}^{reco};#it{p}_{T,jet}^{reco};#DeltaR^{reco}",binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data(), binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());
    histStruct.hMcdPt_vs_ptJet_vs_deltaR->Sumw2();
    histStruct.hMcdPt_vs_ptJet_vs_deltaR->SetStats(0);

    // 2 ------ Fill histograms from TFile data
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
    
    // defining variables for accessing PARTICLE level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    int MCPjetNConst;
    // defining variables for accessing DETECTOR level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;
    int MCDjetNConst, MCDhfMatchedFrom, MCDhfSelectedAs;
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
    tree->SetBranchAddress("fHfMatchedFrom",&MCDhfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&MCDhfSelectedAs);
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // Generator level selection cuts
        double MCPdeltaR = MCPaxisDistance;
        // Generator level selection cuts
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax); //, remove upper bound for particle level
        bool genHfPtRange = ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)); // remove entirely for particle level
        bool genDeltaRRange = ((MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < MCPDeltaRcut));
        bool genAcceptance = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut);
        bool genLevelRange = genJetPtRange && genAcceptance && genDeltaRRange && genHfPtRange; // --> this is new!
        //bool genLevelRange = genJetPtRange && genAcceptance && genHfPtRange;

        // Reconstruction level selection cuts
        double MCDdeltaR = MCDaxisDistance;
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb);
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        bool recoDeltaRRange = ((MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < MCDDeltaRcut));
        bool recoAcceptance = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut);
        bool recoLevelRange = recoJetPtRange && recoHfPtRange && recoDeltaRRange && recoAcceptance;
        //bool recoLevelRange = recoJetPtRange && recoHfPtRange && recoAcceptance;
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isRealD0 = isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs);
        
        if (genLevelRange) {

            // Matched entries must be of real D0s, not reflections or combinatorial background
            if (MCDhfmatch && isRealD0) {
                // Fill particle level entry
                histStruct.hMcpPt[0]->Fill(MCPhfPt);// Fill inclusive histogram
                if (MCPhfprompt) {
                    histStruct.hMcpPt[1]->Fill(MCPhfPt);// fill prompt efficiency histogram
                    histStruct.hMcpPt_vs_ptJet_vs_deltaR->Fill(MCPhfPt, MCPjetPt, MCPdeltaR);
                } else{
                    histStruct.hMcpPt[2]->Fill(MCPhfPt);// fill non-prompt efficiency histogram
                }
                // Fill detector level entry
                if (recoLevelRange && passBDTcut) {
                    // Fill inclusive histogram
                    histStruct.hMcdPt[0]->Fill(MCDhfPt);
                    // fill prompt efficiency histogram
                    if (MCDhfprompt) {
                        histStruct.hMcdPt[1]->Fill(MCDhfPt);
                        histStruct.hMcdPt_vs_ptJet_vs_deltaR->Fill(MCDhfPt, MCDjetPt, MCDdeltaR);
                    } else{
                        // fill non-prompt efficiency histogram
                        histStruct.hMcdPt[2]->Fill(MCDhfPt);
                    }
                } else if (recoLevelRange && !passBDTcut) {
                    // histStruct.hBDTBackgroundScore->Fill(MCDhfMlScore0);
                }
            } else {// fill particle level even if not matched
                histStruct.hMcpPt[0]->Fill(MCPhfPt);// Fill inclusive histogram
                if (MCPhfprompt) {
                    histStruct.hMcpPt[1]->Fill(MCPhfPt);// fill prompt efficiency histogram
                    histStruct.hMcpPt_vs_ptJet_vs_deltaR->Fill(MCPhfPt, MCPjetPt, MCPdeltaR);
                } else{
                    histStruct.hMcpPt[2]->Fill(MCPhfPt);// fill non-prompt efficiency histogram
                }
            }   
        }

    }

    // 3 ----- Calculate efficiency distributions
    const char* efficiencyNames[] = {"inclusive", "prompt", "nonprompt"};
    // Loop through histogram cases: inclusive = 0, prompt only = 1, non-prompt only = 2
    //TCanvas* cPtEffRun2 = new TCanvas("cPtEffRun2", "pT,D distributions for efficiency run 2 style calculation",,1920,1080);
    //cPtEffRun2->Divide(2,2);
    for (int iEff = 0; iEff < histCaseNum; iEff++) {
        // Obtain MC pT distributions
        TH1D* mcpHist = histStruct.hMcpPt[iEff];
        TH1D* mcdHist = histStruct.hMcdPt[iEff];

        // Otain efficiency by dividing tempHist by parameter -> mcd/mcp
        TH1D* efficiencyHist = static_cast<TH1D*>(mcdHist->Clone(Form("efficiency_%s", efficiencyNames[iEff]))); // static needed to explicitly cast the cloned histogram to TH1D. ensures that the function returns a TH1D*, matching the return type.
        efficiencyHist->Divide(mcpHist);
        efficiencyHist->SetMarkerStyle(kFullCircle);
        efficiencyHist->SetMarkerColor(30+iEff*10); // 30 = not so bright green
        efficiencyHist->SetLineColor(30+iEff*10);
        efficiencyHist->SetTitle(";#it{p}_{T,D}^{truth};Efficiency#times Acceptance");
        efficiencyHist->SetStats(0);
        
        // Store efficiency distribution in struct
        histStruct.hSelEff_run2style.push_back(efficiencyHist);

    }

    // Dependence on pT,jet
    // Project 3D hMcpPt and hMcdPt in pT,HF vs pT,jet over full DeltaR space
    TH2D* hMcpPt2d_ptjet = (TH2D*) histStruct.hMcpPt_vs_ptJet_vs_deltaR->Project3D("xy"); // pT,HF vs pT,jet
    TH2D* hMcdPt2d_ptjet = (TH2D*) histStruct.hMcdPt_vs_ptJet_vs_deltaR->Project3D("xy");
    for (size_t iJetPt = 0; iJetPt < binning.ptjetBinEdges_particle.size() - 1; iJetPt++) {
        // Project the 2D into a 1D pT,HF distribution on the corresponding pT,jet bin
        int bin = iJetPt+1;
        TH1D* hMcpPt1d = hMcpPt2d_ptjet->ProjectionY(Form("hMcpPt1d_ptjet_%d", bin), bin, bin);
        TH1D* hMcdPt1d = hMcdPt2d_ptjet->ProjectionY(Form("hMcdPt1d_ptjet_%d", bin), bin, bin);

        // Calculate pT,HF dependant efficiency for each pT,jet interval
        histStruct.hSelEff_run2style_per_jetpt.emplace_back((TH1D*)hMcdPt1d->Clone(Form("hRun2style_per_jetpt_%d",bin)));
        histStruct.hSelEff_run2style_per_jetpt.back()->Divide(hMcpPt1d);
    }
    
    // Dependence on DeltaR
    // Project 3D hMcpPt and hMcdPt in pT,HF vs pT,jet over full DeltaR space
    TH2D* hMcpPt2d_deltaR = (TH2D*) histStruct.hMcpPt_vs_ptJet_vs_deltaR->Project3D("xz"); // pT,HF vs pT,jet
    TH2D* hMcdPt2d_deltaR = (TH2D*) histStruct.hMcdPt_vs_ptJet_vs_deltaR->Project3D("xz");
    for (size_t iDeltaR = 0; iDeltaR < binning.deltaRBinEdges_particle.size() - 1; iDeltaR++) {
        // Project the 2D into a 1D pT,HF distribution on the corresponding pT,jet bin
        int bin = iDeltaR+1;
        TH1D* hMcpPt1d = hMcpPt2d_deltaR->ProjectionY(Form("hMcpPt1d_deltaR_%d", bin), bin, bin);
        TH1D* hMcdPt1d = hMcdPt2d_deltaR->ProjectionY(Form("hMcdPt1d_deltaR_%d", bin), bin, bin);

        // Calculate pT,HF dependant efficiency for each pT,jet interval
        histStruct.hSelEff_run2style_per_deltaR.emplace_back((TH1D*)hMcdPt1d->Clone(Form("hRun2style_per_deltaR_%d",bin)));
        histStruct.hSelEff_run2style_per_deltaR.back()->Divide(hMcpPt1d);
    }

    return histStruct;
}

void foldingObjects(TFile* fClosureInput, EfficiencyData& histStruct, const BinningStruct& binning) { // modified testing version
    // 1 ----- Define histograms
    // Create 2D histograms for prompt D0s: raw and folded
    histStruct.hYieldTruth.first = new TH2D("hYieldTruthPrompt", "Particle level data prompt yield distribution;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    histStruct.hYieldTruth.first->Sumw2();
    histStruct.hYieldMeasured.first = new TH2D("hYieldMeasuredPrompt", "Detector level data prompt yield distribution;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{det} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    // Create 2D histograms for non-prompt D0s: raw and folded
    histStruct.hYieldTruth.second = new TH2D("hYieldTruthNonPrompt", "Particle level data non-prompt yield distribution;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    histStruct.hYieldTruth.second->Sumw2();
    histStruct.hYieldMeasured.second = new TH2D("hYieldMeasuredNonPrompt", "Detector level data non-prompt yield distribution;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{det} (GeV/#it{c})",
                                                binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
                                                binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());

    // Create RooUnfoldResponse objects for prompt and non-prompt D0s
    histStruct.response.first = RooUnfoldResponse(histStruct.hYieldMeasured.first, histStruct.hYieldTruth.first); // prompt D0s
    histStruct.response.second = RooUnfoldResponse(histStruct.hYieldMeasured.second, histStruct.hYieldTruth.second); // non-prompt D0s
    // Create projections of response matrix object for prompt and non-prompt D^{0}
    histStruct.responseProjections.first.push_back(new TH2D("responseProjectionsPtJetPrompt", "Prompt D^{0}'s reponse matrix #it{p}_{T, ch. jet} projection;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, ch. jet}^{gen} (GeV/#it{c})", binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data()));
    histStruct.responseProjections.first.push_back(new TH2D("responseProjectionsPtHFPrompt", "Prompt D^{0}'s reponse matrix p_{T,D^{0}} projection;#it{p}_{T, D^{0}}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data(), binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data()));
    histStruct.responseProjections.second.push_back(new TH2D("responseProjectionsPtJetNonPrompt", "Non-prompt D^{0}'s reponse matrix #it{p}_{T, ch. jet} projection;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, ch. jet}^{gen} (GeV/#it{c})", binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data()));
    histStruct.responseProjections.second.push_back(new TH2D("responseProjectionsPtHFNonPrompt", "Non-prompt D^{0}'s reponse matrix p_{T,D^{0}} projection;#it{p}_{T, D^{0}}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data(), binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data()));
    
    // Kinematic efficiency histograms
    histStruct.hKEffResponseParticle.first = new TH2D("hKEffResponseParticlePrompt", "Truth prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    histStruct.hKEffResponseParticle.second = new TH2D("hKEffResponseParticleNonPrompt", "Truth non-prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    histStruct.hKEffTruthTotalParticle.first = new TH2D("hKEffTruthTotalParticlePrompt", "Truth prompt D^{0}'s within total particle range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})",
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    histStruct.hKEffTruthTotalParticle.second = new TH2D("hKEffTruthTotalParticleNonPrompt", "Truth non-prompt D^{0}'s within total particle range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})",
                                                binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                                binning.ptHFBinEdges_particle.size() - 1, binning.ptHFBinEdges_particle.data());
    histStruct.hKEffResponseDetector.first = new TH2D("hKEffResponseDetectorPrompt", "Reconstructed prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    histStruct.hKEffResponseDetector.second = new TH2D("hKEffResponseDetectorNonPrompt", "Reconstructed non-prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    histStruct.hKEffRecoTotalDetector.first = new TH2D("hKEffRecoTotalDetectorPrompt", "Reconstructed prompt D^{0}'s within total detector range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}",
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    histStruct.hKEffRecoTotalDetector.second = new TH2D("hKEffRecoTotalDetectorNonPrompt", "Reconstructed non-prompt D^{0}'s within total detector range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}",
            binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), 
            binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
    
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
    TTree* tree = (TTree*)fClosureInput->Get("CorrectionTree");
    // Check for correct access
    if (!tree) {
        cout << "Error opening matched closure tree.\n";
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
    int MCDhfMatchedFrom, MCDhfSelectedAs, MCDjetNConst;
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
    MCPhfMass = 1.86484; // D0 rest mass in GeV/c^2 (hard-coded)
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

    TH1D* hEffWeight; // Histogram for efficiency weighting of response matrix
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double MCPDeltaR = MCPaxisDistance;
        double MCDDeltaR = MCDaxisDistance;
        // Generator level selection cuts
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax); //, remove upper bound for particle level
        bool genHfPtRange = ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)); // remove entirely for particle level
        bool genDeltaRRange = ((MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < MCPDeltaRcut));
        bool genAcceptance = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut);
        bool genLevelRange = genJetPtRange && genAcceptance && genDeltaRRange && genHfPtRange; // --> this is new!
        //bool genLevelRange = genJetPtRange && genAcceptance && genHfPtRange;

        // Reconstruction level selection cuts
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = (MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut);
        bool recoDeltaRRange = (MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < MCDDeltaRcut);
        bool recoAcceptance = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut);
        bool recoLevelRange = recoJetPtRange && recoHfPtRange && recoDeltaRRange && recoAcceptance;
        //bool recoLevelRange = recoJetPtRange && recoHfPtRange && recoAcceptance;
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isRealD0 = isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs);

        // 1 --- Selection efficiencies
        if (genLevelRange) {

            // Matched entries must be of real D0s, not reflections or combinatorial background
            if (MCDhfmatch && isRealD0) {
                // Fill particle level entry
                if (MCPhfprompt) {
                    histStruct.hYieldTruth.first->Fill(MCPjetPt, MCPhfPt); // prompt particle level denominator
                } else{
                    histStruct.hYieldTruth.second->Fill(MCPjetPt, MCPhfPt); // non-prompt particle level denonimator
                }
                // Fill detector level entry
                if (recoLevelRange && passBDTcut) {
                    if (MCDhfprompt) {
                        // prompt detector level numerator
                        histStruct.hYieldMeasured.first->Fill(MCDjetPt, MCDhfPt);
                    } else{
                        // non-prompt detector level numerator
                        histStruct.hYieldMeasured.second->Fill(MCDjetPt, MCDhfPt);
                    }
                } else if (recoLevelRange && !passBDTcut) {
                    //histStruct.hBDTBackgroundScore->Fill(MCDhfMlScore0);
                }
            } else {// fill particle level even if not matched
                if (MCPhfprompt) {
                    histStruct.hYieldTruth.first->Fill(MCPjetPt, MCPhfPt); // prompt particle level denominator
                } else{
                    histStruct.hYieldTruth.second->Fill(MCPjetPt, MCPhfPt); // non-prompt particle level denonimator
                }
            }   
        }

        // 2 --- Response matrix: fill histograms considering jet pT and detector acceptance (response range)
        if (MCDhfmatch && isRealD0) {
            if (genLevelRange && recoLevelRange) {

                // fill 2D yields histograms
                if (MCPhfprompt) {

                    // Get efficiency estimate to weight the response matrix
                    hEffWeight = histStruct.hSelEff_run2style[1];
                    // Get efficiency estimate to weight the response matrix: jet pT shape is influenced by D0 pT efficiency
                    int bin = hEffWeight->FindBin(MCDhfPt);
                    double estimatedEfficiency = hEffWeight->GetBinContent(bin);
                    if (estimatedEfficiency == 0) {
                        std::cout << "Warning: Prompt efficiency is zero for pT,HF = " << MCDhfPt << " GeV/c with bin " << bin << " of efficiency_prompt run 2 histogram. Setting it to 1. How to properly deal with these entries?" << std::endl;
                        estimatedEfficiency = 1; // Avoid division by zero
                    }

                    // Fill 4D RooUnfoldResponse object
                    histStruct.response.first.Fill(MCDjetPt, MCDhfPt, MCPjetPt, MCPhfPt, 1 / estimatedEfficiency); // jet pT shape is influenced by D0 pT efficiency
                    histStruct.responseProjections.first[0]->Fill(MCDjetPt, MCPjetPt,1 / estimatedEfficiency); // pT,jet projection, prompt D^{0}
                    histStruct.responseProjections.first[1]->Fill(MCDhfPt, MCPhfPt, 1 / estimatedEfficiency); // pT,HF projection, prompt D^{0}
                } else{

                    // Get efficiency estimate to weight the response matrix
                    hEffWeight = histStruct.hSelEff_run2style[2];
                    // Get efficiency estimate to weight the response matrix: jet pT shape is influenced by D0 pT efficiency
                    int bin = hEffWeight->FindBin(MCDhfPt);
                    double estimatedEfficiency = hEffWeight->GetBinContent(bin);
                    if (estimatedEfficiency == 0) {
                        //std::cout << "Warning: Prompt efficiency is zero for pT = " << MCDhfPt << " with bin " << bin << ". Setting it to 1." << std::endl;
                        estimatedEfficiency = 1; // Avoid division by zero
                    }

                    // Fill 4D RooUnfoldResponse object
                    histStruct.response.second.Fill(MCDjetPt, MCDhfPt, MCPjetPt, MCPhfPt, 1 / estimatedEfficiency); // jet pT shape is influenced by D0 pT efficiency
                    histStruct.responseProjections.second[0]->Fill(MCDjetPt, MCPjetPt,1 / estimatedEfficiency); // pT,jet projection, non-prompt D^{0}
                    histStruct.responseProjections.second[1]->Fill(MCDhfPt, MCPhfPt, 1 / estimatedEfficiency); // pT,HF projection, non-prompt D^{0}
                }
            }
        }

        // 3 --- Kinematic efficiencies
        if (MCDhfmatch && isRealD0) {
            // Particle level kinematic efficiency
            if (genLevelRange) {
                // prompt D0s
                if (MCPhfprompt) {
                    // fill prompt 2D yield: total particle range (denominator)
                    histStruct.hKEffTruthTotalParticle.first->Fill(MCPjetPt, MCPhfPt);
                    if (recoLevelRange) {
                        // fill prompt 2D yield: response range
                        histStruct.hKEffResponseParticle.first->Fill(MCPjetPt, MCPhfPt);
                    }
                } else {
                    // non-prompt D0s
                    histStruct.hKEffTruthTotalParticle.second->Fill(MCPjetPt, MCPhfPt);
                    if (recoLevelRange) {
                        // fill prompt 2D yield: particle level matched
                        histStruct.hKEffResponseParticle.second->Fill(MCPjetPt, MCPhfPt);
                    }
                } 
            }

            // Detector level kinematic efficiency
            if (recoLevelRange) {
                if (MCPhfprompt) {
                    // fill prompt 2D yield: total detector range
                    histStruct.hKEffRecoTotalDetector.first->Fill(MCDjetPt, MCDhfPt);
                    if (genLevelRange) {
                        // fill prompt 2D yield: response range
                        histStruct.hKEffResponseDetector.first->Fill(MCDjetPt, MCDhfPt);
                    }
                } else {
                    histStruct.hKEffRecoTotalDetector.second->Fill(MCDjetPt, MCDhfPt);
                    if (genLevelRange) {
                        // fill prompt 2D yield: detector level matched
                        histStruct.hKEffResponseDetector.second->Fill(MCDjetPt, MCDhfPt);
                    }
                }
            }
        }

    } // end of loop through TTree entries

    // 3 ----- Calculate kinematic efficiencies
    // Truth kinematic efficiency
    histStruct.hKEffResponseParticle_Over_TotalParticle.first = (TH2D*)histStruct.hKEffResponseParticle.first->Clone("hKEffTruthPrompt");
    histStruct.hKEffResponseParticle_Over_TotalParticle.second = (TH2D*)histStruct.hKEffResponseParticle.second->Clone("hKEffTruthNonPrompt");
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->Sumw2(); // necessary for correct error propagation
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->Sumw2();
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->Divide(histStruct.hKEffTruthTotalParticle.first);
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->Divide(histStruct.hKEffTruthTotalParticle.second);
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->SetTitle("#varepsilon_{kinematic} of particle level prompt D0s;#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->SetTitle("#varepsilon_{kinematic} of particle level non-prompt D0s;#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");
    // Reco kinematic efficiency
    histStruct.hKEffResponseDetector_Over_TotalDetector.first = (TH2D*)histStruct.hKEffResponseDetector.first->Clone("hKEffRecoPrompt");
    histStruct.hKEffResponseDetector_Over_TotalDetector.second = (TH2D*)histStruct.hKEffResponseDetector.second->Clone("hKEffRecoNonPrompt");
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->Sumw2();
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->Sumw2();
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->Divide(histStruct.hKEffRecoTotalDetector.first);
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->Divide(histStruct.hKEffRecoTotalDetector.second);
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->SetTitle("#varepsilon_{kinematic} of detector level prompt D0s;#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->SetTitle("#varepsilon_{kinematic} of detector level non-prompt D0s;#it{p}_{T, ch. jet};#it{p}_{T, D^{0}}");

}

void run3StyleEfficiency(TFile* fClosureInput, EfficiencyData& histStruct, const BinningStruct& binning) {
    
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    
    // 2 ----- Convert the pT,D denominator distribution into an appropriate detector level quantity
    // Correct the particle level input histogram: remove entries from input distribution whose matched entry is outside of detector level range encoding
    histStruct.hYieldTruthCorrected.first = (TH2D*)histStruct.hYieldTruth.first->Clone("hYieldTruthCorrectedPrompt");
    histStruct.hYieldTruthCorrected.second = (TH2D*)histStruct.hYieldTruth.second->Clone("hYieldTruthCorrectedNonPrompt");
    histStruct.hYieldTruthCorrected.first->Sumw2();
    histStruct.hYieldTruthCorrected.second->Sumw2();
    histStruct.hYieldTruthCorrected.first->Multiply(histStruct.hKEffResponseParticle_Over_TotalParticle.first);
    histStruct.hYieldTruthCorrected.second->Multiply(histStruct.hKEffResponseParticle_Over_TotalParticle.second);
    // Fold the corrected distribution into the response matrix
    histStruct.hYieldTruthCorrected.first = manualFolding(histStruct.response.first, histStruct.hYieldTruthCorrected.first, histStruct.hYieldMeasured.first);
    histStruct.hYieldTruthCorrected.second = manualFolding(histStruct.response.second, histStruct.hYieldTruthCorrected.second, histStruct.hYieldMeasured.second);
    // Correct the detector level output histogram: add entries to the output distribution which would be present at outside ranges of the response matrix at particle level
    histStruct.hYieldTruthCorrected.first->Divide(histStruct.hKEffResponseDetector_Over_TotalDetector.first);
    histStruct.hYieldTruthCorrected.second->Divide(histStruct.hKEffResponseDetector_Over_TotalDetector.second);

    // 3 ----- Calculate the run 3 style efficiency
    // Obtain pT projection histograms
    int minBin, maxBin;
    minBin = histStruct.hYieldTruth.first->GetXaxis()->FindBin(jetptMin);
    maxBin = histStruct.hYieldTruth.first->GetNbinsX() + 1;
    //maxBin = histStruct.hYieldTruth.first->GetXaxis()->FindBin(jetptMax - 1e-6); // // Tiny epsilon to stay within range: This includes the nearest bin center for jetptMax, but if jetptMax lies between two bins, it might give slightly unintuitive results
    histStruct.hHfPtYieldTruth.first = histStruct.hYieldTruth.first->ProjectionY("hHfPtYieldTruthPrompt", minBin, maxBin);
    histStruct.hHfPtYieldTruth.first->SetTitle("Denominator prompt D^{0} p_{T} distribution before correction; #it{p}_{T, D^{0}}^{gen}; Counts");
    minBin = histStruct.hYieldTruth.second->GetXaxis()->FindBin(jetptMin);
    maxBin = histStruct.hYieldTruth.second->GetNbinsX() + 1;
    //maxBin = histStruct.hYieldTruth.second->GetXaxis()->FindBin(jetptMax - 1e-6);
    histStruct.hHfPtYieldTruth.second = histStruct.hYieldTruth.second->ProjectionY("hHfPtYieldTruthNonPrompt", minBin, maxBin);
    histStruct.hHfPtYieldTruth.second->SetTitle("Denominator non-prompt D^{0} p_{T} distribution before correction; #it{p}_{T, D^{0}}^{gen}; Counts");
    minBin = histStruct.hYieldTruthCorrected.first->GetXaxis()->FindBin(jetptMin);
    //maxBin = histStruct.hYieldTruthCorrected.first->GetXaxis()->FindBin(jetptMax - 1e-6);
    maxBin = histStruct.hYieldTruthCorrected.first->GetNbinsX() + 1;
    histStruct.hHfPtYieldTruthCorrected.first = histStruct.hYieldTruthCorrected.first->ProjectionY("hHfPtYieldTruthCorrectedPrompt", minBin, maxBin);
    histStruct.hHfPtYieldTruthCorrected.first->SetTitle("Denominator prompt D^{0} p_{T} distribution after correction; #it{p}_{T, D^{0}}^{gen}; Counts");
    minBin = histStruct.hYieldTruthCorrected.second->GetXaxis()->FindBin(jetptMin);
    maxBin = histStruct.hYieldTruthCorrected.second->GetNbinsX() + 1;
    //maxBin = histStruct.hYieldTruthCorrected.second->GetXaxis()->FindBin(jetptMax - 1e-6);
    histStruct.hHfPtYieldTruthCorrected.second = histStruct.hYieldTruthCorrected.second->ProjectionY("hHfPtYieldTruthCorrectedNonPrompt", minBin, maxBin);
    histStruct.hHfPtYieldTruthCorrected.second->SetTitle("Denominator non-prompt D^{0} p_{T} distribution after correction; #it{p}_{T, D^{0}}^{gen}; Counts");
    // Calculate the efficiency distributions with the pT projections
    histStruct.hSelectionEfficiency.first = histStruct.hYieldMeasured.first->ProjectionY("hSelectionEfficiencyPrompt", minBin, maxBin);
    histStruct.hSelectionEfficiency.second = histStruct.hYieldMeasured.second->ProjectionY("hSelectionEfficiencyNonPrompt", minBin, maxBin);
    histStruct.hSelectionEfficiency.first->SetTitle("Efficiency prompt D^{0} p_{T} distribution; #it{p}_{T, D^{0}}^{reco}; Efficiency#times Acceptance");
    histStruct.hSelectionEfficiency.second->SetTitle("Efficiency non-prompt D^{0} p_{T} distribution; #it{p}_{T, D^{0}}^{reco}; Efficiency#times Acceptance");
    histStruct.hSelectionEfficiency.first->Sumw2();
    histStruct.hSelectionEfficiency.second->Sumw2();
    histStruct.hSelectionEfficiency.first->Divide(histStruct.hHfPtYieldTruthCorrected.first);
    histStruct.hSelectionEfficiency.second->Divide(histStruct.hHfPtYieldTruthCorrected.second);

}


TH2D* performEfficiencyCorrection(TH3D* hBackgroundSubtracted, EfficiencyData& histStruct) {
    cleanNaNs(hBackgroundSubtracted);
    // 1 - Clone background subtracted distribution
    histStruct.hEfficiencyDataBeforeCorrection.first = (TH3D*)hBackgroundSubtracted->Clone("hEfficiencyDataBeforeCorrection");
    histStruct.hEfficiencyCorrectedData.first = (TH3D*)hBackgroundSubtracted->Clone("h3DEfficiencyCorrected");
    histStruct.hEfficiencyCorrectedData.first->SetTitle("Background subtracted, efficiency corrected");
    int xBins = histStruct.hEfficiencyCorrectedData.first->GetXaxis()->GetNbins();
    int yBins = histStruct.hEfficiencyCorrectedData.first->GetYaxis()->GetNbins();
    int zBins = histStruct.hEfficiencyCorrectedData.first->GetZaxis()->GetNbins();
    for (int xBin = 1; xBin <= xBins; xBin++) {
        for (int yBin = 1; yBin <= yBins; yBin++) {
            for (int zBin = 1; zBin <= zBins; zBin++) { // pT,D bins will be corrected by 1/prompt_efficiency
                
                // Get current content and error
                double content = histStruct.hEfficiencyCorrectedData.first->GetBinContent(xBin, yBin, zBin);
                double error = histStruct.hEfficiencyCorrectedData.first->GetBinError(xBin, yBin, zBin);

                // Get the pT,D bin center from the Z axis of the 3D histogram
                double ptDcenter = histStruct.hEfficiencyCorrectedData.first->GetZaxis()->GetBinCenter(zBin);
                
                // Find corresponding bin in efficiency histogram
                int effBin = histStruct.hSelectionEfficiency.first->GetXaxis()->FindBin(ptDcenter);
                double eff = histStruct.hSelectionEfficiency.first->GetBinContent(effBin);

                // Avoid divide by zero or nonsense values
                if (eff > 0) { // eff > 0
                    double correction = 1. / eff;
                    histStruct.hEfficiencyCorrectedData.first->SetBinContent(xBin, yBin, zBin, content * correction);
                    histStruct.hEfficiencyCorrectedData.first->SetBinError(xBin, yBin, zBin, error * correction);
                } else {
                    // Optionally warn about zero efficiency
                    std::cerr << "Zero or invalid efficiency at pT,HF = " << ptDcenter << " (eff bin = " << effBin << ")" << std::endl;
                }
            }
        }
    }

    // 2 - Project the z axis (pT,D) in order to obtain a 2D distribution of pT,jet (x axis) vs DeltaR (y axis)
    TH2D* h2D = (TH2D*)histStruct.hEfficiencyDataBeforeCorrection.first->Project3D("yx");
    histStruct.hEfficiencyDataBeforeCorrection.second = (TH2D*)h2D->Clone("h2DEfficiencyDataBeforeCorrection");
    delete h2D;
    h2D = (TH2D*)histStruct.hEfficiencyCorrectedData.first->Project3D("yx");
    histStruct.hEfficiencyCorrectedData.second = (TH2D*)h2D->Clone("h2DEfficiencyCorrected");
    delete h2D;
    histStruct.hEfficiencyCorrectedData.second->SetTitle("Background subtracted, efficiency corrected");

    std::cout << "Efficiency (run 3 style) corrected." << std::endl;
    return histStruct.hEfficiencyCorrectedData.second;
}

void plotHistograms(const EfficiencyData& histStruct, const BinningStruct& binning) {

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    //
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    
    TCanvas* cRun2StyleEfficiency = new TCanvas("cRun2StyleEfficiency", "Run 2 style efficiency", 1920,1080);
    histStruct.hSelEff_run2style[1]->SetMarkerStyle(kOpenCircle);
    histStruct.hSelEff_run2style[1]->SetMarkerColor(30+1*10);
    histStruct.hSelEff_run2style[1]->SetLineColor(30+1*10);
    histStruct.hSelEff_run2style[2]->SetMarkerStyle(kOpenCircle);
    histStruct.hSelEff_run2style[2]->SetMarkerColor(30+2*10);
    histStruct.hSelEff_run2style[2]->SetLineColor(30+2*10);
    histStruct.hSelEff_run2style[1]->SetTitle("Run 2 style efficiency;");
    histStruct.hSelEff_run2style[1]->Draw();
    histStruct.hSelEff_run2style[2]->Draw("same");
    TLegend* legendRun2 = new TLegend(0.5, 0.4, 0.7, 0.5);
    legendRun2->AddEntry(histStruct.hSelEff_run2style[1], "Prompt D^{0}", "LP");
    legendRun2->AddEntry(histStruct.hSelEff_run2style[2], "Non-prompt D^{0}", "LP");
    legendRun2->Draw();

    TCanvas* cResponse = new TCanvas("cResponse", "Response matrix projections", 1920,1080);
    cResponse->Divide(2,2);
    cResponse->cd(1);
    histStruct.responseProjections.first[0]->Draw("colz");
    cResponse->cd(2);
    histStruct.responseProjections.second[0]->Draw("colz");
    cResponse->cd(3);
    histStruct.responseProjections.first[1]->Draw("colz");
    cResponse->cd(4);
    histStruct.responseProjections.second[1]->Draw("colz");

    TCanvas* cKinematicEfficiencies = new TCanvas("cKinematicEfficiencies", "Kinematic efficiencies", 1920,1080);
    cKinematicEfficiencies->Divide(2,2);
    cKinematicEfficiencies->cd(1);
    histStruct.hKEffResponseParticle_Over_TotalParticle.first->Draw("text");
    cKinematicEfficiencies->cd(2);
    histStruct.hKEffResponseParticle_Over_TotalParticle.second->Draw("text");
    cKinematicEfficiencies->cd(3);
    histStruct.hKEffResponseDetector_Over_TotalDetector.first->Draw("text");
    cKinematicEfficiencies->cd(4);
    histStruct.hKEffResponseDetector_Over_TotalDetector.second->Draw("text");

    TCanvas* cRun3StyleEfficiency = new TCanvas("cRun3StyleEfficiency", "Run 3 style efficiency", 1920,1080);
    histStruct.hSelectionEfficiency.first->SetMarkerStyle(kFullCircle);
    histStruct.hSelectionEfficiency.first->SetMarkerColor(30+1*10);
    histStruct.hSelectionEfficiency.first->SetLineColor(30+1*10);
    histStruct.hSelectionEfficiency.second->SetMarkerStyle(kFullCircle);
    histStruct.hSelectionEfficiency.second->SetMarkerColor(30+2*10);
    histStruct.hSelectionEfficiency.second->SetLineColor(30+2*10);
    histStruct.hSelectionEfficiency.first->SetTitle("Run 3 style efficiency;");
    histStruct.hSelectionEfficiency.first->Draw();
    histStruct.hSelectionEfficiency.second->Draw("same");
    TLegend* legendRun3 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legendRun3->AddEntry(histStruct.hSelectionEfficiency.first, "Prompt D^{0}", "LP");
    legendRun3->AddEntry(histStruct.hSelectionEfficiency.second, "Non-prompt D^{0}", "LP");

    TCanvas* cBeforeAfterScaling = new TCanvas("cBeforeAfterScaling","Before and after scaling",1920,1080);
    cBeforeAfterScaling->cd();
    TH1D* h1DEfficiencyDataBeforeCorrection = histStruct.hEfficiencyDataBeforeCorrection.second->ProjectionY("h1DEfficiencyDataBeforeCorrection");
    h1DEfficiencyDataBeforeCorrection->SetLineColor(kRed+2);
    TH1D* h1DEfficiencyCorrectedData = histStruct.hEfficiencyCorrectedData.second->ProjectionY("h1DEfficiencyCorrectedData");
    h1DEfficiencyCorrectedData->SetLineColor(kBlue+2);
    TLegend* legBefAft = new TLegend(0.6,0.6,0.8,0.75);
    legBefAft->AddEntry(h1DEfficiencyDataBeforeCorrection,"Before correction","lp");
    legBefAft->AddEntry(h1DEfficiencyCorrectedData,"After correction","lp");
    legBefAft->Draw();

    TCanvas* cDenominators = new TCanvas("cDenominators","Denominator distributions used for efficiency estimation",1920,1080);
    gStyle->SetOptStat(1111);
    cDenominators->Divide(3,3);
    cDenominators->cd(1);
    histStruct.hHfPtYieldTruth.first->SetTitle("Prompt denominator before folding");
    histStruct.hHfPtYieldTruth.first->Draw();
    cDenominators->cd(2);
    histStruct.hHfPtYieldTruthCorrected.first->SetTitle("Prompt denominator after folding");
    histStruct.hHfPtYieldTruthCorrected.first->Draw();
    cDenominators->cd(4);
    histStruct.hHfPtYieldTruth.second->SetTitle("Non-prompt denominator before folding");
    histStruct.hHfPtYieldTruth.second->Draw();
    cDenominators->cd(5);
    histStruct.hHfPtYieldTruthCorrected.second->SetTitle("Non-prompt denominator after folding");
    histStruct.hHfPtYieldTruthCorrected.second->Draw();
    cDenominators->cd(7);
    histStruct.hMcpPt[1]->SetTitle("Prompt denominator from run 2 style efficiency");
    histStruct.hMcpPt[1]->Draw();
    cDenominators->cd(8);
    histStruct.hMcdPt[1]->SetTitle("Prompt numerator from run 2 style efficiency");
    histStruct.hMcdPt[1]->Draw();
    cDenominators->cd(9);
    int minBin, maxBin;
    minBin = histStruct.hYieldTruth.first->GetXaxis()->FindBin(jetptMin);
    //maxBin = histStruct.hYieldTruth.first->GetXaxis()->FindBin(jetptMax - 1e-6);
    maxBin = histStruct.hYieldTruth.first->GetNbinsX() + 1;
    TH1D* hRun3StyleNumerator = histStruct.hYieldMeasured.first->ProjectionY("hRun3StyleNumerator", minBin, maxBin);
    hRun3StyleNumerator->SetTitle("Prompt numerator from run 3 style efficiency");
    hRun3StyleNumerator->Sumw2();
    hRun3StyleNumerator->Draw();
    cDenominators->cd(6);
    TH1D* hNumeratorRun2Run3Ratio = (TH1D*) histStruct.hMcdPt[1]->Clone("hNumeratorRun2Run3Ratio");
    hNumeratorRun2Run3Ratio->SetTitle("(run 2)/(run 3) numerators ratio both at detector level");
    hNumeratorRun2Run3Ratio->Divide(hRun3StyleNumerator);
    hNumeratorRun2Run3Ratio->Draw();
    cDenominators->cd(3);
    TH1D* hDenominatorRun2Run3Ratio = (TH1D*) histStruct.hMcpPt[1]->Clone("hDenominatorRun2Run3Ratio");
    hDenominatorRun2Run3Ratio->SetTitle("(run 2)/(run 3) denominators ratio both at particle level");
    hDenominatorRun2Run3Ratio->Divide(histStruct.hHfPtYieldTruth.first);
    hDenominatorRun2Run3Ratio->Draw();

    // Efficiency per jet pT interval
    TCanvas* cEfficiencyPerJetPt = new TCanvas("cEfficiencyPerJetPt","Efficiency per jet pT interval",1800,1000);
    cEfficiencyPerJetPt->cd();
    TLegend* leg_per_jetpt = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (int iJetPt = 0; iJetPt < histStruct.hSelEff_run2style_per_jetpt.size() - 1; iJetPt++) {
        //
        histStruct.hSelEff_run2style_per_jetpt[iJetPt]->SetLineColor(kBlack+iJetPt);
        histStruct.hSelEff_run2style_per_jetpt[iJetPt]->SetMarkerStyle(kFullCircle);
        double iJetptMin = binning.ptjetBinEdges_particle[iJetPt];
        double iJetptMax = binning.ptjetBinEdges_particle[iJetPt+1];
        histStruct.hSelEff_run2style_per_jetpt[iJetPt]->GetYaxis()->SetTitle("Efficiency#times Acceptance");
        histStruct.hSelEff_run2style_per_jetpt[iJetPt]->GetYaxis()->SetRangeUser(0.,1.0);
        histStruct.hSelEff_run2style_per_jetpt[iJetPt]->Draw(iJetPt == 0 ? "" : "same");
        leg_per_jetpt->AddEntry(histStruct.hSelEff_run2style_per_jetpt[iJetPt],Form("p_{T,jet}#in[%.0f-%.0f] GeV/c",iJetptMin,iJetptMax),"ple");
    }
    latex->DrawLatex(0.15, 0.85, Form("Projected in #DeltaR #in [%.0f, %.2f]", binning.deltaRBinEdges_particle[0], binning.deltaRBinEdges_particle[binning.deltaRBinEdges_particle.size()-1]));
    leg_per_jetpt->Draw();

    // Efficiency per deltaR interval
    TCanvas* cEfficiencyPerDeltaR = new TCanvas("cEfficiencyPerDeltaR","Efficiency per jet pT interval",1800,1000);
    cEfficiencyPerDeltaR->cd();
    TLegend* leg_per_deltaR = new TLegend(0.7, 0.6, 0.9, 0.9);
    for (int iDeltaR = 0; iDeltaR < histStruct.hSelEff_run2style_per_deltaR.size() - 1; iDeltaR++) {
        //
        histStruct.hSelEff_run2style_per_deltaR[iDeltaR]->SetLineColor(kBlack+iDeltaR);
        histStruct.hSelEff_run2style_per_deltaR[iDeltaR]->SetMarkerStyle(kFullCircle);
        double iDeltaRMin = binning.deltaRBinEdges_particle[iDeltaR];
        double iDeltaRMax = binning.deltaRBinEdges_particle[iDeltaR+1];
        histStruct.hSelEff_run2style_per_deltaR[iDeltaR]->GetYaxis()->SetTitle("Efficiency#times Acceptance");
        histStruct.hSelEff_run2style_per_deltaR[iDeltaR]->GetYaxis()->SetRangeUser(0.,1.0);
        histStruct.hSelEff_run2style_per_deltaR[iDeltaR]->Draw(iDeltaR == 0 ? "" : "same");
        leg_per_deltaR->AddEntry(histStruct.hSelEff_run2style_per_deltaR[iDeltaR],Form("#DeltaR#in[%.2f-%.2f]",iDeltaRMin,iDeltaRMax),"ple");
    }
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.0f, %.0f] GeV/c", jetptMin, jetptMax));
    leg_per_deltaR->Draw();

    //
    // Storing in a single pdf file
    //
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../Images/5-ClosureTest/Second/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    cRun2StyleEfficiency->Print(imagePath + Form("closureTest2_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cResponse->Print(imagePath + Form("closureTest2_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cKinematicEfficiencies->Print(imagePath + Form("closureTest2_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cDenominators->Print(imagePath + Form("closureTest2_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cRun3StyleEfficiency->Print(imagePath + Form("closureTest2_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cEfficiencyPerJetPt->Print(imagePath + Form("closureTest2_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cEfficiencyPerDeltaR->Print(imagePath + Form("closureTest2_efficiency_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
}

EfficiencyData EfficiencyClosure(TFile* fClosureInput, TH3D* hBackgroundSubtracted, const BinningStruct& binning) {

    std::cout << "Calculating efficiencies..." << std::endl;

    // Calculate run 2 style efficiency (to be used in run 3 style calculation): inclusive, prompt-only, non-prompt only
    EfficiencyData histStruct = run2StyleEfficiency(fClosureInput, binning);

    // Calculate particle and detector level kinematic efficiencies
    foldingObjects(fClosureInput, histStruct, binning);

    // Calculate run 3 style efficiency: 0 = prompt, 1 = non-prompt
    run3StyleEfficiency(fClosureInput, histStruct, binning);

    // Perform efficiency correction to background subtracted data
    TH2D* hEfficiencyCorrectedData = performEfficiencyCorrection(hBackgroundSubtracted, histStruct);
    
    // Plot selection/kinematic efficiencies and response matrix projections
    plotHistograms(histStruct, binning);

    std::cout << "Efficiencies calculated and correction applied." << std::endl;

    return histStruct;
}
