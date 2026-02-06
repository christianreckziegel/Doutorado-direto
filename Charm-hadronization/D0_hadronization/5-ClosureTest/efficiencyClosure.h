/*
 * Macro for performing selection efficiency estimation and correction to second closure test
 * 
 * 
 * 
 * Author: Christian Reckziegel
**/

using namespace std;

// Estimate the selection efficiency from detector level MC data (run 2 style)
std::vector<TH1D*> calculateSelectionEfficiencyRun2Matched(TFile* fClosureInputMatched, const BinningStruct& binningStruct, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    // The output object of efficiencies: inclusive, prompt only, non-prompt only
    std::vector<TH1D*> hSelEff_run2style;

    // 1 ----- Create pT,D histograms for efficiency calculation (3 cases: inclusive = 0, prompt only = 1, non-prompt only = 2)
    std::vector<TH1D*> hMcpPt;
    std::vector<TH1D*> hMcdPt;
    int histCaseNum = 3;
    for (size_t i = 0; i < histCaseNum; ++i) {
        hMcpPt.push_back(new TH1D(Form("mcp_pt_%zu",i), ";#it{p}_{T,D}^{truth};dN/d#it{p}_{T,D}^{truth}", binningStruct.ptDBinEdges_detector.size() - 1, binningStruct.ptDBinEdges_detector.data()));
        hMcpPt[i]->SetMarkerColor(kBlack);
        hMcpPt[i]->SetLineColor(kBlack);
        hMcpPt[i]->SetMarkerStyle(kOpenCircle);
        hMcpPt[i]->Sumw2();
        hMcpPt[i]->SetStats(0);
        hMcdPt.push_back(new TH1D(Form("mcd_pt_%zu",i), ";#it{p}_{T,D}^{reco};dN/d#it{p}_{T,D}^{reco}", binningStruct.ptDBinEdges_detector.size() - 1, binningStruct.ptDBinEdges_detector.data()));
        hMcdPt[i]->SetMarkerColor(kBlue);
        hMcdPt[i]->SetLineColor(kBlue);
        hMcdPt[i]->SetMarkerStyle(kFullCircle);
        hMcdPt[i]->Sumw2();
        hMcdPt[i]->SetStats(0);
    }

    // 2 ------ Fill histograms from TFile data
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on jet
    const double yCut = 0.8; // on D0
    double jetptMin = binningStruct.ptjetBinEdges_detector[0];
    double jetptMax = binningStruct.ptjetBinEdges_detector[binningStruct.ptjetBinEdges_detector.size() - 1];
    double MCPetaCut = etaCut;
    double MCPyCut = yCut;
    double MCDetaCut = etaCut;
    double MCDyCut = yCut;
    const double deltaRcut = binningStruct.deltaRBinEdges_detector[binningStruct.deltaRBinEdges_detector.size() - 1];
    double MCPDeltaRcut = deltaRcut;
    double MCDDeltaRcut = deltaRcut;
    const double MCPHfPtMincut = binningStruct.ptDBinEdges_particle[0]; // on particle level D0
    const double MCDHfPtMincut = binningStruct.ptDBinEdges_particle[0]; // on detector level D0
    const double MCPHfPtMaxcut = binningStruct.ptDBinEdges_particle[binningStruct.ptDBinEdges_particle.size() - 1]; // on particle level D0
    const double MCDHfPtMaxcut = binningStruct.ptDBinEdges_detector[binningStruct.ptDBinEdges_detector.size() - 1]; // on detector level D0
    TTree* tree = (TTree*)fClosureInputMatched->Get("CorrectionTree");
    if (!tree) {
        cout << "Error opening tree.\n";
    }

    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    //int MCPjetnconst
    int MCPhfmatch;
    float MCPjetnconst;
    bool MCPhfprompt;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    int MCDjetnconst;
    int MCDhfmatch;
    //float MCDjetnconst;
    bool MCDhfprompt;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;

    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetnconst);
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
    tree->SetBranchAddress("fJetNConst",&MCDjetnconst);
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
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));
        
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);
        bool genHfPtRange = (MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut);
        bool genDeltaRRange = (MCPDeltaR >= binningStruct.deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut);
        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && genJetPtRange && genHfPtRange && genDeltaRRange; // currently used
        bool recoJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);
        bool recoHfPtRange = (MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut);
        bool recoDeltaRRange = (MCPDeltaR >= binningStruct.deltaRBinEdges_detector[0]) && (MCPDeltaR < MCPDeltaRcut);
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && recoJetPtRange && recoHfPtRange && recoDeltaRRange; // currently used

        // Get the threshold for this pT range
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;

        // Fill particle level histograms
        if (genLevelRange) {
            
            // Fill inclusive histogram
            hMcpPt[0]->Fill(MCPhfPt);
            // fill prompt efficiency histogram
            if (MCPhfprompt) {
                hMcpPt[1]->Fill(MCPhfPt);
            } else{
                // fill non-prompt efficiency histogram
                hMcpPt[2]->Fill(MCPhfPt);
            }
            
        }

        // Fill detector level histograms
        if (recoLevelRange && passBDTcut) {
            
            // Fill inclusive histogram
            hMcdPt[0]->Fill(MCDhfPt);
            // fill prompt efficiency histogram
            if (MCDhfprompt) {
                hMcdPt[1]->Fill(MCDhfPt);
            } else{
                // fill non-prompt efficiency histogram
                hMcdPt[2]->Fill(MCDhfPt);
            }
            
        }

    }

    // 3 ----- Calculate efficiency distributions
    const char* efficiencyNames[] = {"inclusive", "prompt", "nonprompt"};
    // Loop through histogram cases: inclusive = 0, prompt only = 1, non-prompt only = 2
    for (int iEff = 0; iEff < histCaseNum; iEff++) {
        // Obtain MC pT distributions
        TH1D* mcdHist = hMcdPt[iEff];
        TH1D* mcpHist = hMcpPt[iEff];

        // Otain efficiency by dividing tempHist by parameter -> mcd/mcp
        TH1D* efficiencyHist = static_cast<TH1D*>(mcdHist->Clone(Form("efficiency_%s", efficiencyNames[iEff]))); // static needed to explicitly cast the cloned histogram to TH1D. ensures that the function returns a TH1D*, matching the return type.
        efficiencyHist->Divide(mcpHist);
        efficiencyHist->SetMarkerStyle(kFullCircle);
        efficiencyHist->SetMarkerColor(30+iEff*10); // 30 = not so bright green
        efficiencyHist->SetLineColor(30+iEff*10);
        efficiencyHist->SetTitle(";#it{p}_{T,D}^{truth};Efficiency#times Acceptance");
        efficiencyHist->SetStats(0);
        
        // Store efficiency distribution in struct
        hSelEff_run2style.push_back(efficiencyHist);

    }

    // Plot the efficiencies
    TCanvas* cEffRun2 = new TCanvas("cEffRun2","Run 2 style efficiencies", 1800, 1000);
    cEffRun2->cd();
    hSelEff_run2style[0]->Draw();
    hSelEff_run2style[1]->Draw("same");
    hSelEff_run2style[2]->Draw("same");
    TLegend* legendEff = new TLegend(0.65,0.55,0.8,0.7);
    legendEff->AddEntry(hSelEff_run2style[0],"Inclusive", "lpe");
    legendEff->AddEntry(hSelEff_run2style[1],"Prompt D^{0}", "lpe");
    legendEff->AddEntry(hSelEff_run2style[2],"Non-prompt D^{0}", "lpe");
    legendEff->Draw();
    double statBoxPos = gPad->GetUxmax();
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->SetTextAlign(31); // align right
    latex->DrawLatex(0.85,0.42,Form("|#eta_{jet}| < %.1f, %.1f < #it{p}_{T,jet} < %.1f GeV/c", etaCut, jetptMin, jetptMax));

    return hSelEff_run2style;
}

// Estimate the selection efficiency from detector level MC data (run 2 style)
std::vector<TH1D*> calculateSelectionEfficiencyRun2NonMatched(TFile* fClosureInputNonMatched, const BinningStruct& binningStruct, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    // The output object of efficiencies: inclusive, prompt only, non-prompt only
    std::vector<TH1D*> hSelEff_run2style;

    // 1 ----- Create pT,D histograms for efficiency calculation (3 cases: inclusive = 0, prompt only = 1, non-prompt only = 2)
    std::vector<TH1D*> hMcpPt;
    std::vector<TH1D*> hMcdPt;
    int histCaseNum = 3;
    for (size_t i = 0; i < histCaseNum; ++i) {
        hMcpPt.push_back(new TH1D(Form("mcp_pt_%zu",i), ";#it{p}_{T,D}^{truth};dN/d#it{p}_{T,D}^{truth}", binningStruct.ptDBinEdges_detector.size() - 1, binningStruct.ptDBinEdges_detector.data()));
        hMcpPt[i]->SetMarkerColor(kBlack);
        hMcpPt[i]->SetLineColor(kBlack);
        hMcpPt[i]->SetMarkerStyle(kOpenCircle);
        hMcpPt[i]->Sumw2();
        hMcpPt[i]->SetStats(0);
        hMcdPt.push_back(new TH1D(Form("mcd_pt_%zu",i), ";#it{p}_{T,D}^{reco};dN/d#it{p}_{T,D}^{reco}", binningStruct.ptDBinEdges_detector.size() - 1, binningStruct.ptDBinEdges_detector.data()));
        hMcdPt[i]->SetMarkerColor(kBlue);
        hMcdPt[i]->SetLineColor(kBlue);
        hMcdPt[i]->SetMarkerStyle(kFullCircle);
        hMcdPt[i]->Sumw2();
        hMcdPt[i]->SetStats(0);
    }

    // 2 ------ Fill histograms from TFile data
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on jet
    const double yCut = 0.8; // on D0
    double jetptMin = binningStruct.ptjetBinEdges_detector[0];
    double jetptMax = binningStruct.ptjetBinEdges_detector[binningStruct.ptjetBinEdges_detector.size() - 1];
    double MCPetaCut = etaCut;
    double MCPyCut = yCut;
    double MCDetaCut = etaCut;
    double MCDyCut = yCut;
    const double deltaRcut = binningStruct.deltaRBinEdges_detector[binningStruct.deltaRBinEdges_detector.size() - 1];
    double MCPDeltaRcut = deltaRcut;
    double MCDDeltaRcut = deltaRcut;
    const double MCPHfPtMincut = binningStruct.ptDBinEdges_particle[0]; // on particle level D0
    const double MCDHfPtMincut = binningStruct.ptDBinEdges_particle[0]; // on detector level D0
    const double MCPHfPtMaxcut = binningStruct.ptDBinEdges_particle[binningStruct.ptDBinEdges_particle.size() - 1]; // on particle level D0
    const double MCDHfPtMaxcut = binningStruct.ptDBinEdges_detector[binningStruct.ptDBinEdges_detector.size() - 1]; // on detector level D0
    // Particle level data
    TTree* tree = (TTree*)fClosureInputNonMatched->Get("CorrectionTree/O2mcpjetdisttable");
    if (!tree) {
        cout << "Error opening particle level correction tree.\n";
    }
    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    //int MCPjetnconst
    int MCPhfmatch;
    float MCPjetnconst;
    bool MCPhfprompt;
    

    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetnconst);
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    MCPhfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&MCPhfY);
    tree->SetBranchAddress("fMcHfPrompt",&MCPhfprompt);
    //tree->SetBranchAddress("fMCHfMatch",&MCPhfmatch);
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);
        bool genHfPtRange = (MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut);
        bool genDeltaRRange = (MCPDeltaR >= binningStruct.deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut);
        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && genJetPtRange && genHfPtRange && genDeltaRRange; // currently used

        // Fill particle level histograms
        if (genLevelRange) {
            
            // Fill inclusive histogram
            hMcpPt[0]->Fill(MCPhfPt);
            // fill prompt efficiency histogram
            if (MCPhfprompt) {
                hMcpPt[1]->Fill(MCPhfPt);
            } else{
                // fill non-prompt efficiency histogram
                hMcpPt[2]->Fill(MCPhfPt);
            }
            
        }

    }
    // Detector level data
    tree = (TTree*)fClosureInputNonMatched->Get("CorrectionTree/O2mcdjetdisttable");
    if (!tree) {
        cout << "Error opening detector level correction tree.\n";
    }
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    int MCDjetnconst;
    bool MCDhfmatch;
    //float MCDjetnconst;
    bool MCDhfprompt;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fJetNConst",&MCDjetnconst);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);
    tree->SetBranchAddress("fHfMatch",&MCDhfmatch);
    tree->SetBranchAddress("fHfMlScore0",&MCDhfMlScore0);
    tree->SetBranchAddress("fHfMlScore1",&MCDhfMlScore1);
    tree->SetBranchAddress("fHfMlScore2",&MCDhfMlScore2);
    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // only compute matched detector level candidates, but compute all particle level ones
        if (!MCDhfmatch) {
            continue;
        }

        // calculating delta R
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));
        
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = (MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut);
        bool recoDeltaRRange = (MCDDeltaR >= binningStruct.deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut);
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && recoJetPtRange && recoHfPtRange && recoDeltaRRange; // currently used

        // Get the threshold for this pT range
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;
        //std::cout << "Matched detector level entry." << std::endl;
        if (recoLevelRange) {
            //std::cout << "Passed recoLevelRange flag." << std::endl;
        }
        if (passBDTcut) {
            //std::cout << "Passed passBDTcut flag." << std::endl;
        }
        
        // Fill detector level histograms
        if (recoLevelRange && passBDTcut) {
            //std::cout << "Passed recoLevelRange and passBDTcut flags." << std::endl;
            // Fill inclusive histogram
            hMcdPt[0]->Fill(MCDhfPt);
            // fill prompt efficiency histogram
            if (MCDhfprompt) {
                hMcdPt[1]->Fill(MCDhfPt);
            } else{
                // fill non-prompt efficiency histogram
                hMcdPt[2]->Fill(MCDhfPt);
            }
            
        }

    }

    // 3 ----- Calculate efficiency distributions
    const char* efficiencyNames[] = {"inclusive", "prompt", "nonprompt"};
    // Loop through histogram cases: inclusive = 0, prompt only = 1, non-prompt only = 2
    //TCanvas* cPtEffRun2 = new TCanvas("cPtEffRun2", "pT,D distributions for efficiency run 2 style calculation",1800, 1000);
    //cPtEffRun2->Divide(2,2);
    for (int iEff = 0; iEff < histCaseNum; iEff++) {
        // Obtain MC pT distributions
        TH1D* mcpHist = hMcpPt[iEff];
        TH1D* mcdHist = hMcdPt[iEff];
        // cPtEffRun2->cd(iEff+1);
        // hMcpPt[iEff]->Draw();
        // hMcdPt[iEff]->Draw("same");
        // TLegend* legendPt = new TLegend(0.65,0.55,0.8,0.7);
        // legendPt->AddEntry(hMcpPt[iEff],"Particle level", "lpe");
        // legendPt->AddEntry(hMcdPt[iEff],"Detector level", "lpe");
        // legendPt->Draw();

        // Otain efficiency by dividing tempHist by parameter -> mcd/mcp
        TH1D* efficiencyHist = static_cast<TH1D*>(mcdHist->Clone(Form("efficiency_%s", efficiencyNames[iEff]))); // static needed to explicitly cast the cloned histogram to TH1D. ensures that the function returns a TH1D*, matching the return type.
        efficiencyHist->Divide(mcpHist);
        efficiencyHist->SetMarkerStyle(kFullCircle);
        efficiencyHist->SetMarkerColor(30+iEff*10); // 30 = not so bright green
        efficiencyHist->SetLineColor(30+iEff*10);
        efficiencyHist->SetTitle(";#it{p}_{T,D}^{truth};Efficiency#times Acceptance");
        efficiencyHist->SetStats(0);
        
        // Store efficiency distribution in struct
        hSelEff_run2style.push_back(efficiencyHist);

    }

    // Plot the efficiencies
    TCanvas* cEffRun2 = new TCanvas("cEffRun2","Run 2 style efficiencies",1800, 1000);
    cEffRun2->cd();
    hSelEff_run2style[0]->Draw();
    hSelEff_run2style[1]->Draw("same");
    hSelEff_run2style[2]->Draw("same");
    TLegend* legendEff = new TLegend(0.65,0.55,0.8,0.7);
    legendEff->AddEntry(hSelEff_run2style[0],"Inclusive", "lpe");
    legendEff->AddEntry(hSelEff_run2style[1],"Prompt D^{0}", "lpe");
    legendEff->AddEntry(hSelEff_run2style[2],"Non-prompt D^{0}", "lpe");
    legendEff->Draw();
    double statBoxPos = gPad->GetUxmax();
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->SetTextAlign(31); // align right
    latex->DrawLatex(0.85,0.42,Form("|#eta_{jet}| < %.1f, %.1f < #it{p}_{T,jet} < %.1f GeV/c", etaCut, jetptMin, jetptMax));

    return hSelEff_run2style;
}

struct EfficiencyData {
    // Yield 2D histograms: pT,jet vs pT,D: first = prompt D0 data distribution, second = non-prompt D0 distribution
    std::pair<TH2D*, TH2D*> hYieldTruth; // all particle level entries (denominator): will be kinematically corrected and folded
    std::pair<TH2D*, TH2D*> hYieldMeasured; // detector level entries (numerator): went over smearing effects and passed the selection cuts

    // Response matrices
    std::pair<RooUnfoldResponse, RooUnfoldResponse> response; // first = prompt D0s, second = non-prompt D0s

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

// Calculate particle and detector level kinematic efficiencies
EfficiencyData calculateKinematicEfficiency(TFile* fClosureInputMatched, const std::vector<TH1D*>& hSelEff_run2style, const BinningStruct& binningStruct, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    //
    std::pair<TH2D*, TH2D*> kinematicEfficiency; // first = particle level, second = detector level

    // Create struct to store data
    EfficiencyData histStruct;

    // 1 ----- Define histograms
    // Create 2D histograms for prompt D0s: raw and folded
    histStruct.hYieldTruth.first = new TH2D("hYieldTruthPrompt", "Particle level data prompt yield distribution;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binningStruct.ptjetBinEdges_particle.size() - 1, binningStruct.ptjetBinEdges_particle.data(), 
                                                binningStruct.ptDBinEdges_particle.size() - 1, binningStruct.ptDBinEdges_particle.data());
    histStruct.hYieldTruth.first->Sumw2();
    histStruct.hYieldMeasured.first = new TH2D("hYieldMeasuredPrompt", "Detector level data prompt yield distribution;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{det} (GeV/#it{c})", 
                                                binningStruct.ptjetBinEdges_detector.size() - 1, binningStruct.ptjetBinEdges_detector.data(), 
                                                binningStruct.ptDBinEdges_detector.size() - 1, binningStruct.ptDBinEdges_detector.data());
    // Create 2D histograms for non-prompt D0s: raw and folded
    histStruct.hYieldTruth.second = new TH2D("hYieldTruthNonPrompt", "Particle level data non-prompt yield distribution;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binningStruct.ptjetBinEdges_particle.size() - 1, binningStruct.ptjetBinEdges_particle.data(), 
                                                binningStruct.ptDBinEdges_particle.size() - 1, binningStruct.ptDBinEdges_particle.data());
    histStruct.hYieldTruth.second->Sumw2();
    histStruct.hYieldMeasured.second = new TH2D("hYieldMeasuredNonPrompt", "Detector level data non-prompt yield distribution;#it{p}_{T, ch. jet}^{det} (GeV/#it{c});#it{p}_{T, D^{0}}^{det} (GeV/#it{c})",
                                                binningStruct.ptjetBinEdges_detector.size() - 1, binningStruct.ptjetBinEdges_detector.data(), 
                                                binningStruct.ptDBinEdges_detector.size() - 1, binningStruct.ptDBinEdges_detector.data());

    // Create RooUnfoldResponse objects for prompt and non-prompt D0s
    histStruct.response.first = RooUnfoldResponse(histStruct.hYieldMeasured.first, histStruct.hYieldTruth.first); // prompt D0s
    histStruct.response.second = RooUnfoldResponse(histStruct.hYieldMeasured.second, histStruct.hYieldTruth.second); // non-prompt D0s
    // Kinematic efficiency histograms
    histStruct.hKEffResponseParticle.first = new TH2D("hKEffResponseParticlePrompt", "Truth prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binningStruct.ptjetBinEdges_particle.size() - 1, binningStruct.ptjetBinEdges_particle.data(), 
                                                binningStruct.ptDBinEdges_particle.size() - 1, binningStruct.ptDBinEdges_particle.data());
    histStruct.hKEffResponseParticle.second = new TH2D("hKEffResponseParticleNonPrompt", "Truth non-prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})", 
                                                binningStruct.ptjetBinEdges_particle.size() - 1, binningStruct.ptjetBinEdges_particle.data(), 
                                                binningStruct.ptDBinEdges_particle.size() - 1, binningStruct.ptDBinEdges_particle.data());
    histStruct.hKEffTruthTotalParticle.first = new TH2D("hKEffTruthTotalParticlePrompt", "Truth prompt D^{0}'s within total particle range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})",
                                                binningStruct.ptjetBinEdges_particle.size() - 1, binningStruct.ptjetBinEdges_particle.data(), 
                                                binningStruct.ptDBinEdges_particle.size() - 1, binningStruct.ptDBinEdges_particle.data());
    histStruct.hKEffTruthTotalParticle.second = new TH2D("hKEffTruthTotalParticleNonPrompt", "Truth non-prompt D^{0}'s within total particle range;#it{p}_{T, ch. jet}^{gen} (GeV/#it{c});#it{p}_{T, D^{0}}^{gen} (GeV/#it{c})",
                                                binningStruct.ptjetBinEdges_particle.size() - 1, binningStruct.ptjetBinEdges_particle.data(), 
                                                binningStruct.ptDBinEdges_particle.size() - 1, binningStruct.ptDBinEdges_particle.data());
    histStruct.hKEffResponseDetector.first = new TH2D("hKEffResponseDetectorPrompt", "Reconstructed prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            binningStruct.ptjetBinEdges_detector.size() - 1, binningStruct.ptjetBinEdges_detector.data(), 
            binningStruct.ptDBinEdges_detector.size() - 1, binningStruct.ptDBinEdges_detector.data());
    histStruct.hKEffResponseDetector.second = new TH2D("hKEffResponseDetectorNonPrompt", "Reconstructed non-prompt D^{0}'s within response range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}", 
            binningStruct.ptjetBinEdges_detector.size() - 1, binningStruct.ptjetBinEdges_detector.data(), 
            binningStruct.ptDBinEdges_detector.size() - 1, binningStruct.ptDBinEdges_detector.data());
    histStruct.hKEffRecoTotalDetector.first = new TH2D("hKEffRecoTotalDetectorPrompt", "Reconstructed prompt D^{0}'s within total detector range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}",
            binningStruct.ptjetBinEdges_detector.size() - 1, binningStruct.ptjetBinEdges_detector.data(), 
            binningStruct.ptDBinEdges_detector.size() - 1, binningStruct.ptDBinEdges_detector.data());
    histStruct.hKEffRecoTotalDetector.second = new TH2D("hKEffRecoTotalDetectorNonPrompt", "Reconstructed non-prompt D^{0}'s within total detector range;#it{p}_{T, ch. jet}^{reco};#it{p}_{T, D^{0}}^{reco}",
            binningStruct.ptjetBinEdges_detector.size() - 1, binningStruct.ptjetBinEdges_detector.data(), 
            binningStruct.ptDBinEdges_detector.size() - 1, binningStruct.ptDBinEdges_detector.data());
    
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
    TTree* tree = (TTree*)fClosureInputMatched->Get("CorrectionTree");
    // Check for correct access
    if (!tree) {
        cout << "Error opening matched closure tree.\n";
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
    tree->SetBranchAddress("fHfMlScore0",&MCDhfMlScore0);
    tree->SetBranchAddress("fHfMlScore1",&MCDhfMlScore1);
    tree->SetBranchAddress("fHfMlScore2",&MCDhfMlScore2);
    tree->SetBranchAddress("fHfMatchedFrom",&MCDhfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&MCDhfSelectedAs);

    TH1D* hEffWeight; // Histogram for efficiency weighting of response matrix
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

        //bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= 0.) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        //bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= 0.) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        //including underflow and overflow
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);
        bool genHfPtRange = (MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut);
        bool genDeltaRRange = (MCPDeltaR >= binningStruct.deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut);
        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && genJetPtRange && genHfPtRange && genDeltaRRange; // currently used
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = (MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut);
        bool recoDeltaRRange = (MCDDeltaR >= binningStruct.deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut);
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && recoJetPtRange && recoHfPtRange && recoDeltaRRange; // currently used

        // Get the threshold for this pT range: TODO - do NOT erase this, BDT cuts will be included in next hyperloop wagon
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;

        // Fill histograms considering jet pT and detector acceptance (response range)
        if (genLevelRange && recoLevelRange) { // genLevelRange && recoLevelRange && passBDTcut ?
            // fill 2D yields histograms
            if (MCPhfprompt) {
                // Get efficiency estimate to weight the response matrix
                hEffWeight = hSelEff_run2style[1];
                // Get efficiency estimate to weight the response matrix: jet pT shape is influenced by D0 pT efficiency
                int bin = hEffWeight->FindBin(MCDhfPt);
                double estimatedEfficiency = hEffWeight->GetBinContent(bin);
                if (estimatedEfficiency == 0) {
                    std::cout << "Warning: Prompt efficiency is zero for pT = " << MCDhfPt << " with bin " << bin << ". Setting it to 1." << std::endl;
                    estimatedEfficiency = 1; // Avoid division by zero
                }
                // Fill 4D RooUnfoldResponse object
                histStruct.response.first.Fill(MCDjetPt, MCDhfPt, MCPjetPt, MCPhfPt, 1 / estimatedEfficiency); // jet pT shape is influenced by D0 pT efficiency
                //histStruct.responseProjections.first[0]->Fill(MCDjetPt, MCPjetPt,1 / estimatedEfficiency); // pT,jet projection, prompt D0s
                //histStruct.responseProjections.first[1]->Fill(MCDhfPt, MCPhfPt, 1 / estimatedEfficiency); // pT,D projection, prompt D0s
            } else{
                // Get efficiency estimate to weight the response matrix
                hEffWeight = hSelEff_run2style[2];
                // Get efficiency estimate to weight the response matrix: jet pT shape is influenced by D0 pT efficiency
                int bin = hEffWeight->FindBin(MCDhfPt);
                double estimatedEfficiency = hEffWeight->GetBinContent(bin);
                if (estimatedEfficiency == 0) {
                    std::cout << "Warning: Non-prompt efficiency is zero for pT = " << MCDhfPt << " with bin " << bin << ". Setting it to 1." << std::endl;
                    estimatedEfficiency = 1; // Avoid division by zero
                }
                // Fill 4D RooUnfoldResponse object
                histStruct.response.second.Fill(MCDjetPt, MCDhfPt, MCPjetPt, MCPhfPt, 1 / estimatedEfficiency); // jet pT shape is influenced by D0 pT efficiency
                //histStruct.responseProjections.second[0]->Fill(MCDjetPt, MCPjetPt,1 / estimatedEfficiency); // pT,jet projection, non-prompt D0s
                //histStruct.responseProjections.second[1]->Fill(MCDhfPt, MCPhfPt, 1 / estimatedEfficiency); // pT,D projection, non-prompt D0s
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

    return histStruct;
}

// Estimate the selection efficiency from detector level MC data (run 2 style)
std::vector<TH1D*> calculateSelectionEfficiencyRun3(TFile* fClosureInputNonMatched, const std::vector<TH1D*>& hSelEff_run2style, EfficiencyData& histStruct, const BinningStruct& binningStruct, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    // Efficiency = detector / particle (converted) level

    // 1 ----- Obtain the numerator and denominator histograms: first = prompt, second = non-prompt
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on jet
    const double yCut = 0.8; // on D0
    double jetptMin = binningStruct.ptjetBinEdges_detector[0];
    double jetptMax = binningStruct.ptjetBinEdges_detector[binningStruct.ptjetBinEdges_detector.size() - 1];
    double MCPetaCut = etaCut;
    double MCPyCut = yCut;
    double MCDetaCut = etaCut;
    double MCDyCut = yCut;
    const double deltaRcut = binningStruct.deltaRBinEdges_detector[binningStruct.deltaRBinEdges_detector.size() - 1];
    double MCPDeltaRcut = deltaRcut;
    double MCDDeltaRcut = deltaRcut;
    const double MCPHfPtMincut = binningStruct.ptDBinEdges_particle[0]; // on particle level D0
    const double MCDHfPtMincut = binningStruct.ptDBinEdges_particle[0]; // on detector level D0
    const double MCPHfPtMaxcut = binningStruct.ptDBinEdges_particle[binningStruct.ptDBinEdges_particle.size() - 1]; // on particle level D0
    const double MCDHfPtMaxcut = binningStruct.ptDBinEdges_detector[binningStruct.ptDBinEdges_detector.size() - 1]; // on detector level D0
    // Particle level data
    TTree* tree = (TTree*)fClosureInputNonMatched->Get("CorrectionTree/O2mcpjetdisttable");
    if (!tree) {
        cout << "Error opening particle level correction tree.\n";
    }
    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    //int MCPjetnconst
    int MCPhfmatch;
    float MCPjetnconst;
    bool MCPhfprompt;
    

    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetnconst);
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    MCPhfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&MCPhfY);
    tree->SetBranchAddress("fMcHfPrompt",&MCPhfprompt);
    //tree->SetBranchAddress("fMCHfMatch",&MCPhfmatch);
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);
        bool genHfPtRange = (MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut);
        bool genDeltaRRange = (MCPDeltaR >= binningStruct.deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut);
        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && genJetPtRange && genHfPtRange && genDeltaRRange; // currently used

        // Fill histograms considering jet pT and detector acceptance
        if (genLevelRange) {
            
            // fill prompt efficiency histogram
            if (MCPhfprompt) {
                // fill prompt 2D yield: particle level matched
                histStruct.hYieldTruth.first->Fill(MCPjetPt, MCPhfPt);
            } else{
                // fill non-prompt 2D yield: particle level matched
                histStruct.hYieldTruth.second->Fill(MCPjetPt, MCPhfPt);
            }
            
        }

    }
    // Detector level data
    tree = (TTree*)fClosureInputNonMatched->Get("CorrectionTree/O2mcdjetdisttable");
    if (!tree) {
        cout << "Error opening detector level correction tree.\n";
    }
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    int MCDhfMatchedFrom, MCDhfSelectedAs, MCDjetnconst;
    bool MCDhfmatch;
    //float MCDjetnconst;
    bool MCDhfprompt;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fJetNConst",&MCDjetnconst);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);
    tree->SetBranchAddress("fHfMatch",&MCDhfmatch);
    tree->SetBranchAddress("fHfMlScore0",&MCDhfMlScore0);
    tree->SetBranchAddress("fHfMlScore1",&MCDhfMlScore1);
    tree->SetBranchAddress("fHfMlScore2",&MCDhfMlScore2);
    tree->SetBranchAddress("fHfMatchedFrom",&MCDhfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&MCDhfSelectedAs);
    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = (MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut);
        bool recoDeltaRRange = (MCDDeltaR >= binningStruct.deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut);
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && recoJetPtRange && recoHfPtRange && recoDeltaRRange; // currently used
        bool isReflection = (MCDhfMatchedFrom != MCDhfSelectedAs) ? true : false;

        // only compute matched detector level candidates, but compute all particle level ones
        if (!MCDhfmatch || isReflection) {
            continue;
        }

        // Get the threshold for this pT range
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;
        
        // Fill detector level histograms
        if (recoLevelRange && passBDTcut) { // default = recoLevelRange && passBDTcut
            
            // fill prompt efficiency histogram
            if (MCDhfprompt) {
                // fill prompt 2D yield: detector level matched
                histStruct.hYieldMeasured.first->Fill(MCDjetPt, MCDhfPt);
            } else{
                // fill non-prompt 2D yield: detector level matched
                histStruct.hYieldMeasured.second->Fill(MCDjetPt, MCDhfPt);
            }
            
        }

    }

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
    // Calculate the efficiency distributions with the pT projections
    histStruct.hSelectionEfficiency.first = histStruct.hYieldMeasured.first->ProjectionY("hSelectionEfficiencyPrompt", minBin, maxBin);
    histStruct.hSelectionEfficiency.second = histStruct.hYieldMeasured.second->ProjectionY("hSelectionEfficiencyNonPrompt", minBin, maxBin);
    histStruct.hSelectionEfficiency.first->SetTitle("Efficiency prompt D^{0} p_{T} distribution; #it{p}_{T, D^{0}}^{reco}; Efficiency#times Acceptance");
    histStruct.hSelectionEfficiency.second->SetTitle("Efficiency non-prompt D^{0} p_{T} distribution; #it{p}_{T, D^{0}}^{reco}; Efficiency#times Acceptance");
    histStruct.hSelectionEfficiency.first->Sumw2();
    histStruct.hSelectionEfficiency.second->Sumw2();
    histStruct.hSelectionEfficiency.first->Divide(histStruct.hHfPtYieldTruthCorrected.first);
    histStruct.hSelectionEfficiency.second->Divide(histStruct.hHfPtYieldTruthCorrected.second);


    std::vector<TH1D*> hSelEff_run3style = {nullptr, nullptr};
    hSelEff_run3style[0] = histStruct.hSelectionEfficiency.first;
    hSelEff_run3style[1] = histStruct.hSelectionEfficiency.second;

    // Plot run 3 style efficiencies
    TCanvas* cEffRun3 = new TCanvas("cEffRun3", "Run 3 style efficiency", 1800, 1000);
    cEffRun3->cd();
    hSelEff_run3style[0]->SetLineColor(kBlue+1);
    hSelEff_run3style[0]->Draw();
    hSelEff_run3style[1]->SetLineColor(kRed+1);
    hSelEff_run3style[1]->Draw("same");
    TLegend* legendEffRun3 = new TLegend(0.65,0.55,0.8,0.7);
    legendEffRun3->AddEntry(hSelEff_run3style[0],"Prompt D^{0}", "lpe");
    legendEffRun3->AddEntry(hSelEff_run3style[1],"Non-prompt D^{0}", "lpe");
    legendEffRun3->Draw();

    return hSelEff_run3style;
}

TH2D* performEfficiencyCorrection(TH3D* hBackgroundSubtracted, EfficiencyData& histStruct) {
     // 1 - Clone background subtracted distribution
    histStruct.hEfficiencyCorrected.first = (TH3D*)hBackgroundSubtracted->Clone("h3DEfficiencyCorrected");
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
                if (eff > 0) { // eff > 0
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
    
    // Output object
    TH2D* hEfficiencyCorrected = histStruct.hEfficiencyCorrected.second;

    TCanvas* cEffCorrected = new TCanvas("cEffCorrected", "Background subtracted, efficiency (run 3 style) corrected", 1800, 1000);
    cEffCorrected->cd();
    hEfficiencyCorrected->Draw("colz");

    std::cout << "Efficiency (run 3 style) corrected." << std::endl;
    return hEfficiencyCorrected;
}

EfficiencyData EfficiencyClosure(TFile* fClosureInputNonMatched, TFile* fClosureInputMatched, TH3D* hBackgroundSubtracted, const BinningStruct& binningStruct, const std::vector<std::pair<double, double>>& bdtPtCuts) {

    // Calculate run 2 style efficiency (to be used in run 3 style calculation): inclusive, prompt-only, non-prompt only
    //std::vector<TH1D*> hSelEff_run2style = calculateSelectionEfficiencyRun2Matched(fClosureInputMatched, binningStruct, bdtPtCuts);
    std::vector<TH1D*> hSelEff_run2style = calculateSelectionEfficiencyRun2NonMatched(fClosureInputNonMatched, binningStruct, bdtPtCuts);

    // Calculate particle and detector level kinematic efficiencies
    EfficiencyData histStruct = calculateKinematicEfficiency(fClosureInputMatched, hSelEff_run2style, binningStruct, bdtPtCuts);

    // Calculate run 3 style efficiency: 0 = prompt, 1 = non-prompt
    std::vector<TH1D*> hSelEff_run3style = calculateSelectionEfficiencyRun3(fClosureInputNonMatched, hSelEff_run2style, histStruct, binningStruct, bdtPtCuts);

    // Perform efficiency correction to background subtracted data
    TH2D* hEfficiencyCorrected = performEfficiencyCorrection(hBackgroundSubtracted, histStruct);
    
    std::pair<std::vector<TH1D*>, TH2D*> effOutputs = std::make_pair(hSelEff_run3style, hEfficiencyCorrected);

    return histStruct;
}
