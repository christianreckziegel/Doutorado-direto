/*
 *
 *
 * Macro for plotting HF pT dependent efficiency 
 * and applying correction to BackgroundSubtraction.C
 * and SignalExtraction.C resulting distributions.
 * 
 * 
 * 
 * 
**/

#include <set>

using namespace std;

// calculate number of background subtracted histograms inside file
int HistogramCounter(TFile* file) {
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

struct EfficiencyData {
    std::vector<TH2D*> hMcpPt;              // particle level pT,D distribution: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH2D*> hMcdPt;              // detector level pT,D distribution: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH2D*> hEfficiencies;       // efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hBackSubCorrected;   // background subtracted (sideband method) DeltaR distributions corrected by efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2
    TH1D* hBackSubCorrected_allpt;          // resulting summed all pT,D bins distribution from sidebad method
    
    // Investigation histograms
    TH1D* hBDTBackgroundScore;              // rejected entries that didn't pass the BDT background score cuts
};

// Module to create histograms including interest variable
EfficiencyData createHistograms(const std::vector<double>& ptjetBinEdges_particle, const std::vector<double>& ptDBinEdges_particle,
                                 const std::vector<double>& ptjetBinEdges_detector, const std::vector<double>& ptDBinEdges_detector) {

    // Create struct to store data
    EfficiencyData histStruct;

    // Create 3 histogram cases: inclusive = 0, prompt only = 1, non-prompt only = 2
    int histCaseNum = 3;
    for (size_t i = 0; i < histCaseNum; ++i) {
        histStruct.hMcpPt.push_back(new TH2D(Form("mcp_pt_%zu",i), ";p_{T,jet}^{truth};p_{T,D}^{truth}", ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(),ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data()));
        histStruct.hMcpPt[i]->SetMarkerColor(kBlack);
        histStruct.hMcpPt[i]->SetLineColor(kBlack);
        histStruct.hMcpPt[i]->SetMarkerStyle(kOpenCircle);
        histStruct.hMcpPt[i]->Sumw2();
        histStruct.hMcpPt[i]->SetStats(0);
        histStruct.hMcdPt.push_back(new TH2D(Form("mcd_pt_%zu",i), ";p_{T,jet}^{reco};p_{T,D}^{reco}", ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(),ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data()));
        histStruct.hMcdPt[i]->SetMarkerColor(kBlue);
        histStruct.hMcdPt[i]->SetLineColor(kBlue);
        histStruct.hMcdPt[i]->SetMarkerStyle(kFullCircle);
        histStruct.hMcdPt[i]->Sumw2();
        histStruct.hMcdPt[i]->SetStats(0);
    }

    // Creating investigation histogram
    histStruct.hBDTBackgroundScore = new TH1D("hBDTBackgroundScore", "Entries that didn't pass the cuts;BDT background score;Counts", 100, 0, 1);

    cout << "Histograms created.\n";

    return histStruct;
}

// Find the appropriate cut value based on pT
double GetBkgProbabilityCut(double pT, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    for (size_t i = 0; i < bdtPtCuts.size() - 1; ++i) {
        if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) {
            return bdtPtCuts[i].second;
        }
    }
    return 1.0; // Default: accept all if out of range
}
// Module to fill histograms from TFile data
void fillHistograms(TFile* fSimulated, const EfficiencyData& histStruct, double& jetptMin, double& jetptMax, double& hfptMin, double& hfptMax, const std::vector<std::pair<double, double>>& bdtPtCuts, const double& backProbabilityCut) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on jet
    const double yCut = 0.8; // on D0
    const double deltaRcut = 0.4;
    const double MCPHfPtMincut = hfptMin; // on particle level D0
    const double MCDHfPtMincut = hfptMin; // on detector level D0
    const double MCPHfPtMaxcut = hfptMax; // on particle level D0
    const double MCDHfPtMaxcut = hfptMax; // on detector level D0

    //
    // MC generator level tree and histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)fSimulated->Get("DF_merged/O2mcpjetdisttable");
    
    // Check for correct access
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing the TTree
    float axisDistance, jetPt, jetEta, jetPhi;
    float hfPt, hfEta, hfPhi, hfMass, hfY;
    float jetnconst_float;
    bool hfprompt, hfmatch;

    tree->SetBranchAddress("fMCJetHfDist",&axisDistance);
    tree->SetBranchAddress("fMCJetPt",&jetPt);
    tree->SetBranchAddress("fMCJetEta",&jetEta);
    tree->SetBranchAddress("fMCJetPhi",&jetPhi);
    tree->SetBranchAddress("fMCJetNConst",&jetnconst_float);
    tree->SetBranchAddress("fMCHfPt",&hfPt);
    tree->SetBranchAddress("fMCHfEta",&hfEta);
    tree->SetBranchAddress("fMCHfPhi",&hfPhi);
    hfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMCHfY",&hfY);
    tree->SetBranchAddress("fMCHfPrompt",&hfprompt);
    tree->SetBranchAddress("fMCHfMatch",&hfmatch);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));
        bool genLevelRange = (abs(hfEta) < etaCut) && (abs(hfY) < yCut) && (jetPt >= jetptMin) && (jetPt < jetptMax) && ((deltaR >= 0.) && (deltaR < deltaRcut)) && ((hfPt >= MCPHfPtMincut) && (hfPt < MCPHfPtMaxcut));
        
        // Fill histograms considering jet pT and detector acceptance
        if (genLevelRange) {
        //if ((abs(hfEta) < etaCut) && (jetPt >= jetptMin) && (jetPt < jetptMax)) {
            
            // Fill inclusive histogram
            histStruct.hMcpPt[0]->Fill(jetPt, hfPt);
            // fill prompt efficiency histogram
            if (hfprompt) {
                histStruct.hMcpPt[1]->Fill(jetPt, hfPt);
            } else{
                // fill non-prompt efficiency histogram
                histStruct.hMcpPt[2]->Fill(jetPt, hfPt);
            }
            
        }
        
        
    }
    cout << "Generator level histograms filled.\n";

    //
    // MC detector level tree and histograms
    //
    // Accessing TTree
    tree = (TTree*)fSimulated->Get("DF_merged/O2mcdjetdisttable");
    // Check for correct access
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining ML score variables for accessing the TTree
    float hfMlScore0, hfMlScore1, hfMlScore2;
    int jetnconst_int;
    tree->SetBranchAddress("fJetHfDist",&axisDistance);
    tree->SetBranchAddress("fJetPt",&jetPt);
    tree->SetBranchAddress("fJetEta",&jetEta);
    tree->SetBranchAddress("fJetPhi",&jetPhi);
    tree->SetBranchAddress("fJetNConst",&jetnconst_int);
    tree->SetBranchAddress("fHfPt",&hfPt);
    tree->SetBranchAddress("fHfEta",&hfEta);
    tree->SetBranchAddress("fHfPhi",&hfPhi);
    tree->SetBranchAddress("fHfMass",&hfMass);
    tree->SetBranchAddress("fHfY",&hfY);
    tree->SetBranchAddress("fHfPrompt",&hfprompt);
    tree->SetBranchAddress("fHfMatch",&hfmatch);
    tree->SetBranchAddress("fHfMlScore0",&hfMlScore0);
    tree->SetBranchAddress("fHfMlScore1",&hfMlScore1);
    tree->SetBranchAddress("fHfMlScore2",&hfMlScore2);

    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // only compute matched detector level candidates, but compute all particle level ones
        if (!hfmatch) {
            continue;
        }
        
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));
        bool recoLevelRange = (abs(hfEta) < etaCut) && (abs(hfY) < yCut) && (jetPt >= jetptMin) && (jetPt < jetptMax) && ((deltaR >= 0.) && (deltaR < deltaRcut)) && ((hfPt >= MCDHfPtMincut) && (hfPt < MCDHfPtMaxcut));

        // Fill histograms considering jet pT and detector acceptance
        if (recoLevelRange) {
            
            // Get the threshold for this pT range
            double maxBkgProb = backProbabilityCut * GetBkgProbabilityCut(hfPt, bdtPtCuts);

            // Fill histogram only if the BDT cut is passed
            if (hfMlScore0 < maxBkgProb) {
                // Fill inclusive histogram
                histStruct.hMcdPt[0]->Fill(jetPt, hfPt);
                // fill prompt efficiency histogram
                if (hfprompt) {
                    histStruct.hMcdPt[1]->Fill(jetPt, hfPt);
                } else{
                    // fill non-prompt efficiency histogram
                    histStruct.hMcdPt[2]->Fill(jetPt, hfPt);
                }
            } else {
                histStruct.hBDTBackgroundScore->Fill(hfMlScore0);
            }
                        
        }
        
        
    }
    cout << "Detector level histograms filled.\n";
}

void CalculateEfficiency(EfficiencyData& histStruct){
    const char* efficiencyNames[] = {"inclusive", "prompt", "nonprompt"};

    // Loop through histogram cases: inclusive = 0, prompt only = 1, non-prompt only = 2
    int histCaseNum = 3;
    for (int iEff = 0; iEff < histCaseNum; iEff++) {
        // Obtain MC pT distributions
        TH2D* mcdHist = histStruct.hMcdPt[iEff];
        TH2D* mcpHist = histStruct.hMcpPt[iEff];

        // Otain efficiency by dividing tempHist by parameter -> mcd/mcp
        TH2D* efficiencyHist = static_cast<TH2D*>(mcdHist->Clone(Form("efficiency_%s", efficiencyNames[iEff]))); // static needed to explicitly cast the cloned histogram to TH1D. ensures that the function returns a TH1D*, matching the return type.
        efficiencyHist->Divide(mcpHist);
        efficiencyHist->SetMarkerStyle(kFullCircle);
        efficiencyHist->SetMarkerColor(30+iEff*10); // 30 = not so bright green
        efficiencyHist->SetLineColor(30+iEff*10);
        efficiencyHist->SetTitle("Efficiency#times Acceptance;p_{T,jet}^{truth};p_{T,D}^{truth}");
        efficiencyHist->SetStats(0);
        
        // Store efficiency distribution in struct
        histStruct.hEfficiencies.push_back(efficiencyHist);

    }

    cout << "Efficiency calculated.\n";
}

void PerformCorrections(EfficiencyData& histStruct, TFile* fBackSub) {
    int effOption = 1; // choice of efficiency for correction: 0 = inclusive, 1 = prompt, 2 = non-prompt
    double efficiency = 0;
    int ptDBin;
    TH1D* hBackSub_temp = (TH1D*)fBackSub->Get("h_back_subtracted_0");

    // Initialize all pT,D histogram before adding them
    histStruct.hBackSubCorrected_allpt = (TH1D*)hBackSub_temp->Clone("hBackSubCorrected_allpt");
    histStruct.hBackSubCorrected_allpt->Reset();
    
    // Loop through distributions from all intervals
    int numHistos = HistogramCounter(fBackSub);
    for (size_t iHisto = 0; iHisto < numHistos; iHisto++) {
        // access each distribution
        hBackSub_temp = (TH1D*)fBackSub->Get(Form("h_back_subtracted_%zu",iHisto));

        // obtain efficiency to the correspondent pT interval (bin)
        int firstJetBin = 1;
        int lastJetBin = histStruct.hEfficiencies[effOption]->GetNbinsX();
        TH1D* efficiencyHist = histStruct.hEfficiencies[effOption]->ProjectionY(Form("efficiency_projection_%zu",iHisto), firstJetBin, lastJetBin);
        //bin = histStruct.hEfficiencies[effOption]->GetBin(iHisto);
        ptDBin = iHisto + 1; // bin number starts from 1
        efficiency = efficiencyHist->GetBinContent(ptDBin);
        
        // apply efficiency correction to distribution
        if (efficiency > 0) {
            hBackSub_temp->Scale(1.0 / efficiency); // Scale the histogram by the efficiency value
        } else {
            std::cerr << "Warning: Efficiency is zero or invalid for bin " << ptDBin << "\n";
        }
        
        // modifying visual characteristics
        hBackSub_temp->SetStats(0);
        hBackSub_temp->GetXaxis()->SetTitleSize(0.05);
        hBackSub_temp->GetYaxis()->SetTitleSize(0.05);

        // store corrected distribution to struct
        histStruct.hBackSubCorrected.push_back(hBackSub_temp);

        // add to final distribution of all pT,D's
        histStruct.hBackSubCorrected_allpt->Add(histStruct.hBackSubCorrected[iHisto]);
    }
    
    cout << "Histograms corrected.\n";
}

void PlotHistograms(const EfficiencyData& histStruct, double jetptMin, double jetptMax) {
    cout << "Plotting histograms...\n";

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    std::vector<TLegend*> legendPt;
    legendPt.push_back(new TLegend(0.65,0.49,0.8,0.62));
    legendPt[0]->AddEntry(histStruct.hMcpPt[0],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    legendPt[0]->AddEntry(histStruct.hMcdPt[0],"Reconstructed", "lpe");
    legendPt.push_back(new TLegend(0.65,0.49,0.8,0.62));
    legendPt[1]->AddEntry(histStruct.hMcpPt[1],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    legendPt[1]->AddEntry(histStruct.hMcdPt[1],"Reconstructed", "lpe");
    legendPt.push_back(new TLegend(0.65,0.49,0.8,0.62));
    legendPt[2]->AddEntry(histStruct.hMcpPt[2],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    legendPt[2]->AddEntry(histStruct.hMcdPt[2],"Reconstructed", "lpe");
    legendPt.push_back(new TLegend(0.65,0.47,0.85,0.62));
    legendPt[3]->AddEntry(histStruct.hMcpPt[0],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    legendPt[3]->AddEntry(histStruct.hMcdPt[0],"Inclusive reconstructed", "lpe");
    legendPt[3]->AddEntry(histStruct.hMcdPt[1],"Prompt D^{0} reconstructed", "lpe");
    legendPt[3]->AddEntry(histStruct.hMcdPt[2],"Non-prompt D^{0} reconstructed", "lpe");

    
    //
    // 1 - plot pT,D of generator and reconstruction level
    //
    TCanvas* cEff_1 = new TCanvas("cEff_1","pT,D spectra");
    cEff_1->Divide(2,2);
    cEff_1->SetCanvasSize(1800,1000);
    cEff_1->cd(1);
    histStruct.hMcpPt[0]->Draw();
    histStruct.hMcdPt[0]->Draw("same");
    legendPt[0]->Draw();
    double statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Inclusive");
    cEff_1->cd(2);
    histStruct.hMcpPt[1]->Draw();
    histStruct.hMcdPt[1]->Draw("same");
    legendPt[1]->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Prompt D^{0}'s");
    cEff_1->cd(3);
    histStruct.hMcpPt[2]->Draw();
    histStruct.hMcdPt[2]->Draw("same");
    legendPt[2]->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Non-prompt D^{0}'s");
    cEff_1->cd(4);
    histStruct.hMcdPt[0]->SetMarkerColor(30+0*10); // 30 = not so bright green
    histStruct.hMcdPt[0]->SetLineColor(30+0*10); // 30 = not so bright green
    histStruct.hMcdPt[1]->SetMarkerColor(30+1*10); // 30 = not so bright green
    histStruct.hMcdPt[1]->SetLineColor(30+1*10); // 30 = not so bright green
    histStruct.hMcdPt[2]->SetMarkerColor(30+2*10); // 30 = not so bright green
    histStruct.hMcdPt[2]->SetLineColor(30+2*10); // 30 = not so bright green
    histStruct.hMcdPt[0]->Draw();
    histStruct.hMcdPt[1]->Draw("same");
    histStruct.hMcdPt[2]->Draw("same");
    legendPt[3]->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Reconstructed: inclusive, prompt,\n non-prompt");
    
    //
    // 2 - plot three kinds of efficiencies: [0] = inclusive, [1] = prompt D0, [2] = non-prompt D0
    //
    // Loop through efficiency cases
    TCanvas* cEff_2 = new TCanvas("cEff_2","Prompt/non-prompt fficiency plots");
    cEff_2->SetCanvasSize(1800,1000);
    cEff_2->cd();
    TLegend* legendEff = new TLegend(0.65,0.45,0.8,0.58);
    legendEff->AddEntry(histStruct.hEfficiencies[0],"Inclusive", "lpe");
    legendEff->AddEntry(histStruct.hEfficiencies[1],"Prompt D^{0}", "lpe");
    legendEff->AddEntry(histStruct.hEfficiencies[2],"Non-prompt D^{0}", "lpe");
    for (size_t iEff = 0; iEff < histStruct.hEfficiencies.size(); iEff++) {
        // if the current is the first plot in the canvas
        if (iEff == 0) {
            histStruct.hEfficiencies[iEff]->SetMinimum(0);
            histStruct.hEfficiencies[iEff]->SetMaximum(0.8);
            histStruct.hEfficiencies[iEff]->Draw();
        } else {
            histStruct.hEfficiencies[iEff]->Draw("same");
        }
        
        statBoxPos = gPad->GetUxmax();
        latex->DrawLatex(statBoxPos-0.35, 0.61, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    }
    legendEff->Draw();
    cEff_2->Update();

    

    //
    // 3 - plot efficiency corrected histograms
    //
    TCanvas* cCorrectedBackSub = new TCanvas("cCorrectedBackSub","Efficiency corrected histograms");
    int nHistos = histStruct.hBackSubCorrected.size();
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));
    cCorrectedBackSub->Divide(nCols,nRows);
    
    for (size_t iHisto = 0; iHisto < histStruct.hBackSubCorrected.size(); iHisto++) {
        //
        cCorrectedBackSub->cd(iHisto+1);
        histStruct.hBackSubCorrected[iHisto]->Draw();
        statBoxPos = gPad->GetUxmax();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
        latex->DrawLatex(statBoxPos-0.35, 0.75, "Efficiency corrected");
    }
    cCorrectedBackSub->SetCanvasSize(1800,1000);

    //
    // 4 - plot final corrected 
    //
    TCanvas* cCorrectedBackSub_allpt = new TCanvas("cCorrectedBackSub_allpt","Final efficiency corrected ");
    cCorrectedBackSub_allpt->cd();
    histStruct.hBackSubCorrected_allpt->SetTitle(";#DeltaR_{D^{0}};#frac{dN}{d(#DeltaR_{D^{0}})}");
    histStruct.hBackSubCorrected_allpt->SetStats(0);
    histStruct.hBackSubCorrected_allpt->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.75, "Efficiency corrected");
    cCorrectedBackSub_allpt->SetCanvasSize(1800,1000);


    //
    //
    //
    TCanvas* cInvestigation = new TCanvas("cInvestigation","Investigation");
    cInvestigation->SetCanvasSize(1800,1000);
    cInvestigation->cd();
    histStruct.hBDTBackgroundScore->Draw();
    cInvestigation->SetCanvasSize(1800,1000);
    //
    // Storing images
    //
    TString imagePath = "../Images/2-Efficiency/Run2Style/2D/";
    cEff_1->Update();
    cEff_1->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cEff_1->SaveAs(imagePath + "Efficiency_dNdpT_run2_style.png");
    cEff_2->Update();
    cEff_2->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cEff_2->SaveAs(imagePath + "Efficiency_acceptance_run2_style.png");
    cCorrectedBackSub->Update();
    cCorrectedBackSub->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cCorrectedBackSub->SaveAs(imagePath + "Efficiency_corrected_pT_bins_run2_style.png");
    cCorrectedBackSub_allpt->Update();
    cCorrectedBackSub_allpt->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cCorrectedBackSub_allpt->SaveAs(imagePath + "Efficiency_corrected_run2_style.png");

    //
    // Storing in a single pdf file
    //
    cEff_1->Print(Form(imagePath + "pT_efficiency_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cEff_2->Print(Form(imagePath + "pT_efficiency_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cCorrectedBackSub->Print(Form(imagePath + "pT_efficiency_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cCorrectedBackSub_allpt->Print(Form(imagePath + "pT_efficiency_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

}

void SaveData(const EfficiencyData& histStruct, double& jetptMin, double& jetptMax, const double& backProbabilityCut){
    // Open output file
    //TFile* outFile = new TFile(Form("run2_style_efficiency_%d_to_%d_jetpt_%.1f.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax),backProbabilityCut),"recreate");
    TFile* outFile = new TFile(Form("run2_style_2D_efficiency_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    // Loop over signal extracted histograms
    for (size_t iHisto = 0; iHisto < histStruct.hBackSubCorrected.size(); iHisto++) {
        // store each histogram in file
        histStruct.hBackSubCorrected[iHisto]->Write();
    }
    // Loop over efficiency histograms
    for (size_t iHisto = 0; iHisto < histStruct.hEfficiencies.size(); iHisto++) {
        // store each histogram in file
        histStruct.hEfficiencies[iHisto]->Write();
    }
    // Final all pT,D distribution
    histStruct.hBackSubCorrected_allpt->Write();

    outFile->Close();
    delete outFile;
    
    cout << "Data stored.\n";
}

void Efficiency_run2_style_2D(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;
    // jet pT cuts
    std::vector<double> ptjetBinEdges_particle = {5., 7., 15., 30., 50.}; // TODO: use 5., 7., 15., 30., 50. for final version
    std::vector<double> ptjetBinEdges_detector = {5., 7., 15., 30., 50.};
    double jetptMin = ptjetBinEdges_particle[0]; // GeV
    double jetptMax = ptjetBinEdges_particle[ptjetBinEdges_particle.size() - 1]; // GeV
    // deltaR histogram
    int deltaRbins = 10; // deltaRbins = numberOfPoints, default=100 bins for [0. 1.0]
    const std::vector<double> deltaRBinEdges = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15,0.175, 0.2, 0.25, 0.3, 0.35, 0.4}; // default = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4} chosen by Nima
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
    // pT,D histograms
    std::vector<double> ptDBinEdges_particle = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.}; // default pT,D = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.}
    std::vector<double> ptDBinEdges_detector = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.}; // default pT,D = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.}
    double hfptMin = ptDBinEdges_particle[0];
    double hfptMax = ptDBinEdges_particle[ptDBinEdges_particle.size() - 1];
    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.03 (from first of pair)
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}
    };
    double backProbabilityCut = 1.0; // default value
    
    EfficiencyData histStruct = createHistograms(ptjetBinEdges_particle, ptDBinEdges_particle, ptjetBinEdges_detector, ptDBinEdges_detector); // pT histograms

    // opening files
    TFile* fSimulated = new TFile("../SimulatedData/Hyperloop_output/McEfficiency/New_with_reflections/Merged_AO2D_HF_LHC24d3a_All.root","read"); //../SimulatedData/Hyperloop_output/AO2D_merged.root
    TFile* fBackSub = new TFile(Form("../1-SignalTreatment/SideBand/backSub_%d_to_%d_jetpt_with_reflections.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fSimulated || fSimulated->IsZombie()) {
        std::cerr << "Error: Unable to open simulated data ROOT file." << std::endl;
    }
    if (!fBackSub || fBackSub->IsZombie()) {
        std::cerr << "Error: Unable to open background subtracted data ROOT file." << std::endl;
    }
    
    // Fill histograms
    fillHistograms(fSimulated, histStruct, jetptMin, jetptMax, hfptMin, hfptMax, bdtPtCuts, backProbabilityCut);

    // Calculate efficiency distribution
    CalculateEfficiency(histStruct);

    // Apply efficiency corrections to distributions
    PerformCorrections(histStruct, fBackSub);

    // Plot the efficiency histogram and further corrected histograms
    PlotHistograms(histStruct, jetptMin, jetptMax);

    // Save corrected distributions to file
    SaveData(histStruct, jetptMin, jetptMax, backProbabilityCut);

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    Efficiency_run2_style_2D();
    return 0;
}
