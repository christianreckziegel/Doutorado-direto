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
    std::vector<TH1D*> hMcpPt;              // particle level pT,D distribution: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hMcdPt;              // detector level pT,D distribution: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hEfficiencies;       // efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hBackSubCorrected;   // background subtracted (sideband method) DeltaR distributions corrected by efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2
    TH1D* hBackSubCorrected_allpt;          // resulting summed all pT,D bins distribution from sidebad method
    std::vector<TH1D*> hSigExtCorrected;    // background subtracted (signal extraction method) DeltaR distributions corrected by efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2

    // Investigation histograms
    TH1D* hBDTBackgroundScore;              // rejected entries that didn't pass the BDT background score cuts
};

// Module to create histograms including interest variable
EfficiencyData createHistograms(const std::vector<double>& binEdges) {

    // Create struct to store data
    EfficiencyData histStruct;

    // Create 3 histogram cases: inclusive = 0, prompt only = 1, non-prompt only = 2
    int histCaseNum = 3;
    for (size_t i = 0; i < histCaseNum; ++i) {
        histStruct.hMcpPt.push_back(new TH1D(Form("mcp_pt_%zu",i), ";#it{p}_{T,D}^{truth};dN/d#it{p}_{T,D}^{truth}", binEdges.size() - 1, binEdges.data()));
        histStruct.hMcpPt[i]->SetMarkerColor(kBlack);
        histStruct.hMcpPt[i]->SetLineColor(kBlack);
        histStruct.hMcpPt[i]->SetMarkerStyle(kOpenCircle);
        histStruct.hMcpPt[i]->Sumw2();
        histStruct.hMcpPt[i]->SetStats(0);
        histStruct.hMcdPt.push_back(new TH1D(Form("mcd_pt_%zu",i), ";#it{p}_{T,D}^{reco};dN/d#it{p}_{T,D}^{reco}", binEdges.size() - 1, binEdges.data()));
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
void fillHistograms(TFile* fSimulated, const EfficiencyData& histStruct, double& jetptMin, double& jetptMax, double& hfptMin, double& hfptMax, const std::vector<double>& deltaRBinEdges, const std::vector<std::pair<double, double>>& bdtPtCuts, const double& backProbabilityCut) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on jet
    const double yCut = 0.8; // on D0
    const double deltaRcut = deltaRBinEdges[deltaRBinEdges.size() - 1];
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

    tree->SetBranchAddress("fMcJetHfDist",&axisDistance);
    tree->SetBranchAddress("fMcJetPt",&jetPt);
    tree->SetBranchAddress("fMcJetEta",&jetEta);
    tree->SetBranchAddress("fMcJetPhi",&jetPhi);
    tree->SetBranchAddress("fMcJetNConst",&jetnconst_float);
    tree->SetBranchAddress("fMcHfPt",&hfPt);
    tree->SetBranchAddress("fMcHfEta",&hfEta);
    tree->SetBranchAddress("fMcHfPhi",&hfPhi);
    hfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&hfY);
    tree->SetBranchAddress("fMcHfPrompt",&hfprompt);
    tree->SetBranchAddress("fMcHfMatch",&hfmatch);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));
        bool genLevelRange = (abs(hfEta) < etaCut) && (abs(hfY) < yCut) && (jetPt >= jetptMin) && (jetPt < jetptMax) && ((deltaR >= deltaRBinEdges[0]) && (deltaR < deltaRcut)) && ((hfPt >= MCPHfPtMincut) && (hfPt < MCPHfPtMaxcut));
        
        // Fill histograms considering jet pT and detector acceptance
        if (genLevelRange) {
        //if ((abs(hfEta) < etaCut) && (jetPt >= jetptMin) && (jetPt < jetptMax)) {
            
            // Fill inclusive histogram
            histStruct.hMcpPt[0]->Fill(hfPt);
            // fill prompt efficiency histogram
            if (hfprompt) {
                histStruct.hMcpPt[1]->Fill(hfPt);
            } else{
                // fill non-prompt efficiency histogram
                histStruct.hMcpPt[2]->Fill(hfPt);
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
    int jetnconst_int, hfMatchedFrom, hfSelectedAs;
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
    tree->SetBranchAddress("fHfMatchedFrom",&hfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&hfSelectedAs);

    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // only compute matched detector level candidates, but compute all particle level ones
        bool isReflection = (hfMatchedFrom != hfSelectedAs) ? true : false;
        if (!hfmatch || isReflection) {
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
                histStruct.hMcdPt[0]->Fill(hfPt);
                // fill prompt efficiency histogram
                if (hfprompt) {
                    histStruct.hMcdPt[1]->Fill(hfPt);
                } else{
                    // fill non-prompt efficiency histogram
                    histStruct.hMcdPt[2]->Fill(hfPt);
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
        TH1D* mcdHist = histStruct.hMcdPt[iEff];
        TH1D* mcpHist = histStruct.hMcpPt[iEff];

        // Otain efficiency by dividing tempHist by parameter -> mcd/mcp
        TH1D* efficiencyHist = static_cast<TH1D*>(mcdHist->Clone(Form("efficiency_%s", efficiencyNames[iEff]))); // static needed to explicitly cast the cloned histogram to TH1D. ensures that the function returns a TH1D*, matching the return type.
        efficiencyHist->Divide(mcpHist);
        efficiencyHist->SetMarkerStyle(kFullCircle);
        efficiencyHist->SetMarkerColor(30+iEff*10); // 30 = not so bright green
        efficiencyHist->SetLineColor(30+iEff*10);
        efficiencyHist->SetTitle(";#it{p}_{T,D}^{truth};Efficiency#times Acceptance");
        efficiencyHist->SetStats(0);
        
        // Store efficiency distribution in struct
        histStruct.hEfficiencies.push_back(efficiencyHist);

    }

    cout << "Efficiency calculated.\n";
}

void PerformCorrections(EfficiencyData& histStruct, TFile* fBackSub) {
    int effOption = 1; // choice of efficiency for correction: 0 = inclusive, 1 = prompt, 2 = non-prompt
    double efficiency = 0;
    int bin;
    TH1D* hBackSub_temp = (TH1D*)fBackSub->Get("h_back_subtracted_0");

    // Initialize all pT,D histogram before adding them
    histStruct.hBackSubCorrected_allpt = (TH1D*)hBackSub_temp->Clone("hBackSubCorrected_allpt");
    histStruct.hBackSubCorrected_allpt->Reset();
    
    // Loop through distributions from all intervals
    int numHistos = HistogramCounter(fBackSub);
    for (size_t iHisto = 0; iHisto < numHistos; iHisto++) {
        // access each distribution
        hBackSub_temp = (TH1D*)fBackSub->Get(Form("h_back_subtracted_%zu",iHisto));

        
        // Get pT,D corresponding range from the title of the histogram
        TString title = hBackSub_temp->GetTitle(); // Typical title format: "1 < #it{p}_{T, D^{0}} < 2 GeV/#it{c}"
        // Keep only digits, dots, spaces, and "<"
        std::string clean;
        for (char c : std::string(title.Data())) {
            if ((c >= '0' && c <= '9') || c == '.' || c == '<' || c == ' ') {
                clean += c;
            }
        }
        // Now parse
        double ptBinMin = 0, ptBinMax = 0;
        int nParsed = sscanf(clean.c_str(), "%lf < %*s < %lf", &ptBinMin, &ptBinMax);
        //std::cout << "Parsed pT,D range: " << ptBinMin << " to " << ptBinMax << " GeV/c\n";
        if (nParsed != 2) {
            std::cerr << "Failed to parse pT,D range from histogram title: " << title << std::endl;
            ptBinMin = 0;
            ptBinMax = 0;
        }
        
        // Obtain efficiency to the correspondent pT interval (bin)
        int bin = histStruct.hEfficiencies[effOption]->FindBin((ptBinMin + ptBinMax) / 2.0); // Find the bin corresponding to the average pT,D value in the range
        efficiency = histStruct.hEfficiencies[effOption]->GetBinContent(bin);
        if (efficiency == 0) {
            //std::cout << "Warning: Prompt efficiency is zero for pT = " << MCDhfPt << " with bin " << bin << ". Setting it to 1." << std::endl;
            efficiency = 1; // Avoid division by zero
        }
        
        // apply efficiency correction to distribution
        hBackSub_temp->Scale(1 / efficiency); // Scale the histogram by the efficiency value
        
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
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Inclusive");
    cEff_1->cd(2);
    histStruct.hMcpPt[1]->Draw();
    histStruct.hMcdPt[1]->Draw("same");
    legendPt[1]->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Prompt D^{0}'s");
    cEff_1->cd(3);
    histStruct.hMcpPt[2]->Draw();
    histStruct.hMcdPt[2]->Draw("same");
    legendPt[2]->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
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
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
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
        latex->DrawLatex(statBoxPos-0.35, 0.61, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
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
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
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
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
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
    TString imagePath = "../Images/2-Efficiency/Run2Style/";
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
    TFile* outFile = new TFile(Form("run2_style_efficiency_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
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

std::vector<double> LoadBinning(TFile* fInput, const char* pathInFile) {
    auto* vec = (TVectorD*)fInput->Get(pathInFile);
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", pathInFile));
    }
    return std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
}

void Efficiency_run2_style(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.03 (from first of pair)
    //std::vector<std::pair<double, double>> bdtPtCuts = {
    //    {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}
    //};
    // Dataset: JE_HF_LHC24g5_All_D0
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.4}, {8, 0.6}, {10, 0.8}, {12, 0.8}, {16, 1.0}, {50, 1.0}
    };
    double backProbabilityCut = 1.0; // default value
    
    // Initial values for opening the file
    TFile* fAxes = new TFile(Form("../1-SignalTreatment/SideBand/full_merged_ranges_back_sub.root"),"read");
    if (!fAxes || fAxes->IsZombie()) {
        std::cerr << "Error: Unable to open simulated data ROOT file." << std::endl;
    }

    // Load pT,jet bin edges
    std::vector<double> ptjetBinEdges = LoadBinning(fAxes, "axes/ptjetBinEdges_detector");
    double jetptMin = ptjetBinEdges[0]; // GeV
    double jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 1]; // GeV
    // Load ΔR bin edges
    std::vector<double> deltaRBinEdges = LoadBinning(fAxes, "axes/deltaRBinEdges_detector");
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
    // Load pT,D bin edges
    std::vector<double> ptDBinEdges = LoadBinning(fAxes, "axes/ptDBinEdges_detector");
    double hfptMin = ptDBinEdges[0]; //ptDBinEdges[0] - should start from 0 or from the lowest pT,D value?
    double hfptMax = ptDBinEdges[ptDBinEdges.size() - 1];

    // Opening files
    //TFile* fSimulated = new TFile("../SimulatedData/Hyperloop_output/McEfficiency/New_with_reflections/Merged_AO2D_HF_LHC24d3a_All.root","read"); // previous used file
    TFile* fSimulated = new TFile("../SimulatedData/Hyperloop_output/Train_runs/410602_Eff/AO2D_mergedDFs.root","read");
    TFile* fBackSub = new TFile(Form("../1-SignalTreatment/SideBand/backSub_%d_to_%d_jetpt_with_reflections.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fSimulated || fSimulated->IsZombie()) {
        std::cerr << "Error: Unable to open simulated data ROOT file." << std::endl;
    }
    if (!fBackSub || fBackSub->IsZombie()) {
        std::cerr << "Error: Unable to open background subtracted data ROOT file." << std::endl;
    }

    EfficiencyData histStruct = createHistograms(ptDBinEdges); // pT histograms

    // Fill histograms
    fillHistograms(fSimulated, histStruct, jetptMin, jetptMax, hfptMin, hfptMax, deltaRBinEdges, bdtPtCuts, backProbabilityCut);

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
    Efficiency_run2_style();
    return 0;
}
