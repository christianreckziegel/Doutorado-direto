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


using namespace std;

// calculate number of objects inside file
int HistogramCounter(TFile* file) {
    // Get the list of keys (i.e., objects) in the file
    TList* keys = file->GetListOfKeys();

    int numHistograms = 0;

    // Loop over all keys in the file
    for (int i = 0; i < keys->GetSize(); ++i) {
        TKey* key = (TKey*)keys->At(i);
        TObject* obj = file->Get(key->GetName());

        // Check if the object is a histogram
        if (obj->IsA()->InheritsFrom(TH1::Class())) {
            numHistograms++;
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
    std::vector<TH1D*> hMcpPt;
    std::vector<TH1D*> hMcdPt;
    std::vector<TH1D*> hEfficiencies;       // inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hBackSubCorrected;
    TH1D* hBackSubCorrected_allpt;          // all pT.D bins histograms summed into one
    std::vector<TH1D*> hSigExtCorrected;
};

// Module to create histograms including interest variable
EfficiencyData createHistograms(const std::vector<double>& binEdges) {

    // Create struct to store data
    EfficiencyData histStruct;

    // Create 3 histogram cases: inclusive = 0, prompt only = 1, non-prompt only = 2
    int histCaseNum = 3;
    for (size_t i = 0; i < histCaseNum; ++i) {
        histStruct.hMcpPt.push_back(new TH1D(Form("mcp_pt_%zu",i), ";p_{T,D}^{truth};dN/dp_{T,D}^{truth}", binEdges.size() - 1, binEdges.data()));
        histStruct.hMcpPt[i]->SetMarkerColor(kBlack);
        histStruct.hMcpPt[i]->SetLineColor(kBlack);
        histStruct.hMcpPt[i]->SetMarkerStyle(kOpenCircle);
        histStruct.hMcpPt[i]->Sumw2();
        histStruct.hMcpPt[i]->SetStats(0);
        histStruct.hMcdPt.push_back(new TH1D(Form("mcd_pt_%zu",i), ";p_{T,D}^{reco};dN/dp_{T,D}^{reco}", binEdges.size() - 1, binEdges.data()));
        histStruct.hMcdPt[i]->SetMarkerColor(kBlue);
        histStruct.hMcdPt[i]->SetLineColor(kBlue);
        histStruct.hMcdPt[i]->SetMarkerStyle(kFullCircle);
        histStruct.hMcdPt[i]->Sumw2();
        histStruct.hMcdPt[i]->SetStats(0);
    }

    cout << "Histograms created.\n";

    return histStruct;
}

// Module to fill histograms from TFile data
void fillHistograms(TFile* fSimulated, const EfficiencyData& histStruct, double jetptMin, double jetptMax) {

    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on jet
    const double yCut = 0.8; // on D0
    const double deltaRcut = 0.4;

    //
    // MC generator level tree and histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)fSimulated->Get("O2mcpjetdisttable");

    // Check for correct access
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing the TTree
    float axisDistance, jetPt, jetEta, jetPhi;
    float hfPt, hfEta, hfPhi, hfMass, hfY;
    bool hfprompt, hfmatch;

    tree->SetBranchAddress("fJetHfDist",&axisDistance);
    tree->SetBranchAddress("fJetPt",&jetPt);
    tree->SetBranchAddress("fJetEta",&jetEta);
    tree->SetBranchAddress("fJetPhi",&jetPhi);
    tree->SetBranchAddress("fHfPt",&hfPt);
    tree->SetBranchAddress("fHfEta",&hfEta);
    tree->SetBranchAddress("fHfPhi",&hfPhi);
    hfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fHfY",&hfY);
    tree->SetBranchAddress("fHfPrompt",&hfprompt);
    tree->SetBranchAddress("fHfMatch",&hfmatch);


    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));

        // Fill histograms considering jet pT and detector acceptance
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && (jetPt >= jetptMin) && (jetPt < jetptMax) && (deltaR >= 0.) && (deltaR < deltaRcut)) {
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
    tree = (TTree*)fSimulated->Get("O2mcdjetdisttable");
    // Check for correct access
    if (!tree) {
        cout << "Error opening tree.\n";
    }

    tree->SetBranchAddress("fJetHfDist",&axisDistance);
    tree->SetBranchAddress("fJetPt",&jetPt);
    tree->SetBranchAddress("fJetEta",&jetEta);
    tree->SetBranchAddress("fJetPhi",&jetPhi);
    tree->SetBranchAddress("fHfPt",&hfPt);
    tree->SetBranchAddress("fHfEta",&hfEta);
    tree->SetBranchAddress("fHfPhi",&hfPhi);
    tree->SetBranchAddress("fHfMass",&hfMass);
    tree->SetBranchAddress("fHfY",&hfY);
    tree->SetBranchAddress("fHfPrompt",&hfprompt);
    tree->SetBranchAddress("fHfMatch",&hfmatch);

    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // only compute matched detector level candidates, but compute all particle level ones
        if (!hfmatch) {
            continue;
        }
        
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));

        // Fill histograms considering jet pT and detector acceptance
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && (jetPt > jetptMin) && (jetPt < jetptMax)  && (deltaR >= 0.) && (deltaR < deltaRcut)) {
            // Fill inclusive histogram
            histStruct.hMcdPt[0]->Fill(hfPt);
            // fill prompt efficiency histogram
            if (hfprompt) {
                histStruct.hMcdPt[1]->Fill(hfPt);
            } else{
                // fill non-prompt efficiency histogram
                histStruct.hMcdPt[2]->Fill(hfPt);
            }
            
        }
        
        
    }
    cout << "Detector level histograms filled.\n";
}

void CalculateEfficiency(EfficiencyData& histStruct, double& jetptMin, double& jetptMax){
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
        efficiencyHist->SetTitle(Form("%0.f < p_{T,jet} < %0.f;p_{T,D}^{truth};Efficiency#times Acceptance",jetptMin,jetptMax));
        efficiencyHist->SetStats(0);
        
        // Store efficiency distribution in struct
        histStruct.hEfficiencies.push_back(efficiencyHist);

    }

    cout << "Efficiency calculated.\n";
}

void PerformCorrections(EfficiencyData& histStruct, TFile* fBackSub, double& jetptMin, double& jetptMax) {
    int effOption = 1; // choice of efficiency for correction: 0 = inclusive, 1 = prompt, 2 = non-prompt
    double efficiency = 0;
    int bin;
    TH1D* hBackSub_temp = (TH1D*)fBackSub->Get(Form("JetPtRange_%.0f_%.0f/h_back_subtracted_0", jetptMin, jetptMax));

    // Initialize all pT,D histogram before adding them
    histStruct.hBackSubCorrected_allpt = (TH1D*)hBackSub_temp->Clone("hBackSubCorrected_allpt");
    histStruct.hBackSubCorrected_allpt->SetTitle(Form("Corrected all p_{T,D} interval for %0.f < p_{T,jet} < %0.f", jetptMin, jetptMax));
    histStruct.hBackSubCorrected_allpt->Reset();

    // Loop through distributions from all intervals
    int numHistos = HistogramCounter(fBackSub);
    for (size_t iHisto = 0; iHisto < numHistos; iHisto++) {
        // access each distribution
        hBackSub_temp = (TH1D*)fBackSub->Get(Form("JetPtRange_%.0f_%.0f/h_back_subtracted_%zu", jetptMin, jetptMax, iHisto));

        // obtain efficiency to the correspondent pT interval (bin)
        bin = histStruct.hEfficiencies[effOption]->GetBin(iHisto);
        efficiency = histStruct.hEfficiencies[effOption]->GetBinContent(bin);
        
        // apply efficiency correction to distribution
        hBackSub_temp->Scale(1/efficiency);
        
        // modifying visual characteristics
        hBackSub_temp->SetStats(0);
        hBackSub_temp->GetXaxis()->SetTitleSize(0.05);
        hBackSub_temp->GetYaxis()->SetTitleSize(0.05);

        // store corrected distribution to struct
        histStruct.hBackSubCorrected.push_back(hBackSub_temp);

        // add to final distribution of all pT,D's
        histStruct.hBackSubCorrected_allpt->Add(histStruct.hBackSubCorrected[iHisto]);
    }
    
    cout << "Efficiency corrections applied.\n";
}

void PlotHistograms(const EfficiencyData& histStruct, double jetptMin, double jetptMax) {
    cout << "Plotting histograms...\n";

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    TLegend* legendPt = new TLegend(0.65,0.49,0.8,0.62);
    legendPt->AddEntry(histStruct.hMcpPt[0],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    legendPt->AddEntry(histStruct.hMcdPt[0],"Reconstructed", "lpe");

    
    //
    // 1 - plot pT,D of generator and reconstruction level
    //
    TCanvas* cEff_1 = new TCanvas("cEff_1","Efficiency plots");
    cEff_1->SetCanvasSize(1800,1000);
    cEff_1->cd();
    histStruct.hMcpPt[0]->Draw();
    histStruct.hMcdPt[0]->Draw("same");
    legendPt->Draw();
    double statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    
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
            histStruct.hEfficiencies[iEff]->SetMaximum(0.57);
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
    cCorrectedBackSub->Divide(3,static_cast<int>(histStruct.hBackSubCorrected.size() / 3));
    
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
    // Storing images
    //
    TString imagePath = "../Images/2-Efficiency/";
    cEff_1->Update();
    cEff_1->SaveAs(imagePath + "Efficiency_dNdpT.png");
    cEff_2->Update();
    cEff_2->SaveAs(imagePath + "Efficiency_acceptance.png");
    cCorrectedBackSub->Update();
    cCorrectedBackSub->SaveAs(imagePath + "Efficiency_corrected_pT_bins.png");
    cCorrectedBackSub_allpt->Update();
    cCorrectedBackSub_allpt->SaveAs(imagePath + "Efficiency_corrected.png");

    //
    // Storing in a single pdf file
    //
    cEff_1->Print(Form("pT_efficiency_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cEff_2->Print(Form("pT_efficiency_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cCorrectedBackSub->Print(Form("pT_efficiency_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cCorrectedBackSub_allpt->Print(Form("pT_efficiency_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

}
void UpdateFinalDistribution(std::vector<TH2D*>& hptD0_vs_ptJet, EfficiencyData& histStruct, double& jetptMin) {
    
    // Loop through all efficiency cases: 0 = inclusive, 1 = prompt, 2 = non-prompt
    for (int iEff = 0; iEff < 3; iEff++) {
        // Get the bin index for the current jet pt range in the final 2D histogram
        int binY = hptD0_vs_ptJet[iEff]->GetYaxis()->FindBin(jetptMin);
        
        // Loop over all deltaR bins and update the final 2D histogram
        for (int binX = 1; binX <= histStruct.hEfficiencies[iEff]->GetNbinsX(); binX++) {
            if (!std::isnan(histStruct.hEfficiencies[iEff]->GetBinContent(binX))) {
                double effContent = histStruct.hEfficiencies[iEff]->GetBinContent(binX);
                double effError = histStruct.hEfficiencies[iEff]->GetBinError(binX);
                if (effContent < 0.) {
                    std::cout << "WARNING: negative entry value on 2D final histogram.\n";
                    std::cout << "bin (" << binX << "," << binY << ") content = " << effContent << std::endl << std::endl;
                }
                
                
                // Update the final 2D histogram with the content and error
                hptD0_vs_ptJet[iEff]->SetBinContent(binX, binY, effContent);
                hptD0_vs_ptJet[iEff]->SetBinError(binX, binY, effError);
            }
            
            
        }
    }
    
    
}

void SaveData(TFile* fOut, EfficiencyData& histStruct, double& jetptMin, double& jetptMax){
    
    // Create a directory name for the current pT range
    const std::string dirName = Form("JetPtRange_%.0f_%.0f", jetptMin, jetptMax);

    fOut->cd();
    TDirectory* dir = fOut->GetDirectory(dirName.c_str());

    // Create the directory in case it doesn't exist yet
    if (!dir) {
        dir = fOut->mkdir(dirName.c_str());
    }

    // Enter directory
    dir->cd();

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
    // Final all pT,D summed distribution
    histStruct.hBackSubCorrected_allpt->Write();

    // Return to the root directory
    fOut->cd();
    
    std::cout << "Iteration data stored on ROOT output file.\n\n";
    
}

void ClearData(EfficiencyData& histStruct) {

    // Delete data from data container struct

    for (auto& hist : histStruct.hMcpPt) {
        delete hist;
    }
    histStruct.hMcpPt.clear();

    for (auto& hist : histStruct.hMcdPt) {
        delete hist;
    }
    histStruct.hMcdPt.clear();

    for (auto& hist : histStruct.hEfficiencies) {
        delete hist;
    }
    histStruct.hEfficiencies.clear();

    for (auto& hist : histStruct.hBackSubCorrected) {
        delete hist;
    }
    histStruct.hBackSubCorrected.clear();

    if (histStruct.hBackSubCorrected_allpt != nullptr) {
        delete histStruct.hBackSubCorrected_allpt;
        histStruct.hBackSubCorrected_allpt = nullptr;
    }

    std::cout << "Data cleared for next iteration.\n\n";
}

void AnalyzeJetPtRange(TFile* fSimulated, TFile* fBackSub, std::vector<TH2D*>& hptD0_vs_ptJet, const std::vector<double>& ptDBinEdges, const std::vector<double>& deltaRBinEdges, double& jetptMin, double& jetptMax, TFile* fOut) {
    
    // Create data container
    EfficiencyData histStruct = createHistograms(ptDBinEdges); // pT histograms

    // Fill histograms
    fillHistograms(fSimulated, histStruct, jetptMin, jetptMax);

    // Calculate efficiency distribution
    CalculateEfficiency(histStruct, jetptMin, jetptMax);

    // Apply efficiency corrections to distributions
    PerformCorrections(histStruct, fBackSub, jetptMin, jetptMax);

    // Plot the efficiency histogram and further corrected histograms
    PlotHistograms(histStruct, jetptMin, jetptMax);

    // Save range analysis to 2D DeltaR vs pT,jet distribution
    UpdateFinalDistribution(hptD0_vs_ptJet, histStruct, jetptMin);

    // Save corrected distributions to file
    SaveData(fOut, histStruct, jetptMin, jetptMax);

    // Clean current iteration data in order to avoid memory leaks
    ClearData(histStruct);
}

void EfficiencyEstimation_2D(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // jet pT cuts
    std::vector<double> ptjetBinEdges = {5., 7., 15., 30.}; // GeV/c
    // deltaR histogram
    std::vector<double> deltaRBinEdges = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4}; // chosen by Nima
    // pT,D histograms
    std::vector<double> ptDBinEdges = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.}; // pT,D in GeV/c

    // Opening files
    TFile* fSimulated = new TFile("../SimulatedData/Hyperloop_output/McEfficiency/AO2D_merged_All.root","read"); //../SimulatedData/Hyperloop_output/AO2D_merged.root
    TFile* fBackSub = new TFile(Form("../1-SignalTreatment/SideBand/backSub_%d_to_%d_jetpt_2D.root",static_cast<int>(ptjetBinEdges[0]),static_cast<int>(ptjetBinEdges[ptjetBinEdges.size() - 1])),"read");
    TFile* fSigExt = new TFile(Form("../1-SignalTreatment/SignalExtraction/sigExt_%d_to_%d_jetpt.root",static_cast<int>(ptjetBinEdges[0]),static_cast<int>(ptjetBinEdges[ptjetBinEdges.size() - 1])),"read");
    if (!fSimulated || fSimulated->IsZombie()) {
        std::cerr << "Error: Unable to open simulated data ROOT file." << std::endl;
    }
    if (!fBackSub || fBackSub->IsZombie()) {
        std::cerr << "Error: Unable to open background subtracted data ROOT file." << std::endl;
    }
    if (!fSigExt || fSigExt->IsZombie()) {
        std::cerr << "Error: Unable to open signal extracted ROOT file." << std::endl;
    }
    TFile* fOut = new TFile(Form("backSubEfficiency_%d_to_%d_jetpt_2D.root",static_cast<int>(ptjetBinEdges[0]),static_cast<int>(ptjetBinEdges[ptjetBinEdges.size() - 1])),"recreate");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "Error: Unable to create the output ROOT file." << std::endl;
    }

    // Create final histograms with multiple pT,jet ranges
    std::vector<TH2D*> hptD0_vs_ptJet; // choice of efficiency for correction: 0 = inclusive, 1 = prompt, 2 = non-prompt
    hptD0_vs_ptJet.emplace_back(new TH2D("hptD0_vs_ptJet_inclusive","Inclusive efficiency;p_{T,D}^{part};p_{T,jet}_{part}", ptDBinEdges.size() -1, ptDBinEdges.data(), ptjetBinEdges.size() - 1, ptjetBinEdges.data()));
    hptD0_vs_ptJet.emplace_back(new TH2D("hptD0_vs_ptJet_prompt","Prompt efficiency;p_{T,D}^{part};p_{T,jet}_{part}", ptDBinEdges.size() -1, ptDBinEdges.data(), ptjetBinEdges.size() - 1, ptjetBinEdges.data()));
    hptD0_vs_ptJet.emplace_back(new TH2D("hptD0_vs_ptJet_nonprompt","Non-prompt efficiency;p_{T,D}^{part};p_{T,jet}_{part}", ptDBinEdges.size() -1, ptDBinEdges.data(), ptjetBinEdges.size() - 1, ptjetBinEdges.data()));

    // Loop over all the pT,jet ranges
    for (size_t iPtRange = 0; iPtRange < ptjetBinEdges.size() - 1; iPtRange++) {
        AnalyzeJetPtRange(fSimulated, fBackSub, hptD0_vs_ptJet, ptDBinEdges, deltaRBinEdges, ptjetBinEdges[iPtRange], ptjetBinEdges[iPtRange + 1], fOut);
    }

    // Store final histogram to output file
    hptD0_vs_ptJet[0]->Write();
    hptD0_vs_ptJet[1]->Write();
    hptD0_vs_ptJet[2]->Write();
    fOut->Close();
    delete fOut;

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    EfficiencyEstimation_2D();
    return 0;
}
