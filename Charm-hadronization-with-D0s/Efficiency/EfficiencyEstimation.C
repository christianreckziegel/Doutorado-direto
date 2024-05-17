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
    std::vector<TH1D*> hBackSubCorrected;
    std::vector<TH1D*> hSigExtCorrected;
};

// Module to create TH2D histograms including interest variable
EfficiencyData createHistograms(const std::vector<const char*>& names, const std::vector<const char*>& titles,
                                    int xbins, double xmin, double xmax) {
    EfficiencyData histStruct;
    for (size_t i = 0; i < names.size(); ++i) {
        histStruct.hMcpPt.push_back(new TH1D(Form("mcp_%s",names[i]), titles[i], xbins, xmin, xmax));
        histStruct.hMcpPt[i]->GetXaxis()->SetTitle("p_{T,D}^{truth}");
        histStruct.hMcpPt[i]->GetYaxis()->SetTitle("dN/dp_{T,D}^{truth}");
        histStruct.hMcpPt[i]->SetMarkerColor(kBlack);
        histStruct.hMcpPt[i]->SetLineColor(kBlack);
        histStruct.hMcpPt[i]->SetMarkerStyle(kOpenCircle);
        histStruct.hMcpPt[i]->Sumw2();
        histStruct.hMcdPt.push_back(new TH1D(Form("mcd_%s",names[i]), titles[i], xbins, xmin, xmax));
        histStruct.hMcdPt[i]->GetXaxis()->SetTitle("p_{T,D}^{reco}");
        histStruct.hMcdPt[i]->GetYaxis()->SetTitle("dN/dp_{T,D}^{reco}");
        histStruct.hMcpPt[i]->SetMarkerColor(kBlue);
        histStruct.hMcpPt[i]->SetLineColor(kBlue);
        histStruct.hMcpPt[i]->SetMarkerStyle(kFullCircle);
        histStruct.hMcdPt[i]->Sumw2();
        
    }

    // Efficiency binning corresponds to the steps in the topological cuts
    Double_t binEdges[] = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};
    histStruct.hMcpPt.push_back(new TH1D("hMcpPt", ";p_{T,D}^{reco};dN/dp_{T,D}^{truth}", sizeof(binEdges)/sizeof(binEdges[0])-1, binEdges));
    histStruct.hMcpPt[histStruct.hMcpPt.size()-1]->GetXaxis()->SetTitle("p_{T,D}^{truth}");
    histStruct.hMcpPt[histStruct.hMcpPt.size()-1]->GetYaxis()->SetTitle("dN/dp_{T,D}^{truth}");
    histStruct.hMcpPt[histStruct.hMcpPt.size()-1]->SetMarkerColor(kBlack);
    histStruct.hMcpPt[histStruct.hMcpPt.size()-1]->SetLineColor(kBlack);
    histStruct.hMcpPt[histStruct.hMcpPt.size()-1]->SetMarkerStyle(kOpenCircle);
    histStruct.hMcpPt[histStruct.hMcpPt.size()-1]->Sumw2();
    histStruct.hMcdPt.push_back(new TH1D("hMcdPt", ";p_{T,D}^{reco};dN/dp_{T,D}^{reco}", sizeof(binEdges)/sizeof(binEdges[0])-1, binEdges));
    histStruct.hMcdPt[histStruct.hMcpPt.size()-1]->GetXaxis()->SetTitle("p_{T,D}^{reco}");
    histStruct.hMcdPt[histStruct.hMcpPt.size()-1]->GetYaxis()->SetTitle("dN/dp_{T,D}^{reco}");
    histStruct.hMcdPt[histStruct.hMcdPt.size()-1]->SetMarkerColor(kBlue);
    histStruct.hMcdPt[histStruct.hMcdPt.size()-1]->SetLineColor(kBlue);
    histStruct.hMcdPt[histStruct.hMcdPt.size()-1]->SetMarkerStyle(kFullCircle);
    histStruct.hMcdPt[histStruct.hMcdPt.size()-1]->Sumw2();

    return histStruct;
}

// Module to fill histograms from TFile data
void fillHistograms(TFile* file, const EfficiencyData& histStruct, double jetptMin, double jetptMax) {
    //
    // MC generator level tree and histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)file->Get("DF_2267291774843000/O2mcpjetdisttable");

    // Check for correct access
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing the TTree
    float axisDistance, jetPt, jetEta, jetPhi;
    float hfPt, hfEta, hfPhi, hfMass, hfY;
    bool hfmatch;

    tree->SetBranchAddress("fJetHfDist",&axisDistance);
    tree->SetBranchAddress("fJetPt",&jetPt);
    tree->SetBranchAddress("fJetEta",&jetEta);
    tree->SetBranchAddress("fJetPhi",&jetPhi);
    tree->SetBranchAddress("fHfPt",&hfPt);
    tree->SetBranchAddress("fHfEta",&hfEta);
    tree->SetBranchAddress("fHfPhi",&hfPhi);
    hfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fHfY",&hfY);
    tree->SetBranchAddress("fHfMatch",&hfmatch);


    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // only compute matched candidates
        if (!hfmatch) {
            continue;
        }
        
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));

        // Fill each histogram with their respective pT intervals, consider jet pT and detector acceptance
        if ((abs(hfEta) < 0.9) && jetPt > jetptMin && jetPt < jetptMax) {
            if ((hfPt > 3 && hfPt < 4) && (histStruct.hMcpPt.size() > 0)) {
                histStruct.hMcpPt[0]->Fill(hfPt);
            } else if ((hfPt > 4 && hfPt < 5) && (histStruct.hMcpPt.size() > 1)) {
                histStruct.hMcpPt[1]->Fill(hfPt);
            } else if ((hfPt > 5 && hfPt < 6) && (histStruct.hMcpPt.size() > 2)) {
                histStruct.hMcpPt[2]->Fill(hfPt);
            } else if ((hfPt > 6 && hfPt < 7) && (histStruct.hMcpPt.size() > 3)) {
                histStruct.hMcpPt[3]->Fill(hfPt);
            } else if ((hfPt > 7 && hfPt < 8) && (histStruct.hMcpPt.size() > 4)) {
                histStruct.hMcpPt[4]->Fill(hfPt);
            } else if ((hfPt > 8 && hfPt < 10) && (histStruct.hMcpPt.size() > 5)) {
                histStruct.hMcpPt[5]->Fill(hfPt);
            } else if ((hfPt > 10 && hfPt < 12) && (histStruct.hMcpPt.size() > 6)) {
                histStruct.hMcpPt[6]->Fill(hfPt);
            } else if ((hfPt > 12 && hfPt < 15) && (histStruct.hMcpPt.size() > 7)) {
                histStruct.hMcpPt[7]->Fill(hfPt);
                cout << "hfPt = " << hfPt << endl;
            } else if ((hfPt > 15 && hfPt < 30) && (histStruct.hMcpPt.size() > 8)) {
                histStruct.hMcpPt[8]->Fill(hfPt);
                cout << "hfPt = " << hfPt << endl;
            }
            histStruct.hMcpPt[histStruct.hMcdPt.size()-1]->Fill(hfPt);
        }
        
        
    }
    cout << "Generator level histograms filled.\n";

    //
    // MC detector level tree and histograms
    //
    // Accessing TTree
    tree = (TTree*)file->Get("DF_2267291774843000/O2mcdjetdisttable");
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
    tree->SetBranchAddress("fHfMatch",&hfmatch);

    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // only compute matched candidates
        if (!hfmatch) {
            continue;
        }
        
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));

        // Fill each histogram with their respective pT intervals, consider jet pT and detector acceptance
        if ((abs(hfEta) < 0.9) && jetPt > jetptMin && jetPt < jetptMax) {
            if ((hfPt > 3 && hfPt < 4) && (histStruct.hMcdPt.size() > 0)) {
                histStruct.hMcdPt[0]->Fill(hfPt);
            } else if ((hfPt > 4 && hfPt < 5) && (histStruct.hMcdPt.size() > 1)) {
                histStruct.hMcdPt[1]->Fill(hfPt);
            } else if ((hfPt > 5 && hfPt < 6) && (histStruct.hMcdPt.size() > 2)) {
                histStruct.hMcdPt[2]->Fill(hfPt);
            } else if ((hfPt > 6 && hfPt < 7) && (histStruct.hMcdPt.size() > 3)) {
                histStruct.hMcdPt[3]->Fill(hfPt);
            } else if ((hfPt > 7 && hfPt < 8) && (histStruct.hMcdPt.size() > 4)) {
                histStruct.hMcdPt[4]->Fill(hfPt);
            } else if ((hfPt > 8 && hfPt < 10) && (histStruct.hMcdPt.size() > 5)) {
                histStruct.hMcdPt[5]->Fill(hfPt);
            } else if ((hfPt > 10 && hfPt < 12) && (histStruct.hMcdPt.size() > 6)) {
                histStruct.hMcdPt[6]->Fill(hfPt);
            } else if ((hfPt > 12 && hfPt < 15) && (histStruct.hMcdPt.size() > 7)) {
                histStruct.hMcdPt[7]->Fill(hfPt);
            } else if ((hfPt > 15 && hfPt < 30) && (histStruct.hMcdPt.size() > 8)) {
                histStruct.hMcdPt[8]->Fill(hfPt);
            }
            histStruct.hMcdPt[histStruct.hMcdPt.size()-1]->Fill(hfPt);
            //cout << "hfPt = " << hfPt << endl;

        }
        
        
    }
    cout << "Detector level histograms filled.\n";
}

TH1D* CalculateEfficiency(const EfficiencyData& histStruct){
    // Clone detector level histogram for future division
    TH1D* tempHist = static_cast<TH1D*>(histStruct.hMcdPt[histStruct.hMcdPt.size()-1]->Clone()); // static needed to explicitly cast the cloned histogram to TH1D. ensures that the function returns a TH1D*, matching the return type.

    // Is it just a ratio histogram?
    tempHist->Divide(histStruct.hMcpPt[histStruct.hMcdPt.size()-1]);
    tempHist->SetMarkerStyle(kFullCircle);
    tempHist->SetMarkerColor(30); // 30 = not so bright green
    tempHist->SetLineColor(30);
    tempHist->SetTitle(";p_{T,D}^{truth};Efficiency");
    
    return tempHist;
}

void PerformCorrections(EfficiencyData& histStruct, TH1D* hEfficiency, TFile* fBackSub) {
    double efficiency = 0;
    int bin;
    TH1D* hBackSub_temp;

    // loop through distributions from all intervals
    int numHistos = HistogramCounter(fBackSub);
    for (size_t iHisto = 0; iHisto < numHistos; iHisto++) {
        // access each distribution
        hBackSub_temp = (TH1D*)fBackSub->Get(Form("h_back_subtracted_%zu",iHisto));

        // obtain efficiency to the correspondent pT interval (bin)
        bin = hEfficiency->GetBin(iHisto);
        efficiency = hEfficiency->GetBinContent(bin);
        
        // apply efficiency correction to distribution
        hBackSub_temp->Scale(efficiency);

        // store corrected distribution to struct
        histStruct.hBackSubCorrected.push_back(hBackSub_temp);
    }
    
    
}

void PlotHistograms(const EfficiencyData& histStruct, TH1D* hEfficiency, double jetptMin, double jetptMax) {

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    TLegend* legend = new TLegend(0.65,0.49,0.8,0.62);
    legend->AddEntry(histStruct.hMcpPt[histStruct.hMcdPt.size()-1],"Truth", "lpe");
    legend->AddEntry(histStruct.hMcdPt[histStruct.hMcdPt.size()-1],"Reconstructed", "lpe");

    TCanvas* cEff = new TCanvas("cEff","Efficiency plots");
    cEff->Divide(2,2);
    // obtaining 1D invariant mass histograms from the projection
    /*for (size_t iHist = 0; iHist < histStruct.hMcpPt.size(); iHist++) {
        //
        histStruct.hMcpPt[iHist]->Draw();
    }*/
    cEff->cd(1);
    histStruct.hMcpPt[histStruct.hMcdPt.size()-1]->Draw();
    histStruct.hMcdPt[histStruct.hMcdPt.size()-1]->Draw("same");
    legend->Draw();
    double statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    
    // Plot efficiency distribution
    cEff->cd(2);
    hEfficiency->Draw();

    cEff->SetCanvasSize(1800,1000);
    cEff->Print(Form("pT_efficiency_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));

    // Plot efficiency corrected histograms
    TCanvas* cCorrectedBackSub = new TCanvas("cCorrectedBackSub","Efficiency corrected histograms");
    cCorrectedBackSub->Divide(3,static_cast<int>(histStruct.hBackSubCorrected.size() / 3));
    
    for (size_t iHisto = 0; iHisto < histStruct.hBackSubCorrected.size(); iHisto++) {
        //
        cCorrectedBackSub->cd(iHisto+1);
        histStruct.hBackSubCorrected[iHisto]->Draw();
    }
    cCorrectedBackSub->SetCanvasSize(1800,1000);
    cCorrectedBackSub->Print(Form("pT_efficiency_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
}

void SaveData(const EfficiencyData& histStruct, double jetptMin, double jetptMax){
    // Open output file
    TFile* outFile = new TFile(Form("backSubCorrected_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    // Loop over signal extracted histograms
    for (size_t iHisto = 0; iHisto < histStruct.hBackSubCorrected.size(); iHisto++) {
        // store each histogram in file
        histStruct.hBackSubCorrected[iHisto]->Write();
    }
    outFile->Close();
    delete outFile;
    
    
}

void EfficiencyEstimation(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;
    // jet pT cuts
    double jetptMin = 5; // GeV
    double jetptMax = 30; // GeV
    // deltaR histogram
    int deltaRbins = 10; // deltaRbins = numberOfPoints, default=100 bins for [0. 1.0]
    double minDeltaR = 0.;
    double maxDeltaR = 0.4;
    // mass histogram
    int massBins = 100; 
    double minMass = 1.67;
    double maxMass = 2.1;
    // pT histograms
    int ptBins = 100;
    double minPt = 0.;
    double maxPt = 30.;

    
    std::vector<const char*> names = {"histPt1", "histPt2", "histPt3", 
                                      "histPt4", "histPt5", "histPt6",
                                      "histPt7", "histPt8", "histPt9"}; // Names of histograms
    std::vector<const char*> titles = {"3 < p_{T,D} < 4 GeV/c", "4 < p_{T,D} < 5 GeV/c", "5 < p_{T,D} < 6 GeV/c",
                                       "6 < p_{T,D} < 7 GeV/c", "7 < p_{T,D} < 8 GeV/c", "8 < p_{T,D} < 10 GeV/c",
                                       "10 < p_{T,D} < 12 GeV/c", "12 < p_{T,D} < 15 GeV/c", "15 < p_{T,D} < 30 GeV/c"}; // Titles of histograms
    EfficiencyData histStruct = createHistograms(names, titles, 
                                                 ptBins, minPt, maxPt); // pT histograms

    // opening files
    TFile* fSimulated = new TFile("../SimulatedData/AnalysisResults_trees_original.root","read");
    TFile* fBackSub = new TFile(Form("../SignalTreatment/backSub_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fSigExt = new TFile(Form("../SignalTreatment/sigExt_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fSimulated || fSimulated->IsZombie()) {
        std::cerr << "Error: Unable to open simulated data ROOT file." << std::endl;
    }
    if (!fBackSub || fBackSub->IsZombie()) {
        std::cerr << "Error: Unable to open background subtracted data ROOT file." << std::endl;
    }
    if (!fSigExt || fSigExt->IsZombie()) {
        std::cerr << "Error: Unable to open signal extracted ROOT file." << std::endl;
    }
    
    

    // Fill histograms
    fillHistograms(fSimulated, histStruct, jetptMin, jetptMax);

    // Calculate efficiency distribution
    TH1D* hEfficiency = CalculateEfficiency(histStruct);

    // Apply efficiency corrections to distributions
    PerformCorrections(histStruct, hEfficiency, fBackSub);

    // Plot the efficiency histogram and further corrected histograms
    PlotHistograms(histStruct, hEfficiency, jetptMin, jetptMax);

    // Save corrected distributions to file
    SaveData(histStruct, jetptMin, jetptMax);


    

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    EfficiencyEstimation();
    return 0;
}
