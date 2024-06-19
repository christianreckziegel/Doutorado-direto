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

struct FeedDownData {
    std::vector<RooUnfoldResponse> response; // inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH3D*> hPowheg; // deltaR vs pT,D vs pT,jet; [0] = generator level, [1] = treated but not yet detector level
    std::vector<TH2D*> hAllptDPowheg; // deltaR vs pT,D; [0] = generator level, [1] = folded detector level
    std::vector<TH1D*> hEfficiencies; // inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hBackSubCorrected; // prompt efficiency corrected Delta R distributions
    std::vector<TH1D*> hSBFeedDown; // non-prompt subtracted Delta R distributions
    // remove the above objects
    std::vector<RooUnfoldResponse> responseNonPrompt;
    std::vector<TH1D*> hTrueNonPrompt;
    std::vector<TH1D*> hDetectorNonPrompt;
    std::vector<TH1D*> hSigExtCorrected;
};

// Module to create TH2D histograms including interest variable
FeedDownData createHistograms(const size_t& xNumBinEdges, const double& xBinEdges[], const size_t& yNumBinEdges, const double& yBinEdges[], const double& jetptMin, const double& jetptMax) {
    // Create struct to store data
    FeedDownData dataContainer;

    // Create 2D histogram
    int ptjetBins = 25;
    dataContainer.hPowheg[0] = new TH3D("h_deltaR_vs_pt", ";#Delta R;p_{T,D}^{gen};p_{T,jet}^{ch}", xNumBinEdges-1, xBinEdges, yNumBinEdges-1, yBinEdges, ptjetBins, jetptMin, jetptMax);
    dataContainer.hPowheg[0]->SetMarkerColor(kGreen);
    dataContainer.hPowheg[0]->SetLineColor(kGreen);
    dataContainer.hPowheg[0]->SetMarkerStyle(kCircle);
    dataContainer.hPowheg[0]->Sumw2();
    dataContainer.hPowheg[0]->SetStats(0);

    cout << "Histograms created.\n";

    return dataContainer;
}

// Module to fill histograms from POWHEG+PYTHIA TFile data (all data in file is non-prompt and on particle level)
void fillHistograms(TFile* fPowheg, const FeedDownData& dataContainer, double jetptMin, double jetptMax) {
    // Defining cuts
    double jetRadius = 0.4;
    double etaCut = 0.9 - jetRadius; // on jet
    double yCut = 0.8; // on D0
    //
    // MC generator level tree and histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)fPowheg->Get("tree_D0");

    // Check for correct access
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    
    // defining variables for accessing data on TTree
    float pt_cand, eta_cand, phi_cand, y_cand;
    float pt_jet, eta_jet, phi_jet, delta_r_jet;

    tree->SetBranchAddress("pt_cand",pt_cand);
    tree->SetBranchAddress("eta_cand",eta_cand);
    tree->SetBranchAddress("phi_cand",phi_cand);
    tree->SetBranchAddress("y_cand",y_cand);
    tree->SetBranchAddress("pt_jet",pt_jet);
    tree->SetBranchAddress("eta_jet",eta_jet);
    tree->SetBranchAddress("phi_jet",phi_jet);
    tree->SetBranchAddress("delta_r_jet",delta_r_jet);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        //double deltaR = sqrt(pow(eta_jet-eta_cand,2) + pow(DeltaPhi(phi_jet,phi_cand),2));

        // Fill 2D histogram considering jet pT and detector acceptance
        if ((abs(eta_cand) < etaCut) && (abs(y_cand) < yCut) && (pt_jet > jetptMin) && (pt_jet < jetptMax)) {
            
            dataContainer.hPowheg[0]->Fill(delta_r_jet, pt_cand, pt_jet);
            
        }
        
    }
    cout << "Generator level (POWHEG+PYTHIA) histograms filled.\n";

}

// Module to create response matrices
void createResponseMatrix(FeedDownData& dataContainer, size_t& xNumBinEdges, const double& xBinEdges[], size_t& yNumBinEdges, const double& yBinEdges[]) {
    // Create lambda function for repetitive operation
    auto addResponse = [&](std::vector<RooUnfoldResponse>& container) {
        // 4D response matrix for Delta R and pT,jet
        container.emplace_back(xNumBinEdges, xBinEdges, yNumBinEdges, yBinEdges, xNumBinEdges, xBinEdges, yNumBinEdges, yBinEdges); //(DeltaR_detector, pTjet_detector, DeltaR_particle, pTjet_particle)
    };

    // Create response matrix of inclusive, prompt and non-prompt for overall pT,D
    addResponse(dataContainer.response); // Inclusive
    addResponse(dataContainer.response); // Prompt
    addResponse(dataContainer.response); // Non-prompt
    
    std::cout << "Response matrices created.\n";
}

// Module to build response matrix
void fillResponseMatrix(TFile* fSimulatedO2, const FeedDownData& dataContainer, double jetptMin, double jetptMax){
    // Defining cuts
    double jetRadius = 0.4;
    double etaCut = 0.9 - jetRadius; // on jet
    double yCut = 0.8; // on D0

    // Accessing TTree
    TTree* tree = (TTree*)fSimulatedO2->Get("O2matchtable");
    // Check for correct access
    if (!tree) {
        cout << "Error opening tree.\n";
    }

    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt, MCPhfmatch;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;

    // particle level branches
    tree->SetBranchAddress("fMCJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMCJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMCJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMCJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMCHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMCHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMCHfPhi",&MCPhfPhi);
    MCPhfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMCHfY",&hfY);
    tree->SetBranchAddress("fMCHfPrompt",&MCPhfprompt);
    tree->SetBranchAddress("fMCHfMatch",&MCPhfmatch);
    // detector level branches
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&hfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);
    tree->SetBranchAddress("fHfMatch",&MCDhfmatch);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        // Fill histograms considering jet pT and detector acceptance
        if ((abs(MCPhfEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && (MCPjetPt > MCPjetptMin) && (MCPjetPt < MCPjetptMax)) {
            // fill inclusive histogram
            dataContainer.response[0].Fill(MCDDeltaR, MCDjetPt, MCPDeltaR, MCPjetPt); // Fill(measured, true)

            if (MCPhfprompt) {
                // fill prompt efficiency histogram
                dataContainer.response[1].Fill(MCDDeltaR, MCDjetPt, MCPDeltaR, MCPjetPt); // Fill(measured, true)
            } else{
                // fill non-prompt efficiency histogram
                dataContainer.response[2].Fill(MCDDeltaR, MCDjetPt, MCPDeltaR, MCPjetPt); // Fill(measured, true)
            }
        
        }
        
    }

    cout << "Response matrices filled.\n";
}


// Module for folding particle level data from POWHEG simulation
void smearGeneratorData(FeedDownData& dataContainer, double& luminosity, TFile* fEfficiency, std::vector<const char*>& titles) {
    //
    // 0th step: clone 3D histogram in order to save the smeared at the end (while keeping the original)
    //
    TH3D* hTreatedPowheg = static_cast<TH3D*>(dataContainer.hPowheg[0]->Clone("h_subtraction_Bs"));
    hTreatedPowheg->SetTitle(";#frac{1}{L_{int}}#Delta R^{b #rightarrow D^{0}};#frac{#epsilon_{non-prompt}}{prompt}#frac{1}{L_{int}}p_{T,D}^{b #rightarrow D^{0}};p_{T,jet}^{ch}");

    //
    // 1st step: scale by integrated luminosity
    // (skip this step for now)
    //hTreatedPowheg->Scale(1/luminosity);

    //
    // 2nd step: scale by efficiency ratio
    //
    TH1D* hEffPrompt;
    TH1D* hEffNonPrompt;
    hEffPrompt = (TH1D*)fEfficiency->Get("efficiency_prompt");
    hEffNonPrompt = (TH1D*)fEfficiency->Get("efficiency_nonprompt");
    // manually scaling over each bin of pT,D dimension of 2D histogram (y-axis)
    int nBinsX = hTreatedPowheg->GetNbinsX();
    int nBinsY = hTreatedPowheg->GetNbinsY();
    int nBinsZ = hTreatedPowheg->GetNbinsZ();
    // the order of nested loop dictates which is the scaled axis
    for (int xBin = 1; xBin <= nBinsX; xBin++) {
        for (int zBin = 1; zBin < nBinsZ; zBin++) {
            // Inner loop over y-axis bins: this is scaled axis (pT,D bins)
            for (int yBin = 1; yBin <= nBinsY; yBin++) {
                double binContent = hTreatedPowheg->GetBinContent(xBin, yBin, zBin);
                double binError = hTreatedPowheg->GetBinError(xBin, yBin, zBin);

                // Obtain efficiencies ratio for specific bin
                double effPrompt = hEffPrompt->GetBinContent(yBin);
                double effNonPrompt = hEffNonPrompt->GetBinContent(yBin);

                // Scale the y-axis content and error by non-prompt efficiency/prompt efficiency
                hTreatedPowheg->SetBinContent(xBin, yBin, zBin, binContent * effNonPrompt / effPrompt);
                hTreatedPowheg->SetBinError(xBin, yBin, zBin, binError * effNonPrompt / effPrompt);
            }
        }
        
        
        
    }
    // Save modified 3D histogram
    dataContainer.hPowheg.emplace_back(hTreatedPowheg);
    
    
    //
    // 3rd step: obtain the Delta R vs pT,jet projection for all pT,D bins
    //
    TH2D* hAllptDPow = hTreatedPowheg->ProjectionXZ("h_true_nonprompt", 1, nBinsY);
    hAllptDPow->SetTitle(";#frac{1}{L_{int}}#Delta R^{b #rightarrow D^{0}};p_{T,jet}^{ch}");
    dataContainer.hAllptDPowheg.emplace_back(hAllptDPow);

    //
    // 4th step: fold Delta R vs pT,jet distribution using detector response matrices of non-prompt D0 jets
    // (this can be done by ApplyToTruth() method from a RooUnfoldResponse object or by matrix multiplication - compare methods in the future)
    // by AppluTruth() method
    int responseOption = 2; // inclusive = 0, prompt only = 1, non-prompt only = 2
    TH2D* hFoldedNonPrompt = (TH2D*)dataContainer.response[responseOption].ApplyToTruth(hAllptDPow);
    dataContainer.hAllptDPowheg.emplace_back(hFoldedNonPrompt);
    
    // by matrix multiplication
    // since TH2D objects can't be used for 4D object, each bin needs to be treated directly

    
    std::cout << "Generator data smeared.\n";
}

// Module to subtract non-prompt D0 jets from prompt efficiency corrected distribution
void feedDown(FeedDownData& dataContainer, const double jetptMin, const double jetptMax, std::vector<const char*>& titles) {
    // Temporary histogram for data treatment and storage
    TH1D* h_back_subtracted;
    
    // Loop through Delta R distributions
    for (size_t iHisto = 0; iHisto < count; iHisto++) {
        // Cloning prompt efficiency corrected hitogram from file
        h_back_subtracted = (TH1D*)fEfficiency->Get(Form("h_back_subtracted_%zu",iHisto));

        // Subtracting detector non-prompt distribution from prompt efficiency corrected
        h_back_subtracted->Add(dataContainer.hDetectorNonPrompt[iHisto],-1.0);
        dataContainer.hSBFeedDown.emplace_back(h_back_subtracted);
        dataContainer.hSBFeedDown[iHisto]->SetTitle(titles[iHisto]);
        dataContainer.hSBFeedDown[iHisto]->SetLineColor(kGrey);
        dataContainer.hSBFeedDown[iHisto]->SetMarkerColor(kGrey);
        dataContainer.hSBFeedDown[iHisto]->SetMarkerStyle(kOpenCircle);
        dataContainer.hSBFeedDown[iHisto]->SetStats(0);
    }
    
    std::cout << "Feed-down subtracted from experimental data.\n";
}


void plotHistograms(const FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax) {
    cout << "Plotting histograms...\n";

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    //
    // 2D histograms
    //
    TCanvas* cPowheg = new TCanvas("cPowheg","POWHEG data");
    cPowheg->SetCanvasSize(1800,1000);
    cPowheg->Divide(2,2);
    cPowheg->cd(1);
    dataContainer.hPowheg[0]->Draw("colz");
    cPowheg->cd(2);
    dataContainer.hPowheg[1]->Draw("colz");

    //
    // Response matrix representations in 2D histograms
    //
    TCanvas* cResponse = new TCanvas("cResponse","Response matrices for all pT,D bins");
    cResponse->SetCanvasSize(1800,1000);
    cResponse->Divide(3,static_cast<int>(dataContainer.responseNonPrompt.size() / 3));
    for (size_t iBin = 0; iBin < dataContainer.responseNonPrompt.size(); iBin++) {
        cResponse->cd(iBin+1);
        TH2D* hResponse = dataContainer.responseNonPrompt[iBin].Hresponse();
        hResponse->SetTitle(Form("%s;#Delta R_{reco}^{b #rightarrow D^{0}};#Delta R_{truth}^{b #rightarrow D^{0}}",titles[iBin]));
        hResponse->Draw("colz");
    }
    

    //
    // True prompt and non-prompt Delta R distribution for each pT,D bin
    //
    TLegend* legendNonPrompt = new TLegend(0.65,0.49,0.8,0.62);
    legendNonPrompt->AddEntry(dataContainer.hTrueNonPrompt[0],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    legendNonPrompt->AddEntry(dataContainer.hDetectorNonPrompt[0],"Reconstructed", "lpe");
    TCanvas* cNonPrompt = new TCanvas("cNonPrompt","Non-prompt Delta R plots");
    cNonPrompt->SetCanvasSize(1800,1000);
    cNonPrompt->Divide(3,static_cast<int>(dataContainer.hTrueNonPrompt.size() / 3));
    
    for (size_t iBin = 0; iBin < dataContainer.hTrueNonPrompt.size(); iBin++) {
        cNonPrompt->cd(iBin+1);
        dataContainer.hTrueNonPrompt[iBin]->Draw();
        dataContainer.hDetectorNonPrompt[iBin]->Draw("same");
        legendNonPrompt->Draw();
        double statBoxPos = gPad->GetUxmax();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
        latex->DrawLatex(statBoxPos-0.35, 0.75, "Folded non-prompt D0 jets");
    }

    //
    // Feed-down subtracted Delta R distributions
    //
    TCanvas* cSBFeedDown = new TCanvas("cSBFeedDown","Non-prompt subtracted Delta R plots");
    cSBFeedDown->SetCanvasSize(1800,1000);
    cSBFeedDown->Divide(3,static_cast<int>(dataContainer.hTrueNonPrompt.size() / 3));

    for (size_t iBin = 0; iBin < dataContainer.hTrueNonPrompt.size(); iBin++) {
        cSBFeedDown->cd(iBin+1);
        dataContainer.hSBFeedDown[iBin]->Draw();
        double statBoxPos = gPad->GetUxmax();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
        latex->DrawLatex(statBoxPos-0.35, 0.75, "Feed-down subtracted");
    }


    cPowheg->Print(Form("pT_feeddown_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cNonPrompt->Print(Form("pT_feeddown_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponse->Print(Form("pT_feeddown_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cSBFeedDown->Print(Form("pT_feeddown_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

}

void saveData(const FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax){
    // Open output file
    TFile* outFile = new TFile(Form("backSubFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    // Loop over signal extracted histograms
    for (size_t iHisto = 0; iHisto < dataContainer.hSBFeedDown.size(); iHisto++) {
        // store each histogram in file
        dataContainer.hSBFeedDown[iHisto]->Write();
    }
    outFile->Close();
    delete outFile;
    
    cout << "Data stored in file" << Form("backSubFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)) << endl;
}

void FeedDownSubtraction(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Luminosity (for now arbitrary)
    double luminosity = 10000;

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;
    // jet pT cuts
    double jetptMin = 5; // GeV
    double jetptMax = 30; // GeV
    // deltaR histogram
    int deltaRbins = 10000; // deltaRbins = numberOfPoints, default=10 bins for [0. 0.4]
    double minDeltaR = 0.;
    double maxDeltaR = 0.4;
    std::vector<double> deltaRBinEdges = {0.,0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.4}; // TODO: investigate structure before 0.005
    // mass histogram
    int massBins = 100; 
    double minMass = 1.67;
    double maxMass = 2.1;
    // pT,D histograms
    int ptBins = 100;
    double minPt = 0.;
    double maxPt = 30.;
    double ptDBinEdges[] = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};
    double ptjetBinEdges[] = {5., 10., 15., 20., 25., 30.};

    // opening files
    TFile* fSimulatedO2 = new TFile("../SimulatedData/Hyperloop_output/McChargedMatched/AO2D_merged_All.root","read");
    TFile* fPowheg = new TFile("../SimulatedData/POWHEG/trees_powheg_fd_central.root.root","read");
    TFile* fEfficiency = new TFile(Form()"../2-Efficiency/backSubEfficiency_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fBackSub = new TFile(Form("../1-SignalTreatment/backSub_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fSigExt = new TFile(Form("../1-SignalTreatment/sigExt_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fSimulated || fSimulated->IsZombie()) {
        std::cerr << "Error: Unable to open simulated data ROOT file." << std::endl;
    }
    if (!fBackSub || fBackSub->IsZombie()) {
        std::cerr << "Error: Unable to open background subtracted data ROOT file." << std::endl;
    }
    if (!fSigExt || fSigExt->IsZombie()) {
        std::cerr << "Error: Unable to open signal extracted ROOT file." << std::endl;
    }
    
    std::vector<const char*> names = {"histPt1", "histPt2", "histPt3", 
                                      "histPt4", "histPt5", "histPt6",
                                      "histPt7", "histPt8", "histPt9"};                                                     // Names of histograms
    std::vector<const char*> titles = {"3 < p_{T,D} < 4 GeV/c", "4 < p_{T,D} < 5 GeV/c", "5 < p_{T,D} < 6 GeV/c",
                                       "6 < p_{T,D} < 7 GeV/c", "7 < p_{T,D} < 8 GeV/c", "8 < p_{T,D} < 10 GeV/c",
                                       "10 < p_{T,D} < 12 GeV/c", "12 < p_{T,D} < 15 GeV/c", "15 < p_{T,D} < 30 GeV/c"};    // Titles of histograms
    FeedDownData dataContainer = createHistograms(sizeof(deltaRBinEdges)/sizeof(deltaRBinEdges[0]), deltaRBinEdges,         // Delta R
                                                  sizeof(ptDBinEdges)/sizeof(ptDBinEdges[0]), ptDBinEdges,                  // pT,D
                                                  jetptMin, jetptMax);                                                      //pT,jet

    // Fill histograms with POWHEG simulation data
    fillHistograms(fPowheg, dataContainer, jetptMin, jetptMax);

    // Create response matrices for all pT,D bins considered
    createResponseMatrix(dataContainer, 
                         sizeof(deltaRBinEdges)/sizeof(deltaRBinEdges[0]), deltaRBinEdges, 
                         sizeof(ptDBinEdges)/sizeof(ptDBinEdges[0]), ptDBinEdges);

    // Fill response matrix
    fillResponseMatrix(fSimulatedO2, dataContainer, jetptMin, jetptMax);

    smearGeneratorData(dataContainer, luminosity, fEfficiency, titles);

    // Subtract non-prompt distribution from prompt efficiency corrected ones
    feedDown(dataContainer, jetptMin, jetptMax, titles);

    // Plot the efficiency histogram and further corrected histograms
    plotHistograms(dataContainer, jetptMin, jetptMax);

    // Save corrected distributions to file
    saveData(dataContainer, jetptMin, jetptMax);


    

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    FeedDownSubtraction();
    return 0;
}
