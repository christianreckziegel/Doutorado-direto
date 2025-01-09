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

struct UnfoldData {
    TH2D* hMeasured;                                        // 2D representation of response measured data
    TH2D* hMeasuredTotalRange;                              // 2D total measured range data
    TH2D* hDivMeasuredRange;                                // 2D division of (inside range) / (total range) for measured data
    TH2D* hTruth;                                           // 2D representation of response truth data
    TH2D* hTruthTotalRange;                                 // 2D total truth range data
    TH2D* hDivTruthRange;                                   // 2D division of (inside range) / (total range) for truth data

    RooUnfoldResponse response;                             // response matrix for non-prompt D0s only, for overall pT,D; using RooUnfoldResponse object, method 1
    std::vector<TH1D*> hEfficiencies;                       // inclusive = 0, prompt only = 1, non-prompt only = 2

};

// Obtain the bin edges of a histogram (useful for asymmetrical bin sizes)
std::vector<double> getBinEdges(const TAxis* axis) {
    
    if (!axis) {
        std::cerr << "Error: null axis pointer.\n";
    }
    
    int nBins = axis->GetNbins();

    // vector for storing bin edges
    std::vector<double> binEdges(nBins + 1);

    for (int iBin = 0; iBin <= nBins; iBin++) {
        binEdges[iBin] = axis->GetBinLowEdge(iBin + 1);
    }

    return binEdges;
}

// Module to create TH2D histograms including interest variable
UnfoldData createHistograms(const std::vector<double>& xBinEdges_particle, const std::vector<double>& yBinEdges_particle, const std::vector<double>& zBinEdges_particle,
                              const std::vector<double>& xBinEdges_detector, const std::vector<double>& yBinEdges_detector, const std::vector<double>& zBinEdges_detector) {
                              //const double& jetptMin, const double& jetptMax) {
    // Create struct to store data
    UnfoldData dataContainer;

    int xNumBinEdges = xBinEdges_particle.size();
    int yNumBinEdges = yBinEdges_particle.size();
    int zNumBinEdges = zBinEdges_particle.size();
    

    //
    // Matching histograms for folding process
    //
    // 2D truth and particle level data
    dataContainer.hTruth = new TH2D("hTruth2D", "Truth;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges_particle.data(), zNumBinEdges-1, zBinEdges_particle.data());
    dataContainer.hTruthTotalRange = new TH2D("hTruth2D_totalRange", "Truth;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges_particle.data(), zNumBinEdges-1, zBinEdges_particle.data());
    
    // 2D measured and detector level data
    xNumBinEdges = xBinEdges_detector.size();
    yNumBinEdges = yBinEdges_detector.size();
    zNumBinEdges = zBinEdges_detector.size();
    dataContainer.hMeasured = new TH2D("hMeasured2D", "Inside response range measured;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges_detector.data(), zNumBinEdges-1, zBinEdges_detector.data());
    dataContainer.hMeasuredTotalRange = new TH2D("hMeasured2D_totalRange", "Total measured range;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges_detector.data(), zNumBinEdges-1, zBinEdges_detector.data());


    cout << "Matching histograms created.\n";



    return dataContainer;
}

/**
 * @brief Module to fill histograms from O2 matching task simulation.
 * 
 * This helper function fill two kinds of histograms:
 * - POWHEG+PYTHIA simulated data on particle level with non-prompt D0 
 * - O2 simulated matched data between particle and detector level
 *
 * @param fData The ROOT file with POWHEG+PYTHIA non-prompt D0 simulated data.
 * @param fSimulatedO2 The ROOT file with matched O2 simulated data.
 * @param dataContainer Container with already instanciated histograms.
 * @param jetptMin Min jet pT cut.
 * @param jetptMax Max jet pT cut.
 *
 * @note uses 0-based indexing for calculation but 1-based indexing for ROOT histograms.
 *
 * @see createHistograms() [Instanciate histograms.]
 */
void fillHistograms(TFile* fData, TFile* fSimulatedO2, UnfoldData& dataContainer, double jetptMin, double jetptMax) {
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = 0.4; // on particle level delta R
    const double MCDDeltaRcut = 0.4; // on detector level delta R
    const double MCPHfPtMincut = 3.; // on particle level
    const double MCDHfPtMincut = 3.; // on detector level
    const double MCPHfPtMaxcut = 30.; // on particle level
    const double MCDHfPtMaxcut = 30.; // on detector level

    //
    // Data to be unfolded tree and histograms
    //
    // Accessing TTree
    TTree* tree;
    int nEntries;

    //
    // O2 matching task measured and truth histograms
    //
    // Accessing TTree
    tree = (TTree*)fSimulatedO2->Get("DF_2263915935550653/O2matchtable");
    // Check for correct access
    if (!tree) {
        cout << "Error opening O2 matching tree.\n";
    }

    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
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
    tree->SetBranchAddress("fMCHfY",&MCPhfY);
    tree->SetBranchAddress("fMCHfPrompt",&MCPhfprompt);
    // detector level branches
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);

    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));
        
        

        // Fill histograms considering jet pT and detector acceptance for PROMPT particles, inside response range (truth and measured levels)
        if ((abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= 0.) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)) && MCPhfprompt
            && (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= 0.) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut)) && MCDhfprompt) {
            
            // Filling measured 2D histogram
            dataContainer.hMeasured->Fill(MCDDeltaR, MCDjetPt);

            // Filling truth 2D histogram
            dataContainer.hTruth->Fill(MCPDeltaR, MCPjetPt);
        }

        // Fill reference output total measured range histograms for non-prompt particles
        if ((abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= 0.) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut)) && MCDhfprompt) {
            dataContainer.hMeasuredTotalRange->Fill(MCDDeltaR, MCDjetPt);
        }

    }

    cout << "Response matched histograms filled.\n";

}

/**
 * @brief Module to build 2D response matrix out of flattened 1D input data.
 * 
 * This helper function fill two kinds of histograms:
 * - POWHEG+PYTHIA simulated data on particle level with non-prompt D0 
 * - O2 simulated matched data between particle and detector level
 * 
 * @param dataContainer Container with already instanciated histograms.
 * @param fSimulatedO2 The ROOT file with matched O2 simulated data.
 * @param fEfficiency The ROOT file containing efficiency values.
 * @param jetptMin Min jet pT cut.
 * @param jetptMax Max jet pT cut.
 *
 * @note uses 0-based indexing for calculation but 1-based indexing for ROOT histograms.
 *
 * @see createHistograms() [Instanciate histograms.]
 */
void buildResponseMatrix(UnfoldData& dataContainer, TFile* fSimulatedO2, TFile* fEfficiency, double jetptMin, double jetptMax) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = 0.4; // on particle level delta R
    const double MCDDeltaRcut = 0.4; // on detector level delta R
    const double MCPHfPtMincut = 3.; // on particle level
    const double MCDHfPtMincut = 3.; // on detector level
    const double MCPHfPtMaxcut = 30.; // on particle level
    const double MCDHfPtMaxcut = 30.; // on detector level

    // method 1: create 2D response matrix of prompt flattened Delta R and pT,jet, for overall pT,D
    //(DeltaR_detector, pTjet_detector, DeltaR_particle, pTjet_particle) -> flattened to (detector, particle)
    dataContainer.response = RooUnfoldResponse(dataContainer.hMeasured, dataContainer.hTruth);

    std::cout << "Response matrix created.\n";

    // Access prompt efficiency
    TH1D* hEffPrompt = (TH1D*)fEfficiency->Get("efficiency_prompt");

    //__________________________________________-
    //
    // O2 matching task measured and truth histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)fSimulatedO2->Get("DF_2263915935550653/O2matchtable");
    // Check for correct access
    if (!tree) {
        cout << "Error opening O2 matching tree.\n";
    }

    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
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
    tree->SetBranchAddress("fMCHfY",&MCPhfY);
    tree->SetBranchAddress("fMCHfPrompt",&MCPhfprompt);
    // detector level branches
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        // Fill histograms considering jet pT and detector acceptance
        if ((abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= 0.) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)) && MCPhfprompt
            && (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= 0.) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut)) && MCDhfprompt) {
            
            // Find the bin corresponding to the given pT value
            int bin = hEffPrompt->FindBin(MCDhfPt);
            // Get the efficiency value from the bin content
            double efficiency_prompt = hEffPrompt->GetBinContent(bin);

            // Fill 4D RooUnfoldResponse object
            dataContainer.response.Fill(MCDDeltaR, MCDjetPt, MCPDeltaR, MCPjetPt);
            //dataContainer.response.Fill(MCDDeltaR, MCDjetPt, MCPDeltaR, MCPjetPt, 1./efficiency_prompt); // jet pT shape is influenced by D0 pT efficiency
        }
        
    }
    
    std::cout << "Response matrix filled.\n";
}


/**
 * @brief Remove entries in region not treated by response matrix from the unfolding input data.
 * 
 *
 * @param dataContainer Container with total and intersection range data and the input data to be corrected.
 *
 * @return Corrected measured data input 2D distribution
 *
 * @note 
 *
 * @see smearGeneratorData() [Correct truth POWHEG data before folding.]
 */
TH2D* removeOutsideData(UnfoldData& dataContainer) {

    // Copy original data distribution structure
    //TH2D* hDataCorrected = (TH2D*)dataContainer.hAllptDPowheg[0]->Clone("hPowhegRangeCorrected");

    // Start with histogram of intersection area only (inside response ranges)
    dataContainer.hDivMeasuredRange = dynamic_cast<TH2D*>(dataContainer.hMeasured->Clone("hInsideOverTotalDivision"));

    // Get kinematic efficiency 2D histogram diving intersection over total range
    dataContainer.hDivMeasuredRange->Divide(dataContainer.hMeasuredTotalRange);
    dataContainer.hDivMeasuredRange->SetTitle("#varepsilon_{kinematic} = Inside response range / Total measured range");

    // Correct POWHEG data multiplying each bin by corresponding kinematic efficiency
    //hDataCorrected->Multiply(dataContainer.hDivMeasuredRange);


    return dataContainer.hDivMeasuredRange;
    //return hDataCorrected;

}



void plotHistograms(const UnfoldData& dataContainer, const double& jetptMin, const double& jetptMax) {
    cout << "Plotting histograms...\n";

    gStyle->SetPalette(kRainbow);

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    //
    // Measured kinematic efficiency
    //
    TCanvas* cMeasuredKinEff = new TCanvas("cMeasuredKinEff","Measured data kinematic efficiency");
    cMeasuredKinEff->SetCanvasSize(1800,1000);
    //cMeasuredKinEff->Divide(2,2);
    //cMeasuredKinEff->cd(1);
    //dataContainer.hMeasured->Draw("colz");
    //cMeasuredKinEff->cd(2);
    //dataContainer.hMeasuredTotalRange->Draw("colz");
    //cMeasuredKinEff->cd(3);
    cMeasuredKinEff->cd();
    dataContainer.hDivMeasuredRange->Draw("text");

    //
    // Response matrix representation in 2D histogram
    //
    TCanvas* cResponse = new TCanvas("cResponse","Response matrix");
    cResponse->SetCanvasSize(1800,1000);
    cResponse->cd();
    const TH2* hResponse2D = dataContainer.response.Hresponse();
    TH2D* hResponse2DClone = static_cast<TH2D*>(hResponse2D->Clone("hResponse2DClone"));
    hResponse2DClone->SetTitle("2D representation of RooUnfoldResponse;2D Reconstructed;2D Truth");
    hResponse2DClone->Draw("colz");

    //
    // Storing images
    //
    TString imagePath = "../Images/4-Unfolding/";
    cResponse->Update();
    cResponse->SaveAs(imagePath + "UN_response_matrix.png");
    cMeasuredKinEff->Update();
    cMeasuredKinEff->SaveAs(imagePath + "UN_measured_kin_efficiency.png");
    

}

void saveData(const UnfoldData& dataContainer, const double& jetptMin, const double& jetptMax){
    // Open output file
    TFile* outFile = new TFile(Form("backSubUnfold_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");

    // store each histogram in file
    //dataContainer.hSBUnfolded->Write();
    
    outFile->Close();
    delete outFile;
    
    cout << "Data stored in file" << Form("backSubUnfold_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)) << endl;
}

void Unfolding(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Luminosity (for now arbitrary)
    double luminosity = 1000000;
    double luminosity_powheg = 0;

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;
    // jet pT cuts
    std::vector<double> ptjetBinEdges_particle = {5., 7., 15., 30.};
    std::vector<double> ptjetBinEdges_detector = {5., 7., 15., 30.};
    double jetptMin = ptjetBinEdges_particle[0]; // GeV
    double jetptMax = ptjetBinEdges_particle[ptjetBinEdges_particle.size() - 1]; // GeV
    // deltaR histogram
    int deltaRbins = 10000; // deltaRbins = numberOfPoints, default=10 bins for [0. 0.4]
    std::vector<double> deltaRBinEdges_particle = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4}; // TODO: investigate structure before 0.005: 0.,0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.4
    std::vector<double> deltaRBinEdges_detector = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4};
    double minDeltaR = deltaRBinEdges_particle[0];
    double maxDeltaR = deltaRBinEdges_particle[deltaRBinEdges_particle.size() - 1];
    // mass histogram
    int massBins = 100; 
    double minMass = 1.67;
    double maxMass = 2.1;
    // pT,D histograms
    int ptBins = 100;
    std::vector<double> ptDBinEdges_particle = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};
    std::vector<double> ptDBinEdges_detector = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};
    double minPtD = ptDBinEdges_particle[0];
    double maxPtD = ptDBinEdges_particle[ptDBinEdges_particle.size() - 1];

    // Opening files
    TFile* fData = new TFile("../ExperimentalData/Hyperloop_output/AO2D.root","read");
    TFile* fSimulatedO2 = new TFile("../SimulatedData/Hyperloop_output/McChargedMatched/HF_LHC24d3a_All/AO2D.root","read");
    TFile* fEfficiency = new TFile(Form("../2-Efficiency/backSubEfficiency_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fBackSub = new TFile(Form("../1-SignalTreatment/SideBand/backSub_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fSigExt = new TFile(Form("../1-SignalTreatment/SignalExtraction/sigExt_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fData || fData->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT file." << std::endl;
    }
    if (!fSimulatedO2 || fSimulatedO2->IsZombie()) {
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
    //UnfoldData dataContainer = createHistograms(deltaRBinEdges, ptDBinEdges, jetptMin, jetptMax);
    UnfoldData dataContainer = createHistograms(deltaRBinEdges_particle, ptDBinEdges_particle, ptjetBinEdges_particle,
                                                  deltaRBinEdges_detector, ptDBinEdges_detector, ptjetBinEdges_detector);

    // Fill histograms with POWHEG simulation data
    fillHistograms(fData, fSimulatedO2, dataContainer, jetptMin, jetptMax);

    // Create response matrices for all pT,D bins considered
    buildResponseMatrix(dataContainer, fSimulatedO2, fEfficiency, jetptMin, jetptMax);

    TH2D* hremoveOutsideData = removeOutsideData(dataContainer);

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
    Unfolding();
    return 0;
}
