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

struct ClosureTestData {
    TTree* inputTree; // 20% of entries
    TTree* correctionTree; // 80% of entries
    
    
};



// Module to create TH2D histograms including interest variable
ClosureTestData createHistograms(const std::vector<double>& ptjetBinEdges_particle, const std::vector<double>& deltaRBinEdges_particle, const std::vector<double>& ptDBinEdges_particle,
                              const std::vector<double>& ptjetBinEdges_detector, const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& ptDBinEdges_detector,
                            int& iterationNumber) {
                              //const double& jetptMin, const double& jetptMax) {
    // Create struct to store data
    ClosureTestData dataContainer;
    

    std::cout << "Histograms created." << std::endl;
    return dataContainer;
}

// Get the optimal BDT score cut for the corresponding pT,D of the entry
double GetBkgProbabilityCut(double pT, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    for (size_t i = 0; i < bdtPtCuts.size() - 1; ++i) {
        if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) {
            return bdtPtCuts[i].second;
        }
    }
    return 1.0; // Default: accept all if out of range
}
void BuildDataSamples(TFile* fSimulatedMCMatched, ClosureTestData& dataContainer) {    

    // Accessing TTree
    TTree* originalTree = (TTree*)fSimulatedMCMatched->Get("DF_merged/O2matchtable");
    // Check for correct access
    if (!originalTree) {
        std::cout << "Error opening O2 matching tree.\n";
        return;
    }
    int nEntries = originalTree->GetEntries();
    
    // Create the output file *first*
    TFile* outFile = new TFile("mc_split.root", "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cout << "Error creating output file.\n";
        return;
    }
    // Create two output trees
    outFile->cd();
    dataContainer.inputTree = originalTree->CloneTree(0); // 20% of entries
    dataContainer.correctionTree = originalTree->CloneTree(0); // 80% of entries

    // Shuffle entry indices
    std::vector<int> indices(nEntries);
    std::iota(indices.begin(), indices.end(), 0);  // Fill with 0, 1, 2, ..., nEntries-1
    std::shuffle(indices.begin(), indices.end(), std::mt19937(std::random_device{}())); // Randomize the order

    // Split 20% input / 80% correction
    int splitIndex = int(nEntries * 0.2);

    for (int i = 0; i < splitIndex; ++i) {
        originalTree->GetEntry(indices[i]);
        // ROOT trees internally keep track of branch addresses
        dataContainer.inputTree->Fill();
    }

    for (int i = splitIndex; i < nEntries; ++i) {
        originalTree->GetEntry(indices[i]);
        // ROOT trees internally keep track of branch addresses
        dataContainer.correctionTree->Fill();
    }

    dataContainer.inputTree->Write("InputTree");
    dataContainer.correctionTree->Write("CorrectionTree");
    outFile->Close();
    fSimulatedMCMatched->Close();

    std::cout << "Tree split complete: " << splitIndex << " in InputTree, " << (nEntries - splitIndex) << " in CorrectionTree." << std::endl;
}

// Closure of the unfolding procedure
void FirstClosureTest() {
    // Create 2D (detector level, prompt D0's, matched to particle level) input distribution pT,jet vs DeltaR

    // Create 2D (matched particle level, prompt D0's, matched to the previous detector level distribution) input distribution pT,jet vs DeltaR

    // Create response matrix with correction sample (without efficiency scaling)

    // Unfold the detector level distribution (with particle and detector level kinematic efficiency corrections)

    // Compare the unfolded distribution with the particle level distribution

}

void plotHistograms(const ClosureTestData& dataContainer, const double& jetptMin, const double& jetptMax) {
    cout << "Plotting histograms...\n";

    gStyle->SetPalette(kRainbow);

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    
    
    //
    // Storing images
    //
    TString imagePath = "../Images/5-ClosureTest/";
    //cKinEff->Update();
    //cKinEff->SaveAs(imagePath + "Unfolding_kin_efficiencies.png");    
    
    //
    // Storing in a single pdf file
    //
    //cKinEff->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    //cResponse->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    //cUnfoldedIter->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

}

void saveData(const ClosureTestData& dataContainer, const double& jetptMin, const double& jetptMax){
    // Open output file
    TFile* outFile = new TFile(Form("closure_test_results_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");

    // store each histogram in file
    //dataContainer.hSBUnfolded->Write();
    
    outFile->Close();
    delete outFile;
    
    cout << "Data stored in file" << Form("closure_test_results_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)) << endl;
}

void ClosureTest(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // Number of unfolding procedure iterations
    int iterationNumber = 8;

    // jet pT cuts
    std::vector<double> ptjetBinEdges_particle = {5., 7., 15., 30., 50.};
    std::vector<double> ptjetBinEdges_detector = {5., 7., 15., 30., 50.};
    double jetptMin = ptjetBinEdges_particle[0]; // GeV
    double jetptMax = ptjetBinEdges_particle[ptjetBinEdges_particle.size() - 1]; // GeV

    // deltaR histogram
    std::vector<double> deltaRBinEdges_particle = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15,0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5}; // default = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4} chosen by Nima
    std::vector<double> deltaRBinEdges_detector = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15,0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5}; // default = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4} chosen by Nima
    double minDeltaR = deltaRBinEdges_particle[0];
    double maxDeltaR = deltaRBinEdges_particle[deltaRBinEdges_particle.size() - 1];
    
    // pT,D histograms
    std::vector<double> ptDBinEdges_particle = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.}; // default pT,D = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.}
    std::vector<double> ptDBinEdges_detector = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.}; // default pT,D = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.}
    double hfptMin = ptDBinEdges_particle[0];
    double hfptMax = ptDBinEdges_particle[ptDBinEdges_particle.size() - 1];

    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.03 (from first of pair)
    //std::vector<std::pair<double, double>> bdtPtCuts = {
    //    {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}
    //};
    // Dataset: JE_HF_LHC24g5_All_D0
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.4}, {8, 0.6}, {10, 0.8}, {12, 0.8}, {16, 1.0}, {50, 1.0}
    };

    // Opening files
    TFile* fSimulatedMCNonMatched = new TFile("../SimulatedData/Hyperloop_output/Train_runs/410602_Eff/AO2D_mergedDFs.root","read");
    TFile* fSimulatedMCMatched = new TFile("../SimulatedData/Hyperloop_output/Train_runs/410603_Match/AO2D_mergedDFs.root","read");
    TFile* fEfficiency = new TFile(Form("../2-Efficiency/selection_efficiency_run3style_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fData = new TFile("../ExperimentalData/Hyperloop_output/HF_LHC23_pass4_Thin_small_2P3PDstar_DATA_newMLModel/AnalysisResults.root","read");
    TFile* fFeedDown = new TFile(Form("../3-Feed-Down/outputFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fSimulatedMCMatched || fSimulatedMCMatched->IsZombie()) {
        std::cerr << "Error: Unable to open O2 MC matched ROOT file." << std::endl;
    }
    if (!fEfficiency || fEfficiency->IsZombie()) {
        std::cerr << "Error: Unable to open estimated selection efficiency ROOT file." << std::endl;
    }
    if (!fData || fData->IsZombie()) {
        std::cerr << "Error: Unable to open AnalysisResults.root data ROOT file." << std::endl;
    }
    if (!fFeedDown || fFeedDown->IsZombie()) {
        std::cerr << "Error: Unable to open AnalysisResults.root data ROOT file." << std::endl;
    }
    
    // ClosureTestData dataContainer = createHistograms(deltaRBinEdges, ptDBinEdges, jetptMin, jetptMax);
    ClosureTestData dataContainer = createHistograms(ptjetBinEdges_particle, deltaRBinEdges_particle, ptDBinEdges_particle,
                                                     ptjetBinEdges_detector, deltaRBinEdges_detector, ptDBinEdges_detector,
                                                     iterationNumber);

    // Build MC data samples
    // a) input sample (20%): matched data with detector level to test and particle level to compare with
    // b) correction sample (80%): build all MC level correction objects (efficiencies and response matrices)
    BuildDataSamples(fSimulatedMCMatched, dataContainer);

    // 1 - Unfolding closure test
    FirstClosureTest();

    // 2 - Unfolding closure test + sideband subtraction correction steps

    // 3 - Unfolding closure test + sideband subtraction correction steps + feed-down subtraction steps

    // Plot the efficiency histogram and further corrected histograms
    //plotHistograms(dataContainer, jetptMin, jetptMax);

    // Save corrected distributions to file
    //saveData(dataContainer, jetptMin, jetptMax);


    

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    ClosureTest();
    return 0;
}
