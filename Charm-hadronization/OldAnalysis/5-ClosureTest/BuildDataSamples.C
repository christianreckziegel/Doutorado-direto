void BuildDataSamples() {    

    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    TFile* fSimulatedMCMatched = new TFile("../SimulatedData/Hyperloop_output/Train_runs/410603_Match/AO2D_mergedDFs.root","read");

    // Accessing TTree
    TTree* originalTree = (TTree*)fSimulatedMCMatched->Get("DF_merged/O2matchtable");
    // Check for correct access
    if (!originalTree) {
        std::cout << "Error opening O2 matching tree.\n";
        return;
    }
    int nEntries = originalTree->GetEntries();
    std::cout << "Total entries in original tree: " << nEntries << std::endl;
    
    // Create the output file *first*
    TFile* outFile = new TFile("mc_split.root", "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cout << "Error creating output file.\n";
        return;
    }
    // Create two output trees
    outFile->cd();
    TTree* inputTree = originalTree->CloneTree(0); // 20% of entries
    TTree* correctionTree = originalTree->CloneTree(0); // 80% of entries

    // Shuffle entry indices
    std::vector<int> indices(nEntries);
    std::iota(indices.begin(), indices.end(), 0);  // Fill with 0, 1, 2, ..., nEntries-1
    std::shuffle(indices.begin(), indices.end(), std::mt19937(std::random_device{}())); // Randomize the order

    // Split 20% input / 80% correction
    int splitIndex = int(nEntries * 0.2);
    std::cout << "Splitting at index: " << splitIndex << " (20%)" << std::endl;

    const int progressInterval = std::max(1, splitIndex / 100);  // Print every 1%
    std::cout << "Filling InputTree..." << std::endl;
    for (int i = 0; i < splitIndex; ++i) {
        if (i % progressInterval == 0) {
            std::cout << "  " << (100 * i / splitIndex) << "% done" << std::endl;
            std::cout << "Processing entry " << i << " of " << splitIndex << std::endl;
        }
        //std::cout << "Processing entry " << i << " of " << splitIndex << std::endl;
        originalTree->GetEntry(indices[i]);
        // ROOT trees internally keep track of branch addresses
        inputTree->Fill();
    }

    std::cout << "Filling CorrectionTree..." << std::endl;
    const int correctionSize = nEntries - splitIndex;
    const int progressIntervalCorr = std::max(1, correctionSize / 100);
    for (int i = splitIndex; i < nEntries; ++i) {
        if ((i - splitIndex) % progressIntervalCorr == 0) {
            std::cout << "  " << (100 * (i - splitIndex) / correctionSize) << "% done" << std::endl;
        }
        //std::cout << "Processing entry " << i << " of " << nEntries << std::endl;
        originalTree->GetEntry(indices[i]);
        // ROOT trees internally keep track of branch addresses
        correctionTree->Fill();
    }

    inputTree->Write("InputTree");
    correctionTree->Write("CorrectionTree");
    outFile->Close();
    fSimulatedMCMatched->Close();

    std::cout << "Tree split complete: " << splitIndex << " in InputTree, " << (nEntries - splitIndex) << " in CorrectionTree." << std::endl;

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 
}

int main() {
    BuildDataSamples();
    return 0;
}