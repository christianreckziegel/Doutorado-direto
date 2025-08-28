
void BuildMatchedData() {
    //
    TFile* inFile = TFile::Open("../SimulatedData/Hyperloop_output/Train_runs/410603_Match/AO2D_mergedDFs.root", "READ");
    TTree* originalTree = (TTree*)inFile->Get("DF_merged/O2matchtable");
    if (!originalTree) {
        std::cout << "Error: Couldn't find tree!\n";
        return;
    }

    // Create output file FIRST
    TFile* outFile = new TFile("mc_closure_input_matched_data.root", "RECREATE");
    outFile->cd();

    // Clone original structure to disk-bound tree
    TTree* workingTree = originalTree->CloneTree(0); // empty clone
    int splitFlag;
    TBranch* b_split = workingTree->Branch("SplitGroup", &splitFlag, "SplitGroup/I");

    TRandom3 rand(0);
    Long64_t nEntries = originalTree->GetEntries();

    const Long64_t progressStep = nEntries / 100; // Print every 1%
    if (progressStep == 0) std::cout << "Warning: very few entries, no progress bar.\n";

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (progressStep > 0 && i % progressStep == 0) {
            int percent = static_cast<int>(100.0 * i / nEntries);
            std::cout << "\rFilling workingTree: " << percent << "% done..." << std::flush;
        }
        originalTree->GetEntry(i);
        splitFlag = (rand.Rndm() < 0.2) ? 0 : 1;
        workingTree->Fill();  // Now safely writes to disk
    }

    // Create filtered trees
    TTree* inputTree = workingTree->CopyTree("SplitGroup==0");
    TTree* correctionTree = workingTree->CopyTree("SplitGroup==1");

    inputTree->Write("InputTree");
    correctionTree->Write("CorrectionTree");

    outFile->Close();
    inFile->Close();

    std::cout << "\nFast tree split done with disk-backed tree.\n";
}

void BuildEffData() {
    //
    TFile* inFile = TFile::Open("../SimulatedData/Hyperloop_output/Train_runs/410602_Eff/AO2D_mergedDFs.root", "READ");
    
    TTree* detTree = (TTree*)inFile->Get("DF_merged/O2mcdjetdisttable");
    TTree* partTree = (TTree*)inFile->Get("DF_merged/O2mcpjetdisttable");

    if (!detTree || !partTree) {
        std::cout << "Error: Couldn't find one or both trees!\n";
        return;
    }

    // Create output file FIRST
    TFile* outFile = new TFile("mc_closure_input_non-matched_data.root", "RECREATE");
    // Create two directories: InputTree and CorrectionTree
    TDirectory* inputDir = outFile->mkdir("InputTree");
    TDirectory* corrDir  = outFile->mkdir("CorrectionTree");

    // === Helper lambda to split and write trees === //
    auto processTree = [&](TTree* originalTree, const char* treeName) {
        outFile->cd();

        TTree* workingTree = originalTree->CloneTree(0);
        int splitFlag;
        TBranch* b_split = workingTree->Branch("SplitGroup", &splitFlag, "SplitGroup/I");

        TRandom3 rand(0);  // Random seed can be fixed or dynamic
        Long64_t nEntries = originalTree->GetEntries();
        const Long64_t progressStep = nEntries / 100;

        for (Long64_t i = 0; i < nEntries; ++i) {
            if (progressStep > 0 && i % progressStep == 0) {
                int percent = static_cast<int>(100.0 * i / nEntries);
                std::cout << "\rProcessing " << treeName << ": " << percent << "% done..." << std::flush;
            }
            originalTree->GetEntry(i);
            splitFlag = (rand.Rndm() < 0.2) ? 0 : 1;
            workingTree->Fill();
        }

        std::cout << "\rProcessing " << treeName << ": 100% done.\n";

        TTree* inputTree = workingTree->CopyTree("SplitGroup==0");
        TTree* correctionTree = workingTree->CopyTree("SplitGroup==1");

        // Write input tree into InputTree/ directory
        inputDir->cd();
        inputTree->SetName(treeName);
        inputTree->Write();
        // Write correction tree into CorrectionTree/ directory
        corrDir->cd();
        correctionTree->SetName(treeName);
        correctionTree->Write();
        // Go back to root dir of output file for next call
        outFile->cd();

        // TString inputPath = TString::Format("%s/InputTree", treeName);
        // TString corrPath  = TString::Format("%s/CorrectionTree", treeName);
        // inputTree->Write(inputPath);
        // correctionTree->Write(corrPath);

    };

    // === Process both trees === //
    processTree(detTree, "O2mcdjetdisttable");
    processTree(partTree, "O2mcpjetdisttable");

    outFile->Close();
    inFile->Close();

    std::cout << "\nAll trees split and saved.\n";
}

void FastBuildDataSamples() {
    BuildMatchedData();

    BuildEffData();    
}

int main() {
    FastBuildDataSamples();
    return 0;
}