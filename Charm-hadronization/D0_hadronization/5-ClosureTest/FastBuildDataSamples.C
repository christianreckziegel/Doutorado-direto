#include "../commonUtilities.h"

void BuildMatchedData() {
    
    // Binning file
    TFile* fBinning = new TFile(Form("../1-SignalTreatment/BDTOptimization/binningInfo_2023.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    // Load binning from reflections file
    BinningStruct binning = retrieveBinningFromFile(fBinning);

    //
    //TFile* inFile = TFile::Open("../" + binning.inputMC.second + "/AO2D_mergedDFs.root", "READ");
    // 2022 MC anchored data -> Train_669231
    // 2023 MC anchored data -> Train_671273
    TFile* inFile = TFile::Open("../Data/MonteCarlo/Train_671273/AO2D_mergedDFs.root", "READ");
    TTree* originalTree = (TTree*)inFile->Get("DF_merged/O2matchtable");
    if (!originalTree) {
        std::cout << "Error: Couldn't find tree!\n";
        return;
    }

    // Create output file FIRST
    TFile* outFile = new TFile("mc_closure_input_data_2023.root", "RECREATE");
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
    std::cout << "New iput tree:" << std::endl;
    inputTree->Print();

    inputTree->Write("InputTree");
    correctionTree->Write("CorrectionTree");

    outFile->Close();
    inFile->Close();

    std::cout << "\nFast tree split done with disk-backed tree. Data stored to file" << outFile->GetName() << std::endl;
}


void FastBuildDataSamples() {
    BuildMatchedData(); 
}

int main() {
    FastBuildDataSamples();
    return 0;
}