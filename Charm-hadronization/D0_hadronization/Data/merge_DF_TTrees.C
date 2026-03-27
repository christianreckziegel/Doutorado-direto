#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TChain.h>
#include <iostream>
#include <vector>
#include <map>

void merge_DF_TTrees(const char* inputFileName = "AO2D.root", const char* outputFileName = "AO2D_mergedDFs.root") {
    std::cout << "Starting step 1..." << std::endl;
    // Open input file
    TFile *inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFileName << std::endl;
        return;
    }
    std::cout << "Opened input file: " << inputFileName << std::endl;

    // Create output file
    TFile *outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << outputFileName << std::endl;
        inputFile->Close();
        return;
    }
    std::cout << "Created output file: " << outputFileName << std::endl;

    // Map to store TChains by tree name
    std::map<std::string, TChain*> treeChains;

    // Get list of directories (DF_*)
    inputFile->cd();
    TIter nextDir(inputFile->GetListOfKeys());
    TKey *key;
    bool firstDirectory = true;
    while ((key = (TKey*)nextDir())) {
        if (std::string(key->GetClassName()) != "TDirectoryFile") continue;
        TDirectory *dir = (TDirectory*)key->ReadObj();
        std::string dirName = dir->GetName();
        //std::cout << "Processing directory: " << dirName << std::endl;

        // Get list of trees in the directory
        TIter nextTree(dir->GetListOfKeys());
        TKey *treeKey;
        while ((treeKey = (TKey*)nextTree())) {
            if (std::string(treeKey->GetClassName()) != "TTree") continue;
            std::string treeName = treeKey->GetName();
            //std::cout << "  Found tree: " << treeName << " in " << dirName << std::endl;

            if (firstDirectory) {
                std::string treePath = dirName + "/" + treeName;
                //std::cout << "    First tree path: " << treePath << std::endl;
                firstDirectory = false;  // Stop after printing the first TTree path
            }

            // Add tree to corresponding TChain
            if (treeChains.find(treeName) == treeChains.end()) {
                treeChains[treeName] = new TChain(treeName.c_str());
            }
            std::string treePath = dirName + "/" + treeName;
            if (dir->GetListOfKeys()->Contains(treeName.c_str())) {
                treeChains[treeName]->Add((std::string(inputFileName) + "/" + treePath).c_str());
                //std::cout << "    Added to chain: " << treePath << std::endl;
            } else {
                std::cout << "    Tree " << treeName << " not found in " << dirName << std::endl;
            }
        }
    }

    // Write merged trees to output file
    outputFile->mkdir("DF_merged")->cd();
    for (auto &pair : treeChains) {
        //std::cout << "Merging tree: " << pair.first << std::endl;
        TTree *mergedTree = pair.second->CloneTree(-1);
        if (mergedTree) {
            mergedTree->Write();
        } else {
            std::cerr << "Error: mergedTree is null, skipping write" << std::endl;
        }
        delete pair.second;  // Free memory
        //std::cout << "  Written merged tree: " << pair.first << std::endl;
    }

    // Clean up
    outputFile->Close();
    inputFile->Close();
    std::cout << "Merging completed! Output saved in " << outputFileName << std::endl;

    // Clear the map (optional, for safety)
    treeChains.clear(); 
}

int main() {
    merge_DF_TTrees();
    return 0;
}