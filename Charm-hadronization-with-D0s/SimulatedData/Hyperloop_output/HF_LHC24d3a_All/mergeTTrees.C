/*
 *
 *
 * Macro for merging TTrees from different directories 
 * into single TTree in single directory.
 * 
 * 
 * 
 * 
**/

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TKey.h>
#include <iostream>

using namespace std;

void mergeTTrees() {
    // File names
    const char* inputFileName = "AO2D.root";
    const char* outputFileName = "Merged_a.root";
    
    cout << "Files opened.\n";

    // Open input file
    TFile* inputFile = TFile::Open(inputFileName);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening file: " << inputFileName << std::endl;
        return;
    }

    // Create TChains for each TTree type
    TChain mcdjetdistChain("O2mcdjetdisttable");
    TChain mcpjetdistChain("O2mcpjetdisttable");

    // Get list of keys in the file
    TList* keys = inputFile->GetListOfKeys();
    TIter next(keys);
    TKey* key;

    cout << "Looping through directories in input file...\n";

    // Loop over all keys in file
    while((key = (TKey*)next())) { // use of an assignment as a condition requires parentheses
        // Get the directory
        TDirectoryFile* dir = (TDirectoryFile*)key->ReadObj();
        if (dir) {
            // Add TTrees to the respective TChains
            TTree* mcdjetdistTree = (TTree*)dir->Get("O2mcdjetdisttable");
            if (mcdjetdistTree) {
                TString path = TString::Format("file://%s/%s/O2mcdjetdisttable", inputFileName, dir->GetName());
                mcdjetdistChain.Add(path);
            } else {
                std::cerr << "Error: TTree O2mcdjetdisttable not found in directory " << dir->GetName() << std::endl;
            }
            TTree* mcpjetdistTree = (TTree*)dir->Get("O2mcpjetdisttable");
            if (mcpjetdistTree) {
                TString path = TString::Format("file://%s/%s/O2mcpjetdisttable", inputFileName, dir->GetName());
                mcpjetdistChain.Add(path);
            } else {
                std::cerr << "Error: TTree O2mcpjetdisttable not found in directory " << dir->GetName() << std::endl;
            }
        } else {
            std::cerr << "Error: Directory not found for key " << key->GetName() << std::endl;
        }
        
    } // end of loop

    // Create the output file
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating file: " << outputFileName << std::endl;
        return;
    }
    cout << "Output file created.\n";

    // Merge the TTrees in each TChain into a single TTree and write to the output file
    TTree* mcdjetdistTreeMerged = mcdjetdistChain.CloneTree();
    mcdjetdistTreeMerged->Write();
    //mcdjetdistTreeMerged->Write("", TObject::kOverwrite); // Overwrite if exists

    TTree* mcpjetdistTreeMerged = mcpjetdistChain.CloneTree();
    mcpjetdistTreeMerged->Write();
    //mcdjetdistTreeMerged->Write("", TObject::kOverwrite); // Overwrite if exists

    cout << "TTrees written to output file.\n";

    // Clean up
    delete inputFile;
    outputFile->Close();
    delete outputFile;
}

int main(){    
    // Merge the TTrees
    mergeTTrees();

    return 0;
}
