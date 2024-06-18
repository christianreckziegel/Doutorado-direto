/*
 *
 *
 * Macro for merging TFiles from different directories 
 * into single TFile with same TTree structure.
 * 
 * 
 * 
 * 
**/

#include <TFile.h>
#include <TTree.h>

using namespace std;

void mergeFiles() {
    // Name of files in path
    const char* inputFile1 = "HF_LHC24d3a_All/Merged_a.root";
    const char* inputFile2 = "HF_LHC24d3b_All/Merged_b.root";
    const char* outputFile = "AO2D_merged_All.root";
    
    // Open input files
    TFile* file1 = TFile::Open(inputFile1);
    TFile* file2 = TFile::Open(inputFile2);
    cout << "Input files opened.\n";cout << "Files opened.\n";

    // Create output file
    TFile* outFile = new TFile(outputFile, "RECREATE");
    cout << "Output files opened.\n";

    // Copy contents of the first file
    TTree* tree1 = (TTree*)file1->Get("O2mcdjetdisttable");
    if (tree1) tree1->CloneTree()->Write();

    TTree* tree2 = (TTree*)file1->Get("O2mcpjetdisttable");
    if (tree2) tree2->CloneTree()->Write();

    // Copy contents of the second file
    TTree* tree3 = (TTree*)file2->Get("O2mcdjetdisttable");
    if (tree3) tree3->CloneTree()->Write();

    TTree* tree4 = (TTree*)file2->Get("O2mcpjetdisttable");
    if (tree4) tree4->CloneTree()->Write();

    cout << "TTrees written in " << outputFile << endl;

    // Close all files
    outFile->Close();
    file1->Close();
    file2->Close();

    delete outFile;
    delete file1;
    delete file2;
}

int main() {
    // Merging files
    mergeFiles();

    return 0;
}
