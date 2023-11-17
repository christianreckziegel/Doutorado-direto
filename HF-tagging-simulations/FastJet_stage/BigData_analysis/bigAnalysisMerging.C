/**
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * Compile it with:
 * g++ -g bigAnalysisMerging.C -o bigAnalysisMerging `root-config --cflags --libs`
 * 
*/

#include "TFile.h" // for saving file.
#include "TTree.h" // for accessing data on TTrees
#include "TH1F.h" // for saving data on histograms
#include "TString.h" // for appropriate organization of I/O read and saving
#include <TSystem.h>  // Include this line for gSystem usage
#include <iostream> // for printing with cerr

using namespace std;

int main(int argc, char* argv[]){

    // initial hard coded variables
    TString sNumEvents = argv[1];// 1k, 10k, 100k, 1m
    int numOfRuns = atoi(argv[2]);

    // defining merging histograms
    TH1F* hMergedLeadPt = new TH1F("hMergedLeadPt","Leading jet p_{T};p_{T} (GeV);counts",1000,0,50);
    TH1F* hMergedJetsPt = new TH1F("hMergedJetsPt","Inclusive jets p_{T};p_{T} (GeV);counts",1000,0,40);
    TH1I* hMergedNumJets = new TH1I("hMergedNumJets","Number of jets found per event;# of jets;counts",100,0,100);
    TH1I* hMergedNumMCPart = new TH1I("hMergedNumMCPart","Number of MC particles in 0 jets event;# of particles;counts",1000,0,1000);
    TH1F* hMergedPairDist = new TH1F("hMergedPairDist","Production distance between K^{-} and #pi^{+};d (mm);counts",100,-1,1);
    TH1F* hMergedPartJet_Dist = new TH1F("hMergedPartJet_Dist","Distance between decay candidate and jet it contains;#DeltaR (a.u.);counts",2000,0,10);
    TH1F* hMergedKaonProdDist = new TH1F("hMergedKaonProdDist","K^{-} decay length;L (#mum);counts",10000,0,200000);
    TH1F* hMergedCandInvMass = new TH1F("hMergedCandInvMass","K^{-}#pi^{+} pair candidate invariant mass;m(K^{-}#pi^{+}) MeV;counts",3000,0,3000);
    TH1F* hMergedEnergyDiff = new TH1F("hMergedEnergyDiff","#DeltaE = |E_{jet} - E_{decay cand}|;#DeltaE (GeV);counts",2000,0,1000);

    TFile* file1;
    TH1F* hLeadPt;
    TH1F* hJetsPt;
    TH1I* hNumJets;
    TH1I* hNumMCPart;
    TH1F* hPairDist;
    TH1F* hPartJet_Dist;
    TH1F* hKaonProdDist;
    TH1F* hCandInvMass;
    TH1F* hEnergyDiff;

    //looping over files
    for(int iFile = 0; iFile < numOfRuns; iFile++){
        
        // Check if the file exists
        if (gSystem->AccessPathName(Form("AnalysisResults_"+sNumEvents+"_events_option_b%d.root",iFile))) {
            std::cerr << "File " << Form("AnalysisResults_"+sNumEvents+"_events_option_b%d.root",iFile) << " not found. Skipping." << std::endl;
            continue;
        }
        // opening file
        file1 = new TFile(Form("AnalysisResults_"+sNumEvents+"_events_option_b%d.root",iFile));

        // Leading jet pT histogram
        hLeadPt = (TH1F*)file1->Get("hLeadPt");
        // adding histogram to the merged one
        hMergedLeadPt->Add(hLeadPt);

        // histogram
        hJetsPt = (TH1F*)file1->Get("hJetsPt");
        // adding histogram to the merged one
        hMergedJetsPt->Add(hJetsPt);

        // histogram
        hNumJets = (TH1I*)file1->Get("hNumJets");
        // adding histogram to the merged one
        hMergedNumJets->Add(hNumJets);

        // histogram
        hNumMCPart = (TH1I*)file1->Get("hNumMCPart");
        // adding histogram to the merged one
        hMergedNumMCPart->Add(hNumMCPart);

        // histogram
        hPairDist = (TH1F*)file1->Get("hPairDist");
        // adding histogram to the merged one
        hMergedPairDist->Add(hPairDist);

        // histogram
        hPartJet_Dist = (TH1F*)file1->Get("hPartJet_Dist");
        // adding histogram to the merged one
        hMergedPartJet_Dist->Add(hPartJet_Dist);

        // histogram
        hKaonProdDist = (TH1F*)file1->Get("hKaonProdDist");
        // adding histogram to the merged one
        hMergedKaonProdDist->Add(hKaonProdDist);

        // histogram
        hCandInvMass = (TH1F*)file1->Get("hCandInvMass");
        // adding histogram to the merged one
        hMergedCandInvMass->Add(hCandInvMass);

        // histogram
        hEnergyDiff = (TH1F*)file1->Get("hEnergyDiff");
        // adding histogram to the merged one
        hMergedEnergyDiff->Add(hEnergyDiff);

        

        // Clean up
        file1->Close();
        delete file1;
    }

    // storing in resulting big file
    TFile* mergedOutputFile = new TFile("AnalysisResults_"+sNumEvents+"_events_option_b_Merged.root","RECREATE");
    hMergedLeadPt->Write();
    hMergedJetsPt->Write();
    hMergedNumJets->Write();
    hMergedNumMCPart->Write();
    hMergedPairDist->Write();
    hMergedPartJet_Dist->Write();
    hMergedKaonProdDist->Write();
    hMergedCandInvMass->Write();
    hMergedEnergyDiff->Write();
    mergedOutputFile->Close();

    // Clean up
    // deleting iterator histograms
    hLeadPt = nullptr;
    hJetsPt = nullptr;
    hNumJets = nullptr;
    hNumMCPart = nullptr;
    hPairDist = nullptr;
    hPartJet_Dist = nullptr;
    hKaonProdDist = nullptr;
    hCandInvMass = nullptr;
    hEnergyDiff = nullptr;
    //deleting merged histograms and file
    delete hMergedLeadPt;
    delete hMergedJetsPt;
    delete hMergedNumJets;
    delete hMergedNumMCPart;
    delete hMergedPairDist;
    delete hMergedPartJet_Dist;
    delete hMergedKaonProdDist;
    delete hMergedCandInvMass;
    delete hMergedEnergyDiff;
    delete mergedOutputFile;

    return 0;
}