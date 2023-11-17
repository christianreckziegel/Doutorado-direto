/*
 *
 * Analysis for running inside ROOT environment
 * 
 * Data reading options
 * Option a: in case output file from Pythia stored only one TTree with all events (consume more time reading)
 * Option b: one TTree stored for each event, all particles from the same event are on the same folder/path in the output file (consume more memory in disk)
 * 
 * Running options
 * Option 1 (myAnalysis.C): create a macro for ROOT that uses FastJet, load FastJet with
 *  root [0] gSystem->AddIncludePath("-I/home/christian/Softwares/fastjet-install/include");
 *  root [1] gSystem->Load("/home/christian/Softwares/fastjet-install/lib/libfastjet.so");
 * Option 2 (myAnalysis.cc): create a macro for FastJet that uses ROOT and compile with
 *  g++ myAnalysis.cpp -o myAnalysis `/home/christian/Softwares/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` `root-config --cflags --libs` -lRecursiveTools -lEnergyCorrelator -lNsubjettiness
 *
 * This file: reading option a + running option 1
 * 
 */

#include "TFile.h" // for saving file.
#include "TTree.h" // for saving data on TTrees
#include "fastjet/ClusterSequence.hh"

using namespace std;

void myAnalysis(const char* path="/home/christian/cernbox/Analyses/With_Mauro/Pythia_simulation/Outside_O2/Pythia_examples/SimOutput_10k_events.root"){
    // loading FastJet library
    //gSystem->AddIncludePath("-I/home/christian/Softwares/fastjet-install/include");
    //gSystem->Load("/home/christian/Softwares/fastjet-install/lib/libfastjet.so");
    int numberEvents = 10000;
    // Creating storage data
    TCanvas* cLeadPt = new TCanvas("cLeadPt","Canvas for leading jet p_{T}");
    TH1D* hLeadPt = new TH1D("hLeadPt","Leading jet p_{T};p_{T};counts",1000,1,1);

    // defining particles for accessing in the TTree storage
    int pEvent,pPDG, pStatCode, pMom1, pMom2, pDaughter1, pDaughter2;
    float pPx, pPy, pPz, pEnergy, pProdX, pProdY, pProdZ, pProdT;

    // accessing data in simulation file
    TFile inputFile(path, "READ");
    TTree* tOutPart = (TTree*)inputFile.Get("EventTree");
    tOutPart->SetBranchAddress("partEvent",&pEvent);
    tOutPart->SetBranchAddress("partPDG",&pPDG);
    tOutPart->SetBranchAddress("partStatCode",&pStatCode);
    tOutPart->SetBranchAddress("partMom1",&pMom1);
    tOutPart->SetBranchAddress("partMom2",&pMom2);
    tOutPart->SetBranchAddress("partDaughter1",&pDaughter1);
    tOutPart->SetBranchAddress("partDaughter2",&pDaughter2);
    tOutPart->SetBranchAddress("partPx",&pPx);
    tOutPart->SetBranchAddress("partPy",&pPy);
    tOutPart->SetBranchAddress("partPz",&pPz);
    tOutPart->SetBranchAddress("partEnergy",&pEnergy);
    tOutPart->SetBranchAddress("partProdX",&pProdX);
    tOutPart->SetBranchAddress("partProdY",&pProdY);
    tOutPart->SetBranchAddress("partProdZ",&pProdZ);
    tOutPart->SetBranchAddress("partProdT",&pProdT);

    for (int ev = 0; ev < numberEvents; ev++){
        // get an array of particles for each event
        vector<fastjet::PseudoJet> input_particles;
        int index = 0;
        
        for (int entry = 0; entry < tOutPart->GetEntries(); entry++){
            tOutPart->GetEntry(entry);
            // check if particle belongs to current event being investigated (the one from the loop index)
            if (pEvent == ev && pStatCode > 0){// only collect final state particles
                // collect particles data
                fastjet::PseudoJet particle(pPx, pPy, pPz, pEnergy);
                particle.set_user_index(index);
                input_particles.push_back(particle);
                index++;
            }
            
        }

        // do FastJet analysis
        //////////////////////////////////////////////////
        // Reconstruct the jet with all event particles //
        //////////////////////////////////////////////////
        cout << "Reconstructing jets from event " << ev << endl;
        double R = 0.4;
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
        fastjet::ClusterSequence clust_seq(input_particles, jet_def);
        double ptmin = 10.0; // GeV
        vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
        
        // if jets were found
        if(inclusive_jets.size()){
            hLeadPt->Fill(inclusive_jets[0].perp());
        }

    }
    

    
    cLeadPt->cd();
    hLeadPt->Draw();


}


int main(){
    myAnalysis(const char* path="/home/christian/cernbox/Analyses/With_Mauro/Pythia_simulation/Outside_O2/Pythia_examples/SimOutput_10k_events.root");
    return 0;
}