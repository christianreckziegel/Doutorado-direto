#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include "commonFunctions.h"
#include "MyInfo.h"
#include "TFile.h"
#include "TTree.h"
#include <iomanip> // time precision

using namespace fastjet;
using namespace std;

void reconstructCharmJets() {
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    double comEnergy = 13.6; // 5.02 TeV or 13.6 TeV
    int numberOfEvents = 10; // k

    // Define output TFile and TTree
    TFile* outputFile = new TFile(Form("fastjetCharmEvents_%dTeV_%dk.root", static_cast<int>(comEnergy), numberOfEvents), "RECREATE");
    TTree* tOutput = new TTree("tCharmJets","D0 tagged jets from Pythia events");
    float axisDistance, jetPt, jetEta, jetPhi, jetMass;
    float hfPt, hfEta, hfPhi, hfMass, hfY;
    int jetNConst, iEventNumber;
    int hfPrompt; // 1 if the D0 is prompt, 0 if from b decay, -1 if not a D0 particle
    tOutput->Branch("jetHfDist",&axisDistance,"jetHfDist/F");
    tOutput->Branch("jetPt",&jetPt,"jetPt/F");
    tOutput->Branch("jetEta",&jetEta,"jetEta/F");
    tOutput->Branch("jetPhi",&jetPhi,"jetPhi/F");
    tOutput->Branch("jetMass",&jetMass,"jetMass/F");
    tOutput->Branch("jetNConst",&jetNConst,"jetNConst/I");
    tOutput->Branch("hfPt",&hfPt,"hfPt/F");
    tOutput->Branch("hfEta",&hfEta,"hfEta/F");
    tOutput->Branch("hfPhi",&hfPhi,"hfPhi/F");
    tOutput->Branch("hfMass",&hfMass,"hfMass/F");
    tOutput->Branch("hfY",&hfY,"hfY/F");
    tOutput->Branch("hfPrompt",&hfPrompt,"hfPrompt/I");
    tOutput->Branch("iEventNumber",&iEventNumber,"iEventNumber/I");

    // Open the ROOT file containing the generated events
    TFile* inputFile = TFile::Open(Form("../pythiaCharmEvents_%dTeV_%dk.root", static_cast<int>(comEnergy), numberOfEvents));
    if (!inputFile || inputFile->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }
    // Get the TTree from the file
    TTree* tInput = (TTree*)inputFile->Get("tFinalStateParticles");
    if (!tInput) {
        std::cout << "Error getting TTree!" << std::endl;
        return;
    }

    // Set up variables to read from the TTree
    double px, py, pz, energy;
    int eventNumber, pdg;
    int isPrompt; // 1 if prompt, 0 if from decay from b, -1 if not a D0 particle
    tInput->SetBranchAddress("px", &px);
    tInput->SetBranchAddress("py", &py);
    tInput->SetBranchAddress("pz", &pz);
    tInput->SetBranchAddress("energy", &energy);
    tInput->SetBranchAddress("eventNumber", &eventNumber);
    tInput->SetBranchAddress("pdg", &pdg);
    tInput->SetBranchAddress("isPrompt", &isPrompt);

    // Define a vector of PseudoJets for the current event
    vector<PseudoJet> particles;
    double R = 0.4;
    JetDefinition jet_def(antikt_algorithm, R);

    // Loop over the events in the TTree
    int nEntries = tInput->GetEntries();
    for (int iEntry = 0; iEntry < nEntries; iEntry++) {
        tInput->GetEntry(iEntry);
        
        particles.push_back(PseudoJet(px, py, pz, energy));
        // Store the PDG code and promptness info in the user info of the PseudoJet
        particles.back().set_user_info(new MyInfo(pdg, isPrompt, abs(pdg) == 421, eventNumber));
        const MyInfo& testInfo = particles.back().user_info<MyInfo>();
                
        // If this is the last particle of the event, reconstruct jets for the current event and fill the output tree
        if (isLastEventParticle(iEntry, tInput, eventNumber)) {
            
            // Cluster the particles into jets
            ClusterSequence cs(particles, jet_def);
            vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

            // Look for the D0 jet: which jet contains the D0 as a constituent?
            for (auto& jet : jets) {
                bool hasD0 = false;
                for (const auto& constituent : jet.constituents()) {
                    // const MyInfo* info = dynamic_cast<const MyInfo*>(constituent.user_info());
                    const MyInfo& info = constituent.user_info<MyInfo>();
                    if (abs(info.GetPDG()) == 421) { // Check if the constituent is a D0 or D0bar particle
                        hasD0 = true;
                        // Fill the output tree with the jet and D0 information
                        // axisDistance = jet.delta_R(constituent);
                        axisDistance = calculateDeltaR(jet, constituent);
                        jetPt = jet.pt();
                        jetEta = jet.eta();
                        jetPhi = jet.phi();
                        jetMass = jet.m();
                        jetNConst = jet.constituents().size();
                        hfPt = constituent.pt();
                        hfEta = constituent.eta();
                        hfPhi = constituent.phi();
                        hfMass = constituent.m();
                        hfY = constituent.rapidity();
                        hfPrompt = info.GetPromptness();
                        iEventNumber = info.GetPythiaEventNumber();
                        tOutput->Fill();
                        break; // We found the D0 in this jet, no need to check other constituents
                    }
                } // constituents loop
                // Store whether the jet has a D0 constituent in the user info of the jet for later analysis
                jet.set_user_info(new MyInfo(0, hfPrompt, hasD0, iEventNumber));
            } // jets loop

            // Clear the particles vector for the next event
            particles.clear();
        } // reconstruct jets after all particles from the same event were read
    }

    // Save the output tree to a new ROOT file
    outputFile->cd();  // Add this line!
    tOutput->Write();
    std::cout << "Data stored to file " << outputFile->GetName() << "." << std::endl;
    outputFile->Close();
    inputFile->Close();

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed 
         << time_taken/60 << std::setprecision(5); 
    std::cout << " min " << std::endl;
}

int main() {
    reconstructCharmJets();
    return 0;
}