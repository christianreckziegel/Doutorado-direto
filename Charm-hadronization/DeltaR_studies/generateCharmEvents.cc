
/**
 * This code generates charm events using Pythia8, stores the kinematic information of final state particles in a TTree, and saves it to a ROOT file. The settings for the event generation are read from a card file (cardfile.cmd). The code checks for the presence of D0 particles in the events and classifies them as prompt or from b decay before filling the tree.
 * The main steps are:
 * 1. Initialize Pythia and read settings from a card file.
 * 2. Define a TTree to store the kinematic information of final state particles.
 * 3. Loop over a specified number of events, generate each event, and analyze the particles.
 * 4. Check for the presence of D0 particles and classify them as prompt or from b decay.
 * 5. Fill the TTree with the kinematic information of final state particles.
 * 6. Save the TTree to a ROOT file.
 * 
 * To do: only store the D0 events in which these would decay specifically to K- pi+ (or the charge conjugate), since this would bias the final distribution results
 * To do: run it with 5.02 and 13.6 TeV and compare the results
*/

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "commonFunctions.h"
#include <iomanip> // time precision

using namespace Pythia8;

void generateCharmEvents() {
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // --- Initialization ---
    Pythia pythia;
    Event& event = pythia.event;

    // Read in settings
    pythia.readFile("cardfile.cmd"); // via file
    int numberOfEvents = pythia.mode("Main:numberOfEvents") / 1000; // k
    int maxEvents = numberOfEvents * 1000;
    double eCM = pythia.info.eCM(); // Get the center of mass energy from the settings
    double aliceAcceptanceEta = 0.9; // ALICE central barrel acceptance in pseudorapidity

    // Define TTrees and variables
    TTree* tFinalStateParticles = new TTree("tFinalStateParticles","Final state particles from Pythia events");
    double px, py, pz, energy;
    int eventNumber, pdg;
    int isPrompt; // 1 if prompt, 0 if from decay from b, -1 if not a D0 particle
    tFinalStateParticles->Branch("px",&px,"px/D");
    tFinalStateParticles->Branch("py",&py,"py/D");
    tFinalStateParticles->Branch("pz",&pz,"pz/D");
    tFinalStateParticles->Branch("energy",&energy,"energy/D");
    tFinalStateParticles->Branch("eventNumber",&eventNumber,"eventNumber/I");
    tFinalStateParticles->Branch("pdg",&pdg,"pdg/I");
    tFinalStateParticles->Branch("isPrompt",&isPrompt,"isPrompt/I");

    pythia.init();

    // --- The event loop ---
    int iEvent = 0;
    // Generate maxEvents with each one containing at least one D0 particle
    while (iEvent < maxEvents) {
        
        // Generate next evet (error message in case it fails)
        if (!pythia.next()) {
            std::cout << "Error generating event " << iEvent << "!" << std::endl;
        }

        // Analyse event: loop over event particles
        // check if there is a D0 particle in the event
        bool hasD0 = eventHasD0(pythia.event);
        if (!hasD0) {
            continue; // skip to next event if no D0 is found
        }
        
        for (int iParticle = 0; iParticle < pythia.event.size(); iParticle++) {
            //
            // Store final state particle data on tree
            if (event[iParticle].isFinal() && (abs(event[iParticle].eta()) < aliceAcceptanceEta)) { // Only consider final state particles within the ALICE central barrel acceptance (|eta| < 0.9)
                px = event[iParticle].px();
                py = event[iParticle].py();
                pz = event[iParticle].pz();
                energy = event[iParticle].e();
                pdg = event[iParticle].id();
                eventNumber = iEvent;

                // Check if the D0 is prompt or from b decay
                if (abs(pdg) == 421) { // is it a D0 or D0bar particle?
                    isPrompt = hasBHadronAncestor(iParticle, pythia.event) ? 0 : 1; // 0 if from b decay, 1 if prompt
                } else {
                    isPrompt = -1; // Not a D0 particle
                    if (!event[iParticle].isCharged()) {
                        continue; // skip neutral particles that are not D0, since ALICE's TPC only detects charged particle tracks
                    }
                    
                }
                
                // Store the particle data in the tree
                tFinalStateParticles->Fill();
            }
        }

        iEvent++;
    }
    
    // Store the tree in a ROOT file
    TFile* outFile = new TFile(Form("pythiaCharmEvents_%dTeV_%dk.root", static_cast<int>(eCM), numberOfEvents), "RECREATE");
    tFinalStateParticles->Write();
    std::cout << "Data stored to file " << outFile->GetName() << "." << std::endl;
    outFile->Close();

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed 
         << time_taken/60 << std::setprecision(5); 
    std::cout << " min " << std::endl;
}

int main() {
    generateCharmEvents();
    return 0;
}