/**
 * Simulation for running inside Pythia 8.3 environment
 * Purpose: initial HF charm jet tagging
 * Author: Christian Reckziegel
 * 
 * 
 * -> Data storing options
 *   Option a: output file from Pythia stored in only one TTree with all events (consume more time reading)
 *   Option b: one TTree stored for each event, all particles from the same event are on the same folder/path in the output file (consume more memory in disk)
 * 
 * This file: storing option a
 * 
 * main95_option_a.cc
 * 
*/

#include "Pythia8/Pythia.h" // access to Pythia objects.

// ROOT
#include "TFile.h" // for saving file.
#include "TTree.h" // for saving data on TTrees
#include "TDirectory.h" // for organizing events in TTree

using namespace Pythia8;
using namespace std;

// allow simplified notation.
int main() {
  // saving starting clock
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  int maxEvents = 1000; // The number of events to run: 100k.
  // defining particles for TTree storage
  int partEvent;
  int partPDG;
  int partStatCode;
  int partMom1;
  int partMom2;
  int partDaughter1;
  int partDaughter2;
  float partPx;
  float partPy;
  float partPz;
  float partEnergy;
  float partProdX;
  float partProdY;
  float partProdZ;
  float partProdT;
  TFile *fOutPart = TFile::Open("SimOutput_1k_events.root","RECREATE");

  // creating directory structure for runs in output file
  TTree *tOutPart = tOutPart = new TTree("EventTree","Tree containing event particles data from event");
  tOutPart->Branch("partEvent",&partEvent,"partEvent/I");
  tOutPart->Branch("partPDG",&partPDG,"partPDG/I");
  tOutPart->Branch("partStatCode",&partStatCode,"partStatCode/I");
  tOutPart->Branch("partMom1",&partMom1,"partMom1/I");
  tOutPart->Branch("partMom2",&partMom2,"partMom2/I");
  tOutPart->Branch("partDaughter1",&partDaughter1,"partDaughter1/I");
  tOutPart->Branch("partDaughter2",&partDaughter2,"partPDG/I");
  tOutPart->Branch("partPx",&partPx,"partPx/F");
  tOutPart->Branch("partPy",&partPy,"partPy/F");
  tOutPart->Branch("partPz",&partPz,"partPz/F");
  tOutPart->Branch("partEnergy",&partEnergy,"partEnergy/F");
  tOutPart->Branch("partProdX",&partProdX,"partProdX/F");
  tOutPart->Branch("partProdY",&partProdY,"partProdY/F");
  tOutPart->Branch("partProdZ",&partProdZ,"partProdZ/F");
  tOutPart->Branch("partProdT",&partProdT,"partProdT/F");

  // --- Initialization ---
  Pythia pythia;
  // Define Pythia object.
  Event& event = pythia.event; // quick access to current event.
  // Read in settings via file
  pythia.readFile("cardfile.cmnd");

  pythia.init();
  // Initialize
  // --- The event loop ---
  for(int iEvent = 0; iEvent < maxEvents; iEvent++){

    // Generate next event;
    // Produce the next event, returns true on success.
    if(!pythia.next()) {
      // Any error handling goes here.
      cout << "Error generating event " << iEvent << endl;
      continue;
    }
    //cout << "Looping through particles from event " << iEvent << endl;
    for (int i = 0; i < pythia.event.size(); ++i){
      partEvent = iEvent;
      partPDG = pythia.event[i].id();
      partStatCode = pythia.event[i].status();
      partMom1 = pythia.event[i].mother1();
      partMom2 = pythia.event[i].mother2();
      partDaughter1 = pythia.event[i].daughter1();
      partDaughter2 = pythia.event[i].daughter2();
      partPx = pythia.event[i].px();
      partPy = pythia.event[i].py();
      partPz = pythia.event[i].pz();
      partEnergy = pythia.event[i].e();
      partProdX = pythia.event[i].xProd();
      partProdY = pythia.event[i].yProd();
      partProdZ = pythia.event[i].zProd();
      partProdT = pythia.event[i].tProd();
      tOutPart->Fill();
    }
    // Analyse event; fill histograms etc.
    
  } // End event loop.
  
  // --- Calculate final statistics ---
  pythia.stat();

  // Save TTree
  tOutPart->Write();
  fOutPart->Close();
  delete fOutPart;
  
  // ending clock
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time elapsed: " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000 << " seconds.\n";

  return 0;
}
