/**
 * Simulation for running inside Pythia 8.3 environment
 * Purpose: initial HF charm jet tagging
 * Author: Christian Reckziegel
 * 
 * -> Data storing options
 *   Option a: output file from Pythia stored in only one TTree with all events (consume more time reading)
 *   Option b: one TTree stored for each event, all particles from the same event are on the same folder/path in the output file (consume more memory in disk)
 * 
 * This file: storing option b
 * 
 * main95_mainCode.cc
 * 
*/

#include "Pythia8/Pythia.h" // access to Pythia objects.

// ROOT
#include "TFile.h" // for saving file.
#include "TTree.h" // for saving data on TTrees
#include "TString.h" // for appropriate organization of I/O read and saving
#include "TDirectory.h" // for organizing events in TTree
#include "TCanvas.h" // for looking at graphical results on the fly
#include "TH1.h" // for plotting 1D histograms
#include "TH2.h" // for plotting 2D histograms
#include <chrono> // for measuring elapsed time

using namespace Pythia8;
using namespace std;

// allow simplified notation.
int main() {
  // saving starting clock
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  int maxEvents = 100000; // The number of events to run: 100k.
  TString sNumEvents = "100k";
  TString path = "SimOutput_"+sNumEvents+"_events_option_b.root";
  //TH2F* hPairDistPt = new TH2F("hPairDistPt","D^{0} decay products production distance in function of D^{0} p_{T};p_{T} (GeV/c);d (mm);counts",1000,1,1,1000,1,1);
  //TH1F* hPairDist = new TH1F("hPairDist","D^{0} decay products production distance;d (mm);counts",1000,1,1);
  //TH1F* hD0Pt = new TH1F("hD0Pt","D^{0} p_{T};p_{T} (GeV/c);counts",1000,1,1);
  // defining particles for TTree storage
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
  TFile *fOutPart = TFile::Open(path,"RECREATE");
  //TFile *fOutPart = TFile::Open("Ignore.root","RECREATE");
  TDirectory *topDir = fOutPart->mkdir("top");
  topDir->cd();
  // creating directory structure for runs in output file
  TDirectory *runDirs;
  TTree *tOutPart = tOutPart = new TTree("EventTree","Tree containing event particles data from event");
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
    runDirs = topDir->mkdir(Form("Event_%d",iEvent));
	  runDirs->cd();
    

    // Generate next event;
    // Produce the next event, returns true on success.
    if(!pythia.next()) {
      // Any error handling goes here.
      cout << "Error generating event " << iEvent << endl;
      continue;
    }
    //cout << "Looping through particles from event " << iEvent << endl;
    for (int i = 0; i < pythia.event.size(); ++i){
      // if D0 particle has been found
      if(pythia.event[i].id() == 421){
        int daughter1Index = pythia.event[i].daughter1();
        int daughter2Index = pythia.event[i].daughter2();
        //cout << "D0 found!\n";
        // check for kaon daughter
        if(pythia.event[daughter1Index].id() == -321){
          if(pythia.event[daughter2Index].id() == 211){
            //double distance = sqrt(pow(pythia.event[daughter1Index].xProd()-pythia.event[daughter2Index].xProd(),2) + pow(pythia.event[daughter1Index].yProd()-pythia.event[daughter2Index].yProd(),2) + pow(pythia.event[daughter1Index].zProd()-pythia.event[daughter2Index].zProd(),2));
            //hPairDist->Fill(distance);
            //double D0_pt = sqrt(pow(pythia.event[i].px(),2) + pow(pythia.event[i].py(),2));
            //hPairDistPt->Fill(D0_pt, distance);
            //hD0Pt->Fill(D0_pt);
          }
        }else if(pythia.event[daughter1Index].id() == 211){
          if(pythia.event[daughter2Index].id() == -321){
            //double distance = sqrt(pow(pythia.event[daughter1Index].xProd()-pythia.event[daughter2Index].xProd(),2) + pow(pythia.event[daughter1Index].yProd()-pythia.event[daughter2Index].yProd(),2) + pow(pythia.event[daughter1Index].zProd()-pythia.event[daughter2Index].zProd(),2));
            //hPairDist->Fill(distance);
            //double D0_pt = sqrt(pow(pythia.event[i].px(),2) + pow(pythia.event[i].py(),2));
            //hPairDistPt->Fill(D0_pt, distance);
            //hD0Pt->Fill(D0_pt);
          }
        }
      }

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
      tOutPart->Fill(); // fill TTree
    }
    // Analyse event; fill histograms etc.
    tOutPart->Write(); // save TTree
    tOutPart->Reset(); // reset TTree for next event
  } // End event loop.
  
  // --- Calculate final statistics ---
  pythia.stat();

  // Save TTree
  
  fOutPart->Close();
  delete fOutPart;

  //TFile *fView = TFile::Open("Teste.root","RECREATE");
  //hPairDist->Write();
  //hPairDistPt->Write();
  //hD0Pt->Write();
  //fView->Close();
  //delete fView;
  
  // ending clock
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time elapsed: " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000 << " seconds.\n";

  return 0;
}
