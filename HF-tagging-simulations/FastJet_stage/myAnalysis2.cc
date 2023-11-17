/********************************************************************************************************************************
 *
 * Analysis for running inside FastJet environment
 * Purpose: initial HF charm jet tagging
 * Author: Christian Reckziegel
 * 
 * -> Data reading options
 *   Option a: in case output file from Pythia stored only one TTree with all events (consume more time reading)
 *   Option b: one TTree stored for each event, all particles from the same event are on the same folder/path in the output file (consume more memory in disk)
 * -> Running options
 *   Option 1 (myAnalysis.C): create a macro for ROOT that uses FastJet, load FastJet with
 *     root [0] gSystem->AddIncludePath("-I/home/christian/Softwares/fastjet-install/include");
 *     root [1] gSystem->Load("/home/christian/Softwares/fastjet-install/lib/libfastjet.so");
 *   Option 2 (myAnalysis.cc): create a macro for FastJet that uses ROOT and compile with
 *     g++ myAnalysis.cc -o myAnalysis `/home/christian/Softwares/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` `root-config --cflags --libs`
 *     or
 *     g++ myAnalysis.cc -o myAnalysis.exe $(/home/christian/Softwares/root/bin/root-config --cflags --glibs --libs) -ltbb $(/home/christian/Softwares/fastjet-install/bin/fastjet-config --cxxflags --libs)
 * 
 * This file: reading option b + running option 2
 * 
 * Obs.: dedicated function created for processing each event and struct with pointers of vectors used to store data
 * 
 * ************************************************************
 * Strategy for tagging heavy flavour (D0->K-Pi+ | PDG=421 | mass=1864.84 +- 0.17 MeV/c^2) jets
 * 1 - search for K- (PDG=-321) a Pi+ (PDG=211) asking for the PDG
 * 2 - search for geometrical association of the pairs: once a K- is found, look for Pi+ close to it
 * 3 - sum their quadri momentum and obtain invariant mass
 * 4 - if obtained mass is close from D0's mass from literature look for jets geometrically close from the eta-phi position
 * 5 - search for K- and Pi+ identical (same energy?) from obtained in jet constituents
 * 6 - HF jet tagged (identified)
 * ************************************************************
 * 
 * 
 * 
 ********************************************************************************************************************************/



#include "TFile.h" // for saving file.
#include "TTree.h" // for accessing data on TTrees
#include "TH1F.h" // for saving data on histograms
#include "TLorentzVector.h" // for four-vector calculations
#include "fastjet/ClusterSequence.hh"
#include "TString.h" // for appropriate organization of I/O read and saving
#include <chrono> // for measuring elapsed time
#include <map> // for creating mapping between kaons and pions found
#include <cmath> // for calculating absolute values

using namespace std;
using namespace fastjet;

struct myParticle{
    int myEvent, myPDG, myStatCode, myIndex;
    float myPx, myPy, myPz, myEnergy, myProdX, myProdY, myProdZ, myProdT;
};
float distanceCalc(float x1, float y1, float z1, float x2, float y2, float z2){
    return sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2));
}

struct EventData{
    std::vector<fastjet::PseudoJet>* input_particles;
    std::vector<myParticle>* kaonList;
    std::vector<myParticle>* pionList;
};

// Derived class to store user information about jets: to tag as heavy flavour
class JetUserInfo : public fastjet::PseudoJet::UserInfoBase{
    public:
        // Constructor
        JetUserInfo(bool value1): isHeavyFlavour(value1) {}
        // Destructor
        ~JetUserInfo(){};
        //Accessor method to get the heavy flavour status
        bool GetIsHeavyFlavour() const { return isHeavyFlavour; }
        void SetHeavyFlavour(bool value2) { isHeavyFlavour = value2; }

    private:
        // Data member to store heavy flavour status
        bool isHeavyFlavour = false;
};

void processEvent(int ev, TFile* inputFile, TTree* tOutPart, 
                       TH1F* hLeadPt, TH1F* hJetsPt, TH1I* hNumJets,TH1I* hNumMCPart, TH1F* hPairDist, TH1F* hCandInvMass,TH1F* hPartJet_Dist, TH1F* hEnergyDiff){
    // defining particles for accessing in the TTree storage
    int pEvent,pPDG, pStatCode, pMom1, pMom2, pDaughter1, pDaughter2;
    float pPx, pPy, pPz, pEnergy, pProdX, pProdY, pProdZ, pProdT;

    // accessing data from current event
    inputFile->GetObject(Form("top/Event_%d/EventTree",ev),tOutPart);
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

    EventData eventData;
    // Assign vectors to the struct
    // a vector of particles for each event
    eventData.input_particles = new std::vector<fastjet::PseudoJet>();
    // vector for storing decay candidates
    eventData.kaonList = new std::vector<myParticle>();
    eventData.pionList = new std::vector<myParticle>();

    // ... (my existing code to process an event)
    struct myParticle saveParticle;
    
    // search for kaons and pions and save Pseudojets
    for (int iPart = 0; iPart < tOutPart->GetEntries(); iPart++){
        tOutPart->GetEntry(iPart);
        // check if it's a final state particle
        if (pStatCode > 0){
            // collect particles data
            fastjet::PseudoJet particle(pPx, pPy, pPz, pEnergy);
            particle.set_user_index(iPart);
            eventData.input_particles->push_back(particle);
            
            // if particle observed is a K- save data
            if(pPDG == -321){
                saveParticle.myEvent = ev;
                saveParticle.myPDG = pPDG;
                saveParticle.myStatCode = pStatCode;
                saveParticle.myIndex = iPart;
                saveParticle.myPx = pPx;
                saveParticle.myPy = pPy;
                saveParticle.myPz = pPz;
                saveParticle.myEnergy = pEnergy;
                saveParticle.myProdX = pProdX;
                saveParticle.myProdY = pProdY;
                saveParticle.myProdZ = pProdZ;
                saveParticle.myProdT = pProdT;
                eventData.kaonList->push_back(saveParticle);
            } else if(pPDG == 211){// if particle observed is a pi+ save data
                saveParticle.myEvent = ev;
                saveParticle.myPDG = pPDG;
                saveParticle.myStatCode = pStatCode;
                saveParticle.myIndex = iPart;
                saveParticle.myPx = pPx;
                saveParticle.myPy = pPy;
                saveParticle.myPz = pPz;
                saveParticle.myEnergy = pEnergy;
                saveParticle.myProdX = pProdX;
                saveParticle.myProdY = pProdY;
                saveParticle.myProdZ = pProdZ;
                saveParticle.myProdT = pProdT;
                eventData.pionList->push_back(saveParticle);
            }
            
        }
        
        
    }// end particles from same event loop
    

    
    //////////////////////////////////////////////////
    // Reconstruct the jets with all event particles //
    //////////////////////////////////////////////////
    double R = 0.4;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
    fastjet::ClusterSequence clust_seq(*eventData.input_particles, jet_def);
    double ptmin = 5.0; // GeV
    vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
    // if jets were found
    if(inclusive_jets.size()){
        hLeadPt->Fill(inclusive_jets[0].pt());
        hNumJets->Fill(inclusive_jets.size());
    } else{
        hNumMCPart->Fill(eventData.input_particles->size());
    }

    // Create a map of kaon index in vector to pion index in vector
    map<int, int> pairMap;

    // search for close pair and obtain invariant mass
    for(int j = 0; j < eventData.kaonList->size(); j++){
        float distance = 1000; // 0.002
        int pionIndex = -1; // closest pion in stack
        TLorentzVector kaon4Vec(eventData.kaonList->at(j).myPx, eventData.kaonList->at(j).myPy, eventData.kaonList->at(j).myPz, eventData.kaonList->at(j).myEnergy);
        // look for closest pion from current kaon being observed
        for(int k = 0; k < eventData.pionList->size(); k++){
            // calculate distance between kaon and pion production vertex
            if(distanceCalc(eventData.kaonList->at(j).myProdX,eventData.kaonList->at(j).myProdY,eventData.kaonList->at(j).myProdZ,eventData.pionList->at(k).myProdX,eventData.pionList->at(k).myProdY,eventData.pionList->at(k).myProdZ) < distance){
                distance = distanceCalc(eventData.kaonList->at(j).myProdX,eventData.kaonList->at(j).myProdY,eventData.kaonList->at(j).myProdZ,eventData.pionList->at(k).myProdX,eventData.pionList->at(k).myProdY,eventData.pionList->at(k).myProdZ);
                pionIndex = k;
                pairMap[j] = k;
            }
        }
        // save data from closest pion found
        if(pionIndex != -1){
            TLorentzVector pion4Vec(eventData.pionList->at(pionIndex).myPx, eventData.pionList->at(pionIndex).myPy, eventData.pionList->at(pionIndex).myPz, eventData.pionList->at(pionIndex).myEnergy);
            float invMass = sqrt((kaon4Vec+pion4Vec).Mag());
            hPairDist->Fill(distance);
            hCandInvMass->Fill(invMass*1000);
        }

        eventData.kaonList->clear();
        eventData.pionList->clear();

        int jetIndex = -1;
        int constIndex = -1;
        double DeltaR = 1000;
        vector<fastjet::PseudoJet> constituents;
        // search for closest jet to pair found in eta-phi plane
        for(int iJet = 0; iJet < inclusive_jets.size(); iJet++){
            // accessing all jets
            constituents = inclusive_jets[iJet].constituents();
            hJetsPt->Fill(inclusive_jets[iJet].pt());
            // accessing constituents from jet
            for(int iConst = 0; iConst < constituents.size(); iConst++){
                // searching for geometrically closest constituent in eta-phi plane
                if(sqrt(pow(constituents[iConst].eta()-kaon4Vec.Eta(),2) + pow(constituents[iConst].phi()-kaon4Vec.Phi(),2)) < DeltaR){
                    // saving distance, jet and constituent index
                    DeltaR = sqrt(pow(constituents[iConst].eta()-kaon4Vec.Eta(),2) + pow(constituents[iConst].phi()-kaon4Vec.Phi(),2));
                    jetIndex = iJet;
                    constIndex = iConst;
                }
            }
            constituents.clear();
            
        }
        

        
        // if closest constituent from kaon was found
        if (inclusive_jets.size()){
            if(jetIndex >= 0 && constIndex >= 0 && DeltaR < 0.02){
                // mark jet with heavy flavour
                inclusive_jets[jetIndex].set_user_info(new JetUserInfo(true));
                hPartJet_Dist->Fill(DeltaR);
                constituents = inclusive_jets[jetIndex].constituents();
                hEnergyDiff->Fill(abs(constituents[constIndex].e()-kaon4Vec.E()));
            } else{
                inclusive_jets[jetIndex].set_user_info(new JetUserInfo(false));
            }
            // IMPORTANT NOTE: set_user_info(...) takes a pointer as an
            // argument. It will "own" that pointer i.e. will delete it when
            // all the PseudoJet's using it will be deleted.
        }
        constituents.clear();
    }

    // accessing heavy flavour jets
    for(int iJet = 0; iJet < inclusive_jets.size(); iJet++){
        if(inclusive_jets[iJet].has_user_info()){
            bool isHF_jet = inclusive_jets[iJet].user_info<JetUserInfo>().GetIsHeavyFlavour();
            if(isHF_jet){
                //cout << "Jet number " << iJet << " tagged with heavy flavour particle.\n";
            }
        }
        
    }

    // Clear the vectors for reuse in the next event
    eventData.input_particles->clear();
    //eventData.kaonList->clear();
    //eventData.pionList->clear();

    // Free the memory and set pointers to nullptr
    delete eventData.input_particles;
    delete eventData.kaonList;
    delete eventData.pionList;
    eventData.input_particles = nullptr;
    eventData.kaonList = nullptr;
    eventData.pionList = nullptr;

    //return eventData;
}

int main(){
    // saving starting clock
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // initial hard coded variables
    int numberEvents = 100000;
    TString sNumEvents = "100k";
    TString path = "/home/christian/cernbox/Analyses/With_Mauro/Pythia_simulation/Outside_O2/Pythia_examples/SimOutput_"+sNumEvents+"_events_option_b.root";
    cout << "Input file path = " << path << endl;
    // Creating storage data
    TH1F* hLeadPt = new TH1F("hLeadPt","Leading jet p_{T};p_{T} (GeV);counts",1000,1,1);
    TH1F* hJetsPt = new TH1F("hJetsPt","Inclusive jets p_{T};p_{T} (GeV);counts",1000,1,1);
    TH1I* hNumJets = new TH1I("hNumJets","Number of jets found per event;# of jets;counts",1000,1,1);
    TH1I* hNumMCPart = new TH1I("hNumMCPart","Number of MC particles in 0 jets event;# of particles;counts",1000,1,1);
    TH1F* hPairDist = new TH1F("hPairDist","K^{-}#pi^{+} production distance;d (mm);counts",1000,1,1);
    //TH1F* hKaonProdDist = new TH1F("hKaonProdDist","K^{-} production distance from PV;d (mm);counts",1000,1,1);
    TH1F* hCandInvMass = new TH1F("hCandInvMass","K^{-}#pi^{+} pair candidate invariant mass;m (MeV);counts",1000,1,1);
    TH1F* hPartJet_Dist = new TH1F("hPartJet_Dist","Distance between kaon candidate and closest jet;#DeltaR (a.u.);counts",1000,1,1);
    TH1F* hEnergyDiff = new TH1F("hEnergyDiff","#DeltaE = |E_{constituent} - E_{kaon}|;#DeltaE (GeV);counts",1000,1,1);

    // accessing data in simulation file
    TFile* inputFile = TFile::Open(path,"READ");
    TTree* tOutPart;

    cout << "Starting to loop over " << sNumEvents << " events...\n";
    for (int ev = 0; ev < numberEvents; ev++){

        // Call the function to process an event
        processEvent(ev, inputFile, tOutPart, 
                     hLeadPt, hJetsPt, hNumJets,hNumMCPart, hPairDist, hCandInvMass,hPartJet_Dist, hEnergyDiff);
        
        // Access the vectors through pointers
        /*std::vector<fastjet::PseudoJet>* input_particles = eventData.input_particles;
        std::vector<myParticle>* kaonList = eventData.kaonList;
        std::vector<myParticle>* pionList = eventData.pionList;

        // Free the memory and set pointers to nullptr
        delete input_particles;
        delete kaonList;
        delete pionList;
        eventData.input_particles = nullptr;
        eventData.kaonList = nullptr;
        eventData.pionList = nullptr;*/


        // show processing progress
        if(ev == numberEvents*5/100){
            cout << "5\% done\n";
            cout << "|=>                   |\n";
        } else if(ev == numberEvents*10/100){
            cout << "10\% done\n";
            cout << "|==>                  |\n";
        } else if(ev == numberEvents*15/100){
            cout << "15\% done\n";
            cout << "|===>                 |\n";
        } else if(ev == numberEvents*20/100){
            cout << "20\% done\n";
            cout << "|====>                |\n";
        } else if(ev == numberEvents*25/100){
            cout << "25\% done\n";
            cout << "|=====>               |\n";
        } else if(ev == numberEvents*30/100){
            cout << "30\% done\n";
            cout << "|======>              |\n";
        } else if(ev == numberEvents*35/100){
            cout << "35\% done\n";
            cout << "|=======>             |\n";
        } else if(ev == numberEvents*40/100){
            cout << "40\% done\n";
            cout << "|========>            |\n";
        } else if(ev == numberEvents*45/100){
            cout << "45\% done\n";
            cout << "|=========>           |\n";
        } else if(ev == numberEvents*50/100){
            cout << "50\% done\n";
            cout << "|==========>          |\n";
        } else if(ev == numberEvents*55/100){
            cout << "55\% done\n";
            cout << "|===========>         |\n";
        } else if(ev == numberEvents*60/100){
            cout << "60\% done\n";
            cout << "|============>        |\n";
        } else if(ev == numberEvents*65/100){
            cout << "65\% done\n";
            cout << "|=============>       |\n";
        } else if(ev == numberEvents*70/100){
            cout << "70\% done\n";
            cout << "|==============>      |\n";
        } else if(ev == numberEvents*75/100){
            cout << "75\% done\n";
            cout << "|===============>     |\n";
        } else if(ev == numberEvents*80/100){
            cout << "80\% done\n";
            cout << "|================>    |\n";
        } else if(ev == numberEvents*85/100){
            cout << "85\% done\n";
            cout << "|=================>   |\n";
        } else if(ev == numberEvents*90/100){
            cout << "90\% done\n";
            cout << "|==================>  |\n";
        } else if(ev == numberEvents*95/100){
            cout << "95\% done\n";
            cout << "|===================> |\n";
        } else if(ev == numberEvents-1){
            cout << "100\% done\n";
            cout << "|====================>|\n";
            cout << endl;
        }

    }// end events loop
    
    cout << "Events loop ended.\n" 
         << "Closing input file and saving data to output file...\n";

    // Closing input file
    inputFile->Close();

    // creating file for storing results
    TFile f("AnalysisResults_"+sNumEvents+"_events_option_b.root","recreate");
    hLeadPt->Write();
    hJetsPt->Write();
    hNumJets->Write();
    hNumMCPart->Write();
    hPairDist->Write();
    hPartJet_Dist->Write();
    hCandInvMass->Write();
    hEnergyDiff->Write();
    cout << "Done.\n";
    // ending clock
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed: " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000 << " seconds.\n";

    return 0;
}