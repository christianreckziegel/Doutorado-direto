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
 * Obs.: save data about PDG and D0 decay product on JetUserInfo class
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

// Derived class to store user information about jets: to tag as heavy flavour
class JetUserInfo : public fastjet::PseudoJet::UserInfoBase{
    public:
        // TODO: create new constructor
        // Constructor 1
        JetUserInfo(): isHeavyFlavour(false),  isKaon(false), isPion(false){}
        // Constructor 2
        JetUserInfo(bool value1, bool value2, bool value3): isHeavyFlavour(value1),  isKaon(value2), isPion(value3){}
        // Destructor
        virtual ~JetUserInfo(){};
        //Accessor method to get the heavy flavour status
        bool GetIsHeavyFlavour() const { return isHeavyFlavour; }
        void SetHeavyFlavour(bool value4) { isHeavyFlavour = value4; }
        bool GetIsKaon() const { return isKaon; }
        void SetKaon(bool value5) { isKaon = value5; }
        bool GetIsPion() const { return isPion; }
        void SetPion(bool value6) { isPion = value6; }
        float GetXProd() const { return xProd; }
        void SetXProd(float newXProd) {xProd = newXProd; }
        float GetYProd() const { return yProd; }
        void SetYProd(float newYProd) {yProd = newYProd; }
        float GetZProd() const { return zProd; }
        void SetZProd(float newZProd) {zProd = newZProd; }

    private:
        // Data member to store heavy flavour status
        bool isHeavyFlavour = false; // if asked to input_particles/constituents tells if it's a decay candidate, if asked to jets tells if it's a HF jet candidate
        bool isKaon = false;
        bool isPion = false;
        float xProd = -1000;
        float yProd = -1000;
        float zProd = -1000;
};

int main(){
    // saving starting clock
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // initial hard coded variables
    int numberEvents = 1000;
    TString sNumEvents = "1k";
    TString path = "/home/christian/cernbox/Analyses/With_Mauro/Pythia_simulation/Outside_O2/Pythia_examples/SimOutput_"+sNumEvents+"_events_option_b.root";
    cout << "Input file path = " << path << endl;
    // Creating storage data
    TH1F* hLeadPt = new TH1F("hLeadPt","Leading jet p_{T};p_{T} (GeV);counts",1000,1,1);
    TH1F* hJetsPt = new TH1F("hJetsPt","Inclusive jets p_{T};p_{T} (GeV);counts",1000,1,1);
    TH1I* hNumJets = new TH1I("hNumJets","Number of jets found per event;# of jets;counts",1000,1,1);
    TH1I* hNumMCPart = new TH1I("hNumMCPart","Number of MC particles in 0 jets event;# of particles;counts",1000,1,1);
    TH1F* hPairDist = new TH1F("hPairDist","Production distance between K^{-} and #pi^{+};d (mm);counts",1000,1,1);
    TH1F* hKaonProdDist = new TH1F("hKaonProdDist","K^{-} decay length;L (#mum);counts",1000,1,1);
    TH1F* hCandInvMass = new TH1F("hCandInvMass","K^{-}#pi^{+} pair candidate invariant mass;m(K^{-}#pi^{+}) MeV;counts",1000,1,1);
    TH1F* hPartJet_Dist = new TH1F("hPartJet_Dist","Distance between decay candidate and jet it contains;#DeltaR (a.u.);counts",1000,1,1);
    TH1F* hEnergyDiff = new TH1F("hEnergyDiff","#DeltaE = |E_{jet} - E_{decay cand}|;#DeltaE (GeV);counts",1000,1,1);

    // defining particles for accessing in the TTree storage
    int pEvent,pPDG, pStatCode, pMom1, pMom2, pDaughter1, pDaughter2;
    float pPx, pPy, pPz, pEnergy, pProdX, pProdY, pProdZ, pProdT;

    // accessing data in simulation file
    TFile* inputFile = TFile::Open(path,"READ");
    if(!inputFile || inputFile->IsZombie()){
        cout << "Error opening file, finishing the process...\n";
        return 0;
    }

    TTree* tOutPart;
    cout << "Starting to loop over " << sNumEvents << " events...\n";
    for (int ev = 0; ev < numberEvents; ev++){
        
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

        // get a vector of particles for each event
        vector<fastjet::PseudoJet> input_particles;
        //vector<JetUserInfo*> user_info_pointers;
        //JetUserInfo* userInfoHandler;

        // index for accessing final state particle stored
        int finalStateIndex = 0;
        // search for kaons and pions and save Pseudojets
        for (int iPart = 0; iPart < tOutPart->GetEntries(); iPart++){
            tOutPart->GetEntry(iPart); // TODO: check that this index is appropriate for managing the input_particles vector inside the if conditional
            // check if it's a final state particle
            if (pStatCode > 0){
                // collect particles data
                fastjet::PseudoJet particle(pPx, pPy, pPz, pEnergy);
                particle.set_user_index(iPart);
                input_particles.push_back(particle);
                // if particle observed is a K- save data
                if(pPDG == -321){
                    JetUserInfo* inputPartUIKaon = new JetUserInfo(false, true, false); // mark as a kaon - option 2
                    input_particles[finalStateIndex].set_user_info(inputPartUIKaon); // option 2
                    //input_particles[finalStateIndex].set_user_info(new JetUserInfo(false, true, false)); // mark as a kaon - option 1
                    inputPartUIKaon->SetXProd(pProdX);
                    inputPartUIKaon->SetYProd(pProdY);
                    inputPartUIKaon->SetZProd(pProdZ);
                    //double kaonProdDistance = 1000*sqrt(pow(pProdX,2) + pow(pProdY,2) + pow(pProdZ,2)); // in micrometers
                } else if(pPDG == 211){// if particle observed is a pi+ save data
                    JetUserInfo* inputPartUIPion = new JetUserInfo(false, false, true); // mark as a pion
                    input_particles[finalStateIndex].set_user_info(inputPartUIPion); // mark as a pion
                    //input_particles[finalStateIndex].set_user_info(new JetUserInfo(false, false, true)); // mark as a pion
                    inputPartUIPion->SetXProd(pProdX);
                    inputPartUIPion->SetYProd(pProdY);
                    inputPartUIPion->SetZProd(pProdZ);
                    //double pionProdDistance = 1000*sqrt(pow(pProdX,2) + pow(pProdY,2) + pow(pProdZ,2)); // in micrometers
                }
                finalStateIndex++;
            }
            
            
        }// end particles from same event loop

        // IMPORTANT NOTE: set_user_info(...) takes a pointer as an
        // argument. It will "own" that pointer i.e. will delete it when
        // all the PseudoJet's using it will be deleted.

        /////////////////////////////////////////////////////
        // Search for close pair and obtain invariant mass //
        /////////////////////////////////////////////////////

        // Create a map of kaon index in vector to pion index in vector
        map<int, int> pairMap;
        for(int jKaon = 0; jKaon < input_particles.size(); jKaon++){
            // protection for not trying to access non existent data
            if(input_particles[jKaon].has_user_info()){
                // ask if it's a kaon
                if(input_particles[jKaon].user_info<JetUserInfo>().GetIsKaon()){
                    float distance = 0.002;//distance between kaon and pion produced cut of 0.002 mm
                    int pionIndex = -1;// closest pion in stack
                    TLorentzVector kaon4Vec(input_particles[jKaon].px(), input_particles[jKaon].py(), input_particles[jKaon].pz(), input_particles[jKaon].e());
                    // look for closest pion
                    for(int kPion = 0; kPion < input_particles.size(); kPion++){
                        // protection for against trying to access non existent data
                        if(input_particles[kPion].has_user_info()){
                            if(input_particles[kPion].user_info<JetUserInfo>().GetIsPion()){
                                // calculate distance between kaon and pion production vertex
                                float distDelta = distanceCalc(input_particles[jKaon].user_info<JetUserInfo>().GetXProd(),input_particles[jKaon].user_info<JetUserInfo>().GetYProd(),input_particles[jKaon].user_info<JetUserInfo>().GetZProd(), input_particles[kPion].user_info<JetUserInfo>().GetXProd(),input_particles[kPion].user_info<JetUserInfo>().GetYProd(),input_particles[kPion].user_info<JetUserInfo>().GetZProd());
                                if(distDelta < distance){
                                    distance = distDelta;
                                    pionIndex = kPion;
                                    pairMap[jKaon] = kPion;
                                }
                            }
                        }
                        
                    }
                    // closest pion found
                    if(pionIndex > -1){
                        TLorentzVector pion4Vec(input_particles[pionIndex].px(), input_particles[pionIndex].py(), input_particles[pionIndex].pz(), input_particles[pionIndex].e());
                        // kinematical D0 production mass cut
                        const double d0Mass = 1.864; // The known mass of the D0 meson in GeV/c^2
                        //if(true){
                        if((kaon4Vec+pion4Vec).E() >= d0Mass){
                            
                            float invMass = (kaon4Vec+pion4Vec).Mag(); // square root of the invariant
                            hPairDist->Fill(distance);
                            // kaon production distance in micrometers
                            double kaonProdDist = 1000*sqrt(pow(input_particles[jKaon].user_info<JetUserInfo>().GetXProd(),2) + pow(input_particles[jKaon].user_info<JetUserInfo>().GetYProd(),2) + pow(input_particles[jKaon].user_info<JetUserInfo>().GetZProd(),2));
                            
                            // apply topological cut
                            if(kaonProdDist > 0){
                                // tag HF decay candidates
                                //const JetUserInfo* inputPartUserInfoKaon = &(input_particles[jKaon].user_info<JetUserInfo>()); // option 2
                                JetUserInfo* inputPartUserInfoKaon = const_cast<JetUserInfo*>(&(input_particles[jKaon].user_info<JetUserInfo>())); // option 1
                                inputPartUserInfoKaon->SetHeavyFlavour(true); // option 1
                                //const JetUserInfo* inputPartUserInfoPion = &(input_particles[pionIndex].user_info<JetUserInfo>()); // option 2
                                JetUserInfo* inputPartUserInfoPion = const_cast<JetUserInfo*>(&(input_particles[pionIndex].user_info<JetUserInfo>())); // option 1
                                inputPartUserInfoPion->SetHeavyFlavour(true); // option1
                                hKaonProdDist->Fill(kaonProdDist); 
                                hCandInvMass->Fill(invMass*1000);
                            }

                            
                        }
                        
                    }
                }
            }
        }
        
        ///////////////////////////////////////////////////
        // Reconstruct the jets with all event particles //
        ///////////////////////////////////////////////////
        double R = 0.4;
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
        fastjet::ClusterSequence clust_seq(input_particles, jet_def);
        double ptmin = 0.0; // 5 GeV for HF
        vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
        // if jets were found
        if(inclusive_jets.size()){
            hLeadPt->Fill(inclusive_jets[0].pt());
            hNumJets->Fill(inclusive_jets.size());
        } else{
            hNumMCPart->Fill(input_particles.size());
        }
        
        ///////////////////////////////////////////////////
        //             Tag HF jet candidate              //
        ///////////////////////////////////////////////////
        double DeltaR = 1000;
        vector<fastjet::PseudoJet> constituents;
        for(int iJet = 0; iJet < inclusive_jets.size(); iJet++){
            // accessing all jets
            constituents = inclusive_jets[iJet].constituents();
            hJetsPt->Fill(inclusive_jets[iJet].pt());
            // accessing constituents from jet
            for(int iConst = 0; iConst < constituents.size(); iConst++){
                if(constituents[iConst].has_user_info()){
                    // check if it has D0 decay candidates
                    if(constituents[iConst].user_info<JetUserInfo>().GetIsHeavyFlavour()){
                        // tag jet as HF jet candidate
                        //JetUserInfo* inclusiveJUInfo = const_cast<JetUserInfo*>(&(inclusive_jets[iJet].user_info<JetUserInfo>())); // option 1
                        //inclusiveJUInfo->SetHeavyFlavour(true); // option 1
                        inclusive_jets[iJet].set_user_info(new JetUserInfo(true, false, false)); // mark as jet containing HF decay product candidates
                        //inclusiveJUInfo = nullptr;
                        //user_info_pointersJets[jKaon]->SetHeavyFlavour(true);
                        //inclusive_jets[iJet].user_info<JetUserInfo>().SetHeavyFlavour(true);

                        // calculate distance between HF decay candidate and jet axis in eta-phi plane
                        DeltaR = sqrt(pow(constituents[iConst].eta()-inclusive_jets[iJet].eta(),2) + pow(constituents[iConst].phi()-inclusive_jets[iJet].phi(),2));
                        
                        // fill histograms
                        hPartJet_Dist->Fill(DeltaR);
                        hEnergyDiff->Fill(abs(inclusive_jets[iJet].e()-constituents[iConst].e()));

                        // if there are other HF decay candidates, they are ignored (their distance from the jet axis is not calculated)

                        // get out of constituents loop
                        iConst = constituents.size();
                    }
                }
            }// end of constituents loop
            constituents.clear();
            
        }// end of jets loop
        
        // accessing heavy flavour jets
        for(int iJet = 0; iJet < inclusive_jets.size(); iJet++){
            if(inclusive_jets[iJet].has_user_info()){
                bool isHF_jet = inclusive_jets[iJet].user_info<JetUserInfo>().GetIsHeavyFlavour();
                if(isHF_jet){
                    //cout << "Jet number " << iJet << " tagged with heavy flavour particle.\n";
                }
            }
            
        }

        // clearing vectors
        input_particles.clear();
        inclusive_jets.clear();

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
    
    //cLeadPt->cd();
    //hLeadPt->Draw();
    // closing input file
    inputFile->Close();

    // creating file for storing results
    TFile f("AnalysisResults_"+sNumEvents+"_events_option_b.root","recreate");
    hLeadPt->Write();
    hJetsPt->Write();
    hNumJets->Write();
    hNumMCPart->Write();
    hPairDist->Write();
    hPartJet_Dist->Write();
    hKaonProdDist->Write();
    hCandInvMass->Write();
    hEnergyDiff->Write();
    cout << "Done.\n";
    f.Close();
    delete hLeadPt;
    delete hJetsPt;
    delete hNumJets;
    delete hNumMCPart;
    delete hPairDist;
    delete hPartJet_Dist;
    delete hKaonProdDist;
    delete hCandInvMass;
    delete hEnergyDiff;
    // ending clock
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed: " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000 << " seconds.\n";

    return 0;
}