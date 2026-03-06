#include "fastjet/ClusterSequence.hh"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"

bool isLastEventParticle(int treeEntry, TTree* tInput, int& currentEventNumber) {
    int nEntries = tInput->GetEntries();
    int nextEventNumber;
    
    // Save the current event number (already in currentEventNumber from main)
    
    if (treeEntry == nEntries - 1) {
        return true; // This is the last entry in the tree
    }
    
    // Temporarily read the next entry to check its event number
    tInput->SetBranchAddress("eventNumber", &nextEventNumber);
    tInput->GetEntry(treeEntry + 1);
    
    // Restore the branch address to point to the main variable
    tInput->SetBranchAddress("eventNumber", &currentEventNumber);
    
    // IMPORTANT: Restore the current entry!
    tInput->GetEntry(treeEntry);
    
    return (currentEventNumber != nextEventNumber);
}

double DeltaPhi(double phi1, double phi2) {
    // Compute the absolute difference between phi1 and phi2
    double dphi = std::abs(phi1 - phi2); 
    if (dphi > M_PI) {
        // subtract 2pi if the difference if bigger than pi
        dphi = dphi - 2*M_PI;
    }

    return dphi;
}

double calculateDeltaR(fastjet::PseudoJet jet, fastjet::PseudoJet constituent) {
    double deltaEta = jet.eta() - constituent.eta();
    double deltaPhi = DeltaPhi(jet.phi(), constituent.phi());
    return sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
    
}

TH1D* buildDeltaRHistogram(TFile* fInput, const char* histName, const char* histTitle, const std::vector<double>& deltaRBinEdges, 
                           const double& ptjetMin, const double& ptjetMax, const double& hfPtMin, const double& hfPtMax, const double& etaCut, const double& yCut) {

    // Get the TTree from the file
    TTree* tInput = (TTree*)fInput->Get("tCharmJets");
    if (!tInput) {
        std::cout << "Error getting TTree!" << std::endl;
        return nullptr;
    }
    float axisDistance, jetPt, jetEta, jetPhi, jetMass;
    float hfPt, hfEta, hfPhi, hfMass, hfY;
    int jetNConst, iEventNumber;
    int hfPrompt; // 1 if the D0 is prompt, 0 if from b decay, -1 if not a D0 particle
    tInput->SetBranchAddress("jetHfDist",&axisDistance);
    tInput->SetBranchAddress("jetPt",&jetPt);
    tInput->SetBranchAddress("jetEta",&jetEta);
    tInput->SetBranchAddress("jetPhi",&jetPhi);
    tInput->SetBranchAddress("jetMass",&jetMass);
    tInput->SetBranchAddress("jetNConst",&jetNConst);
    tInput->SetBranchAddress("hfPt",&hfPt);
    tInput->SetBranchAddress("hfEta",&hfEta);
    tInput->SetBranchAddress("hfPhi",&hfPhi);
    tInput->SetBranchAddress("hfMass",&hfMass);
    tInput->SetBranchAddress("hfY",&hfY);
    tInput->SetBranchAddress("hfPrompt",&hfPrompt);
    tInput->SetBranchAddress("iEventNumber",&iEventNumber);

    TH1D* hDeltaR = new TH1D(histName, histTitle, deltaRBinEdges.size() - 1, deltaRBinEdges.data());
    hDeltaR->GetXaxis()->SetTitle("#Delta R");
    hDeltaR->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#Delta R}");

    int nEntries = tInput->GetEntries();
    for (int iEntry = 0; iEntry < nEntries; iEntry++) {
        tInput->GetEntry(iEntry);
        bool jetRange = (jetPt >= ptjetMin && jetPt < ptjetMax) ? true : false;
        bool hfPtRange = (hfPt >= hfPtMin && hfPt < hfPtMax) ? true : false;
        if (!jetRange || !hfPtRange || (abs(jetEta) > etaCut) || (abs(hfY) > yCut)) {
            continue; // skip entries that do not satisfy the cuts
        } else {
            hDeltaR->Fill(axisDistance);
        }
    }

    // Normalize histograms by number of jets in the pT range
    hDeltaR->Scale(1.0 / hDeltaR->Integral(), "width"); // Normalize by total number of jets and bin width


    return hDeltaR;
}