/*
 *
 *
 * Macro for performing unfolding on prompt D0s measured 
 * data and applying correction to FeedDownSubtraction.C
 * resulting distributions.
 * 
 * 
 * 
 * 
**/

#include "commonFunctions.h"
#include "sidebandClosure.h"
#include "efficiencyClosure.h"
#include "unfoldingClosure.h"

using namespace std;

// Already defined in sidebandClosure header file: calculate delta phi such that 0 < delta phi < 2*pi
// double DeltaPhi(double phi1, double phi2) {
//     // Compute the absolute difference between phi1 and phi2
//     double dphi = std::abs(phi1 - phi2); 
//     if (dphi > M_PI) {
//         // subtract 2pi if the difference if bigger than pi
//         dphi = dphi - 2*M_PI;
//     }

//     return dphi;
// }

// Already defined in sidebandClosure header file: get the optimal BDT score cut for the corresponding pT,D of the entry
// double GetBkgProbabilityCut(double pT, const std::vector<std::pair<double, double>>& bdtPtCuts) {
//     for (size_t i = 0; i < bdtPtCuts.size() - 1; ++i) {
//         if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) {
//             return bdtPtCuts[i].second;
//         }
//     }
//     return 1.0; // Default: accept all if out of range
// }


struct ClosureTestData2 {

    // Input objects: pT,jet vs DeltaR vs pT,D0
    TH3D* hInputParticle = nullptr;
    TH3D* hInputDetector = nullptr;

    // Correction objects
    RooUnfoldResponse* response;                                        // response matrix
    std::vector<TH2D*> hKineEffParticle = {nullptr, nullptr, nullptr};  // particle level kinematic efficiency -> 0 = numerator, 1 = denominator, 2 = efficiency
    std::vector<TH2D*> hKineEffDetector = {nullptr, nullptr, nullptr};  // detector level kinematic efficiency -> 0 = numerator, 1 = denominator, 2 = efficiency
    
    // Unfolding objects
    std::vector<RooUnfoldBayes*> unfold;                                // unfolding objects, there are iterationNumber unfolding objects
    std::vector<TH2D*> hUnfolded;                                       // unfolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects
};


ClosureTestData2 createAnalysisObjects(TFile* fClosureInput, double& jetptMin, double& jetptMax, double& hfptMin, double& hfptMax, 
                                             const std::vector<double>& ptjetBinEdges_particle, const std::vector<double>& deltaRBinEdges_particle, const std::vector<double>& ptDBinEdges_particle,
                                             const std::vector<double>& ptjetBinEdges_detector, const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& ptDBinEdges_detector, 
                                             const std::vector<std::pair<double, double>>& bdtPtCuts) {
    // Create empty struct to store data
    ClosureTestData2 dataContainer;
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = deltaRBinEdges_particle[deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = deltaRBinEdges_detector[deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double MCPHfPtMincut = hfptMin; // on particle level D0
    const double MCDHfPtMincut = hfptMin; // on detector level D0
    const double MCPHfPtMaxcut = hfptMax; // on particle level D0
    const double MCDHfPtMaxcut = hfptMax; // on detector level D0
    
    // Create 2D (detector level, prompt D0's, matched to particle level) input distribution pT,jet vs DeltaR
    dataContainer.hInputDetector = new TH3D("hInputDetector", "Detector level prompt D^{0} jets distribution;p_{T,jet}^{reco} (GeV);#DeltaR^{reco};p_{T,D^{0}}^{reco};", 
                        ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                        deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data(),
                        ptDBinEdges_detector.size() - 1, ptDBinEdges_detector.data());
    dataContainer.hInputDetector->Sumw2();
    // Create 2D (matched particle level, prompt D0's, matched to the previous detector level distribution) input distribution pT,jet vs DeltaR
    dataContainer.hInputParticle = new TH3D("hInputParticle", "Particle level prompt D^{0} jets distribution;p_{T,jet}^{gen} (GeV);#DeltaR^{gen};p_{T,D^{0}}^{gen}", 
                        ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(),
                        deltaRBinEdges_particle.size() - 1, deltaRBinEdges_particle.data(),
                        ptDBinEdges_particle.size() - 1, ptDBinEdges_particle.data());
    dataContainer.hInputParticle->Sumw2();

    // Create kinematic efficiency histograms
    dataContainer.hKineEffParticle[0] = new TH2D("hKineEffParticleNumerator", "Particle level kinematic efficiency numerator (non-prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                                                  deltaRBinEdges_particle.size() - 1, deltaRBinEdges_particle.data());
    dataContainer.hKineEffParticle[1] = new TH2D("hKineEffParticleDenominator", "Particle level kinematic efficiency denominator (non-prompt D^{0}'s);p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_particle.size() - 1, ptjetBinEdges_particle.data(), 
                                                                                  deltaRBinEdges_particle.size() - 1, deltaRBinEdges_particle.data());
    dataContainer.hKineEffDetector[0] = new TH2D("hKineEffDetectorNumerator", "Detector level kinematic efficiency numerator (non-prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                                                  deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data());
    dataContainer.hKineEffDetector[1] = new TH2D("hKineEffDetectorDenominator", "Detector level kinematic efficiency denominator (non-prompt D^{0}'s);p_{T,jet}^{reco ch} (GeV/#it{c});#DeltaR", 
                                                                                  ptjetBinEdges_detector.size() - 1, ptjetBinEdges_detector.data(), 
                                                                                  deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data());

    // Create response matrix for unfolding
    dataContainer.response = new RooUnfoldResponse(dataContainer.hKineEffDetector[0], dataContainer.hKineEffParticle[0]);


    //_______________________________________________Correction data_______________________________________________
    //
    // Accessing TTree for correction data
    //
    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;
    // Accessing TTree
    TTree* tree = (TTree*)fClosureInput->Get("CorrectionTree");
    // Check for correct access
    if (!tree) {
        cout << "Error opening correction data tree.\n";
    }
    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    MCPhfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&MCPhfY);
    tree->SetBranchAddress("fMcHfPrompt",&MCPhfprompt);
    // detector level branches
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);
    tree->SetBranchAddress("fHfMlScore0",&MCDhfMlScore0);
    tree->SetBranchAddress("fHfMlScore1",&MCDhfMlScore1);
    tree->SetBranchAddress("fHfMlScore2",&MCDhfMlScore2);

    //
    // Filling objects
    //
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // Apply prompt selection (i.e., only c → D0)
        if (!MCDhfprompt) {
            continue;
        }

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        // Get the threshold for this pT range
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;

        // Fill response matrix and kinematic efficiency histograms
        if (genLevelRange && recoLevelRange && passBDTcut) {
            
            // Fill 4D RooUnfoldResponse object (jet pT shape is influenced by D0 pT efficiency)
            dataContainer.response->Fill(MCDjetPt, MCDDeltaR, MCPjetPt, MCPDeltaR);
            
            // Fill kinematic efficiency numerator histograms
            dataContainer.hKineEffParticle[0]->Fill(MCPjetPt, MCPDeltaR);
            dataContainer.hKineEffDetector[0]->Fill(MCDjetPt, MCDDeltaR);
        }
        if (genLevelRange) {
            // Fill kinematic efficiency denominator histogram for full particle level range
            dataContainer.hKineEffParticle[1]->Fill(MCPjetPt, MCPDeltaR);
            
        }
        if (recoLevelRange && passBDTcut) {
            // Fill kinematic efficiency denominator histogram for full detector level range
            dataContainer.hKineEffDetector[1]->Fill(MCDjetPt, MCDDeltaR);
        }
    }

    // Calculate kinematic efficiency
    dataContainer.hKineEffParticle[0]->Sumw2(); // numerator
    dataContainer.hKineEffParticle[1]->Sumw2(); // denominator
    dataContainer.hKineEffParticle[2] = (TH2D*)dataContainer.hKineEffParticle[0]->Clone("hKineEffParticleEfficiency");
    dataContainer.hKineEffParticle[2]->Divide(dataContainer.hKineEffParticle[1]); // A = A / B:  A = A->Divide(B)
    dataContainer.hKineEffDetector[0]->Sumw2(); // numerator
    dataContainer.hKineEffDetector[1]->Sumw2(); // denominator
    dataContainer.hKineEffDetector[2] = (TH2D*)dataContainer.hKineEffDetector[0]->Clone("hKineEffDetectorEfficiency");
    dataContainer.hKineEffDetector[2]->Divide(dataContainer.hKineEffDetector[1]); // A = A / B:  A = A->Divide(B)
    //_______________________________________________Input data_______________________________________________
    //
    // Accessing TTree for correction data
    //
    tree = (TTree*)fClosureInput->Get("InputTree");
    // Check for correct access
    if (!tree) {
        cout << "Error opening correction data tree.\n";
    }
    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    MCPhfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&MCPhfY);
    tree->SetBranchAddress("fMcHfPrompt",&MCPhfprompt);
    // detector level branches
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);
    tree->SetBranchAddress("fHfMlScore0",&MCDhfMlScore0);
    tree->SetBranchAddress("fHfMlScore1",&MCDhfMlScore1);
    tree->SetBranchAddress("fHfMlScore2",&MCDhfMlScore2);
    float MCDhfMatchedFrom, MCDhfSelectedAs;
    tree->SetBranchAddress("fHfMatchedFrom",&MCDhfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&MCDhfSelectedAs);

    //
    // Filling objects
    //
    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // Apply prompt (real+reflection) selection (i.e., only c → D0)
        if (!MCDhfprompt) {
            continue;
        }

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= deltaRBinEdges_detector[0]) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        // Get the threshold for this pT range
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;

        if (genLevelRange) {
            // Fill input distribution for particle level
            dataContainer.hInputParticle->Fill(MCPjetPt, MCPDeltaR, MCPhfPt);
            
        }
        if (recoLevelRange && passBDTcut) {
            // Fill input distribution for detector level
            dataContainer.hInputDetector->Fill(MCDjetPt, MCDDeltaR, MCDhfPt);
        }
    }

    return dataContainer;

}

void unfoldInputDetector(ClosureTestData2& dataContainer, int& iterationNumber) {
    // Correct input with detector level kinematic efficiency
    dataContainer.hInputDetector->Multiply(dataContainer.hKineEffDetector[2]); // A = A * B:  A = A->Multiply(B)

    // Unfold with multiple iterations
    // Unfold in multiple iterations
    for (int iIter = 1; iIter <= iterationNumber; iIter++) {
        dataContainer.unfold.push_back(new RooUnfoldBayes(dataContainer.response, dataContainer.hInputDetector, iIter));
        //dataContainer.unfold[iIter - 1] = new RooUnfoldBayes(dataContainer.response, dataContainer.hInputDetector, iIter);
        dataContainer.hUnfolded.push_back((TH2D*) dataContainer.unfold[iIter - 1]->Hreco(RooUnfold::kCovariance)); // kCovariance = 2 = Errors from the square root of of the covariance matrix given by the unfolding
        dataContainer.hUnfolded[iIter - 1]->SetTitle(Form("Unfolded 2D histogram with %d iterations;p_{T,jet}^{gen} (GeV/#it{c});#DeltaR^{gen}", iIter));
        dataContainer.hUnfolded[iIter - 1]->SetName(Form("hUnfolded_iter%d", iIter));
        dataContainer.hUnfolded[iIter - 1]->Sumw2();

        // Correct distribution with particle level kinematic efficiency (add entries)
        dataContainer.hUnfolded[iIter - 1]->Divide(dataContainer.hKineEffParticle[2]); // A = A / B:  A = A->Divide(B)
        
    }
}
void plotHistograms(const ClosureTestData2& dataContainer, const double& jetptMin, const double& jetptMax) {
    cout << "Plotting histograms...\n";

    gStyle->SetPalette(kRainbow);

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    TCanvas* cFirstClosureTest = new TCanvas("cFirstClosureTest", "First Closure Test: Unfolding");
    cFirstClosureTest->cd();
    TLegend* lUnfoldedIter = new TLegend(0.6,0.57,0.7,0.77);
    std::vector<TH1D*> hUnfoldedProj(dataContainer.hUnfolded.size());
    int secondBin = 2;
    int lastButOneBin = dataContainer.hUnfolded[0]->GetXaxis()->GetNbins() - 1; // penultimate bin
    for (size_t iIter = 0; iIter < dataContainer.hUnfolded.size(); iIter++) {
        // Project Y axis (DeltaR) using X axis (pTjet) range [secondBin, oneBeforeLastBin] (excluding the padding bins)
        hUnfoldedProj[iIter] = dataContainer.hUnfolded[iIter]->ProjectionY(Form("hProjIter_%zu", iIter),secondBin, lastButOneBin); // bins specified in X (pT,jet) dimension
        hUnfoldedProj[iIter]->SetLineColor(kBlack + iIter);
        lUnfoldedIter->AddEntry(hUnfoldedProj[iIter],Form("Iteration %zu", iIter+1), "le");
        if (iIter == 0) {
            hUnfoldedProj[iIter]->SetTitle("Closure test 1: unfolding");
            hUnfoldedProj[iIter]->Draw();
        } else {
            hUnfoldedProj[iIter]->Draw("same");
        }
        
    }
    TH1D* hInputParticleProj = dataContainer.hInputParticle->ProjectionY("hInputParticleProj", secondBin, lastButOneBin);
    hInputParticleProj->SetLineColor(kRed);
    hInputParticleProj->SetLineStyle(2);
    hInputParticleProj->Draw("same");
    lUnfoldedIter->AddEntry(hInputParticleProj, "Input particle level", "le");
    lUnfoldedIter->Draw();
    double reportingJetPtMin = dataContainer.hUnfolded[0]->GetXaxis()->GetBinLowEdge(secondBin);
    double reportingJetPtMax = dataContainer.hUnfolded[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));
    
    //
    // Storing images
    //
    TString imagePath = "../Images/5-ClosureTest/";
    cFirstClosureTest->Update();
    cFirstClosureTest->SaveAs(imagePath + "ClosureTest1_unfolding.png");    
    
    //
    // Storing in a single pdf file
    //
    //cKinEff->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cFirstClosureTest->Print(imagePath + Form("closureTest1_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    //cUnfoldedIter->Print(imagePath + Form("unfolding_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

}

void saveData(const ClosureTestData2& dataContainer, const double& jetptMin, const double& jetptMax){
    // Open output file
    TFile* outFile = new TFile(Form("closure_test_1_results_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");

    // store each histogram in file
    //dataContainer.hSBUnfolded->Write();
    
    outFile->Close();
    delete outFile;
    
    std::cout << "Data stored in file" << Form("closure_test_1_results_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)) << endl;
}

std::vector<double> LoadBinning(TFile* fInput, const char* pathInFile) {
    auto* vec = (TVectorD*)fInput->Get(pathInFile);
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", pathInFile));
    }
    return std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
}

TH2D* CompareClosureTest(TFile* fClosureInputNonMatched, std::vector<TH2D*>& hUnfKinCorrected, const BinningStruct& binningStruct) {
    // 1 ----- Build particle level distribution the same way as the detector level was built
    TH2D* hInputParticle = new TH2D("hInputParticle", "Particle level prompt D^{0} jets distribution; pT,jet (GeV); #DeltaR", 
                        binningStruct.ptjetBinEdges_particle.size() - 1, binningStruct.ptjetBinEdges_particle.data(), 
                        binningStruct.deltaRBinEdges_particle.size() - 1, binningStruct.deltaRBinEdges_particle.data());
    hInputParticle->Sumw2();
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = binningStruct.deltaRBinEdges_particle[binningStruct.deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = binningStruct.deltaRBinEdges_detector[binningStruct.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double jetptMin = binningStruct.ptjetBinEdges_particle[0]; // on both levels jet
    const double jetptMax = binningStruct.ptjetBinEdges_particle[binningStruct.ptjetBinEdges_particle.size() - 1]; // on both levels jet
    const double MCPHfPtMincut = binningStruct.ptDBinEdges_particle[0]; // on particle level D0
    const double MCDHfPtMincut = binningStruct.ptDBinEdges_detector[0]; // on detector level D0
    const double MCPHfPtMaxcut = binningStruct.ptDBinEdges_particle[binningStruct.ptDBinEdges_particle.size() - 1]; // on particle level D0
    const double MCDHfPtMaxcut = binningStruct.ptDBinEdges_detector[binningStruct.ptDBinEdges_detector.size() - 1]; // on detector level D0
    TTree* tree = (TTree*)fClosureInputNonMatched->Get("InputTree/O2mcpjetdisttable");
    // Check for correct access
    if (!tree) {
        cout << "Error opening correction data tree.\n";
    }
    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    float MCPjetNconst;
    bool MCPhfprompt, MCPhfmatch;
    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetNconst);
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    MCPhfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&MCPhfY);
    tree->SetBranchAddress("fMcHfPrompt",&MCPhfprompt);
    tree->SetBranchAddress("fMcHfMatch",&MCPhfmatch);
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // Apply prompt (including reflections) selection (i.e., only c → D0)
        if (!MCPhfprompt) {
            continue;
        }

        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));

        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= binningStruct.deltaRBinEdges_particle[0]) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        
        if (genLevelRange) {
            // Fill input distribution for particle level
            hInputParticle->Fill(MCPjetPt, MCPDeltaR);
            
        }
    }
    // 2 ----- Plot input particle level distribution against unfolded and kinematically corrected distributions
    TCanvas* cClosureTest = new TCanvas("cClosureTest","Second closure test", 1800, 1000);
    cClosureTest->cd();
    TLegend* lUnfoldedIter = new TLegend(0.6,0.57,0.7,0.77);
    std::vector<TH1D*> hUnfoldedProj(hUnfKinCorrected.size());
    int secondBin = 2;
    int lastButOneBin = hUnfKinCorrected[0]->GetXaxis()->GetNbins() - 1; // penultimate bin
    for (size_t iIter = 0; iIter < hUnfKinCorrected.size(); iIter++) {
        // Project Y axis (DeltaR) using X axis (pTjet) range [secondBin, oneBeforeLastBin] (excluding the padding bins)
        hUnfoldedProj[iIter] = hUnfKinCorrected[iIter]->ProjectionY(Form("hProjIter_%zu", iIter),secondBin, lastButOneBin); // bins specified in X (pT,jet) dimension
        hUnfoldedProj[iIter]->SetLineColor(kBlack + iIter);
        lUnfoldedIter->AddEntry(hUnfoldedProj[iIter],Form("Iteration %zu", iIter+1), "le");
        if (iIter == 0) {
            hUnfoldedProj[iIter]->SetTitle("Closure test 2: background subtraction+efficiency+unfolding");
            hUnfoldedProj[iIter]->Draw();
        } else {
            hUnfoldedProj[iIter]->Draw("same");
        }
        
    }
    TH1D* hInputParticleProj = hInputParticle->ProjectionY("hInputParticleProj", secondBin, lastButOneBin);
    hInputParticleProj->SetLineColor(kRed);
    hInputParticleProj->SetLineStyle(2);
    hInputParticleProj->Draw("same");
    lUnfoldedIter->AddEntry(hInputParticleProj, "Input particle level", "le");
    lUnfoldedIter->Draw();
    double reportingJetPtMin = hUnfKinCorrected[0]->GetXaxis()->GetBinLowEdge(secondBin);
    double reportingJetPtMax = hUnfKinCorrected[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));    

    // 3 ----- Plot same last distributions, but normalized
    TCanvas* cClosureTestNorm = new TCanvas("cClosureTestNorm","Second closure test", 1800, 1000);
    cClosureTestNorm->cd();
    TLegend* lUnfoldedIterNorm = new TLegend(0.6,0.57,0.7,0.77);
    std::vector<TH1D*> hUnfoldedProjNorm(hUnfKinCorrected.size());
    secondBin = 2;
    lastButOneBin = hUnfKinCorrected[0]->GetXaxis()->GetNbins() - 1; // penultimate bin
    for (size_t iIter = 0; iIter < hUnfKinCorrected.size(); iIter++) {
        // Project Y axis (DeltaR) using X axis (pTjet) range [secondBin, oneBeforeLastBin] (excluding the padding bins)
        hUnfoldedProjNorm[iIter] = hUnfKinCorrected[iIter]->ProjectionY(Form("hProjIterNorm_%zu", iIter),secondBin, lastButOneBin); // bins specified in X (pT,jet) dimension
        hUnfoldedProjNorm[iIter]->Scale(1. / hUnfoldedProjNorm[iIter]->GetEntries());
        hUnfoldedProjNorm[iIter]->SetLineColor(kBlack + iIter);
        lUnfoldedIterNorm->AddEntry(hUnfoldedProjNorm[iIter],Form("Iteration %zu", iIter+1), "le");
        if (iIter == 0) {
            hUnfoldedProjNorm[iIter]->SetTitle("Closure test 2: background subtraction+efficiency+unfolding (normalized)");
            hUnfoldedProjNorm[iIter]->Draw();
        } else {
            hUnfoldedProjNorm[iIter]->Draw("same");
        }
        
    }
    TH1D* hInputParticleProjNorm = (TH1D*)hInputParticleProj->Clone("hInputParticleProjNorm");
    hInputParticleProjNorm->Scale(1. / hInputParticleProjNorm->GetEntries());
    hInputParticleProjNorm->SetLineColor(kRed);
    hInputParticleProjNorm->SetLineStyle(2);
    hInputParticleProjNorm->Draw("same");
    lUnfoldedIterNorm->AddEntry(hInputParticleProjNorm, "Input particle level", "le");
    lUnfoldedIterNorm->Draw();
    reportingJetPtMin = hUnfKinCorrected[0]->GetXaxis()->GetBinLowEdge(secondBin);
    reportingJetPtMax = hUnfKinCorrected[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
    TLatex* latexNorm = new TLatex();
    latexNorm->SetNDC(); // Set the coordinates to be normalized device coordinates
    latexNorm->SetTextSize(0.03);
    latexNorm->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));

    //
    // Storing images
    //
    TString imagePath = "../Images/5-ClosureTest/Second/";
    cClosureTest->Print(imagePath + Form("ClosureTest2_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cClosureTestNorm->Print(imagePath + Form("ClosureTest2_normalized_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));

    return hInputParticle;
}

// 2 - Sideband subtraction + efficiency correction + unfolding closure test
void SecondClosureTest(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // Number of unfolding procedure iterations
    int iterationNumber = 8;

    TFile* fAxes = new TFile(Form("../1-SignalTreatment/SideBand/full_merged_ranges_back_sub.root"),"read");
    if (!fAxes || fAxes->IsZombie()) {
        std::cerr << "Error: Unable to open simulated data ROOT file." << std::endl;
    }
    // Load pTjet bin edges
    std::vector<double> ptjetBinEdges_detector = LoadBinning(fAxes, "axes/ptjetBinEdges_detector");
    double jetptMin = ptjetBinEdges_detector[0]; // GeV
    double jetptMax = ptjetBinEdges_detector[ptjetBinEdges_detector.size() - 1]; // GeV
    // Load ΔR bin edges
    std::vector<double> deltaRBinEdges_detector = LoadBinning(fAxes, "axes/deltaRBinEdges_detector");
    double minDeltaR = deltaRBinEdges_detector[0];
    double maxDeltaR = deltaRBinEdges_detector[deltaRBinEdges_detector.size() - 1];
    // Load pTD bin edges
    std::vector<double> ptDBinEdges_detector = LoadBinning(fAxes, "axes/ptDBinEdges_detector");
    double hfptMin = ptDBinEdges_detector[0]; //ptDBinEdges[0] - should start from 0 or from the lowest pT,D value?
    double hfptMax = ptDBinEdges_detector[ptDBinEdges_detector.size() - 1];
    fAxes->Close();

    TFile* fEfficiency = new TFile(Form("../2-Efficiency/selection_efficiency_run3style_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fEfficiency || fEfficiency->IsZombie()) {
        std::cerr << "Error: Unable to open estimated selection efficiency ROOT file." << std::endl;
    }
    std::vector<double> ptjetBinEdges_particle = LoadBinning(fEfficiency, "axes/ptjetBinEdges_particle");
    std::vector<double> deltaRBinEdges_particle = LoadBinning(fEfficiency, "axes/deltaRBinEdges_particle");
    std::vector<double> ptDBinEdges_particle = LoadBinning(fEfficiency, "axes/ptDBinEdges_particle");

    // Create struct to hold all binning distributions
    BinningStruct binningStruct;
    binningStruct.ptjetBinEdges_particle = ptjetBinEdges_particle;
    binningStruct.deltaRBinEdges_particle = deltaRBinEdges_particle;
    binningStruct.ptDBinEdges_particle = ptDBinEdges_particle;
    binningStruct.ptjetBinEdges_detector = ptjetBinEdges_detector;
    binningStruct.deltaRBinEdges_detector = deltaRBinEdges_detector;
    binningStruct.ptDBinEdges_detector = ptDBinEdges_detector;

    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.03 (from first of pair)
    //std::vector<std::pair<double, double>> bdtPtCuts = {
    //    {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}
    //};
    // Dataset: JE_HF_LHC24g5_All_D0
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.4}, {8, 0.6}, {10, 0.8}, {12, 0.8}, {16, 1.0}, {50, 1.0}
    };

    // Opening files
    TFile* fSimulatedMCNonMatched = new TFile("../SimulatedData/Hyperloop_output/Train_runs/410602_Eff/AO2D_mergedDFs.root","read");
    TFile* fSimulatedMCMatched = new TFile("../SimulatedData/Hyperloop_output/Train_runs/410603_Match/AO2D_mergedDFs.root","read");
    TFile* fData = new TFile("../ExperimentalData/Hyperloop_output/HF_LHC23_pass4_Thin_small_2P3PDstar_DATA_newMLModel/AnalysisResults.root","read");
    TFile* fFeedDown = new TFile(Form("../3-Feed-Down/outputFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fClosureInputMatched = new TFile("mc_closure_input_matched_data.root","read");
    TFile* fClosureInputNonMatched = new TFile("mc_closure_input_non-matched_data.root","read");
    if (!fSimulatedMCMatched || fSimulatedMCMatched->IsZombie()) {
        std::cerr << "Error: Unable to open O2 MC matched ROOT file." << std::endl;
    }
    if (!fData || fData->IsZombie()) {
        std::cerr << "Error: Unable to open AnalysisResults.root data ROOT file." << std::endl;
    }
    if (!fFeedDown || fFeedDown->IsZombie()) {
        std::cerr << "Error: Unable to open AnalysisResults.root data ROOT file." << std::endl;
    }

    // ----1: Perform side-band subtraction
    // Example: selecting a model
    FitModelType modelToUse = FitModelType::SignalReflectionsOnly;
    TH3D* hBackgroundSubtracted = SidebandClosure(fClosureInputNonMatched, ptjetBinEdges_detector, deltaRBinEdges_detector, ptDBinEdges_detector, bdtPtCuts, modelToUse);
    //TH3D* hBackgroundSubtracted;

    // ----2: Perform efficiency correction
    std::pair<std::vector<TH1D*>, TH2D*> effOutputs = EfficiencyClosure(fClosureInputNonMatched, fClosureInputMatched, hBackgroundSubtracted, binningStruct, bdtPtCuts);
    std::vector<TH1D*> hSelEff_run3style = effOutputs.first;
    TH2D* hEfficiencyCorrected = effOutputs.second;

    // ----3: Perform unfolding
    std::vector<TH2D*> hUnfKinCorrected = UnfoldingClosure(fClosureInputMatched, hSelEff_run3style, hEfficiencyCorrected, binningStruct, bdtPtCuts);

    // ----4: Compare input (MC particle level) with output (background subtracted, efficiency corrected, unfolded) distributions
    TH2D* inputMCP = CompareClosureTest(fClosureInputNonMatched, hUnfKinCorrected, binningStruct);

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    SecondClosureTest();
    return 0;
}
