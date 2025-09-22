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
    latex->DrawLatex(0.55, 0.5, Form("Unfolded integral (MC detector) = %.3f", hUnfoldedProj[hUnfoldedProj.size()-1]->Integral()));
    latex->DrawLatex(0.55, 0.45, Form("Input integral (MC particle) = %.3f", hInputParticleProj->Integral()));
    latex->DrawLatex(0.55, 0.4, Form("Unfolded / Input = %.3f", hUnfoldedProj[hUnfoldedProj.size()-1]->Integral() / hInputParticleProj->Integral()));

    // 3 ----- Plot same last distributions, but normalized
    TCanvas* cClosureTestNorm = new TCanvas("cClosureTestNorm","Second closure test normalized", 1800, 1000);
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

void StoreSecondClosureTest(TH3D* hBackgroundSubtracted, EfficiencyData& efficiencyDatacontainer, UnfoldData& unfoldDataContainer) {
    //
    TFile* fOutput = new TFile("SecondClosureTestResults.root","recreate");
    if (!fOutput || fOutput->IsZombie()) {
        std::cerr << "Error: Unable to create output ROOT file." << std::endl;
    }
    fOutput->cd();
    TDirectory* dirOuputSideband = fOutput->mkdir("1_BackgroundSubtraction");
    dirOuputSideband->cd();
    hBackgroundSubtracted->Write("hBackgroundSubtracted");

    TDirectory* dirOuputEfficiency = fOutput->mkdir("2_Efficiency");
    dirOuputEfficiency->cd();
    efficiencyDatacontainer.hSelectionEfficiency.first->Write("hSelectionEfficiencyPrompt");
    efficiencyDatacontainer.hSelectionEfficiency.second->Write("hSelectionEfficiencyNonPrompt");
    efficiencyDatacontainer.hEfficiencyCorrected.second->Write("hEfficiencyCorrected");
    efficiencyDatacontainer.hKEffResponseParticle_Over_TotalParticle.first->Write("hKinEff_particleLevel_prompt");
    efficiencyDatacontainer.hKEffResponseParticle_Over_TotalParticle.second->Write("hKinEff_particleLevel_prompt");
    efficiencyDatacontainer.hKEffResponseDetector_Over_TotalDetector.first->Write("hKinEff_detectorLevel_prompt");
    efficiencyDatacontainer.hKEffResponseDetector_Over_TotalDetector.second->Write("hKinEff_detectorLevel_nonprompt");
    
    TDirectory* dirOuputUnfolding = fOutput->mkdir("3_Unfolding/Unfolded/");
    dirOuputUnfolding->cd();
    for (size_t i = 0; i < unfoldDataContainer.hUnfoldedKinCorrected.size(); ++i) {
        unfoldDataContainer.hUnfoldedKinCorrected[i]->Write(Form("hUnfoldedKinCorrected_%zu", i));
    }
    unfoldDataContainer.hKineEffParticle[2]->Write("UnfoldKineEffParticle");
    unfoldDataContainer.hKineEffDetector[2]->Write("UnfoldKineEffDetector");
    fOutput->Close();
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
    EfficiencyData efficiencyDatacontainer = EfficiencyClosure(fClosureInputNonMatched, fClosureInputMatched, hBackgroundSubtracted, binningStruct, bdtPtCuts);
    std::vector<TH1D*> hSelEff_run3style = {efficiencyDatacontainer.hSelectionEfficiency.first,efficiencyDatacontainer.hSelectionEfficiency.second};
    TH2D* hEfficiencyCorrected = efficiencyDatacontainer.hEfficiencyCorrected.second;

    // ----3: Perform unfolding
    UnfoldData unfoldDataContainer = UnfoldingClosure(fClosureInputMatched, hSelEff_run3style, hEfficiencyCorrected, binningStruct, bdtPtCuts);

    // ----4: Compare input (MC particle level) with output (background subtracted, efficiency corrected, unfolded) distributions
    TH2D* inputMCP = CompareClosureTest(fClosureInputNonMatched, unfoldDataContainer.hUnfoldedKinCorrected, binningStruct);

    StoreSecondClosureTest(hBackgroundSubtracted, efficiencyDatacontainer, unfoldDataContainer);

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
