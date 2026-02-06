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
        cout << "Error opening input data tree.\n";
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

        // Apply prompt (excluding reflections) selection (i.e., only c → D0)
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
    bool useAllEffData = false;
    if (useAllEffData) {
        tree = (TTree*)fClosureInputNonMatched->Get("CorrectionTree/O2mcpjetdisttable");
        // Check for correct access
        if (!tree) {
            cout << "Error opening correction data tree.\n";
        }
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

            // Apply prompt (excluding reflections) selection (i.e., only c → D0)
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
        hUnfoldedProjNorm[iIter]->Scale(1. / hUnfoldedProjNorm[iIter]->Integral(), "width");
        std::cout << "Unfolded iteration " << iIter << " integral is " << hUnfoldedProjNorm[iIter]->Integral() << std::endl;
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
    hInputParticleProjNorm->Scale(1. / hInputParticleProjNorm->Integral(), "width");
    std::cout << "Input particle level distribution integral is " << hInputParticleProjNorm->Integral() << std::endl;
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

TCanvas* AnalysisEvolutionSteps(TH2D* inputMCP, TH3D* hBackgroundSubtracted, TH2D* hEfficiencyCorrected, TH2D* hKinCorrBeforeUnfolding, TH2D* hUnfolded, TH2D* hUnfoldedKinCorrected, const BinningStruct& binningStruct) {
    //
    TCanvas* cStepsEvolution = new TCanvas("cStepsEvolution", "Analysis steps evolution", 1800, 1000);
    cStepsEvolution->Divide(2,2);

    int secondBin = 2; // default: 2
    int lastButOneBin = inputMCP->GetXaxis()->GetNbins() - 1; // default: inputMCP->GetXaxis()->GetNbins() - 1 (penultimate bin)
    TH1D* hInputParticleProj = inputMCP->ProjectionY("hInputParticleProj_StepsEvolution0", secondBin, lastButOneBin);
    hInputParticleProj->SetLineColor(kRed+2);
    hInputParticleProj->SetLineStyle(2);
    hInputParticleProj->SetLineWidth(2);

    // 1: Background subtracted
    cStepsEvolution->cd(1);
    TH2D* hBackgroundSubtracted2d = (TH2D*)hBackgroundSubtracted->Project3D("yx"); // project (pT,jet; DeltaR; pT,D) = (x; y; z) -> (x; y)
    TH1D* hBackgroundSubtractedProj = (TH1D*)hBackgroundSubtracted2d->ProjectionY("hBackgroundSubtractedProj_StepsEvolution", secondBin, lastButOneBin);
    hBackgroundSubtractedProj->SetLineColor(kBlue+2);
    hBackgroundSubtractedProj->SetLineWidth(2);
    //hBackgroundSubtractedProj->SetTitle("Background subtracted distribution");
    TH1D* hInputParticleProj1 = (TH1D*)hInputParticleProj->Clone("hInputParticleProj_StepsEvolution1");
    hInputParticleProj1->SetTitle("Step 1: background subtraction");
    double max1 = hInputParticleProj1->GetMaximum();
    double max2 = hBackgroundSubtractedProj->GetMaximum();
    double ymax = TMath::Max(max1, max2) * 1.2; // Add 20% margin
    hInputParticleProj1->GetYaxis()->SetRangeUser(0., ymax);
    hInputParticleProj1->Draw();
    hBackgroundSubtractedProj->Draw("same");

    // 2: Efficiency corrected
    cStepsEvolution->cd(2);
    TH1D* hEfficiencyCorrectedProj = (TH1D*)hEfficiencyCorrected->ProjectionY("hEfficiencyCorrectedProj_StepsEvolution", secondBin, lastButOneBin);
    hEfficiencyCorrectedProj->SetLineColor(kGreen+2);
    hEfficiencyCorrectedProj->SetLineWidth(2);
    //hBackgroundSubtractedProj->SetTitle("Background subtracted distribution");
    TH1D* hInputParticleProj2 = (TH1D*)hInputParticleProj->Clone("hInputParticleProj_StepsEvolution2");
    hInputParticleProj2->SetTitle("Step 2: efficiency correction");
    hInputParticleProj2->Draw();
    hEfficiencyCorrectedProj->Draw("same");

    // 3: Unfolded + kinematic efficiency corrected
    cStepsEvolution->cd(3);
    TH1D* hKinCorrBeforeUnfoldingProj = (TH1D*)hKinCorrBeforeUnfolding->ProjectionY("hKinCorrBeforeUnfoldingProj_StepsEvolution", secondBin, lastButOneBin);
    hKinCorrBeforeUnfoldingProj->SetLineColor(kMagenta+2);
    hKinCorrBeforeUnfoldingProj->SetLineWidth(2);
    TH1D* hUnfoldedProj = (TH1D*)hUnfolded->ProjectionY("hUnfoldedProj_StepsEvolution", secondBin, lastButOneBin);
    hUnfoldedProj->SetLineColor(kOrange+2);
    hUnfoldedProj->SetLineWidth(2);
    TH1D* hUnfoldedKinCorrectedProj = (TH1D*)hUnfoldedKinCorrected->ProjectionY("hUnfoldedKinCorrectedProj_StepsEvolution", secondBin, lastButOneBin);
    hUnfoldedKinCorrectedProj->SetLineColor(kCyan+2);
    hUnfoldedKinCorrectedProj->SetLineWidth(2);
    //hBackgroundSubtractedProj->SetTitle("Background subtracted distribution");
    TH1D* hInputParticleProj3 = (TH1D*)hInputParticleProj->Clone("hInputParticleProj_StepsEvolution3");
    hInputParticleProj3->SetTitle("Step 3: unfolding and kinematic efficiency correction");
    hInputParticleProj3->Draw();
    hKinCorrBeforeUnfoldingProj->Draw("same");
    hUnfoldedProj->Draw("same");
    hUnfoldedKinCorrectedProj->Draw("same");

    // 4: altogether
    cStepsEvolution->cd(4);
    //hBackgroundSubtractedProj->SetTitle("Background subtracted distribution");
    hInputParticleProj->SetTitle("All steps combined");
    hInputParticleProj->Draw();
    hBackgroundSubtractedProj->Draw("same");
    hEfficiencyCorrectedProj->Draw("same");
    hKinCorrBeforeUnfoldingProj->Draw("same");
    hUnfoldedProj->Draw("same");
    hUnfoldedKinCorrectedProj->Draw("same");

    TLegend* legend = new TLegend(0.5,0.5,0.85,0.85);
    legend->AddEntry(hInputParticleProj, "Input particle level", "le");
    legend->AddEntry(hBackgroundSubtractedProj, "1 - After background subtraction", "le");
    legend->AddEntry(hEfficiencyCorrectedProj, "2 - After efficiency correction", "le");
    legend->AddEntry(hKinCorrBeforeUnfoldingProj, "3.1 - #varepsilon_{kin} correction before unfolding", "le");
    legend->AddEntry(hUnfoldedProj, "3.2 - After unfolding", "le");
    legend->AddEntry(hUnfoldedKinCorrectedProj, "3.3 - #varepsilon_{kin} correction after unfolding", "le");
    legend->Draw();

    double reportingJetPtMin = hUnfoldedKinCorrected->GetXaxis()->GetBinLowEdge(secondBin);
    double reportingJetPtMax = hUnfoldedKinCorrected->GetXaxis()->GetBinUpEdge(lastButOneBin);

    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));


    
    return cStepsEvolution;
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
    
    // Create parent directory first
    TDirectory* dirParent = fOutput->mkdir("3_Unfolding");
    if (dirParent) {
        dirParent->cd();
        TDirectory* dirOuputUnfolding = dirParent->mkdir("Unfolded");
        if (!dirOuputUnfolding) {
            std::cerr << "Error: Failed to create directory '3_Unfolding/Unfolded/'" << std::endl;
            fOutput->Close();
            return;
        }
        dirOuputUnfolding->cd();
        for (size_t i = 0; i < unfoldDataContainer.hUnfoldedKinCorrected.size(); ++i) {
            unfoldDataContainer.hUnfoldedKinCorrected[i]->Write(Form("hUnfoldedKinCorrected_%zu", i));
        }
        unfoldDataContainer.hKineEffParticle[2]->Write("UnfoldKineEffParticle");
        unfoldDataContainer.hKineEffDetector[2]->Write("UnfoldKineEffDetector");
    }
    
    fOutput->Close();
}

// void InvestigateNonmatchingProblem(UnfoldData& unfoldDataContainer, std::vector<TH2D*>& hUnfKinCorrected, TH2D* inputMCP) {
//     //
//     TCanvas* cProblemRatio = new TCanvas("cProblemRatio","Investigate non-matching of final and inicial distributions", 1800, 1000);
//     cProblemRatio->cd();
//     TLegend* lUnfoldedIter = new TLegend(0.6,0.57,0.7,0.77);
//     std::vector<TH1D*> hUnfoldedProj(hUnfKinCorrected.size());
//     int secondBin = 2;
//     int lastButOneBin = hUnfKinCorrected[0]->GetXaxis()->GetNbins() - 1; // penultimate bin
//     for (size_t iIter = 0; iIter < hUnfKinCorrected.size(); iIter++) {
//         // Project Y axis (DeltaR) using X axis (pTjet) range [secondBin, oneBeforeLastBin] (excluding the padding bins)
//         hUnfoldedProj[iIter] = hUnfKinCorrected[iIter]->ProjectionY(Form("hProjIterInvestigation_%zu", iIter),secondBin, lastButOneBin); // bins specified in X (pT,jet) dimension
//         hUnfoldedProj[iIter]->SetLineColor(kBlack + iIter);
//         lUnfoldedIter->AddEntry(hUnfoldedProj[iIter],Form("Iteration %zu", iIter+1), "le");
//         if (iIter == 0) {
//             hUnfoldedProj[iIter]->SetTitle("Closure test 2: background subtraction+efficiency+unfolding");
//             //hUnfoldedProj[iIter]->Draw();
//         } else {
//             //hUnfoldedProj[iIter]->Draw("same");
//         }
//     }
//     //hUnfoldedProj[hUnfoldedProj.size() - 1]->Draw();
//     TH1D* hInputParticleProj = inputMCP->ProjectionY("hInputParticleProjInvestigation", secondBin, lastButOneBin);
//     hInputParticleProj->SetLineColor(kRed);
//     hInputParticleProj->SetLineStyle(2);
//     //hInputParticleProj->Draw("same");
//     hInputParticleProj->Divide(hUnfoldedProj[hUnfoldedProj.size() - 1]);
//     hInputParticleProj->Draw();
//     lUnfoldedIter->AddEntry(hInputParticleProj, "Input particle level", "le");
//     //lUnfoldedIter->Draw();
//     double reportingJetPtMin = hUnfKinCorrected[0]->GetXaxis()->GetBinLowEdge(secondBin);
//     double reportingJetPtMax = hUnfKinCorrected[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
//     TLatex* latex = new TLatex();
//     latex->SetNDC(); // Set the coordinates to be normalized device coordinates
//     latex->SetTextSize(0.03);
//     latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));
//     // latex->DrawLatex(0.55, 0.5, Form("Unfolded integral (MC detector) = %.3f", hUnfoldedProj[hUnfoldedProj.size()-1]->Integral()));
//     // latex->DrawLatex(0.55, 0.45, Form("Input integral (MC particle) = %.3f", hInputParticleProj->Integral()));
//     // latex->DrawLatex(0.55, 0.4, Form("Unfolded / Input = %.3f", hUnfoldedProj[hUnfoldedProj.size()-1]->Integral() / hInputParticleProj->Integral()));
//     // Investigate response matrices
//     // TCanvas* cResponseMatrices = new TCanvas("cResponseMatrices","Response matrices investigation", 1800, 1000);
//     // cResponseMatrices->Divide(2,2);
//     // for (size_t iIter = 0; iIter < unfoldDataContainer.unfold.size(); ++iIter) {
//     //     cResponseMatrices->cd(iIter + 1);
//     //     TH2D* hResponseMatrix = (TH2D*)unfoldDataContainer.unfold[iIter]->Hresponse();
//     //     hResponseMatrix->SetTitle(Form("Response matrix - iteration %zu", iIter + 1));
//     //     hResponseMatrix->GetXaxis()->SetTitle("Detector level: p_{T,jet} (GeV/c)");
//     //     hResponseMatrix->GetYaxis()->SetTitle("Particle level: p_{T,jet} (GeV/c)");
//     //     hResponseMatrix->Draw("text");
//     // }
// }

void InvestigateNonmatchingProblem(UnfoldData& unfoldDataContainer, std::vector<TH2D*>& hUnfKinCorrected, TH2D* inputMCP) {
    TCanvas* cProblemRatio = new TCanvas("cProblemRatio","Investigate non-matching of final and inicial distributions", 1800, 1000);
    cProblemRatio->cd();
    TLegend* lUnfoldedIter = new TLegend(0.6,0.57,0.7,0.77);
    std::vector<TH1D*> hUnfoldedProj(hUnfKinCorrected.size());
    int secondBin = 2;
    int lastButOneBin = hUnfKinCorrected[0]->GetXaxis()->GetNbins() - 1; // penultimate bin
    
    for (size_t iIter = 0; iIter < hUnfKinCorrected.size(); iIter++) {
        hUnfoldedProj[iIter] = hUnfKinCorrected[iIter]->ProjectionY(Form("hProjIterInvestigation_%zu", iIter),secondBin, lastButOneBin);
        hUnfoldedProj[iIter]->SetLineColor(kBlack + iIter);
        lUnfoldedIter->AddEntry(hUnfoldedProj[iIter],Form("Iteration %zu", iIter+1), "le");
    }
    
    TH1D* hInputParticleProj = inputMCP->ProjectionY("hInputParticleProjInvestigation", secondBin, lastButOneBin);
    hInputParticleProj->SetLineColor(kRed);
    hInputParticleProj->SetLineStyle(2);
    
    // Calculate percentage deviation: (unfolded - particle) / particle * 100%
    TLegend* legDeviationHistos = new TLegend(0.7, 0.6, 0.9, 0.75);
    legDeviationHistos->SetFillStyle(0);
    
    std::vector<TH1D*> hPercentageDeviation(hUnfKinCorrected.size());
    for (size_t iHisto = 0; iHisto < hPercentageDeviation.size(); iHisto++) {
        hPercentageDeviation[iHisto] = (TH1D*)hUnfoldedProj[iHisto]->Clone(Form("hPercentageDeviation_%zu", iHisto));
        hPercentageDeviation[iHisto]->Add(hInputParticleProj, -1.0); // unfolded - particle
        hPercentageDeviation[iHisto]->Divide(hInputParticleProj);    // (unfolded - particle) / particle
        hPercentageDeviation[iHisto]->Scale(100.0);                  // Convert to percentage

        // Style the percentage deviation plot
        hPercentageDeviation[iHisto]->SetLineColor(kBlack + iHisto);
        hPercentageDeviation[iHisto]->SetLineWidth(2);
        hPercentageDeviation[iHisto]->SetMarkerColor(kBlack + iHisto);
        hPercentageDeviation[iHisto]->SetMarkerStyle(20);
        hPercentageDeviation[iHisto]->SetMarkerSize(0.8);
        if (iHisto == hPercentageDeviation.size()-1) {
            hPercentageDeviation[iHisto]->SetLineColor(28);
            hPercentageDeviation[iHisto]->SetMarkerColor(28);
        }
        

        // Set titles and labels for the deviation plot
        hPercentageDeviation[iHisto]->SetTitle(Form("Percentage Deviation from Particle Level, iteration %zu; #DeltaR; Deviation (%%)", iHisto));
        hPercentageDeviation[iHisto]->GetYaxis()->SetTitle("(Unfolded - Particle) / Particle #times 100%");

        // Draw the percentage deviation plot        
        hPercentageDeviation[iHisto]->Draw(iHisto == 0 ? "E1" : "E1 same");

        legDeviationHistos->AddEntry(hPercentageDeviation[iHisto], "Percentage Deviation", "lep");
    }
    
    // Create reference lines for ±10% deviation
    TLine* linePlus10 = new TLine(hPercentageDeviation[hPercentageDeviation.size()-1]->GetXaxis()->GetXmin(), 10.0, 
                                  hPercentageDeviation[hPercentageDeviation.size()-1]->GetXaxis()->GetXmax(), 10.0);
    TLine* lineMinus10 = new TLine(hPercentageDeviation[hPercentageDeviation.size()-1]->GetXaxis()->GetXmin(), -10.0, 
                                   hPercentageDeviation[hPercentageDeviation.size()-1]->GetXaxis()->GetXmax(), -10.0);
    TLine* lineZero = new TLine(hPercentageDeviation[hPercentageDeviation.size()-1]->GetXaxis()->GetXmin(), 0.0, 
                                hPercentageDeviation[hPercentageDeviation.size()-1]->GetXaxis()->GetXmax(), 0.0);
    
    // Style the reference lines
    linePlus10->SetLineColor(kRed);
    linePlus10->SetLineStyle(2);
    linePlus10->SetLineWidth(1);
    
    lineMinus10->SetLineColor(kRed);
    lineMinus10->SetLineStyle(2);
    lineMinus10->SetLineWidth(1);
    
    lineZero->SetLineColor(kBlack);
    lineZero->SetLineStyle(1);
    lineZero->SetLineWidth(1);
    
    // Draw reference lines
    lineZero->Draw("same");
    linePlus10->Draw("same");
    lineMinus10->Draw("same");
    
    // Add legend for deviation plot
    TLegend* legDeviation = new TLegend(0.7, 0.75, 0.9, 0.9);
    legDeviation->SetBorderSize(0);
    legDeviation->SetFillStyle(0);
    //legDeviation->AddEntry(hPercentageDeviation, "Percentage Deviation", "lep");
    legDeviation->AddEntry(linePlus10, "#pm10% Reference", "l");
    legDeviation->AddEntry(lineZero, "Perfect Agreement (0%)", "l");
    legDeviation->Draw();
    legDeviationHistos->Draw();
    
    // Calculate and display statistics
    double meanDeviation = hPercentageDeviation[hPercentageDeviation.size()-1]->GetMean();
    double rmsDeviation = hPercentageDeviation[hPercentageDeviation.size()-1]->GetRMS();
    double maxDeviation = hPercentageDeviation[hPercentageDeviation.size()-1]->GetMaximum();
    double minDeviation = hPercentageDeviation[hPercentageDeviation.size()-1]->GetMinimum();
    
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.02); // default = 0.3
    double latexVerticalDistance = 0.03; // default = 0.04
    double reportingJetPtMin = hUnfKinCorrected[0]->GetXaxis()->GetBinLowEdge(secondBin);
    double reportingJetPtMax = hUnfKinCorrected[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));
    
    // Add deviation statistics
    latex->DrawLatex(0.15, 0.80, Form("Mean deviation (last iteration): %.1f%%", meanDeviation));
    latex->DrawLatex(0.15, 0.8 - latexVerticalDistance, Form("RMS deviation (last iteration): %.1f%%", rmsDeviation));
    latex->DrawLatex(0.15, 0.8 - 2*latexVerticalDistance, Form("Max deviation (last iteration): %.1f%%", maxDeviation));
    latex->DrawLatex(0.15, 0.8 - 3*latexVerticalDistance, Form("Min deviation (last iteration): %.1f%%", minDeviation));
    
    // Count bins outside ±10% range
    int binsOutside10Percent = 0;
    int totalBins = hPercentageDeviation[hPercentageDeviation.size()-1]->GetNbinsX();
    for (int i = 1; i <= totalBins; i++) {
        double deviation = hPercentageDeviation[hPercentageDeviation.size()-1]->GetBinContent(i);
        if (fabs(deviation) > 10.0) {
            binsOutside10Percent++;
        }
    }
    double percentOutside10Percent = (double)binsOutside10Percent / totalBins * 100.0;
    
    latex->DrawLatex(0.15, 0.62, Form("Bins outside #pm10%% (last iteration): %d/%d (%.1f%%)", 
                                     binsOutside10Percent, totalBins, percentOutside10Percent));
    
    // Also keep your original ratio plot if you want both
    // You might want to create a second canvas for the ratio plot
    TCanvas* cOriginalRatio = new TCanvas("cOriginalRatio", "Original Ratio Plot", 1800, 1000);
    cOriginalRatio->cd();
    
    // Your original ratio plotting code here...
    hInputParticleProj->Divide(hUnfoldedProj[hUnfoldedProj.size() - 1]);
    hInputParticleProj->Draw();
    lUnfoldedIter->AddEntry(hInputParticleProj, "Input particle level", "le");
    //lUnfoldedIter->Draw();
    //double reportingJetPtMin = hUnfKinCorrected[0]->GetXaxis()->GetBinLowEdge(secondBin);
    //double reportingJetPtMax = hUnfKinCorrected[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
    TLatex* origlatex = new TLatex();
    origlatex->SetNDC(); // Set the coordinates to be normalized device coordinates
    origlatex->SetTextSize(0.03);
    origlatex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));
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

    // ----5: Compare correction after each step of the analysis chain
    TCanvas* cStepsEvolution = AnalysisEvolutionSteps(inputMCP, hBackgroundSubtracted, hEfficiencyCorrected, unfoldDataContainer.hCorrectedInput, unfoldDataContainer.hUnfolded[iterationNumber - 1], unfoldDataContainer.hUnfoldedKinCorrected[iterationNumber - 1], binningStruct);

    StoreSecondClosureTest(hBackgroundSubtracted, efficiencyDatacontainer, unfoldDataContainer);

    InvestigateNonmatchingProblem(unfoldDataContainer, unfoldDataContainer.hUnfoldedKinCorrected, inputMCP);

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
