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

#include "../commonUtilities.h"
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

TH2D* CompareClosureTest(TFile* fClosureInput, const SidebandClosureResult& sidebandDataContainer, const EfficiencyData& efficiencyDatacontainer, const UnfoldData& unfoldDataContainer, const BinningStruct& binning) {
    // Number of iterations of unfolding considered as the solution
    size_t selectedIteration = 4;

    // Creating input particle level histogram
    TH2D* hDeltaRInputMCP = new TH2D("hDeltaRInputMCP",";p_{T,jet}^{gen ch} (GeV/#it{c});#DeltaR",
                                     binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), 
                                     binning.deltaRBinEdges_particle.size() - 1, binning.deltaRBinEdges_particle.data());
    // Filling histogram
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = binning.deltaRBinEdges_particle[binning.deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double jetptMin = binning.ptjetBinEdges_particle[0]; // on both levels jet
    const double jetptMax = binning.ptjetBinEdges_particle[binning.ptjetBinEdges_particle.size() - 1]; // on both levels jet
    const double MCPHfPtMincut = binning.ptHFBinEdges_particle[0]; // on particle level D0
    const double MCDHfPtMincut = binning.ptHFBinEdges_detector[0]; // on detector level D0
    const double MCPHfPtMaxcut = binning.ptHFBinEdges_particle[binning.ptHFBinEdges_particle.size() - 1]; // on particle level D0
    const double MCDHfPtMaxcut = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1]; // on detector level D0

    // Accessing TTree
    TTree* tree = (TTree*)fClosureInput->Get("InputTree");
    // Check for correct access
    if (!tree) {
        cout << "Error opening O2 matching tree.\n";
    }
    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    int MCPjetNConst;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    int MCDhfMatchedFrom, MCDhfSelectedAs;
    int MCDjetNConst;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;

    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetNConst); // float
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
    tree->SetBranchAddress("fJetNConst",&MCDjetNConst); // int
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);
    tree->SetBranchAddress("fHfMlScore0",&MCDhfMlScore0);
    tree->SetBranchAddress("fHfMlScore1",&MCDhfMlScore1);
    tree->SetBranchAddress("fHfMlScore2",&MCDhfMlScore2);
    tree->SetBranchAddress("fHfMatchedFrom",&MCDhfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&MCDhfSelectedAs);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // Apply prompt-only selection (i.e., only c → D0, not b → B → D0), reflections and non-matched entries are allowed here
        bool isPrompt = MCPhfprompt; // unmatched entries will possess MCDhfprompt = -2
        if (!isPrompt) {
            continue;
        }
        // Reflections should not be included on particle or detector level
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        if (MCDhfmatch && !isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs)) {
            continue;
        }

        // Generator level selection cuts
        double MCPDeltaR = MCPaxisDistance;
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax); // removed upper bound, --> this is new!
        bool genHfPtRange = ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)); // removed entirely, --> this is new!
        bool genDeltaRRange = ((MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < MCPDeltaRcut));
        bool genAcceptance = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut);
        bool genLevelRange = genAcceptance && genJetPtRange && genDeltaRRange && genHfPtRange; // --> this is new!
        // Reconstruction level selection cuts
        double MCDdeltaR = MCDaxisDistance;
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb);
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        bool recoDeltaRRange = ((MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < MCDDeltaRcut));
        bool recoAcceptance = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut);
        bool recoLevelRange = recoAcceptance && recoJetPtRange && recoHfPtRange && recoDeltaRRange;

        // Should the detector level matched must always satisfy the kinematic cuts as applied to data?
        if (MCDhfmatch && !recoLevelRange) { // --> this is new!
            continue;
        }

        // Fill response matrix and kinematic efficiency histograms
        if (genLevelRange) {

            // Find which jet pT bin this entry belongs to
            int iJetPt = -1;
            for (size_t i = 0; i < binning.ptjetBinEdges_particle.size() - 1; ++i) {
                if (MCPjetPt >= binning.ptjetBinEdges_particle[i] && MCPjetPt <  binning.ptjetBinEdges_particle[i+1]) {
                    iJetPt = i;
                    break;
                }
            }
            if (iJetPt < 0) continue; // outside range

            // Find which pT,D bin this entry belongs to
            int iPtHF = -1;
            for (size_t i = 0; i < binning.ptHFBinEdges_particle.size() - 1; ++i) {
                if (MCPhfPt >= binning.ptHFBinEdges_particle[i] && MCPhfPt <  binning.ptHFBinEdges_particle[i+1]) {
                    iPtHF = i;
                    break;
                }
            }
            if (iPtHF < 0) continue; // outside pT,D range

            // Check if this pT,HF bin survived the sideband fit for this jet pT bin
            std::vector<bool> workingFits = sidebandDataContainer.workingFitsPerJetPt[iJetPt];
            if (!workingFits.empty() && !eraseHistogram(workingFits, iPtHF)) { // Also apply eraseHistogram logic to determine the surviving block
                // use Emma's run 2 cuts in case useEmmaYeatsBins is true
                // if (!binning.useEmmaYeatsBins || (MCDhfmatch && passEmmaCut(MCDjetPt, MCDhfPt)) || !MCDhfmatch) { // --> this is new!
                //     // Emma cut should only be checked for the matched detector level entry MCDs?
                //     // Verify Emma's cuts only in case these are matched entries, since the cuts were applied to detector level
                //     // Fill input distribution
                //     hDeltaRInputMCP->Fill(MCPjetPt, MCPDeltaR);
                // }
                if (!binning.useEmmaYeatsBins || passEmmaCut(MCPjetPt, MCPhfPt)) { // --> this is new!
                    // Emma cut should only be checked for the matched detector level entry MCDs?
                    // Verify Emma's cuts only in case these are matched entries, since the cuts were applied to detector level
                    // Fill input distribution
                    hDeltaRInputMCP->Fill(MCPjetPt, MCPDeltaR);
                }
            }
        }
    }
    // Check that eraseHistogram is selecting the correct bins
    for (size_t iJetPt = 0; iJetPt < binning.ptjetBinEdges_particle.size() - 1; ++iJetPt) {

        // Check if this pT,HF bin survived the sideband fit for this jet pT bin
        std::vector<bool> workingFits = sidebandDataContainer.workingFitsPerJetPt[iJetPt];
        std::cout << "Erase histogram check for jets in [" << binning.ptjetBinEdges_particle[iJetPt] << ", " << binning.ptjetBinEdges_particle[iJetPt+1] << "] GeV/c:" << std::endl;
        std::cout << "Bins: ";
        for (size_t iHFPt = 0; iHFPt < workingFits.size(); ++iHFPt) { // binning.ptHFBinEdges_particle
            
            if (workingFits[iHFPt]) {
                std::cout << "True | ";
            } else {
                std::cout << "False | ";
            }
        }
        std::cout << std::endl;
    }
    

    // Start comparison
    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);
    
    // 1- Compare all iterations with the MCP input
    TCanvas* cFullyCorrected1D = new TCanvas("cFullyCorrected1D","Fully corrected deltaR 1D distribution iterations",1800,1000);
    TLegend* lUnfoldedIter = new TLegend(0.6,0.57,0.7,0.77);
    std::vector<TH1D*> hUnfoldedKinCorrectedProj(unfoldDataContainer.hUnfoldedKinCorrected.size());
    int secondBin = hDeltaRInputMCP->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[1] + 1e-6);
    int lastButOneBin = hDeltaRInputMCP->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[binning.ptjetBinEdges_particle.size()-2] - 1e-6);
    TH1D* hDeltaRInputMCPProj = hDeltaRInputMCP->ProjectionY("hDeltaRInputMCPProj",secondBin,lastButOneBin);
    for (size_t iIter = 0; iIter < unfoldDataContainer.hUnfoldedKinCorrected.size(); iIter++) {
        cFullyCorrected1D->cd();
        secondBin = unfoldDataContainer.hUnfoldedKinCorrected[iIter]->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[1] + 1e-6);
        lastButOneBin = unfoldDataContainer.hUnfoldedKinCorrected[iIter]->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[binning.ptjetBinEdges_particle.size()-2] - 1e-6);
        hUnfoldedKinCorrectedProj[iIter] = unfoldDataContainer.hUnfoldedKinCorrected[iIter]->ProjectionY(Form("hUnfoldedKinCorrected_iter%zu",iIter), secondBin, lastButOneBin);
        hUnfoldedKinCorrectedProj[iIter]->SetLineColor(kBlack + iIter);
        lUnfoldedIter->AddEntry(hUnfoldedKinCorrectedProj[iIter],Form("Iteration %zu", iIter+1), "le");
        if (iIter == 0) {
            hUnfoldedKinCorrectedProj[iIter]->GetYaxis()->SetRangeUser(0.,hDeltaRInputMCPProj->GetMaximum() * 1.2);
            hUnfoldedKinCorrectedProj[iIter]->SetTitle("Fully corrected distribution projection");
            hUnfoldedKinCorrectedProj[iIter]->Draw();
        } else {
            hUnfoldedKinCorrectedProj[iIter]->Draw("same");
        }
    }
    hDeltaRInputMCPProj->GetYaxis()->SetTitle("dN");
    hDeltaRInputMCPProj->SetLineStyle(kDashed);
    hDeltaRInputMCPProj->SetLineColor(kRed);
    hDeltaRInputMCPProj->Sumw2();
    hDeltaRInputMCPProj->Draw("same");
    lUnfoldedIter->AddEntry(hDeltaRInputMCPProj,"Input particle level", "le");
    double reportingJetPtMin = unfoldDataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetBinLowEdge(secondBin);
    double reportingJetPtMax = unfoldDataContainer.hUnfoldedKinCorrected[0]->GetXaxis()->GetBinUpEdge(lastButOneBin);
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));

    // 2 - Compare only the 4th iteration with the MCP input
    TCanvas* c4thIterationFullyCorrected = new TCanvas("c4thIterationFullyCorrected","Unfolded (and kin. eff. corrected) 4th iteration fully corrected distribution",1920,1080);
    c4thIterationFullyCorrected->cd();
    hUnfoldedKinCorrectedProj[selectedIteration-1]->GetYaxis()->SetRangeUser(0.,hUnfoldedKinCorrectedProj[selectedIteration-1]->GetMaximum() * 1.2);
    hDeltaRInputMCPProj->Draw();
    hUnfoldedKinCorrectedProj[selectedIteration-1]->Draw("same");
    TLegend* leg4thIterFullyCorrected = new TLegend(0.6,0.57,0.7,0.77);
    leg4thIterFullyCorrected->AddEntry(hUnfoldedKinCorrectedProj[selectedIteration-1],Form("Iteration %zu", selectedIteration), "le");
    leg4thIterFullyCorrected->AddEntry(hDeltaRInputMCPProj,"Input particle level", "le");
    leg4thIterFullyCorrected->Draw();
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));

    // 3 -Compare only the 4th iteration with the MCP input normalized
    TCanvas* c4thIterationFullyCorrectedNorm = new TCanvas("c4thIterationFullyCorrectedNorm","Unfolded (and kin. eff. corrected) 4th iteration fully corrected distribution",1920,1080);
    TH1D* hUnfoldedKinCorrectedProj4thNorm = (TH1D*) hUnfoldedKinCorrectedProj[selectedIteration-1]->Clone("hUnfoldedKinCorrectedProj4thNorm");
    hUnfoldedKinCorrectedProj4thNorm->Scale(1 / hUnfoldedKinCorrectedProj4thNorm->Integral(),"width");
    TH1D* hDeltaRInputMCPProjNorm = (TH1D*) hDeltaRInputMCPProj->Clone("hDeltaRInputMCPProjNorm");
    hDeltaRInputMCPProjNorm->Scale(1 / hDeltaRInputMCPProjNorm->Integral(),"width");
    hDeltaRInputMCPProjNorm->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
    c4thIterationFullyCorrectedNorm->cd();
    hDeltaRInputMCPProjNorm->Draw();
    hUnfoldedKinCorrectedProj4thNorm->Draw("same");
    TLegend* leg4thIterFullyCorrectedNorm = new TLegend(0.6,0.57,0.7,0.77);
    leg4thIterFullyCorrectedNorm->AddEntry(hUnfoldedKinCorrectedProj4thNorm, "Iteration 4", "le");
    leg4thIterFullyCorrectedNorm->AddEntry(hDeltaRInputMCPProjNorm,"Input particle level", "le");
    leg4thIterFullyCorrectedNorm->Draw();
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));

    // 4 - Divide the 4th iteration by the MCP input
    TCanvas* cRatio = new TCanvas("cRatio","Input divided by output",1800,1000);
    TH1D* hRatio = (TH1D*) hDeltaRInputMCPProj->Clone("hRatio");
    hRatio->Sumw2();
    hUnfoldedKinCorrectedProj[selectedIteration-1]->Sumw2();
    hRatio->Divide(hUnfoldedKinCorrectedProj[selectedIteration-1]);
    hRatio->SetTitle("Input distribution divided by fully corrected");
    hRatio->GetYaxis()->SetTitle("#frac{Input reference}{fully corrected}");
    hRatio->SetLineColor(kBlue);
    hRatio->SetLineStyle(kSolid);
    hRatio->GetYaxis()->SetRangeUser(0.75 * std::min(hRatio->GetMinimum(), 1.), 1.1 * std::max(hRatio->GetMaximum(), 1.25));
    hRatio->Draw();
    TLine* lineHorizontal = new TLine(0., 1., hRatio->GetXaxis()->GetBinUpEdge(hRatio->GetXaxis()->GetNbins()), 1.);
    lineHorizontal->SetLineStyle(kDashed);
    lineHorizontal->SetLineColor(kBlack);
    lineHorizontal->Draw("same");
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));

    // 5 - Evolution through the analysis steps
    TLegend* legEvol = new TLegend(0.6,0.57,0.80,0.77);
    TCanvas* cAnalysisEvolution = new TCanvas("cAnalysisEvolution","Evolution through the analysis steps",1920,1080);
    cAnalysisEvolution->Divide(2,2);
    TH2D* hSidebandSub2D = (TH2D*)sidebandDataContainer.hBackgroundSubtracted->Project3D("yx");
    int xBinMin = hSidebandSub2D->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[1] + 1e-6);
    int xBinMax = hSidebandSub2D->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[binning.ptjetBinEdges_particle.size()-2] - 1e-6);
    TH1D* hSidebandSub1D = hSidebandSub2D->ProjectionY("hSidebandSub1D", xBinMin, xBinMax);
    hSidebandSub1D->SetLineColor(kMagenta+1);
    xBinMin = efficiencyDatacontainer.hEfficiencyCorrectedData.second->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[1] + 1e-6);
    xBinMax = efficiencyDatacontainer.hEfficiencyCorrectedData.second->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[binning.ptjetBinEdges_particle.size()-2] - 1e-6);
    TH1D* hEfficienCorr = efficiencyDatacontainer.hEfficiencyCorrectedData.second->ProjectionY("hEfficienCorr", xBinMin, xBinMax);
    hEfficienCorr->SetLineColor(kGreen+2);
    hEfficienCorr->GetYaxis()->SetTitle("dN");
    legEvol->AddEntry(hDeltaRInputMCPProj,"Input particle","lp");
    legEvol->AddEntry(hSidebandSub1D,"Background subtracted","lp");
    legEvol->AddEntry(hEfficienCorr,"Efficiency corrected","lp");
    legEvol->AddEntry(hUnfoldedKinCorrectedProj[selectedIteration-1],"Unfolded fully corrected","lp");
    cAnalysisEvolution->cd(1);
    hEfficienCorr->Draw();
    hSidebandSub1D->Draw("same");
    cAnalysisEvolution->cd(2);
    hEfficienCorr->Draw();
    hUnfoldedKinCorrectedProj[selectedIteration-1]->Draw("same");
    cAnalysisEvolution->cd(3);
    hDeltaRInputMCPProj->Draw();
    hUnfoldedKinCorrectedProj[selectedIteration-1]->Draw("same");
    cAnalysisEvolution->cd(4);
    double maxY = 1.1 * std::max({hDeltaRInputMCPProj->GetMaximum(), hSidebandSub1D->GetMaximum(), hEfficienCorr->GetMaximum(), hUnfoldedKinCorrectedProj[selectedIteration-1]->GetMaximum()});
    hDeltaRInputMCPProj->GetYaxis()->SetRangeUser(0.,maxY);
    hDeltaRInputMCPProj->Draw();
    hSidebandSub1D->Draw("same");
    hEfficienCorr->Draw("same");
    hUnfoldedKinCorrectedProj[selectedIteration-1]->Draw("same");
    legEvol->Draw();
    latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", reportingJetPtMin, reportingJetPtMax));

    // 6 - Evolution through the analysis steps for each pT,jet bin
    std::vector<TCanvas*> cEvolPerJetPt(binning.ptjetBinEdges_particle.size() - 1);
    std::vector<TLegend*> legEvolPerJetPt(binning.ptjetBinEdges_particle.size() - 1);
    std::vector<TH1D*> hSidebandSub1DPerJetPt(binning.ptjetBinEdges_particle.size() - 1);
    std::vector<TH1D*> hEfficienCorrPerJetPt(binning.ptjetBinEdges_particle.size() - 1);
    std::vector<TH1D*> hUnfoldedKinCorrectedProjPerJetPt(binning.ptjetBinEdges_particle.size() - 1);
    std::vector<TH1D*> hDeltaRInputMCPProjPerJetPt(binning.ptjetBinEdges_particle.size() - 1);
    for (size_t iJetBin = 0; iJetBin < binning.ptjetBinEdges_particle.size() - 1; iJetBin++) {

        // Side-band subtracted
        xBinMin = hSidebandSub2D->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[iJetBin] + 1e-6);
        xBinMax = hSidebandSub2D->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[iJetBin+1] - 1e-6);
        hSidebandSub1DPerJetPt[iJetBin] = hSidebandSub2D->ProjectionY(Form("hSidebandSub1D_%zu", iJetBin), xBinMin, xBinMax);
        hSidebandSub1DPerJetPt[iJetBin]->SetLineColor(kMagenta+1);
        hSidebandSub1DPerJetPt[iJetBin]->GetYaxis()->SetTitle("dN");

        // Efficiency corrected
        xBinMin = efficiencyDatacontainer.hEfficiencyCorrectedData.second->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[iJetBin] + 1e-6);
        xBinMax = efficiencyDatacontainer.hEfficiencyCorrectedData.second->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[iJetBin+1] - 1e-6);
        hEfficienCorrPerJetPt[iJetBin] = efficiencyDatacontainer.hEfficiencyCorrectedData.second->ProjectionY(Form("hEfficienCorr_bin%zu", iJetBin), xBinMin, xBinMax);
        // hEfficienCorrPerJetPt[iJetBin] = efficiencyDatacontainer.hEfficiencyCorrectedData.second->ProjectionY(Form("hEfficienCorr_bin%zu", iJetBin), iJetBin+1, iJetBin+1);
        hEfficienCorrPerJetPt[iJetBin]->SetLineColor(kGreen+2);
        hEfficienCorrPerJetPt[iJetBin]->GetYaxis()->SetTitle("dN");

        // Unfolded, kinematic efficiency corrected
        xBinMin = unfoldDataContainer.hUnfoldedKinCorrected[selectedIteration-1]->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[iJetBin] + 1e-6);
        xBinMax = unfoldDataContainer.hUnfoldedKinCorrected[selectedIteration-1]->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[iJetBin+1] - 1e-6);
        hUnfoldedKinCorrectedProjPerJetPt[iJetBin] = unfoldDataContainer.hUnfoldedKinCorrected[selectedIteration-1]->ProjectionY(Form("hUnfoldedKinCorrectedProj_bin%zu", iJetBin), xBinMin, xBinMax);
        hUnfoldedKinCorrectedProjPerJetPt[iJetBin]->SetLineColor(kBlue+2);
        hUnfoldedKinCorrectedProjPerJetPt[iJetBin]->GetYaxis()->SetTitle("dN");

        // Reference input distribution
        xBinMin = hDeltaRInputMCP->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[iJetBin] + 1e-6);
        xBinMax = hDeltaRInputMCP->GetXaxis()->FindBin(binning.ptjetBinEdges_particle[iJetBin+1] - 1e-6);
        hDeltaRInputMCPProjPerJetPt[iJetBin] = hDeltaRInputMCP->ProjectionY(Form("hDeltaRInputMCPProj_bin%zu", iJetBin), xBinMin, xBinMax);
        hDeltaRInputMCPProjPerJetPt[iJetBin]->SetLineColor(kRed);
        hDeltaRInputMCPProjPerJetPt[iJetBin]->SetLineStyle(kDashed);
        hDeltaRInputMCPProjPerJetPt[iJetBin]->Sumw2();
        hDeltaRInputMCPProjPerJetPt[iJetBin]->GetYaxis()->SetTitle("dN");

        // Create canvas and legend
        cEvolPerJetPt[iJetBin] = new TCanvas(Form("cEvolPerJetPt_bin%zu", iJetBin), Form("Evolution through the analysis steps for p_{T,jet} in [%.1f, %.1f] GeV/c", binning.ptjetBinEdges_particle[iJetBin], binning.ptjetBinEdges_particle[iJetBin+1]), 1920, 1080);
        legEvolPerJetPt[iJetBin] = new TLegend(0.6,0.57,0.80,0.77);
        legEvolPerJetPt[iJetBin]->AddEntry(hDeltaRInputMCPProj,"Input particle","lp");
        legEvolPerJetPt[iJetBin]->AddEntry(hSidebandSub1DPerJetPt[iJetBin],"Background subtracted","lp");
        legEvolPerJetPt[iJetBin]->AddEntry(hEfficienCorrPerJetPt[iJetBin],"Efficiency corrected","lp");
        legEvolPerJetPt[iJetBin]->AddEntry(hUnfoldedKinCorrectedProjPerJetPt[iJetBin],"Unfolded fully corrected","lp");
        legEvolPerJetPt[iJetBin]->AddEntry(hDeltaRInputMCPProjPerJetPt[iJetBin],"Input particle","lp");

        // Draw histograms
        cEvolPerJetPt[iJetBin]->cd();
        double maxY = 1.2 * std::max({hDeltaRInputMCPProjPerJetPt[iJetBin]->GetMaximum(), hSidebandSub1DPerJetPt[iJetBin]->GetMaximum(), hEfficienCorrPerJetPt[iJetBin]->GetMaximum(), hUnfoldedKinCorrectedProjPerJetPt[iJetBin]->GetMaximum()});
        hDeltaRInputMCPProjPerJetPt[iJetBin]->GetYaxis()->SetRangeUser(0.,maxY);
        hDeltaRInputMCPProjPerJetPt[iJetBin]->Draw();
        hSidebandSub1DPerJetPt[iJetBin]->Draw("same");
        hEfficienCorrPerJetPt[iJetBin]->Draw("same");
        hUnfoldedKinCorrectedProjPerJetPt[iJetBin]->Draw("same");
        legEvolPerJetPt[iJetBin]->Draw();
        latex->DrawLatex(0.15, 0.85, Form("Projected in p_{T,jet} #in [%.1f, %.1f] GeV/c", binning.ptjetBinEdges_particle[iJetBin], binning.ptjetBinEdges_particle[iJetBin+1]));
    }
    
    // 7 - Ratio of input over fully corrected for each pT,jet bin
    TCanvas* cRatioPerJetPt = new TCanvas("cRatioPerJetPt","Ratio of input over fully corrected for each pT,jet bin",1800,1000); // 1920 x 1080?
    std::vector<TH1D*> hRatioPerJetPt(binning.ptjetBinEdges_particle.size() - 1);
    TLegend* legRatioPerJetPt = new TLegend(0.6,0.57,0.77,0.77);
    double minHistoRange = 1;
    double maxHistoRange = 1;
    // build ratio histograms
    for (size_t iJetBin = 0; iJetBin < binning.ptjetBinEdges_particle.size() - 1; iJetBin++) {
        
        hRatioPerJetPt[iJetBin] = (TH1D*) hDeltaRInputMCPProjPerJetPt[iJetBin]->Clone(Form("hRatioPerJetPt_bin%zu", iJetBin));
        hRatioPerJetPt[iJetBin]->Sumw2();
        hUnfoldedKinCorrectedProjPerJetPt[iJetBin]->Sumw2();
        hRatioPerJetPt[iJetBin]->Divide(hUnfoldedKinCorrectedProjPerJetPt[iJetBin]);
        hRatioPerJetPt[iJetBin]->SetTitle(Form("Input distribution divided by fully corrected for each p_{T,jet} bin  GeV/c"));
        hRatioPerJetPt[iJetBin]->GetYaxis()->SetTitle("#frac{Input reference}{fully corrected}");
        hRatioPerJetPt[iJetBin]->SetLineColor(kBlack + iJetBin);
        hRatioPerJetPt[iJetBin]->SetLineStyle(kSolid);
        if (hRatioPerJetPt[iJetBin]->GetMinimum() < minHistoRange) {
            minHistoRange = hRatioPerJetPt[iJetBin]->GetMinimum();
        }
        if (hRatioPerJetPt[iJetBin]->GetMaximum() > maxHistoRange) {
            maxHistoRange = hRatioPerJetPt[iJetBin]->GetMaximum();
        }
        legRatioPerJetPt->AddEntry(hRatioPerJetPt[iJetBin], Form("p_{T,jet} #in [%.1f, %.1f] GeV/c", binning.ptjetBinEdges_particle[iJetBin], binning.ptjetBinEdges_particle[iJetBin+1]), "le");
    }
    // plot histograms
    std::cout << "minHistoRange: " << minHistoRange << ", maxHistoRange: " << maxHistoRange << std::endl;
    for (size_t iJetBin = 0; iJetBin < binning.ptjetBinEdges_particle.size() - 1; iJetBin++) {
        if (iJetBin == 0) {
            minHistoRange = -0.5;
            maxHistoRange = 2.8;
            hRatioPerJetPt[iJetBin]->GetYaxis()->SetRangeUser(0.9 * minHistoRange, 1.1 * maxHistoRange);
            hRatioPerJetPt[iJetBin]->Draw();
        } else {
            hRatioPerJetPt[iJetBin]->Draw("same");
        }
    }
    lineHorizontal->Draw("same");
    legRatioPerJetPt->Draw();

    //
    // Storing images
    //
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../Images/5-ClosureTest/Second/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    cFullyCorrected1D->Print(imagePath + Form("closureTest2_comparison_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    c4thIterationFullyCorrected->Print(imagePath + Form("closureTest2_comparison_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    c4thIterationFullyCorrectedNorm->Print(imagePath + Form("closureTest2_comparison_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cRatio->Print(imagePath + Form("closureTest2_comparison_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cAnalysisEvolution->Print(imagePath + Form("closureTest2_comparison_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    for (size_t iJetBin = 0; iJetBin < binning.ptjetBinEdges_particle.size() - 1; iJetBin++) {
        cEvolPerJetPt[iJetBin]->Print(imagePath + Form("closureTest2_comparison_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    }
    cRatioPerJetPt->Print(imagePath + Form("closureTest2_comparison_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
    

    return hDeltaRInputMCP;
}

// 2 - Sideband subtraction + efficiency correction + unfolding closure test
void SecondClosureTest(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // Number of unfolding procedure iterations
    int iterationNumber = 4;

    bool useEmmaYeatsBins = false;
    TString sEmmaBins;
    if (useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    // Select data period to retrieve the corresponding binning information and BDT score thresholds
    TString dataPeriod = "2023";
    // Open binning information file with optimal BDT score thresholds
    TFile* fBinning = new TFile("../1-SignalTreatment/BDTOptimization/binningInfo_" + dataPeriod + "_" + sEmmaBins + ".root", "read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    // Force by hand new binning regardless of previous steps
    binning.useEmmaYeatsBins = useEmmaYeatsBins;
    if (binning.useEmmaYeatsBins) {
        // pT,jet cuts
        binning.ptjetBinEdges_detector = {5., 7., 10., 20., 50.};
        binning.ptjetBinEdges_particle = binning.ptjetBinEdges_detector;
        // DeltaR bins
        binning.deltaRBinEdges_detector = {0., 0.01, 0.03, 0.05, 0.12, 0.2};
        binning.deltaRBinEdges_particle = binning.deltaRBinEdges_detector;
    }
    // Choose file period
    binning.dataPeriod = dataPeriod;
    if (binning.dataPeriod == "2023") {
        // DATA
        binning.inputDATA.first = "JE_HF_LHC23_pass4_Thin_2P3PDstar_D0CJ_4_D0_1";
        binning.inputDATA.second = "Data/Experimental/Train_643652";
        // Anchored MC
        binning.inputMC.first = "HF_LHC24h1c_All_D0";
        binning.inputMC.second = "Data/MonteCarlo/Train_671273";
    } else if (binning.dataPeriod == "2022") {
        // DATA
        binning.inputDATA.first = "JE_HF_LHC22o_pass7_minBias_2P3PDstar_D0CJ_4_D0_1";
        binning.inputDATA.second = "Data/Experimental/Train_659513";
        // Anchored MC
        binning.inputMC.first = "HF_LHC24g5_All_D0";
        binning.inputMC.second = "Data/MonteCarlo/Train_669231";
    }
    
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    double minDeltaR = binning.deltaRBinEdges_detector[0];
    double maxDeltaR = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1];
    double hfptMin = binning.ptHFBinEdges_detector[0]; //ptHFBinEdges[0] - should start from 0 or from the lowest pT,D value?
    double hfptMax = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1];


    // Opening files
    TFile* fClosureInput = new TFile("mc_closure_input_data_"+ binning.dataPeriod + ".root","read");
    if (!fClosureInput || fClosureInput->IsZombie()) {
        std::cerr << "Error: Unable to open 1st closure input ROOT file." << std::endl;
    }

    // ----1: Perform side-band subtraction
    // Example: selecting a model
    FitModelType modelToUse = FitModelType::SignalReflectionsOnly;
    SidebandClosureResult sidebandDataContainer = SidebandClosure(fClosureInput, binning, modelToUse);

    // ----2: Perform efficiency correction
    EfficiencyData efficiencyDatacontainer = EfficiencyClosure(fClosureInput, sidebandDataContainer.hBackgroundSubtracted, binning);
    std::vector<TH1D*> hSelEff_run3style = {efficiencyDatacontainer.hSelectionEfficiency.first,efficiencyDatacontainer.hSelectionEfficiency.second};
    TH2D* hEfficiencyCorrectedData = efficiencyDatacontainer.hEfficiencyCorrectedData.second;

    // ----3: Perform unfolding
    UnfoldData unfoldDataContainer = UnfoldingClosure(fClosureInput, hEfficiencyCorrectedData, hSelEff_run3style, binning);

    // ----4: Compare input (MC particle level) with output (background subtracted, efficiency corrected, unfolded) distributions
    TH2D* hDeltaRInputMCP = CompareClosureTest(fClosureInput, sidebandDataContainer, efficiencyDatacontainer, unfoldDataContainer, binning);

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    std::time_t now = std::time(nullptr);
    std::cout << "Finished at: " << std::ctime(&now);
}

int main(){
    SecondClosureTest();
    return 0;
}
