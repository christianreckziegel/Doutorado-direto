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

std::pair<TH1D*, std::vector<TBox*>> createEmmasDistribution(TCanvas*& cEmmaRawPlot) {
    // Define bin edges
    const int nbins = 4;
    double bins[nbins+1] = {0.0, 0.01, 0.03, 0.05, 0.12};
    
    // Create TH1D with the binning
    TH1D* hDeltaR = new TH1D("hDeltaREmma", "Emma's distribution;#DeltaR;#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}", nbins, bins);
    
    // Set bin values and statistical errors
    double values[nbins] = {16.6786, 14.4632, 13.0093, 3.49983};
    double stat_errors[nbins] = {2.27639, 2.03747, 2.33001, 0.703231};
    double syst_errors[nbins] = {1.07718, 0.365236, 0.455006, 0.169006};
    
    for (int i = 0; i < nbins; i++) {
        hDeltaR->SetBinContent(i+1, values[i]);
        hDeltaR->SetBinError(i+1, stat_errors[i]);
    }
    
    // Style the main histogram
    hDeltaR->SetLineColor(kRed);
    hDeltaR->SetMarkerColor(kRed);
    //hDeltaR->SetMarkerStyle(20);
    hDeltaR->SetMarkerSize(1.0);
    hDeltaR->SetFillStyle(0);
    
    // Create systematic uncertainty boxes
    std::vector<TBox*> systBoxes(nbins);
    
    for (int i = 0; i < nbins; i++) {
        double xlow = hDeltaR->GetBinLowEdge(i+1);
        double xhigh = xlow + hDeltaR->GetBinWidth(i+1);
        double yval = values[i];
        double ysyst = syst_errors[i];
        
        systBoxes[i] = new TBox(xlow, yval - ysyst, xhigh, yval + ysyst);
        systBoxes[i]->SetFillColor(kRed);
        systBoxes[i]->SetFillStyle(3001); // 3001 = transparent fill
        systBoxes[i]->SetLineColor(kRed);
        systBoxes[i]->SetLineStyle(2);
    }
    gStyle->SetOptStat(0);
    // Draw everything
    cEmmaRawPlot = new TCanvas("cEmmaRawPlot", "DeltaR Distribution with Systematic Bands", 1800, 1000);
    cEmmaRawPlot->SetGrid();
    
    hDeltaR->Draw();

    // Draw systematic boxes first (so they are in the background)
    for (int i = 0; i < nbins; i++) {
        systBoxes[i]->Draw("l same");
    }
    
    // Then draw the data points with statistical errors
    //hDeltaR->Draw("E1 SAME");
    
    // Customize the plot
    hDeltaR->GetYaxis()->SetRangeUser(0, 20); // Adjust range as needed
    
    // Create legend
    // TLegend* leg = new TLegend(0.65, 0.75, 0.85, 0.85);
    // leg->SetBorderSize(0);
    // leg->SetFillStyle(0);
    // leg->AddEntry(hDeltaR, "Data #pm stat.", "lep");
    // leg->AddEntry(systBoxes[0], "Syst. uncertainty", "f");
    // leg->Draw();

    // Create legend
    TLegend* leg = new TLegend(0.65, 0.75, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hDeltaR, "Data #pm stat.", "lep");
    leg->AddEntry(systBoxes[0], "Syst. uncertainty", "f");
    leg->Draw();

    // Add TLatex text with all the information
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.035);
    
    // Top left corner - Process information
    latex->SetTextAlign(11); // left-bottom aligned
    latex->DrawLatex(0.15, 0.85, "pp #rightarrow D^{0}, pp #rightarrow #bar{D}^{0}");
    latex->DrawLatex(0.15, 0.81, "#sqrt{s} = 5.02 TeV");
    
    // Top right corner - Jet information  
    latex->SetTextAlign(31); // right-bottom aligned
    latex->DrawLatex(0.85, 0.65, "anti-k_{T}, R = 0.4");
    latex->DrawLatex(0.85, 0.61, "p_{T}^{ch jet}: 10 - 20 GeV/c");
    latex->DrawLatex(0.85, 0.57, "|#eta_{jet}| #leq 0.5");
    
    // Bottom left corner - D0 information
    latex->SetTextAlign(11); // left-bottom aligned
    latex->DrawLatex(0.15, 0.25, "p_{T}^{D^{0}}: 5 - 20 GeV/c");
    latex->DrawLatex(0.15, 0.21, "|y_{D^{0}}| #leq 0.8");

    return std::make_pair(hDeltaR, systBoxes);
}

// Get my unfolded data DeltaR distribution
TH1D* getMyDataDistribution(TFile* fMyData, const BinningStruct& binning) {
    TCanvas* cMyDataNormalization = new TCanvas("cMyDataNormalization", "My Data Normalization Check", 1800, 1000);
    cMyDataNormalization->SetGrid();
    cMyDataNormalization->Divide(2,2);
    TH2D* hpTjet_vs_DeltaR = (TH2D*)fMyData->Get("Unfolded/hUnfoldedKinCor_iter4");
    if (!hpTjet_vs_DeltaR) {
        throw std::runtime_error("Could not find TH2D 'Unfolded/hUnfoldedKinCor_iter4' in the provided file.");
    }
    // pT,D⁰ #in [5;20] GeV/c -> this is ensured since the beginning of the analysis on background subtraction
    
    // pT,jet #in [10;20] GeV/c
    
    ptjetBinMin = hpTjet_vs_DeltaR->GetXaxis()->FindBin(10.0); // Find the bin corresponding to 10 GeV/c
    ptjetBinMax = hpTjet_vs_DeltaR->GetXaxis()->FindBin(20.0) - 1; // Find the bin corresponding to 20 GeV/c and take the previous bin as max to be inclusive of the upper edge
    std::cout << "X bins of 2D histogram: " << hpTjet_vs_DeltaR->GetNbinsX() << std::endl;
    std::cout << "Selected pT,jet bin range: [" << hpTjet_vs_DeltaR->GetXaxis()->GetBinLowEdge(ptjetBinMin) << ";" << hpTjet_vs_DeltaR->GetXaxis()->GetBinUpEdge(ptjetBinMax) << "] GeV/c" << std::endl;
    // std::cout << "Y bins of 2D histogram: " << hpTjet_vs_DeltaR->GetNbinsY() << std::endl;
    // std::cout << "Selected pT,D⁰ bin range: [" << hpTjet_vs_DeltaR->GetYaxis()->GetBinLowEdge(deltaRBinMin) << ";" << hpTjet_vs_DeltaR->GetYaxis()->GetBinUpEdge(deltaRBinMax) << "] GeV/c" << std::endl;

    // Project onto DeltaR axis (Y-axis) - only specify X bin range
    TH1D* hDeltaR = hpTjet_vs_DeltaR->ProjectionY("hMyDeltaRDistribution", ptjetBinMin, ptjetBinMax);

    // ΔR #in [0;0.12] - specify Y-axis bin range too!
    int deltaRBinMin = hpTjet_vs_DeltaR->GetYaxis()->FindBin(0.0);
    int deltaRBinMax = hpTjet_vs_DeltaR->GetYaxis()->FindBin(0.12) - 1;

    // Get DeltaR bin edges from detector level binning
    // keep only bins <= 0.12
    std::vector<double> deltaRBinEdges_truncated;

    for (double edge : binning.deltaRBinEdges_detector) {
        if (edge <= 0.12)
            deltaRBinEdges_truncated.push_back(edge);
    }

    
    
    // Check without truncation of the last bin
    bool doTruncation = true;
    TH1D* hDeltaR_final;
    if (doTruncation) {
        // Create a new histogram with only the desired DeltaR range
        hDeltaR_final = new TH1D("hMyDeltaRDistribution_final_v1", ";#DeltaR^{truth};#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}", deltaRBinEdges_truncated.size() - 1,deltaRBinEdges_truncated.data());
        
        // Copy the content and errors
        for (int i = 1; i <= hDeltaR_final->GetNbinsX(); i++) {
            int origBin = deltaRBinMin + i - 1;
            hDeltaR_final->SetBinContent(i, hDeltaR->GetBinContent(origBin));
            hDeltaR_final->SetBinError(i, hDeltaR->GetBinError(origBin));
        }
    } else {
        // If not truncating, just clone the original histogram and rename it
        hDeltaR_final = (TH1D*)hDeltaR->Clone("hMyDeltaRDistribution_final");
        hDeltaR_final->SetTitle(";#DeltaR^{truth};#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
    }
    

    // Clean up the temporary histogram
    delete hDeltaR;

    // Style the histogram
    hDeltaR_final->SetLineColor(kBlue);
    hDeltaR_final->SetMarkerColor(kBlue);
    //hDeltaR_final->SetMarkerStyle(21);
    hDeltaR_final->SetMarkerSize(1.0);
    hDeltaR_final->SetFillStyle(0);

    cMyDataNormalization->cd(1);
    hDeltaR_final->GetYaxis()->SetTitle("dN");
    hDeltaR_final->Draw();
    // NORMALIZE the distribution
    TH1D* hDeltaR_final_v2 = (TH1D*)hDeltaR_final->Clone("hMyDeltaRDistribution_final_v2");
    double integral = hDeltaR_final_v2->Integral(); // Integral with bin width
    hDeltaR_final_v2->Scale(1.0, "width");
    hDeltaR_final_v2->GetYaxis()->SetTitle("#frac{dN}{d#DeltaR}");
    cMyDataNormalization->cd(2);
    hDeltaR_final_v2->Draw();
    TH1D* hDeltaR_final_v3 = (TH1D*)hDeltaR_final_v2->Clone("hMyDeltaRDistribution_final");
    hDeltaR_final_v3->Scale(1.0 / integral);
    hDeltaR_final_v3->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
    cMyDataNormalization->cd(3);
    hDeltaR_final_v3->Draw();

    // count manually the number of jets for each bin
    int numberOfJets = 0;
    std::cout << "My original histogram has " << hDeltaR_final->GetNbinsX() << " bins" << std::endl;
    for (int iBin = 1; iBin <= hDeltaR_final->GetNbinsX(); iBin++) {
        numberOfJets += hDeltaR_final->GetBinContent(iBin);
    }
    std::cout << "Number of jets in my data distribution: " << numberOfJets << std::endl;
    std::cout << "Integral (with bin width) before normalization: " << integral << std::endl;

    // if (integral > 0) {
    //     hDeltaR_final->Scale(1.0 / integral, "width");
    // }
    
    return hDeltaR_final_v3;
}

TH1D* createPythiaDistribution(TFile* fSimulatedMCMatched, const BinningStruct& binning) {
    // Create histogram with Emma's binning
    TH1D* hMyDeltaRPythia = new TH1D("hMyDeltaRPythia", "Pythia distribution;#DeltaR;#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}", binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());

    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level HF
    const double MCDyCut = 0.8; // on detector level HF
    const double MCPDeltaRcut = binning.deltaRBinEdges_particle[binning.deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double MCPHfPtMincut = 5.; // binning.ptHFBinEdges_particle[0]; // on particle level HF
    const double MCDHfPtMincut = 5.; // binning.ptHFBinEdges_detector[0]; // on detector level HF
    const double MCPHfPtMaxcut = 20.; // binning.ptHFBinEdges_particle[binning.ptHFBinEdges_particle.size() - 1]; // on particle level HF
    const double MCDHfPtMaxcut = 20.; // binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1]; // on detector level HF
    const double jetptMin = 10.; // binning.ptjetBinEdges_particle[0];
    const double jetptMax = 20.; // binning.ptjetBinEdges_particle[binning.ptjetBinEdges_particle.size() - 1];

    // Get TTree
    TTree* tree = (TTree*)fSimulatedMCMatched->Get("DF_merged/O2matchtable");
    // Check for correct access
    if (!tree) {
        cout << "Error opening O2 matching tree.\n";
    }

    // Fill histogram
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
    tree->SetBranchAddress("fMcJetNConst",&MCPjetNConst);
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
    tree->SetBranchAddress("fJetNConst",&MCDjetNConst);
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

        // Generator level selection cuts
        double MCPDeltaR = MCPaxisDistance;
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);
        bool genHfPtRange = ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        bool genDeltaRRange = ((MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < MCPDeltaRcut));
        bool genAcceptance = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut);
        bool genLevelRange = genAcceptance && genJetPtRange && genDeltaRRange && genHfPtRange; // --> this is new!
        // Reconstruction level selection cuts
        double MCDDeltaR = MCDaxisDistance;
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb);
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        bool recoDeltaRRange = ((MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < MCDDeltaRcut));
        bool recoAcceptance = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut);
        bool recoLevelRange = recoAcceptance && recoJetPtRange && recoHfPtRange && recoDeltaRRange;

        // Apply non-prompt selection (i.e., only B → D0)
        //bool isReflection = (MCDhfMatchedFrom != MCDhfSelectedAs) ? true : false;
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isRealD0 = isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs); // can't be reflection

        // Fill with particle level, real (non-reflections) prompt D0s (regardless of their reconstruction level matching) that are within the defined cuts
        if (MCPhfprompt && genLevelRange) {
            // Matched entries can't be reflections
            if (MCDhfmatch) {
                if (isRealD0) {
                    if (passEmmaCut(MCPjetPt, MCPhfPt)) {
                        hMyDeltaRPythia->Fill(MCPDeltaR);
                    }
                }
            } else {
                if (passEmmaCut(MCPjetPt, MCPhfPt)) {
                    hMyDeltaRPythia->Fill(MCPDeltaR);
                }
            }
        }
    }

    // Normalize by bin width and number of jets
    hMyDeltaRPythia->Scale(1.0 / hMyDeltaRPythia->Integral(), "width");

    // Drawing style
    hMyDeltaRPythia->SetLineColor(kBlue);
    hMyDeltaRPythia->SetMarkerColor(kBlue);
    hMyDeltaRPythia->SetLineStyle(kDashed);
    //hMyDeltaRPythia->SetMarkerStyle(21);
    hMyDeltaRPythia->SetMarkerSize(1.0);
    hMyDeltaRPythia->SetFillStyle(0);

    return hMyDeltaRPythia;
}

std::pair<std::vector<TH1D*>, std::vector<TBox*>> compareEmmasAndMine(TFile* fMyData, const double& jetptMin, const double& jetptMax, const BinningStruct& binning) {
    TCanvas* cEmmaRawPlot = nullptr;
    // Create Emma's ΔR-STD distribution histogram
    std::pair<TH1D*, std::vector<TBox*>> emmaDistribution = createEmmasDistribution(cEmmaRawPlot);
    TH1D* hEmmasDeltaR = emmaDistribution.first;
    std::vector<TBox*> systBoxes = emmaDistribution.second;

    // Fetch my unfolded ΔR distribution from data
    TH1D* hDataDeltaR = getMyDataDistribution(fMyData, binning);

    // Create my pythia distribution
    TFile* fSimulatedMCMatched = new TFile("../" + binning.inputMC.second + "/AO2D_mergedDFs.root","read");
    TH1D* hMyDeltaRPythia = createPythiaDistribution(fSimulatedMCMatched, binning);

    TCanvas* cCompare = new TCanvas("cCompare", "Comparison of #DeltaR Distributions", 1800, 1000);
    cCompare->cd();
    hEmmasDeltaR->SetTitle("Comparison of #DeltaR_{STD} Distributions;#DeltaR_{STD};#frac{1}{N_{jets}}#frac{dN}{d#DeltaR_{STD}}");
    hDataDeltaR->GetYaxis()->SetRangeUser(0, 1.1 * std::max(hDataDeltaR->GetMaximum(), hEmmasDeltaR->GetMaximum())); // Adjust range as needed
    hDataDeltaR->Draw();
    hMyDeltaRPythia->Draw("same");
    hEmmasDeltaR->Draw("same");
    // hDataDeltaR->Draw("same");
    // Draw systematic boxes first (so they are in the background)
    for (int iBox = 0; iBox < systBoxes.size(); iBox++) {
        systBoxes[iBox]->Draw("l same");
    }

    TLegend* legendCompare = new TLegend(0.6, 0.7, 0.85, 0.85);
    legendCompare->SetBorderSize(0);
    legendCompare->SetFillStyle(0);
    // legendCompare->AddEntry(hEmmasDeltaR, "(Emma Yeats') Run 2, pp #rightarrow #sqrt{s} = 5.02 TeV", "lep");
    legendCompare->AddEntry(hEmmasDeltaR, "Run 2 (published), pp #rightarrow #sqrt{s} = 5.02 TeV", "lep");
    legendCompare->AddEntry(hDataDeltaR, Form("Run 3 - %s, pp #rightarrow #sqrt{s} = 13.6 TeV", binning.dataPeriod.Data()), "lep");
    legendCompare->AddEntry(hMyDeltaRPythia, Form("Anchored MC - %s, pp #rightarrow #sqrt{s} = 13.6 TeV", binning.dataPeriod.Data()), "lep");
    legendCompare->Draw();

    TLatex* latexCompare = new TLatex();
    latexCompare->SetNDC();
    latexCompare->SetTextFont(42);
    latexCompare->SetTextSize(0.035);
    latexCompare->SetTextAlign(11); // left-bottom aligned
    //latexCompare->DrawLatex(0.15, 0.85, "pp #rightarrow D^{0}, pp #rightarrow #bar{D}^{0}");
    //latexCompare->DrawLatex(0.15, 0.81, "#sqrt{s} = 5.02 TeV");
    //latexCompare->DrawLatex(0.15, 0.77, "anti-k_{T}, R = 0.4");
    latexCompare->DrawLatex(0.15, 0.86, "p_{T}^{ch jet}: 10 - 20 GeV/c, |#eta_{jet}| #leq 0.5");
    latexCompare->DrawLatex(0.15, 0.80, "p_{T}^{D^{0}}: 5 - 20 GeV/c, |y_{D^{0}}| #leq 0.8");

    //
    // Storing images
    //
    TString imagePath = "../Images/5-ClosureTest/Fourth/" + binning.dataPeriod + "/";
    cEmmaRawPlot->Print(imagePath + Form("closureTest4_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cCompare->Print(imagePath + Form("closureTest4_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
    std::vector<TH1D*> histogramsToSave = {hEmmasDeltaR, hDataDeltaR};

    return std::make_pair(histogramsToSave, emmaDistribution.second);
    

}
// Compare data final distribution to Emma Yeats' reported distribution
void FourthClosureTest() {
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Load binning from reflections file
    TFile* fBinning = new TFile(Form("../1-SignalTreatment/Reflections/binningInfo.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    // Create struct to hold all binning distributions
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    binning.dataPeriod = "2023";
    // Opening files
    TFile* fMyData = new TFile(Form("../4-Unfolding/EmmaYeatsBins/unfolding_%.0f_to_%.0f_jetpt_" + binning.dataPeriod + ".root", jetptMin, jetptMax), "read");
    if (!fMyData || fMyData->IsZombie()) {
        std::cerr << "Error: Unable to open unfolded data ROOT file." << std::endl;
    }
    // First call to get the corresponding data period, second call to fetch the binning with Emma's binning
    binning = retrieveBinningFromFile(fMyData);

    // Compare Emma Yeats' reported DeltaR distribution to my unfolded data distribution
    std::pair<std::vector<TH1D*>, std::vector<TBox*>> emmaAndMine = compareEmmasAndMine(fMyData, jetptMin, jetptMax, binning);
    TH1D* hEmmasDeltaR = emmaAndMine.first[0];
    TH1D* hMyDeltaR = emmaAndMine.first[1];
    std::vector<TBox*> systematicBands = emmaAndMine.second;

    TFile* fOutput = new TFile(Form("FourthClosureTestResults_%.0f_to_%.0fGeV_jetpt_" + binning.dataPeriod + ".root",jetptMin,jetptMax), "RECREATE");
    fOutput->cd();
    hEmmasDeltaR->Write();
    hMyDeltaR->Write();
    TDirectory* systDir = fOutput->mkdir("SystematicBands");
    std::cout << "Number of systematic bands = " << systematicBands.size() << std::endl;
    for (size_t iBox = 0; iBox < systematicBands.size(); iBox++) {
        TDirectory* binDir = systDir->mkdir(
            Form("SystematicBands_bin%d", (int)iBox)
        );
        binDir->cd();
        systematicBands[iBox]->Write();
    }
    fOutput->Close();

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
    FourthClosureTest();
    return 0;
}
