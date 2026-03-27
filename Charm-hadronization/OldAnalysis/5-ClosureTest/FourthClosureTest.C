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

std::vector<double> LoadBinning(TFile* fInput, const char* pathInFile) {
    auto* vec = (TVectorD*)fInput->Get(pathInFile);
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", pathInFile));
    }
    return std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
}

TH1D* createEmmasDistribution(TCanvas*& cEmmaRawPlot) {
    // Define bin edges
    const int nbins = 4;
    double bins[nbins+1] = {0.0, 0.01, 0.03, 0.05, 0.12};
    
    // Create TH1D with the binning
    TH1D* hDeltaR = new TH1D("hDeltaREmma", 
                             "Emma's distribution;#DeltaR;#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}", 
                             nbins, bins);
    
    // Set bin values and statistical errors
    double values[nbins] = {16.6786, 14.4632, 13.0093, 3.49983};
    double stat_errors[nbins] = {2.27639, 2.03747, 2.33001, 0.703231};
    double syst_errors[nbins] = {1.07718, 0.365236, 0.455006, 0.169006};
    
    for (int i = 0; i < nbins; i++) {
        hDeltaR->SetBinContent(i+1, values[i]);
        hDeltaR->SetBinError(i+1, stat_errors[i]);
    }
    
    // Style the main histogram
    hDeltaR->SetLineColor(kBlack);
    hDeltaR->SetMarkerColor(kBlack);
    hDeltaR->SetMarkerStyle(20);
    hDeltaR->SetMarkerSize(1.0);
    hDeltaR->SetFillStyle(0);
    
    // Create systematic uncertainty boxes
    TBox* systBoxes[nbins];
    
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
    cEmmaRawPlot = new TCanvas("cEmmaRawPlot", "DeltaR Distribution with Systematic Bands", 800, 600);
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

    return hDeltaR;
}


// Get my unfolded data DeltaR distribution
TH1D* getMyDataDistribution(TFile* fMyData) {
    TH2D* hpTjet_vs_DeltaR = (TH2D*)fMyData->Get("Unfolded/hUnfoldedKinCor_iter8");
    if (!hpTjet_vs_DeltaR) {
        throw std::runtime_error("Could not find TH2D 'Unfolded/hUnfoldedKinCor_iter8' in the provided file.");
    }
    // pT,D⁰  #in [5;20] GeV/c
    
    // pT,jet #in [10;20] GeV/c
    int ptjetBinMin = hpTjet_vs_DeltaR->GetXaxis()->FindBin(10.0);
    int ptjetBinMax = hpTjet_vs_DeltaR->GetXaxis()->FindBin(20.0) - 1; // inclusive upper edge

    // Project onto DeltaR axis (Y-axis) - only specify X bin range
    TH1D* hDeltaR = hpTjet_vs_DeltaR->ProjectionY("hMyDeltaRDistribution", ptjetBinMin, ptjetBinMax);

    // ΔR #in [0;0.12] - specify Y-axis bin range too!
    int deltaRBinMin = hpTjet_vs_DeltaR->GetYaxis()->FindBin(0.0);
    int deltaRBinMax = hpTjet_vs_DeltaR->GetYaxis()->FindBin(0.12) - 1;

    // Get DeltaR bin edges from detector level binning
    std::vector<double> deltaRBinEdges_detector = LoadBinning(fMyData, "axes/deltaRBinEdges_detector");

    // Create a new histogram with only the desired DeltaR range
    TH1D* hDeltaR_final = new TH1D("hMyDeltaRDistribution_final", 
                                  "My Data Distribution;#DeltaR;#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}",
                                  deltaRBinEdges_detector.size() - 1,deltaRBinEdges_detector.data());
    
    // Copy the content and errors
    for (int i = 1; i <= hDeltaR_final->GetNbinsX(); i++) {
        int origBin = deltaRBinMin + i - 1;
        hDeltaR_final->SetBinContent(i, hDeltaR->GetBinContent(origBin));
        hDeltaR_final->SetBinError(i, hDeltaR->GetBinError(origBin));
    }
    
    // Clean up the temporary histogram
    delete hDeltaR;

    // Style the histogram
    hDeltaR_final->SetLineColor(kBlue);
    hDeltaR_final->SetMarkerColor(kBlue);
    hDeltaR_final->SetMarkerStyle(21);
    hDeltaR_final->SetMarkerSize(1.0);
    hDeltaR_final->SetFillStyle(0);

    // NORMALIZE the distribution
    double integral = hDeltaR_final->Integral(); // Integral with bin width
    // count manually the number of jets for each bin
    int numberOfJets = 0;
    for (int iBin = 1; iBin <= hDeltaR_final->GetNbinsX(); iBin++) {
        numberOfJets += hDeltaR_final->GetBinContent(iBin);
    }
    std::cout << "Number of jets in my data distribution: " << numberOfJets << std::endl;
    std::cout << "Integral (with bin width) before normalization: " << integral << std::endl;

    if (integral > 0) {
        hDeltaR_final->Scale(1.0 / integral, "width");
    }
    
    return hDeltaR_final;
}

void compareEmmasAndMine(TFile* fMyData, const double& jetptMin, const double& jetptMax) {
    TCanvas* cEmmaRawPlot = nullptr;
    // Create Emma's ΔR-STD distribution histogram
    TH1D* hEmmasDeltaR = createEmmasDistribution(cEmmaRawPlot);

    // Fetch my unfolded ΔR distribution from data
    TH1D* hDataDeltaR = getMyDataDistribution(fMyData);

    TCanvas* cCompare = new TCanvas("cCompare", "Comparison of #DeltaR Distributions", 1800, 1000);
    cCompare->cd();
    hEmmasDeltaR->SetTitle("Comparison of #DeltaR_{STD} Distributions;#DeltaR_{STD};#frac{1}{N_{jets}}#frac{dN}{d#DeltaR_{STD}}");
    hEmmasDeltaR->Draw();
    hDataDeltaR->Draw("same");

    TLegend* legendCompare = new TLegend(0.6, 0.7, 0.85, 0.85);
    legendCompare->SetBorderSize(0);
    legendCompare->SetFillStyle(0);
    legendCompare->AddEntry(hEmmasDeltaR, "(Emma Yeats') Run 2, pp #rightarrow #sqrt{s} = 5.02 TeV", "lep");
    legendCompare->AddEntry(hDataDeltaR, "Run 3, pp #rightarrow #sqrt{s} = 13 TeV", "lep");
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
    TString imagePath = "../Images/5-ClosureTest/Fourth/";
    cEmmaRawPlot->Print(imagePath + Form("closureTest4_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cCompare->Print(imagePath + Form("closureTest4_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
    

}
// Compare data final distribution to Emma Yeats' reported distribution
void FourthClosureTest(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution


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
    TFile* fClosureInputMatched = new TFile("mc_closure_input_matched_data.root","read");
    TFile* fClosureInputNonMatched = new TFile("mc_closure_input_non-matched_data.root","read");
    TFile* fPowheg = new TFile("../SimulatedData/POWHEG/trees_powheg_fd_central.root","read");
    TFile* fMyData = new TFile("../4-Unfolding/unfolding_5_to_30_jetpt.root","read");
    if (!fClosureInputMatched || fClosureInputMatched->IsZombie()) {
        std::cerr << "Error: Unable to open O2 MC matched ROOT file." << std::endl;
    }
    if (!fClosureInputNonMatched || fClosureInputNonMatched->IsZombie()) {
        std::cerr << "Error: Unable to open AnalysisResults.root data ROOT file." << std::endl;
    }
    if (!fPowheg || fPowheg->IsZombie()) {
        std::cerr << "Error: Unable to open POWHEG data ROOT file." << std::endl;
    }
    if (!fMyData || fMyData->IsZombie()) {
        std::cerr << "Error: Unable to open unfolded data ROOT file." << std::endl;
    }

    // Compare Emma Yeats' reported DeltaR distribution to my unfolded data distribution
    compareEmmasAndMine(fMyData, jetptMin, jetptMax);

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    FourthClosureTest();
    return 0;
}
