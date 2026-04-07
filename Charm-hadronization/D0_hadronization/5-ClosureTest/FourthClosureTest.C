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

std::pair<TH1D*, std::vector<TBox*>> createEmmasDistribution(TCanvas*& cEmmaRawPlot) {
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
    hDeltaR->SetLineColor(kRed);
    hDeltaR->SetMarkerColor(kRed);
    hDeltaR->SetMarkerStyle(20);
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
    TH2D* hpTjet_vs_DeltaR = (TH2D*)fMyData->Get("Unfolded/hUnfoldedKinCor_iter8");
    if (!hpTjet_vs_DeltaR) {
        throw std::runtime_error("Could not find TH2D 'Unfolded/hUnfoldedKinCor_iter8' in the provided file.");
    }
    // pT,D⁰ #in [5;20] GeV/c -> this is ensured since the beginning of the analysis on background subtraction
    
    // pT,jet #in [10;20] GeV/c
    int ptjetBinMin = 1 + 1; // 10 GeV/c on 2nd bin
    int ptjetBinMax = hpTjet_vs_DeltaR->GetNbinsX() - 1; // inclusive upper edge, 20 GeV/c is the upper edge of the last bin, so we take the previous bin as max
    std::cout << "(ptjetBinMin,ptjetBinMax) = (" << ptjetBinMin << "," << ptjetBinMax << ")" << std::endl;
    std::cout << "hpTjet_vs_DeltaR->GetNbinsY() = " << hpTjet_vs_DeltaR->GetNbinsY() << std::endl;

    // loop through x bins to find the one corresponding to 10 GeV/c and 20 GeV/c, and print the bin edges for verification
    for (size_t iXbin = 0; iXbin < hpTjet_vs_DeltaR->GetNbinsX(); iXbin++) {
        double binLowEdge = hpTjet_vs_DeltaR->GetXaxis()->GetBinLowEdge(iXbin + 1);
        double binUpEdge = hpTjet_vs_DeltaR->GetXaxis()->GetBinUpEdge(iXbin + 1);
        std::cout << "X bin " << iXbin + 1 << "\t[binLowEdge,binUpEdge]: [" << binLowEdge << ";" << binUpEdge << "] GeV/c" << std::endl;
        if (binLowEdge >= 10 && binUpEdge <= 20) {
            std::cout << "Selected pT,jet bin " << iXbin + 1 << ": [" << binLowEdge << ";" << binUpEdge << "] GeV/c" << std::endl;
        }
    }
    
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
        hDeltaR_final = new TH1D("hMyDeltaRDistribution_final_v1", 
                                    "My Data Distribution;#DeltaR;#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}",
                                    deltaRBinEdges_truncated.size() - 1,deltaRBinEdges_truncated.data());
        
        // Copy the content and errors
        for (int i = 1; i <= hDeltaR_final->GetNbinsX(); i++) {
            int origBin = deltaRBinMin + i - 1;
            hDeltaR_final->SetBinContent(i, hDeltaR->GetBinContent(origBin));
            hDeltaR_final->SetBinError(i, hDeltaR->GetBinError(origBin));
        }
    } else {
        // If not truncating, just clone the original histogram and rename it
        hDeltaR_final = (TH1D*)hDeltaR->Clone("hMyDeltaRDistribution_final");
        hDeltaR_final->SetTitle("My Data Distribution;#DeltaR;#frac{1}{N_{jets}}#frac{dN}{d#DeltaR}");
    }
    

    // Clean up the temporary histogram
    delete hDeltaR;

    // Style the histogram
    hDeltaR_final->SetLineColor(kBlue);
    hDeltaR_final->SetMarkerColor(kBlue);
    hDeltaR_final->SetMarkerStyle(21);
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

std::pair<std::vector<TH1D*>, std::vector<TBox*>> compareEmmasAndMine(TFile* fMyData, const double& jetptMin, const double& jetptMax, const BinningStruct& binning) {
    TCanvas* cEmmaRawPlot = nullptr;
    // Create Emma's ΔR-STD distribution histogram
    std::pair<TH1D*, std::vector<TBox*>> emmaDistribution = createEmmasDistribution(cEmmaRawPlot);
    TH1D* hEmmasDeltaR = emmaDistribution.first;
    std::vector<TBox*> systBoxes = emmaDistribution.second;

    // Fetch my unfolded ΔR distribution from data
    TH1D* hDataDeltaR = getMyDataDistribution(fMyData, binning);

    TCanvas* cCompare = new TCanvas("cCompare", "Comparison of #DeltaR Distributions", 1800, 1000);
    cCompare->cd();
    hEmmasDeltaR->SetTitle("Comparison of #DeltaR_{STD} Distributions;#DeltaR_{STD};#frac{1}{N_{jets}}#frac{dN}{d#DeltaR_{STD}}");
    hDataDeltaR->GetYaxis()->SetRangeUser(0, 20); // Adjust range as needed
    hDataDeltaR->Draw();
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
    std::vector<TH1D*> histogramsToSave = {hEmmasDeltaR, hDataDeltaR};

    return std::make_pair(histogramsToSave, emmaDistribution.second);
    

}
// Compare data final distribution to Emma Yeats' reported distribution
void FourthClosureTest(){
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

    // Opening files
    TFile* fMyData = new TFile(Form("../4-Unfolding/unfolding_%.0f_to_%.0f_jetpt.root", jetptMin, jetptMax), "read");
    if (!fMyData || fMyData->IsZombie()) {
        std::cerr << "Error: Unable to open unfolded data ROOT file." << std::endl;
    }

    // Compare Emma Yeats' reported DeltaR distribution to my unfolded data distribution
    std::pair<std::vector<TH1D*>, std::vector<TBox*>> emmaAndMine = compareEmmasAndMine(fMyData, jetptMin, jetptMax, binning);
    TH1D* hEmmasDeltaR = emmaAndMine.first[0];
    TH1D* hMyDeltaR = emmaAndMine.first[1];
    std::vector<TBox*> systematicBands = emmaAndMine.second;

    TFile* fOutput = new TFile(Form("FourthClosureTestResults_%.0f_to_%.0fGeV_jetpt.root",jetptMin,jetptMax), "RECREATE");
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

    
}

int main(){
    FourthClosureTest();
    return 0;
}
