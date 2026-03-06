#include "../FastJet/commonFunctions.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <iomanip> // time precision

void plotDeltaR() {
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    int numberOfEvents = 10; // k
    
    // Rapidity cuts
    double R = 0.4;
    double etaCut = 0.9 - R; // ALICE central barrel acceptance in pseudorapidity
    double yCut = 0.8; // on detector level D0, to be consistent with ALICE's analysis cuts

    // Open the ROOT file containing the reconstructed jets
    TFile* inputFile5 = TFile::Open(Form("../FastJet/fastjetCharmEvents_5TeV_%dk.root", numberOfEvents));
    if (!inputFile5 || inputFile5->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }
    TFile* inputFile13 = TFile::Open(Form("../FastJet/fastjetCharmEvents_13TeV_%dk.root", numberOfEvents));
    if (!inputFile13 || inputFile13->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }

    // pT,jet cuts
    std::vector<double> ptjetBinEdges;
    // pT,D bins
    std::vector<double> ptDBinEdges;
    // DeltaR bin edges (asymmetric binning)
    std::vector<double> deltaRBinEdges;
    bool useEmmasCuts = true;
    if (useEmmasCuts) {
        // pT,jet cuts
        ptjetBinEdges = {5., 10., 15., 20., 30.};
        // pT,D bins
        ptDBinEdges = {5., 6., 7., 8., 9., 10., 12., 20.};
        // DeltaR bin edges (asymmetric binning)
        deltaRBinEdges = {0., 0.01, 0.03, 0.05, 0.12}; // Emma binning = {0., 0.01, 0.03, 0.05, 0.12, 0.2}
    } else {
        // pT,jet cuts
        ptjetBinEdges = {5., 7., 10., 15., 30., 50.};
        // pT,D bins
        ptDBinEdges = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.};
        // DeltaR bin edges (asymmetric binning)
        deltaRBinEdges = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5};
    }
    
    double hfPtMin = ptDBinEdges[0];
    double hfPtMax = ptDBinEdges[ptDBinEdges.size() - 1];
    double jetptMin = ptjetBinEdges[1];
    double jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 2];
    
    TH1D* hDeltaR5 = buildDeltaRHistogram(inputFile5, "hDeltaR5", Form("Distance between D^{0} and jet axis (%dk events)", numberOfEvents), deltaRBinEdges, jetptMin, jetptMax, hfPtMin, hfPtMax, etaCut, yCut);
    TH1D* hDeltaR13 = buildDeltaRHistogram(inputFile13, "hDeltaR13", Form("Distance between D^{0} and jet axis (%dk events)", numberOfEvents), deltaRBinEdges, jetptMin, jetptMax, hfPtMin, hfPtMax, etaCut, yCut);
    if (!hDeltaR5 || !hDeltaR13) {
        std::cout << "Error: Failed to create histograms!" << std::endl;
        return;
    }

    TCanvas* cDeltaR = new TCanvas("cDeltaR", "Delta R between D0 and jet axis", 800, 600);
    hDeltaR5->SetLineColor(kBlack);
    hDeltaR5->SetMarkerColor(kBlack);
    hDeltaR5->SetMarkerStyle(20);
    hDeltaR13->SetLineColor(kBlue);
    hDeltaR13->SetMarkerColor(kBlue);
    hDeltaR13->SetMarkerStyle(21);
    if (hDeltaR5->GetMaximum() > hDeltaR13->GetMaximum()) {
        hDeltaR5->GetYaxis()->SetRangeUser(0, hDeltaR5->GetMaximum() * 1.5);
        hDeltaR5->Draw();
        hDeltaR13->Draw("SAME");
    } else {
        hDeltaR13->GetYaxis()->SetRangeUser(0, hDeltaR13->GetMaximum() * 1.5);
        hDeltaR13->Draw();
        hDeltaR5->Draw("SAME");
    }

    TLegend* legend = new TLegend(0.7,0.8,0.88,0.88);
    legend->AddEntry(hDeltaR5, "Pythia, pp #rightarrow #sqrt{s} = 5.02 TeV", "lp");
    legend->AddEntry(hDeltaR13, "Pythia, pp #rightarrow #sqrt{s} = 13.6 TeV", "lp");
    legend->Draw();

    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.03);
    latex->DrawLatex(0.15, 0.85, Form("Jet p_{T} in [%.0f, %.0f] GeV/c, |#eta_{jet}| #leq %.1f", jetptMin, jetptMax, etaCut));
    latex->DrawLatex(0.15, 0.80, Form("D^{0} p_{T} in [%.0f, %.0f] GeV/c, |y_{D^{0}}| #leq %.1f", hfPtMin, hfPtMax, yCut));

    // Store data and figure in output file
    TFile* outputFile = new TFile(Form("DeltaR_comparison_%dk.root", numberOfEvents), "RECREATE");
    outputFile->cd();
    hDeltaR5->Write();
    hDeltaR13->Write();
    cDeltaR->Write();
    cDeltaR->Print(Form("DeltaR_comparison_%dk.pdf", numberOfEvents));

    // Clean up - IMPORTANT: Close files but don't delete objects
    outputFile->Close();
    inputFile5->Close();
    inputFile13->Close();

    // Delete canvas (optional, but good practice)
    delete cDeltaR;

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed 
         << time_taken/60 << std::setprecision(5); 
    std::cout << " min " << std::endl;
}

int main() {
    plotDeltaR();
    return 0;
}