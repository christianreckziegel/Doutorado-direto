/**
 * @file create3DBackgroundSubtracted.C
 * @author Christian Reckziegel
 * @brief Macro for merging DeltaR distributions from all pT,jet and pT,D bins in one single TH3 histogram
**/

using namespace std;


void create3DBackgroundSubtracted(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // jet pT cuts
    std::vector<double> ptjetBinEdges = {5., 7., 15., 30., 50.};
    double jetptMin = ptjetBinEdges[0]; // GeV
    double jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 1]; // GeV
    // deltaR histogram
    int deltaRbins = 10000; // deltaRbins = numberOfPoints, default=10 bins for [0. 0.4]
    std::vector<double> deltaRBinEdges = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15,0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5}; // chosen by Nima
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
    // pT,D bins
    std::vector<double> ptDBinEdges = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.};

    TH3D* h3D = new TH3D("h3DBackgroundSubtracted", "Background subtracted; p_{T,jet} (GeV/c); #DeltaR; p_{T,D^{0}} (GeV/c)",
                     ptjetBinEdges.size()-1, ptjetBinEdges.data(),
                     deltaRBinEdges.size()-1, deltaRBinEdges.data(),
                     ptDBinEdges.size()-1, ptDBinEdges.data());
    h3D->Sumw2();

    for (size_t iJetBin = 0; iJetBin < ptjetBinEdges.size() - 1; iJetBin++) {
        TFile* fJetRange = new TFile(Form("backSub_%0.f_to_%0.f_jetpt_with_reflections.root",ptjetBinEdges[iJetBin],ptjetBinEdges[iJetBin+1]),"read");
        if (!fJetRange || fJetRange->IsZombie()) {
            std::cerr << "Error opening file " << Form("backSub_%0.f_to_%0.f_jetpt_with_reflections.root",ptjetBinEdges[iJetBin],ptjetBinEdges[iJetBin+1]) << std::endl;
            continue;
        }

        for (int iHist = 0; iHist < 100; ++iHist) {
            TString histName = Form("h_back_subtracted_%d", iHist);
            TH1D* hDeltaR = (TH1D*)fJetRange->Get(histName);
            if (!hDeltaR) break; // No more histograms
            //std::cout << "Histogram x axis range: " << hDeltaR->GetXaxis()->GetXmin() << " to " << hDeltaR->GetXaxis()->GetXmax() << std::endl;

            // Get bin centers
            double ptJetCenter = 0.5 * (ptjetBinEdges[iJetBin] + ptjetBinEdges[iJetBin+1]);
            TString title = hDeltaR->GetTitle();

            // If the title is missing, try to get it from an older cycle
            if (title.IsNull() || title.IsWhitespace()) {
                TString fallbackName = histName + ";1";
                TH1D* hFallback = (TH1D*)fJetRange->Get(fallbackName);
                if (hFallback) {
                    std::cout << "Using title from fallback for " << histName << std::endl;
                    title = hFallback->GetTitle(); // use title only, keep content from cycle ;2
                } else {
                    std::cerr << "No valid title found for " << histName << std::endl;
                    continue;
                }
            }

            double ptDlow = -1, ptDhigh = -1;
            if (sscanf(title.Data(), "%lf < #it{p}_{T, D^{0}} < %lf GeV/#it{c}", &ptDlow, &ptDhigh) != 2) {
                std::cerr << "Could not parse pT,D range from histogram title: " << title << std::endl;
                std::cout << "Trying to parse title: '" << title << "'" << std::endl;
                std::cout << "hDeltaR histogram -> " << hDeltaR->GetName() << std::endl;
                std::cout << "In file " << fJetRange->GetName() << std::endl;
                continue;
            }
            double ptDCenter = 0.5 * (ptDlow + ptDhigh);


            for (int iBin = 1; iBin <= hDeltaR->GetNbinsX(); ++iBin) {
                double deltaRcenter = hDeltaR->GetBinCenter(iBin);
                double content = hDeltaR->GetBinContent(iBin);
                double error = hDeltaR->GetBinError(iBin);

                // Fill the TH3D using centers
                h3D->Fill(ptJetCenter, deltaRcenter, ptDCenter, content);
                // Optionally, if you want to preserve error propagation:
                int binX = h3D->GetXaxis()->FindBin(ptJetCenter);
                int binY = h3D->GetYaxis()->FindBin(deltaRcenter);
                int binZ = h3D->GetZaxis()->FindBin(ptDCenter);
                h3D->SetBinError(binX, binY, binZ, error); // only if needed
            }
        }

    }
    
    h3D->Draw("colz");

    TFile* fOut = new TFile(Form("full_merged_ranges_back_sub.root"), "RECREATE");
    h3D->Write();
    //fOut->Close();
    std::cout << "3D histogram created and saved to " << fOut->GetName() << " with 3D histogram." << std::endl;
    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    create3DBackgroundSubtracted();
    return 0;
}



