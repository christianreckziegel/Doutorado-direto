/**
 * @file ReflectionsTreatment.C
 * @author Christian Reckziegel
 * @brief Macro for obtaining the fit parameters, integrals and reflection scaling parameter for each invariant mass histogram fit
**/

using namespace std;

//--- bit manipulation ---------------------------------------------------------
#define BIT(n)       (1ULL << (n))
#define SETBIT(n,i)  ((n) |= BIT(i))
#define CLRBIT(n,i)  ((n) &= ~BIT(i))
#define TESTBIT(n,i) ((Bool_t)(((n) & BIT(i)) != 0))
//-----------------------------------------------------------------------------Structs---------------------------------------------------------------------------------------------------------
// Histograms containers for each case (each container has a number of histograms corresponding to the invariant mass intervals)
struct HistogramGroup {
    std::vector<TH2D*> signals;
    std::vector<TH2D*> reflections;
    std::vector<TH2D*> signals_and_reflections;

    std::vector<TH1D*> signals_1d;
    std::vector<TH1D*> reflections_1d;
    std::vector<TH1D*> signals_and_reflections_1d;
};
struct FitsGroup {
    std::vector<TF1*> signals;
    std::vector<TF1*> reflections;
    std::vector<TF1*> signals_and_reflections;
};
struct AnalysisData {

    // Histograms
    HistogramGroup histograms;

    // Fits
    FitsGroup fits;
};
//-----------------------------------------------------------------------------Support functions---------------------------------------------------------------------------------------------------------
double DeltaPhi(double phi1, double phi2) {
    // Compute the absolute difference between phi1 and phi2
    double dphi = std::abs(phi1 - phi2); 
    if (dphi > M_PI) {
        // subtract 2pi if the difference if bigger than pi
        dphi = dphi - 2*M_PI;
    }

    return dphi;
}
//-----------------------------------------------------------------------------Fit functions---------------------------------------------------------------------------------------------------------
// Custom pure signal fit function
Double_t pureSignalFunction(Double_t* x, Double_t* par) {

    Double_t m = x[0];
    Double_t C_1 = par[0];
    Double_t C_2 = par[1];
    Double_t m0 = par[2]; // constrain: same mean value parameter for both single gaussians
    Double_t sigma_1 = par[3];
    Double_t sigma_2 = par[4];

    // Defining the custom function
    Double_t result = C_1 * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma_1), 2)) + C_2 * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma_2), 2));
    return result;
}
// Custom pure reflections fit function
Double_t pureReclectionsFunction(Double_t* x, Double_t* par) {

    Double_t m = x[0];
    Double_t C_1 = par[0];
    Double_t C_2 = par[1];
    Double_t m0_1 = par[2];
    Double_t m0_2 = par[3];
    Double_t sigma_1 = par[4];
    Double_t sigma_2 = par[5];

    // Defining the custom function
    Double_t result = C_1 * TMath::Exp(-TMath::Power((m - m0_1) / (2 * sigma_1), 2)) + C_2 * TMath::Exp(-TMath::Power((m - m0_2) / (2 * sigma_2), 2));
    return result;
}
// Custom signal and reflections fit function
Double_t signalAndReflectionsFunction(Double_t* x, Double_t* par) {
    
    return pureSignalFunction(x,par) + pureReclectionsFunction(x,&par[5]);
}
//-----------------------------------------------------------------------------Main workflow functions---------------------------------------------------------------------------------------------------------
// Module to create TH2D histograms including interest variable
HistogramGroup createHistograms(const std::vector<double>& ptDBinEdges, int xbins, double xmin, double xmax, const std::vector<double>& yBinEdges) {

    // create histograms containers for each case (each container has a number of histograms corresponding to the invariant mass intervals)
    HistogramGroup histograms;

    for (size_t i = 0; i < ptDBinEdges.size() - 1; ++i) {
        // pure reflections histograms
        histograms.reflections.push_back(new TH2D(Form("histMass_reflections%zu",i+1), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        histograms.reflections[i]->Sumw2();

        // pure signals histograms
        histograms.signals.push_back(new TH2D(Form("histMass_signals%zu",i+1), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        histograms.signals[i]->Sumw2();
        
        // signals+reflections histograms: calculate ratios
        histograms.signals_and_reflections.push_back(new TH2D(Form("histMass_signals_and_reflections%zu",i+1), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        histograms.signals_and_reflections[i]->Sumw2();
        
    }
    return histograms;
}

// Module to fill 2D histograms from TTree data
void fillHistograms(TFile* fInputMC, HistogramGroup& histograms, const double& jetptMin, const double& jetptMax, std::vector<double>& ptDBinEdges) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = 0.4; // on particle level delta R

    // Accessing detector level data TTree
    TTree* tree = (TTree*)fInputMC->Get("DF_combined/O2mcdjetdisttable");

    // Assuming histograms and tree data correspond in some way
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing the TTree
    float axisDistance, jetPt, jetEta, jetPhi, jetNConst;
    float hfPt, hfEta, hfPhi, hfMass, hfY, hfMlScore0;
    int hfMatchedFrom, hfSelectedAs;

    tree->SetBranchAddress("fJetHfDist",&axisDistance);
    tree->SetBranchAddress("fJetPt",&jetPt);
    tree->SetBranchAddress("fJetEta",&jetEta);
    tree->SetBranchAddress("fJetPhi",&jetPhi);
    tree->SetBranchAddress("fJetNConst",&jetNConst);
    tree->SetBranchAddress("fHfPt",&hfPt);
    tree->SetBranchAddress("fHfEta",&hfEta);
    tree->SetBranchAddress("fHfPhi",&hfPhi);
    tree->SetBranchAddress("fHfMass",&hfMass);
    tree->SetBranchAddress("fHfY",&hfY);
    tree->SetBranchAddress("fHfMlScore0",&hfMlScore0);
    tree->SetBranchAddress("fHfMatchedFrom",&hfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&hfSelectedAs);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // correct selection BIT shift info: Checks whether BIT(i) is set, regardless of other bits
        if (hfSelectedAs & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as D0
            hfSelectedAs = 1;
        } else if (hfSelectedAs & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as D0bar
            hfSelectedAs = -1;
        }
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));
        
        // Fill each histogram with their respective pT intervals
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= 0.) && (deltaR < DeltaRcut))) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < ptDBinEdges.size() - 1 && !filled; iEdge++) {
                if ((hfPt >= ptDBinEdges[iEdge]) && (hfPt < ptDBinEdges[iEdge + 1])) {

                    // D0 = +1, D0bar = -1, neither = 0
                    if ((hfMatchedFrom != 0) && (hfSelectedAs != 0)) {
                        if (hfMatchedFrom == hfSelectedAs) {
                            // pure signals
                            histograms.signals[iEdge]->Fill(hfMass, deltaR); // 2D case: (hfMass, deltaR);
                        } else {
                            // pure reflections
                            histograms.reflections[iEdge]->Fill(hfMass, deltaR);
                        }
                        
                        // signals and reflections altogether, withOUT "neither = 0" entries
                        histograms.signals_and_reflections[iEdge]->Fill(hfMass, deltaR);
                    }
                    // signals and reflections altogether, wiTH "neither = 0" entries (should I include "neither" entries too? Altough they should not appear)
                    //histograms.signals_and_reflections[iEdge]->Fill(hfMass, deltaR);
                    filled = true; // Exit the loop once the correct histogram is found
                }
                
            }
            
        }

        
        
        
    }
    cout << "Histograms filled.\n";
}

FitsGroup performFits(HistogramGroup& histograms, std::vector<double>& ptDBinEdges, const double& xmin, const double& xmax) {

    FitsGroup fits;

    if ((histograms.signals.size() != histograms.reflections.size()) || (histograms.signals.size() != histograms.signals_and_reflections.size()) || (histograms.reflections.size() != histograms.signals_and_reflections.size())) {
        std::cout << "Warning: different size of histogram vectors (signals, reflections, signal_and_reflections)" << std::endl;
    }
    
    // Helper lambda for projection and fitting
    /*auto performFit = [&](std::vector<TH2D*>& inputHistograms, 
                          std::vector<TH1D*>& outputHistograms, 
                          std::vector<TF1*>& fitResults, 
                          const char* fitName, 
                          TF1* fitFunction) {
        for (size_t iInterval = 0; iInterval < ptDBinEdges.size() - 1; iInterval++) {

            
            TH1D* projected = inputHistograms[iInterval]->ProjectionX(Form("%s_1d_%zu", fitName, iInterval + 1), 1, -1);

            // Initialize parameters from the histogram distribution characteristics
            A1 = 0.7 * projected->GetMaximum(); // Assume first Gaussian dominates
            A2 = 0.3 * projected->GetMaximum(); // Second Gaussian smaller
            mean = projected->GetMean();       // Mean of the distribution
            sigma1 = projected->GetRMS() / 2; // First Gaussian narrower
            sigma2 = projected->GetRMS();    // Second Gaussian wider

            outputHistograms.push_back(projected);
            projected->Fit(fitFunction, "RQ");
            fitResults.push_back(fitFunction);
        }
    };

    // Perform fits for signals
    performFit(histograms.signals, histograms.signals_1d, fits.signals, 
               "signal", new TF1("signalFit", pureSignalFunction, xmin, xmax, 5));

    // Perform fits for reflections
    performFit(histograms.reflections, histograms.reflections_1d, fits.reflections, 
               "reflections", new TF1("reflectionsFit", pureReclectionsFunction, xmin, xmax, 6));

    // Perform combined fits
    performFit(histograms.signals_and_reflections, histograms.signals_and_reflections_1d, fits.signals_and_reflections, 
               "combined", new TF1("combinedFit", signalAndReflectionsFunction, xmin, xmax, 11));*/

    // Signal fits
    for (size_t iInterval = 0; iInterval < ptDBinEdges.size() - 1; iInterval++) {
        TF1* signalFit = new TF1(Form("signalFit_%zu", iInterval), pureSignalFunction, xmin, xmax, 5);

        // Obtain 1D mass histogram projection
        histograms.signals_1d.push_back( histograms.signals[iInterval]->ProjectionX(Form("histMass_signals_1d_%zu", iInterval+1), 1, -1) );

        // Initialize parameters from the histogram distribution characteristics
        double C1, C2, m0, sigma1, sigma2;
        C1 = 0.7 * histograms.signals_1d[iInterval]->GetMaximum(); // Assume first Gaussian dominates
        C2 = 0.3 * histograms.signals_1d[iInterval]->GetMaximum(); // Second Gaussian smaller
        m0 = histograms.signals_1d[iInterval]->GetMean();       // Mean of the distribution
        sigma1 = histograms.signals_1d[iInterval]->GetRMS() / 2; // First Gaussian narrower
        sigma2 = histograms.signals_1d[iInterval]->GetRMS();    // Second Gaussian wider
        signalFit->SetParameters(C1, C2, m0, sigma1, sigma2);
        signalFit->SetParName(0, "C1");
        signalFit->SetParName(1, "C2");
        signalFit->SetParName(2, "m0");
        signalFit->SetParName(3, "sigma1");
        signalFit->SetParName(4, "sigma2");
        //signalFit->SetParLimits(3, 0.5 * sigma1, 2.0 * sigma1); // Example constraint for sigma1

        histograms.signals_1d[iInterval]->Fit(signalFit, "RQ"); // "Q" option performs quiet fit without drawing the fit function
        fits.signals.push_back(signalFit);
    }

    // Reflections fits
    for (size_t iInterval = 0; iInterval < ptDBinEdges.size() - 1; iInterval++) {
        TF1* reflectionsFit = new TF1(Form("reflectionsFit_%zu", iInterval), pureReclectionsFunction, xmin, xmax, 6);

        // Obtain 1D mass histogram projection
        histograms.reflections_1d.push_back( histograms.reflections[iInterval]->ProjectionX(Form("histMass_reflections_1d_%zu", iInterval+1), 1, -1) );

        // Initialize parameters from the histogram distribution characteristics
        double C1, C2, m0_1, m0_2, sigma1, sigma2;
        C1 = 0.6 * histograms.reflections_1d[iInterval]->GetMaximum();   // Dominant reflection
        C2 = 0.4 * histograms.reflections_1d[iInterval]->GetMaximum();   // Secondary reflection
        m0_1 = histograms.reflections_1d[iInterval]->GetMean() - histograms.reflections_1d[iInterval]->GetRMS() / 2; // Offset from mean
        m0_2 = histograms.reflections_1d[iInterval]->GetMean() + histograms.reflections_1d[iInterval]->GetRMS() / 2; // Offset from mean
        sigma1 = histograms.reflections_1d[iInterval]->GetRMS() / 2;  // Narrower
        sigma2 = histograms.reflections_1d[iInterval]->GetRMS();      // Wider
        reflectionsFit->SetParameters(C1, C2, m0_1, m0_2, sigma1, sigma2);

        histograms.reflections_1d[iInterval]->Fit(reflectionsFit, "RQ");
        fits.reflections.push_back(reflectionsFit);
    }

    // Combined fits: true signal + reflections
    for (size_t iInterval = 0; iInterval < ptDBinEdges.size() - 1; iInterval++) {
        TF1* combinedFit = new TF1(Form("combinedFit_%zu", iInterval), signalAndReflectionsFunction, xmin, xmax, 11);

        // Obtain 1D mass histogram projection
        histograms.signals_and_reflections_1d.push_back( histograms.signals_and_reflections[iInterval]->ProjectionX(Form("histMass_signals_and_reflections_1d_%zu", iInterval+1), 1, -1) );

        // Initialize parameters from previous fits
        combinedFit->SetParameters(fits.signals[iInterval]->GetParameter(0), fits.signals[iInterval]->GetParameter(1), fits.signals[iInterval]->GetParameter(2), fits.signals[iInterval]->GetParameter(3), fits.signals[iInterval]->GetParameter(4), 
                                      fits.reflections[iInterval]->GetParameter(0), fits.reflections[iInterval]->GetParameter(1), fits.reflections[iInterval]->GetParameter(2), fits.reflections[iInterval]->GetParameter(3), fits.reflections[iInterval]->GetParameter(4), fits.reflections[iInterval]->GetParameter(5));

        histograms.signals_and_reflections_1d[iInterval]->Fit(combinedFit, "RQ");
        fits.signals_and_reflections.push_back(combinedFit);
    }
    

    return fits;
};

void PlotHistograms(const HistogramGroup& histograms, FitsGroup& fits, const double jetptMin, const double jetptMax) {
    
    // creating 1D mass projection histograms
    TH1D* tempHist;
    //std::vector<TH1D*> histograms_1d;

    // Defining canvases before plotting
    TCanvas* cSignals = new TCanvas("cSignals","Pure signal histograms");
    //cSignals->SetCanvasSize(1800,1000);
    cSignals->Divide(3,static_cast<int>(histograms.signals.size() / 3)); // columns, lines
    TCanvas* cReflections = new TCanvas("cReflections","Pure reflections histograms");
    //cReflections->SetCanvasSize(1800,1000);
    cReflections->Divide(3,static_cast<int>(histograms.signals.size() / 3)); // columns, lines
    TCanvas* cSignalsAndReflections = new TCanvas("cSignalsAndReflections","Signal and reflections histograms");
    //cSignalsAndReflections->SetCanvasSize(1800,1000);
    cSignalsAndReflections->Divide(3,static_cast<int>(histograms.signals.size() / 3)); // columns, lines

    // Creating legend for the three components in the final fit
    /*TLegend* legend = new TLegend(0.6,0.57,0.9,0.77);
    legend->AddEntry(fits.signals[0],"Signal", "l");
    legend->AddEntry(fits.reflections[0],"Reflections", "l");
    legend->AddEntry(fits.signals_and_reflections[0],"Total fit", "l");*/
    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.05);

    // Loop over pT,D intervals/histograms
    for (size_t iHist = 0; iHist < histograms.signals.size(); iHist++) {

        //
        // Pure signals
        //
        tempHist = histograms.signals[iHist]->ProjectionX(Form("h_mass_signals_proj_%zu", iHist));
        cSignals->cd(iHist+1);
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        tempHist->SetMarkerStyle(kDot); // kFullDotMedium
        tempHist->SetMarkerColor(kBlack); // kBlack
        tempHist->SetLineColor(kBlack); // kGray
        tempHist->GetYaxis()->SetTitle("counts");
        tempHist->Draw();
        fits.signals[iHist]->SetLineColor(8); // pastel green
        fits.signals[iHist]->Draw("same");
        //fittings[iHist]->Draw("same");
        //int backgroundHist = histograms.size() + iHist;
        //fittings[backgroundHist]->Draw("same");
        //double m_0 = fittings[iHist]->GetParameter(3); // Get the value of parameter 'm_0'
        //double sigma = fittings[iHist]->GetParameter(4); // Get the value of parameter 'sigma'
        double chi2 = fits.signals[iHist]->GetChisquare();
        double degOfFreedom = fits.signals[iHist]->GetNDF();
        // original position at statBoxPos-0.35, 0.70 with 0.03 of size
        latex->DrawLatex(statBoxPos-0.4, 0.65, Form("m_{0} = %.3f GeV/c^{2}", fits.signals[iHist]->GetParameter(2))); // Display parameter 'm_0' value
        latex->DrawLatex(statBoxPos-0.87, 0.70,"MC derived data:"); // Derived data from LHC24d3a (MC)
        latex->DrawLatex(statBoxPos-0.87, 0.65,"HF_LHC24d3a_All"); // Derived data from LHC24d3a (MC)
        //latex->DrawLatex(statBoxPos-0.3, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax)); // Display jet pT cut applied
        latex->DrawLatex(statBoxPos-0.4, 0.58, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));

        // Drawing 2D histograms
        //cSignals_2d->cd(iHist+1);
        //histograms2d[iHist]->GetYaxis()->SetTitle("#DeltaR");
        //histograms2d[iHist]->Draw("colz");

        //
        // Pure reflections
        //
        tempHist = histograms.reflections[iHist]->ProjectionX(Form("h_mass_reflections_proj_%zu", iHist));
        cReflections->cd(iHist+1);
        statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        tempHist->SetMarkerStyle(kDot);
        tempHist->SetMarkerColor(kBlack);
        tempHist->SetLineColor(kBlack);
        tempHist->GetYaxis()->SetTitle("counts");
        tempHist->Draw();
        chi2 = fits.reflections[iHist]->GetChisquare();
        degOfFreedom = fits.reflections[iHist]->GetNDF();
        fits.reflections[iHist]->SetLineColor(46); // pastel red
        fits.reflections[iHist]->Draw("same");
        latex->DrawLatex(statBoxPos-0.87, 0.40,"MC derived data:"); // Derived data from LHC24d3a (MC)
        latex->DrawLatex(statBoxPos-0.87, 0.35,"HF_LHC24d3a_All"); // Derived data from LHC24d3a (MC)
        latex->DrawLatex(statBoxPos-0.87, 0.30, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));

        //
        // Signal and reflections
        //
        tempHist = histograms.signals_and_reflections[iHist]->ProjectionX(Form("h_mass_signals_and_reflections_proj_%zu", iHist));
        cSignalsAndReflections->cd(iHist+1);
        //statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        tempHist->SetMarkerStyle(kDot);
        tempHist->SetMarkerColor(kBlack);
        tempHist->SetLineColor(kBlack); // kGray
        tempHist->GetYaxis()->SetTitle("counts");
        tempHist->Draw();
        chi2 = fits.signals_and_reflections[iHist]->GetChisquare();
        degOfFreedom = fits.signals_and_reflections[iHist]->GetNDF();
        fits.signals_and_reflections[iHist]->SetLineColor(38); // pastel blue
        latex->DrawLatex(statBoxPos-0.87, 0.70,"MC derived data:"); // Derived data from LHC24d3a (MC)
        latex->DrawLatex(statBoxPos-0.87, 0.65,"HF_LHC24d3a_All"); // Derived data from LHC24d3a (MC)
        latex->DrawLatex(statBoxPos-0.87, 0.60, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));

        // Creating legend for the three components in the final fit
        TLegend* legend = new TLegend(0.6,0.57,0.85,0.77);
        legend->AddEntry(fits.signals[iHist],"Pure signal", "l");
        legend->AddEntry(fits.reflections[iHist],"Reflections", "l");
        legend->AddEntry(fits.signals_and_reflections[iHist],"Total fit", "l");

        fits.signals[iHist]->Draw("same");
        fits.reflections[iHist]->Draw("same");
        fits.signals_and_reflections[iHist]->Draw("same");
        legend->Draw();

    }

    //
    // Storing images
    //
    TString imagePath = "../../Images/1-SignalTreatment/SideBand/Reflections/";
    TString jetPtRange = Form("%.0f_to_%.0fGeV", jetptMin, jetptMax);

    
    cSignals->Update();
    cSignals->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cSignals->SaveAs(imagePath + "SB_reflections_invMass_signals_fits_" + jetPtRange + ".png");
    cReflections->Update();
    cReflections->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cReflections->SaveAs(imagePath + "SB_reflections_invMass_reflections_fits_" + jetPtRange + ".png");
    cSignalsAndReflections->Update();
    cSignalsAndReflections->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cSignalsAndReflections->SaveAs(imagePath + "SB_reflections_invMass_signalsAndReflections_fits_" + jetPtRange + ".png");

    //
    // Storing in a single pdf file
    //
    cSignals->Print(Form("reflections_fits_%s.pdf(",jetPtRange.Data()));
    cReflections->Print(Form("reflections_fits_%s.pdf",jetPtRange.Data()));
    cSignalsAndReflections->Print(Form("reflections_fits_%s.pdf)",jetPtRange.Data()));
    
}

void ReflectionsTreatment(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;

    // pT,jet cuts
    std::vector<double> ptjetBinEdges = {5., 7., 15., 30.};
    const double jetptMin = ptjetBinEdges[0]; // GeV
    const double jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 1]; // GeV
    // DeltaR bins
    std::vector<double> deltaRBinEdges = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4};
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
    // pT,D bins
    std::vector<double> ptDBinEdges = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};

    // mass histogram
    int massBins = 50; // default=100 
    double minMass = 1.72; // use from 1.72, used to use 1.67
    double maxMass = 2.1;    

    // Opening data file
    TFile* fInputMC = new TFile("../../SimulatedData/Hyperloop_output/McEfficiency/New_with_reflections/output_HF_LHC24d3a_All/AO2Ds/Merged_HF_LHC24d3a_All.root","read");
    if (!fInputMC || fInputMC->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT file." << std::endl;
    }

    // Create multiple histograms
    HistogramGroup histograms = createHistograms(ptDBinEdges,                                 // the pT,D edges will determine the number of mass histograms
                                                     massBins, minMass, maxMass,                  // mass histograms binning
                                                     deltaRBinEdges);                             // deltaR histograms with asymmetrical bin widths

    // Fill histograms
    fillHistograms(fInputMC, histograms, jetptMin, jetptMax, ptDBinEdges);

    // Perform fits
    //std::vector<TF1*> fittings = performFit(histograms, parametersVectors, minMass, maxMass);
    FitsGroup fits = performFits(histograms, ptDBinEdges, minMass, maxMass);

    // signal/side-band region parameters
    int signalSigmas = 2; // 2
    int startingBackSigma = 4; // 4
    int backgroundSigmas = 4; // 4

    // Subtract side-band from signal
    //SubtractionResult finalDeltaR = SideBand(histograms, fittings, signalSigmas, startingBackSigma, backgroundSigmas);
    
    // Plot histograms
    PlotHistograms(histograms, fits, jetptMin, jetptMax);

    // Storing final histograms to output file
    //SaveData(fOut, outputStruct,jetptMin,jetptMax);

    // Save range analysis to 2D DeltaR vs pT,jet distribution
    //UpdateFinalDistribution(hDeltaR_vs_ptJet, outputStruct, jetptMin);
    
    // Clean current iteration data in order to avoid memory leaks
    //ClearData(outputStruct, massHistograms);

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    ReflectionsTreatment();
    return 0;
}



