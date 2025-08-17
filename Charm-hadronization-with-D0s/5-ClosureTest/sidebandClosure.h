/*
 * Macro for performing side-band subtraction procedure to second closure test
 * 
 * 
 * 
 * Author: Christian Reckziegel
**/


using namespace std;

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
Double_t backgroundFunction(Double_t* x, Double_t* par) {
    Double_t m = x[0];
    Double_t a = par[0];
    Double_t b = par[1];

    // Defining the custom function
    Double_t result = a * TMath::Power(m, b);
    return result;
}
Double_t signalFunction(Double_t* x, Double_t* par) {
    Double_t m = x[0];
    Double_t A1Signal = par[2]; // Free parameter
    Double_t A1toA2MCSignalRatio = par[3]; // Fixed parameter
    Double_t A2Signal = A1Signal / A1toA2MCSignalRatio; // Constrained to primary/secondary integral obtained from MC data fit
    Double_t m0 = par[4]; // Free parameter
    Double_t sigma1 = par[5]; // Free parameter
    Double_t Sigma1toSigma2MCSignalRatio = par[6]; // Fixed parameter
    Double_t sigma2 = sigma1 / Sigma1toSigma2MCSignalRatio; // Constrained to primary/secondary width obtained from MC data fit

    // Defining the custom function
    Double_t result = A1Signal * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma1), 2)) + A2Signal * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma2), 2));
    return result;
}
Double_t reflectionFunction(Double_t* x, Double_t* par) {
    // m0 and sigma values fixed from fits obtained from MC data fits
    Double_t m = x[0];
    Double_t A1SignalToA1ReflectionMCRatio = par[7]; // Fixed parameter
    // Extract A1Signal from the parameter array (same index as in `signalFunction`)
    Double_t A1Signal = par[2];
    Double_t A1Reflection = A1Signal / A1SignalToA1ReflectionMCRatio; // Constrained to signal/reflection integral obtained from MC data fits
    Double_t A1toA2MCReflectionsRatio = par[8]; // Fixed parameter
    Double_t A2Reflection = A1Reflection / A1toA2MCReflectionsRatio; // Constrained to primary/secondary integral obtained from MC data fit
    Double_t m0_1 = par[9]; // Fixed parameter
    Double_t m0_2 = par[10]; // Fixed parameter
    Double_t sigma_1 = par[11]; // Fixed parameter
    Double_t sigma_2 = par[12]; // Fixed parameter

    // Defining the custom function
    Double_t result = A1Reflection * TMath::Exp(-TMath::Power((m - m0_1) / (2 * sigma_1), 2)) + A2Reflection * TMath::Exp(-TMath::Power((m - m0_2) / (2 * sigma_2), 2));
    return result;
}
Double_t customFitFunction(Double_t* x, Double_t* par) {
    Double_t m = x[0];
    
    // Background component
    Double_t background = backgroundFunction(x, par); 

    // Signal component
    Double_t signalComponent = signalFunction(x, par);  // Signal function uses signal parameters
    
    // Reflection component
    Double_t reflectionComponent = reflectionFunction(x, par);  // Reflection function uses reflection parameters

    // Total fit: sum of signal and reflection components
    Double_t totalFit = background + signalComponent + reflectionComponent;
    return totalFit;
}
// Stand-alone functions for signal and reflection components
Double_t signalOnlyFunction(Double_t* x, Double_t* par) {
    Double_t m = x[0];
    Double_t A1Signal = par[0]; // Free parameter
    Double_t A1toA2MCSignalRatio = par[1]; // Fixed parameter
    Double_t A2Signal = A1Signal / A1toA2MCSignalRatio; // Constrained to primary/secondary integral obtained from MC data fit
    Double_t m0 = par[2]; // Free parameter
    Double_t sigma1 = par[3]; // Free parameter
    Double_t Sigma1toSigma2MCSignalRatio = par[4]; // Fixed parameter
    Double_t sigma2 = sigma1 / Sigma1toSigma2MCSignalRatio; // Constrained to primary/secondary width obtained from MC data fit

    // Defining the custom function
    Double_t result = A1Signal * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma1), 2)) + A2Signal * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma2), 2));
    return result;
}
Double_t reflectionOnlyFunction(Double_t* x, Double_t* par) {
    // m0 and sigma values fixed from fits obtained from MC data fits
    Double_t m = x[0];
    Double_t A1SignalToA1ReflectionMCRatio = par[0]; // Fixed parameter
    // Extract A1Signal from the parameter array (same index as in `signalFunction`)
    Double_t A1Signal = par[1];
    Double_t A1Reflection = A1Signal / A1SignalToA1ReflectionMCRatio; // Constrained to signal/reflection integral obtained from MC data fits
    Double_t A1toA2MCReflectionsRatio = par[2]; // Fixed parameter
    Double_t A2Reflection = A1Reflection / A1toA2MCReflectionsRatio; // Constrained to primary/secondary integral obtained from MC data fit
    Double_t m0_1 = par[3]; // Fixed parameter
    Double_t m0_2 = par[4]; // Fixed parameter
    Double_t sigma_1 = par[5]; // Fixed parameter
    Double_t sigma_2 = par[6]; // Fixed parameter

    // Defining the custom function
    Double_t result = A1Reflection * TMath::Exp(-TMath::Power((m - m0_1) / (2 * sigma_1), 2)) + A2Reflection * TMath::Exp(-TMath::Power((m - m0_2) / (2 * sigma_2), 2));
    return result;
}
// calculate delta phi such that 0 < delta phi < 2*pi
double DeltaPhi(double phi1, double phi2) {
    // Compute the absolute difference between phi1 and phi2
    double dphi = std::abs(phi1 - phi2); 
    if (dphi > M_PI) {
        // subtract 2pi if the difference if bigger than pi
        dphi = dphi - 2*M_PI;
    }

    return dphi;
}

// Get the optimal BDT score cut for the corresponding pT,D of the entry
double GetBkgProbabilityCut(double pT, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    for (size_t i = 0; i < bdtPtCuts.size() - 1; ++i) {
        if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) {
            return bdtPtCuts[i].second;
        }
    }
    return 1.0; // Default: accept all if out of range
}

// Module to perform fits to histograms
struct FitContainer {
    std::vector<TF1*> fitBackgroundOnly;        // background only fits
    std::vector<TF1*> fitSignalOnly;            // signal only fits
    std::vector<TF1*> fitReflectionsOnly;       // reflection only fits
    std::vector<TF1*> fitSignalsAndReflections; // 
    std::vector<TF1*> fitTotal;                 // total fits
};
// Histograms containers for each case
struct HistogramGroup {
    std::vector<TH2D*> signals;
    std::vector<TH2D*> reflections;
    std::vector<TH2D*> signals_and_reflections;

    std::vector<TH1D*> signals_1d;
    std::vector<TH1D*> reflections_1d;
    std::vector<TH1D*> signals_and_reflections_1d;
};
std::pair<FitContainer, TCanvas*> calculateFitTemplates(TFile* fClosureInput, const double& jetptMin, const double& jetptMax, const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& ptDBinEdges_detector,
                                   const std::vector<std::pair<double, double>>& bdtPtCuts, int xbins, double xmin, double xmax) {
    // Create template fits container
    FitContainer fTemplateFits;

    // Create template histograms for template fits from detector level correction data sample
    HistogramGroup histogramTemplates;
    for (size_t i = 0; i < ptDBinEdges_detector.size() - 1; ++i) {
        // So that the title adapts to fractional binning title
        if (std::fmod(ptDBinEdges_detector[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(ptDBinEdges_detector[i+1], 1.0) != 0) {
                // pure reflections histograms
                histogramTemplates.reflections.push_back(new TH2D(Form("histMassTemplate_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
                // pure signals histograms
                histogramTemplates.signals.push_back(new TH2D(Form("histMassTemplate_signals%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
                // signals+reflections histograms: calculate ratios
                histogramTemplates.signals_and_reflections.push_back(new TH2D(Form("histMassTemplate_signals_and_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
            } else {
                // pure reflections histograms
                histogramTemplates.reflections.push_back(new TH2D(Form("histMassTemplate_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
                // pure signals histograms
                histogramTemplates.signals.push_back(new TH2D(Form("histMassTemplate_signals%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
                // signals+reflections histograms: calculate ratios
                histogramTemplates.signals_and_reflections.push_back(new TH2D(Form("histMassTemplate_signals_and_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
            }
        } else {
            if (std::fmod(ptDBinEdges_detector[i+1], 1.0) != 0) {
                // pure reflections histograms
                histogramTemplates.reflections.push_back(new TH2D(Form("histMassTemplate_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
                // pure signals histograms
                histogramTemplates.signals.push_back(new TH2D(Form("histMassTemplate_signals%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
                // signals+reflections histograms: calculate ratios
                histogramTemplates.signals_and_reflections.push_back(new TH2D(Form("histMassTemplate_signals_and_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
            } else {
                // pure reflections histograms
                histogramTemplates.reflections.push_back(new TH2D(Form("histMassTemplate_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
                // pure signals histograms
                histogramTemplates.signals.push_back(new TH2D(Form("histMassTemplate_signals%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
                // signals+reflections histograms: calculate ratios
                histogramTemplates.signals_and_reflections.push_back(new TH2D(Form("histMassTemplate_signals_and_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), xbins, xmin, xmax, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
            }
        }
        // pure reflections histograms
        histogramTemplates.reflections[i]->Sumw2();
        // pure signals histograms
        histogramTemplates.signals[i]->Sumw2();
        // signals+reflections histograms: calculate ratios
        histogramTemplates.signals_and_reflections[i]->Sumw2();
    }

    // Fill histograms with appropriate cuts
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = deltaRBinEdges_detector[deltaRBinEdges_detector.size() - 1]; // on particle level delta R
    // Accessing detector level data TTree
    TTree* tree = (TTree*)fClosureInput->Get("CorrectionTree");
    // Assuming histograms and tree data correspond in some way
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing the TTree
    float axisDistance, jetPt, jetEta, jetPhi;
    float hfPt, hfEta, hfPhi, hfMass, hfY, hfMlScore0;
    int hfMatchedFrom, hfSelectedAs, jetNConst;
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
        /*if (hfSelectedAs & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as D0
            hfSelectedAs = 1;
        } else if (hfSelectedAs & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as D0bar
            hfSelectedAs = -1;
        }*/
        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));
        // Fill each histogram with their respective pT intervals
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= deltaRBinEdges_detector[0]) && (deltaR < DeltaRcut))) {
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < ptDBinEdges_detector.size() - 1 && !filled; iEdge++) {
                if ((hfPt >= ptDBinEdges_detector[iEdge]) && (hfPt < ptDBinEdges_detector[iEdge + 1])) {
                    // D0 = +1, D0bar = -1, neither = 0
                    if ((hfMatchedFrom != 0) && (hfSelectedAs != 0)) {
                        if (hfMatchedFrom == hfSelectedAs) {
                            // Get the threshold for this pT range
                            double maxBkgProb = GetBkgProbabilityCut(hfPt, bdtPtCuts);
                            // Fill histogram only if the cut is passed
                            if (hfMlScore0 < maxBkgProb) {
                                // pure signals
                                histogramTemplates.signals[iEdge]->Fill(hfMass, deltaR); // 2D case: (hfMass, deltaR);
                            }
                        } else {
                            // Get the threshold for this pT range
                            double maxBkgProb = GetBkgProbabilityCut(hfPt, bdtPtCuts);
                            // Fill histogram only if the cut is passed
                            if (hfMlScore0 < maxBkgProb) {
                                // pure reflections
                                histogramTemplates.reflections[iEdge]->Fill(hfMass, deltaR); // 2D case: (hfMass, deltaR);
                            }
                        }
                        // Get the threshold for this pT range
                        double maxBkgProb = GetBkgProbabilityCut(hfPt, bdtPtCuts);
                        // Fill histogram only if the cut is passed
                        if (hfMlScore0 < maxBkgProb) {
                            // signals and reflections altogether, withOUT "neither = 0" entries
                            histogramTemplates.signals_and_reflections[iEdge]->Fill(hfMass, deltaR); // 2D case: (hfMass, deltaR);
                        }
                    }
                    // signals and reflections altogether, wiTH "neither = 0" entries (should I include "neither" entries too? Altough they should not appear)
                    //histogramTemplates.signals_and_reflections[iEdge]->Fill(hfMass, deltaR);
                    filled = true; // Exit the loop once the correct histogram is found
                }
            }
        }
    }

    // Fit signal and reflection models to the histograms
    // Signal fits
    for (size_t iInterval = 0; iInterval < ptDBinEdges_detector.size() - 1; iInterval++) {
        TF1* signalFit = new TF1(Form("signalFit_%zu_%.0f_to_%.0fGeV", iInterval, jetptMin, jetptMax), pureSignalFunction, xmin, xmax, 5);

        // Obtain 1D mass histogram projection
        histogramTemplates.signals_1d.push_back( histogramTemplates.signals[iInterval]->ProjectionX(Form("histMass_signals_1d_%zu_%.0f_to_%.0fGeV", iInterval+1, jetptMin, jetptMax), 1, -1) );

        // Initialize parameters from the histogram distribution characteristics
        double C1, C2, m0, sigma1, sigma2;
        C1 = 0.7 * histogramTemplates.signals_1d[iInterval]->GetMaximum(); // Assume first Gaussian dominates
        C2 = 0.3 * histogramTemplates.signals_1d[iInterval]->GetMaximum(); // Second Gaussian smaller
        m0 = histogramTemplates.signals_1d[iInterval]->GetMean();       // Mean of the distribution
        sigma1 = histogramTemplates.signals_1d[iInterval]->GetRMS() / 2; // First Gaussian narrower
        sigma2 = histogramTemplates.signals_1d[iInterval]->GetRMS();    // Second Gaussian wider
        signalFit->SetParameters(C1, C2, m0, sigma1, sigma2);
        signalFit->SetParName(0, "C1");
        signalFit->SetParName(1, "C2");
        signalFit->SetParName(2, "m0");
        signalFit->SetParName(3, "sigma1");
        signalFit->SetParName(4, "sigma2");
        //signalFit->SetParLimits(3, 0.5 * sigma1, 2.0 * sigma1); // Example constraint for sigma1

        histogramTemplates.signals_1d[iInterval]->Fit(signalFit, "RQ"); // "Q" option performs quiet fit without drawing the fit function
        fTemplateFits.fitSignalOnly.push_back(signalFit);
        //fTemplateFits.fitSignalOnly[iInterval]->Print("V");
    }

    // Reflections fits
    for (size_t iInterval = 0; iInterval < ptDBinEdges_detector.size() - 1; iInterval++) {
        TF1* reflectionsFit = new TF1(Form("reflectionsFit_%zu_%.0f_to_%.0fGeV", iInterval, jetptMin, jetptMax), pureReclectionsFunction, xmin, xmax, 6);

        // Obtain 1D mass histogram projection
        histogramTemplates.reflections_1d.push_back( histogramTemplates.reflections[iInterval]->ProjectionX(Form("histMass_reflections_1d_%zu_%.0f_to_%.0fGeV", iInterval+1, jetptMin, jetptMax), 1, -1) );

        // Initialize parameters from the histogram distribution characteristics
        double C1, C2, m0_1, m0_2, sigma1, sigma2;
        C1 = 0.6 * histogramTemplates.reflections_1d[iInterval]->GetMaximum();   // Dominant reflection
        C2 = 0.4 * histogramTemplates.reflections_1d[iInterval]->GetMaximum();   // Secondary reflection
        m0_1 = histogramTemplates.reflections_1d[iInterval]->GetMean() - histogramTemplates.reflections_1d[iInterval]->GetRMS() / 2; // Offset from mean
        m0_2 = histogramTemplates.reflections_1d[iInterval]->GetMean() + histogramTemplates.reflections_1d[iInterval]->GetRMS() / 2; // Offset from mean
        sigma1 = histogramTemplates.reflections_1d[iInterval]->GetRMS() / 2;  // Narrower
        sigma2 = histogramTemplates.reflections_1d[iInterval]->GetRMS();      // Wider
        reflectionsFit->SetParameters(C1, C2, m0_1, m0_2, sigma1, sigma2);

        histogramTemplates.reflections_1d[iInterval]->Fit(reflectionsFit, "RQ");
        fTemplateFits.fitReflectionsOnly.push_back(reflectionsFit);
        fTemplateFits.fitReflectionsOnly[iInterval]->Print("V");
        std::cout << "A1/A2 = " << fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(0) / fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(1) << std::endl;
    }

    // Combined fits: true signal + reflections
    for (size_t iInterval = 0; iInterval < ptDBinEdges_detector.size() - 1; iInterval++) {
        TF1* combinedFit = new TF1(Form("combinedFit_%zu_%.0f_to_%.0fGeV", iInterval, jetptMin, jetptMax), signalAndReflectionsFunction, xmin, xmax, 11);

        // Obtain 1D mass histogram projection
        histogramTemplates.signals_and_reflections_1d.push_back( histogramTemplates.signals_and_reflections[iInterval]->ProjectionX(Form("histMass_signals_and_reflections_1d_%zu_%.0f_to_%.0fGeV", iInterval+1, jetptMin, jetptMax), 1, -1) );

        // Initialize parameters from previous fits
        combinedFit->SetParameters(fTemplateFits.fitSignalOnly[iInterval]->GetParameter(0), fTemplateFits.fitSignalOnly[iInterval]->GetParameter(1), fTemplateFits.fitSignalOnly[iInterval]->GetParameter(2), fTemplateFits.fitSignalOnly[iInterval]->GetParameter(3), fTemplateFits.fitSignalOnly[iInterval]->GetParameter(4), 
                                      fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(0), fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(1), fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(2), fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(3), fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(4), fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(5));

        histogramTemplates.signals_and_reflections_1d[iInterval]->Fit(combinedFit, "RQ");
        fTemplateFits.fitSignalsAndReflections.push_back(combinedFit);
        //fTemplateFits.fitSignalsAndReflections[iInterval]->Print("V");
    }

    // Plot templates histograms and fits
    // Create a canvas with a pad for each pT,D interval
    TCanvas* cTemplateFits = new TCanvas("cTemplateFits",Form("Template histograms and fits %.0f < pT,jet < %0.f GeV/c", jetptMin, jetptMax));
    int nHistos = histogramTemplates.signals_and_reflections_1d.size();
    // Start with a square layout (or close to it)
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));
    cTemplateFits->Divide(nCols,nRows);
    for (size_t iInterval = 0; iInterval < ptDBinEdges_detector.size() - 1; iInterval++) {
        cTemplateFits->cd(iInterval + 1);
        // Draw the 1D histograms
        histogramTemplates.signals_and_reflections_1d[iInterval]->SetMarkerStyle(kDot); //kFullDotMedium
        histogramTemplates.signals_and_reflections_1d[iInterval]->SetMarkerColor(kBlack);
        histogramTemplates.signals_and_reflections_1d[iInterval]->SetLineColor(kBlack);
        histogramTemplates.signals_and_reflections_1d[iInterval]->GetYaxis()->SetTitle("counts");
        histogramTemplates.signals_and_reflections_1d[iInterval]->SetTitle(Form("%.1f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c}", ptDBinEdges_detector[iInterval], ptDBinEdges_detector[iInterval + 1]));
        histogramTemplates.signals_and_reflections_1d[iInterval]->Draw();
        // Draw the fits
        fTemplateFits.fitSignalsAndReflections[iInterval]->SetLineColor(kRed);
        fTemplateFits.fitSignalsAndReflections[iInterval]->SetLineWidth(2);
        fTemplateFits.fitSignalsAndReflections[iInterval]->Draw("same");
        // Draw the signal fit
        fTemplateFits.fitSignalOnly[iInterval]->SetLineColor(kBlue);
        fTemplateFits.fitSignalOnly[iInterval]->SetLineWidth(2);
        fTemplateFits.fitSignalOnly[iInterval]->Draw("same");
        // Draw the reflections fit
        fTemplateFits.fitReflectionsOnly[iInterval]->SetLineColor(kGreen);
        fTemplateFits.fitReflectionsOnly[iInterval]->SetLineWidth(2);
        fTemplateFits.fitReflectionsOnly[iInterval]->Draw("same");
        // Add a legend
        TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
        legend->AddEntry(histogramTemplates.signals_and_reflections_1d[iInterval], "Data", "p");
        legend->AddEntry(fTemplateFits.fitSignalsAndReflections[iInterval], "Signal + Reflections Fit", "l");
        legend->AddEntry(fTemplateFits.fitSignalOnly[iInterval], "Signal Fit", "l");
        legend->AddEntry(fTemplateFits.fitReflectionsOnly[iInterval], "Reflections Fit", "l");
        legend->Draw();
    }
    std::pair<FitContainer, TCanvas*> fTemplateFitsAndCanvas = std::make_pair(fTemplateFits, cTemplateFits);
    return fTemplateFitsAndCanvas;
};
// Save images to pdf file
void storeImages(std::vector<TH2D*>& hInvariantMass2D, std::vector<TH1D*>& hInvariantMass1D, const double& jetptMin, const double& jetptMax, TCanvas* cTemplateFits, TCanvas* c1d_fit) {
    
    int nHistos = hInvariantMass2D.size();
    // Start with a square layout (or close to it)
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));
    
    // Creating canvases
    TCanvas* cInvariantMass2D = new TCanvas("cInvariantMass2D","DeltaR as a function of invariant mass plots");
    cInvariantMass2D->SetCanvasSize(1800,1000);
    cInvariantMass2D->Divide(nCols,nRows); // columns, lines
    TCanvas* cInvariantMass1D = new TCanvas("cInvariantMass1D","Invariant mass plots");
    cInvariantMass1D->SetCanvasSize(1800,1000);
    cInvariantMass1D->Divide(nCols,nRows); // columns, lines

    // Loop through all histograms (and fitting functions in the future)
    for(size_t iHisto = 0; iHisto < hInvariantMass2D.size() - 1; ++iHisto) {
        
        // Skip not filled histograms
        if (hInvariantMass2D[iHisto]->GetEntries() == 0) {
            continue;
        }
        
        // Drawing 2D histograms
        cInvariantMass2D->cd(iHisto+1);
        hInvariantMass2D[iHisto]->Draw("colz");

        // Drawing 1D histograms
        cInvariantMass1D->cd(iHisto+1);
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        hInvariantMass1D[iHisto]->SetMarkerStyle(kDot); //kFullDotMedium
        hInvariantMass1D[iHisto]->SetMarkerColor(kBlack);
        hInvariantMass1D[iHisto]->SetLineColor(kBlack);
        hInvariantMass1D[iHisto]->GetYaxis()->SetTitle("counts");
        //hInvariantMass1D->SetMinimum(0);
        hInvariantMass1D[iHisto]->Draw();
    }

    //
    // Storing images in a single pdf file
    //
    TString imagePath = "../Images/5-ClosureTest/Second/";
    cTemplateFits->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cInvariantMass2D->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cInvariantMass1D->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    c1d_fit->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

    delete cInvariantMass2D;
    delete cInvariantMass1D;

}

TH2D* AnalyzeJetPtRange(TFile* fClosureInput, const double& jetptMin, const double& jetptMax, 
                        const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& ptDBinEdges_detector,
                        const std::vector<std::pair<double, double>>& bdtPtCuts) {

    // mass histogram
    int massBins = 50; // default=100 
    double minMass = 1.72; // use from 1.72, used to use 1.67
    double maxMass = 2.06; // old = 2.1, new = 2.05 + 0.01
    // Mass and std reference values for fits
    double m_0_reference = 1.86484; // D0 mass in GeV/c^2
    double sigma_reference = 0.012;

    // ----- Create and fill invariant mass distributions with detector level data
    std::vector<TH2D*> hInvariantMass2D;
    for (size_t i = 0; i < ptDBinEdges_detector.size() - 1; ++i) {
        // So that the title adapts to fractional binning title
        if (std::fmod(ptDBinEdges_detector[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(ptDBinEdges_detector[i+1], 1.0) != 0) {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), massBins, minMass, maxMass, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
        
            } else {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), massBins, minMass, maxMass, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
        
            }
        } else {
            if (std::fmod(ptDBinEdges_detector[i+1], 1.0) != 0) {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), massBins, minMass, maxMass, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
        
            } else {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), massBins, minMass, maxMass, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
        
            }
        }
        hInvariantMass2D[i]->Sumw2();
    }
    // ----- Fill data
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCDyCut = 0.8; // on detector level D0
    const double MCDDeltaRcut = deltaRBinEdges_detector[deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double MCDjetptMin = jetptMin;
    const double MCDjetptMax = jetptMax;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;
    // variables for reflection contribution selection
    int MCDhfMatchedFrom, MCDhfSelectedAs;
    // Accessing TTree
    TTree* tree = (TTree*)fClosureInput->Get("InputTree");
    // Check for correct access
    if (!tree) {
        cout << "Error opening correction data tree.\n";
    }
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
    tree->SetBranchAddress("fHfMatchedFrom",&MCDhfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&MCDhfSelectedAs);
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        bool isReflection = (MCDhfMatchedFrom != MCDhfSelectedAs) ? true : false;

        // Apply prompt (real+reflections) selection
        if (!MCDhfprompt) {
            continue;
        }
        

        // calculating delta R
        double MCDdeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));
        
        // Fill each histogram with their respective pT intervals
        if ((abs(MCDhfEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= MCDjetptMin) && (MCDjetPt < MCDjetptMax)) && ((MCDdeltaR >= deltaRBinEdges_detector[0]) && (MCDdeltaR < MCDDeltaRcut))) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < ptDBinEdges_detector.size() - 1 && !filled; iEdge++) {
                if ((MCDhfPt >= ptDBinEdges_detector[iEdge]) && (MCDhfPt < ptDBinEdges_detector[iEdge + 1])) {
                    // Get the threshold for this pT range
                    double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);

                    // Fill histogram only if the cut is passed
                    if (MCDhfMlScore0 < maxBkgProb) {
                        hInvariantMass2D[iEdge]->Fill(MCDhfMass, MCDdeltaR);
                    }
                    filled = true; // Exit the loop once the correct histogram is found (alternative: break)
                }
                
            }
            
        }

        
        
        
    }
    // ----- Obtain 1D projections: invariant mass distributions
    std::vector<TH1D*> hInvariantMass1D;
    for (size_t iHisto = 0; iHisto < hInvariantMass2D.size() - 1; iHisto++) {

        hInvariantMass1D.push_back(hInvariantMass2D[iHisto]->ProjectionX(Form("h_mass_proj_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax)));
    }
    
    // ----- Obtain template fits from MC particle level data
    std::pair<FitContainer, TCanvas*> fTemplateFitsAndCanvas = calculateFitTemplates(fClosureInput, jetptMin, jetptMax, deltaRBinEdges_detector, ptDBinEdges_detector, bdtPtCuts, massBins, minMass, maxMass);
    FitContainer fTemplateFits = fTemplateFitsAndCanvas.first;

    // ----- Fit detector level data
    TCanvas* c1d_fit = new TCanvas("c1d_fit", "Input sample 1D histograms with fit");
    c1d_fit->SetCanvasSize(1800,1000);
    int nHistos = hInvariantMass2D.size();
    // Start with a square layout (or close to it)
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));
    c1d_fit->Divide(nCols,nRows); // columns, lines
    FitContainer fDataFits;
    for (size_t iHisto = 0; iHisto < hInvariantMass2D.size()-1; iHisto++) {

        
        c1d_fit->cd(iHisto+1);

        double A1Signal = fTemplateFits.fitSignalOnly[iHisto]->GetParameter(0); // A1Signal
        double A2Signal = fTemplateFits.fitSignalOnly[iHisto]->GetParameter(1); // A2Signal
        double m0Signal = fTemplateFits.fitSignalOnly[iHisto]->GetParameter(2); // m0
        double sigma1Signal = fTemplateFits.fitSignalOnly[iHisto]->GetParameter(3); // sigma1
        double sigma2Signal = fTemplateFits.fitSignalOnly[iHisto]->GetParameter(4); // sigma2
        //std::cout << "A1Signal = " << A1Signal << ", A2Signal = " << A2Signal << ", m0Signal = " << m0Signal << ", sigma1Signal = " << sigma1Signal << ", sigma2Signal = " << sigma2Signal << std::endl;

        double A1Reflection = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(0); // A1Reflection
        double A2Reflection = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(1); // A2Reflection
        double m0_1Reflections = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(2); // m0_1
        double m0_2Reflections = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(3); // m0_2
        double sigma1Reflections = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(4); // sigma1
        double sigma2Reflections = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(5); // sigma2
        //std::cout << "A1Reflection = " << A1Reflection << ", A2Reflection = " << A2Reflection << ", m0_1Reflections = " << m0_1Reflections << ", m0_2Reflections = " << m0_2Reflections << ", sigma1Reflections = " << sigma1Reflections << ", sigma2Reflections = " << sigma2Reflections << std::endl;

        // -> signal constrains
        double A1toA2MCSignalRatio = A1Signal / A2Signal;
        double Sigma1toSigma2MCSignalRatio = sigma1Signal / sigma2Signal;

        // -> reflection constrains
        double A1toA2MCReflectionsRatio = A1Reflection / A2Reflection;
        double Sigma1toSigma2MCReflectionsRatio = sigma1Reflections / sigma2Reflections;
        double A1SignalToA1ReflectionMCRatio = A1Signal / A1Reflection;

        // ----Perform fits----
        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        fDataFits.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), customFitFunction, minMass, maxMass, 13)); // 12 parameters = 8 fixed + 5 free
        // Set initial values and fix parameters
        Double_t params[13] = {200000, -10.0, 1000, 1.5, m_0_reference, 0.02, 1.2, 2.0, 1.3, 1.83, 1.85, 0.02, 0.03};
        fDataFits.fitTotal[iHisto]->SetParameters(params);
        // Set parameter names
        fDataFits.fitTotal[iHisto]->SetParName(0, "Background A");
        fDataFits.fitTotal[iHisto]->SetParName(1, "Background B");
        fDataFits.fitTotal[iHisto]->SetParName(2, "A1 Signal");
        fDataFits.fitTotal[iHisto]->SetParName(3, "A1/A2 Signal Ratio");
        fDataFits.fitTotal[iHisto]->SetParName(4, "Signal Mean m0");
        fDataFits.fitTotal[iHisto]->SetParName(5, "Signal Width Sigma1");
        fDataFits.fitTotal[iHisto]->SetParName(6, "Sigma Ratio (Sig1/Sig2)");
        fDataFits.fitTotal[iHisto]->SetParName(7, "A1Signal/A1Reflection Ratio");
        fDataFits.fitTotal[iHisto]->SetParName(8, "A1/A2 Reflection Ratio");
        fDataFits.fitTotal[iHisto]->SetParName(9, "Reflection Mean m0_1");
        fDataFits.fitTotal[iHisto]->SetParName(10, "Reflection Mean m0_2");
        fDataFits.fitTotal[iHisto]->SetParName(11, "Reflection Width Sigma_1");
        fDataFits.fitTotal[iHisto]->SetParName(12, "Reflection Width Sigma_2");
        // Fix the parameters from MC fits
        fDataFits.fitTotal[iHisto]->FixParameter(3, A1toA2MCSignalRatio);   // A1toA2MCSignalRatio
        fDataFits.fitTotal[iHisto]->FixParameter(6, Sigma1toSigma2MCSignalRatio);   // Sigma1toSigma2MCSignalRatio
        fDataFits.fitTotal[iHisto]->FixParameter(7, A1SignalToA1ReflectionMCRatio);   // A1SignalToA1ReflectionMCRatio
        fDataFits.fitTotal[iHisto]->FixParameter(8, A1toA2MCReflectionsRatio);   // A1toA2MCReflectionsRatio
        fDataFits.fitTotal[iHisto]->FixParameter(9, m0_1Reflections);  // m0_1
        fDataFits.fitTotal[iHisto]->FixParameter(10, m0_2Reflections); // m0_2
        fDataFits.fitTotal[iHisto]->FixParameter(11, sigma1Reflections); // sigma_1
        fDataFits.fitTotal[iHisto]->FixParameter(12, sigma2Reflections); // sigma_2
        // Apply range limits to the parameters
        fDataFits.fitTotal[iHisto]->SetParLimits(4, 0.95 * m_0_reference, 1.05 * m_0_reference);
        fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.35 * sigma_reference, 3.0 * sigma_reference);
        // Choose line color
        fDataFits.fitTotal[iHisto]->SetLineColor(kBlack);
        // Perform fit with "Q" (quiet) option: no drawing of the fit function
        hInvariantMass1D[iHisto]->Fit(fDataFits.fitTotal[iHisto], "Q");

        // Getting total parameter values
        double a_par = fDataFits.fitTotal[iHisto]->GetParameter(0); // Get the value of parameter 'a'
        double b_par = fDataFits.fitTotal[iHisto]->GetParameter(1); // Get the value of parameter 'b'
        fDataFits.fitBackgroundOnly.push_back(new TF1(Form("backgroundOnlyFit_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax), backgroundFunction, minMass, maxMass, 2));
        fDataFits.fitBackgroundOnly[iHisto]->SetParameters(a_par,b_par); // fDataFits.size()-1 = the latest added to the vector
        fDataFits.fitBackgroundOnly[iHisto]->SetLineStyle(kDashed);
        fDataFits.fitBackgroundOnly[iHisto]->SetLineColor(kRed+1);

        // Extract signal-related parameters from the total fit
        double a1Signal = fDataFits.fitTotal[iHisto]->GetParameter(2); // A1 signal
        A1toA2MCSignalRatio = fDataFits.fitTotal[iHisto]->GetParameter(3); // A2 signal
        m0Signal = fDataFits.fitTotal[iHisto]->GetParameter(4); // Signal mean m0
        double sigmaSignal = fDataFits.fitTotal[iHisto]->GetParameter(5); // Signal width sigma
        Sigma1toSigma2MCSignalRatio = fDataFits.fitTotal[iHisto]->GetParameter(6); // Sigma1/Sigma2 signal
        // Perfoming fit with acquired parameters from total fit
        fDataFits.fitSignalOnly.push_back(new TF1(Form("signalOnlyFit_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax), signalOnlyFunction, minMass, maxMass, 5));
        fDataFits.fitSignalOnly[iHisto]->SetParameters(a1Signal, A1toA2MCSignalRatio, m0Signal, sigmaSignal, Sigma1toSigma2MCSignalRatio); // fDataFits.size()-1 = the latest added to the vector
        fDataFits.fitSignalOnly[iHisto]->SetLineStyle(kDashed);
        fDataFits.fitSignalOnly[iHisto]->SetLineColor(kBlue+1);

        // Extract reflection-related parameters from the total fit
        double A1SignalToA1ReflectionMCRatios = fDataFits.fitTotal[iHisto]->GetParameter(7); // A1 signal / A1 reflection
        A1Signal = fDataFits.fitTotal[iHisto]->GetParameter(2); // A1 signal
        A1toA2MCReflectionsRatio = fDataFits.fitTotal[iHisto]->GetParameter(8); // A1 reflection / A2 reflection
        double m0_1Reflection = fDataFits.fitTotal[iHisto]->GetParameter(9); // Reflection mean m0_1
        double m0_2Reflection = fDataFits.fitTotal[iHisto]->GetParameter(10); // Reflection mean m0_2
        double sigma1Reflection = fDataFits.fitTotal[iHisto]->GetParameter(11); // Reflection sigma_1
        double sigma2Reflection = fDataFits.fitTotal[iHisto]->GetParameter(12); // Reflection sigma_2
        //std:cout << "A1SignalToA1ReflectionMCRatios = " << A1SignalToA1ReflectionMCRatios << ", A1Signal = " << A1Signal << ", A1toA2MCReflectionsRatio = " << A1toA2MCReflectionsRatio << ", m0_1Reflection = " << m0_1Reflection << ", m0_2Reflection = " << m0_2Reflection << ", sigma1Reflection = " << sigma1Reflection << ", sigma2Reflection = " << sigma2Reflection << std::endl;

        // Perfoming fit with acquired parameters from total fit
        fDataFits.fitReflectionsOnly.push_back(new TF1(Form("reflectionOnlyFit_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax), reflectionOnlyFunction, minMass, maxMass, 7));
        fDataFits.fitReflectionsOnly[iHisto]->SetParameters(A1SignalToA1ReflectionMCRatios, A1Signal, A1toA2MCReflectionsRatio, m0_1Reflection, m0_2Reflection, sigma1Reflection, sigma2Reflection); // fDataFits.size()-1 = the latest added to the vector
        fDataFits.fitReflectionsOnly[iHisto]->SetLineStyle(kDashed);
        fDataFits.fitReflectionsOnly[iHisto]->SetLineColor(kGreen+3);

        
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        hInvariantMass1D[iHisto]->Draw();
        fDataFits.fitTotal[iHisto]->Draw("same");
        fDataFits.fitSignalOnly[iHisto]->Draw("same");
        fDataFits.fitBackgroundOnly[iHisto]->Draw("same");
        fDataFits.fitReflectionsOnly[iHisto]->Draw("same");

    }

    // ----- Perform side-band subtraction to fitted detector level data

    // ----- Create 2D distribution of DeltaR vs. pT,D0
    TH2D* hDeltaR_vs_PtD;

    // Store images
    storeImages(hInvariantMass2D, hInvariantMass1D, jetptMin, jetptMax, fTemplateFitsAndCanvas.second, c1d_fit);

    // Clean vector of pointers
    for (auto hist : hInvariantMass2D) {
        delete hist;
    }
    hInvariantMass2D.clear();
    for (auto hist : hInvariantMass1D) {
        delete hist;
    }
    hInvariantMass1D.clear();

    // Return resulting ditribution
    return hDeltaR_vs_PtD;
}

TH3D* create3DBackgroundSubtracted(const std::vector<TH2D*>& outputHistograms) {
    //

    TH3D* hBackgroundSubtracted;

    return hBackgroundSubtracted;
}

TH3D* SidebandClosure(TFile* fClosureInput, const std::vector<double>& ptjetBinEdges_detector, const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& ptDBinEdges_detector, const std::vector<std::pair<double, double>>& bdtPtCuts) {

    // One TH2D for each pT,jet range computed
    std::vector<TH2D*> outputHistograms;
    for (size_t iPtJetRange = 0; iPtJetRange < ptjetBinEdges_detector.size() - 1; iPtJetRange++) {
        TH2D* hDeltaR_vs_PtD = AnalyzeJetPtRange(fClosureInput, ptjetBinEdges_detector[iPtJetRange], ptjetBinEdges_detector[iPtJetRange + 1], deltaRBinEdges_detector, ptDBinEdges_detector, bdtPtCuts);
        outputHistograms.push_back(hDeltaR_vs_PtD);
    }

    TH3D* hBackgroundSubtracted = create3DBackgroundSubtracted(outputHistograms);
    return hBackgroundSubtracted;
}
