/**
 * @file ReflectionsTreatment.C
 * @author Christian Reckziegel
 * @brief Macro for obtaining the fit parameters, integrals and reflection scaling parameter for each invariant mass histogram fit
**/

using namespace std;

//--- bit manipulation --------------------------------------------------------- (TODO: need of direct bit manipulation removed after latest PR)
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

    // Individual components
    std::vector<std::pair<TF1*,TF1*>> individualSignals;
    std::vector<std::pair<TF1*,TF1*>> individualReflections;
    std::vector<std::pair<TF1*,TF1*>> individualSignals_and_reflections;
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
// Custom single gaussian fit function
Double_t singleGaussianFunction(Double_t* x, Double_t* par) {

    Double_t m = x[0];
    Double_t C_1 = par[0];
    Double_t m0 = par[1]; // constrain: same mean value parameter for both single gaussians
    Double_t sigma_1 = par[2];

    // Defining the custom function
    Double_t result = C_1 * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma_1), 2));
    return result;
}
// Custom pure signal fit function: signal components (amplitude,mean,sigma) = (0,2,3) and (1,2,4)
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
// Custom pure reflections fit function: reflection components (amplitude,mean,sigma) = (0,2,4) and (1,3,5)
Double_t pureReflectionsFunction(Double_t* x, Double_t* par) {

    Double_t m = x[0];
    Double_t C_1 = par[0];      // total parameter index: 5
    Double_t C_2 = par[1];      // total parameter index: 6
    Double_t m0_1 = par[2];     // total parameter index: 7
    Double_t m0_2 = par[3];     // total parameter index: 8
    Double_t sigma_1 = par[4];  // total parameter index: 9
    Double_t sigma_2 = par[5];  // total parameter index: 10

    // Defining the custom function
    Double_t result = C_1 * TMath::Exp(-TMath::Power((m - m0_1) / (2 * sigma_1), 2)) + C_2 * TMath::Exp(-TMath::Power((m - m0_2) / (2 * sigma_2), 2));
    return result;
}
// Custom signal and reflections fit function
Double_t signalAndReflectionsFunction(Double_t* x, Double_t* par) {
    // Gaussian components (amplitude,mean,sigma):
    // signal primary: (0,2,3)
    // signal secondary: (1,2,4)
    // reflection primary: reflection(0,2,4) =  total(5,7,9)
    // reflection secondary: reflection(1,3,5) =  total(6,8,10)
    
    return pureSignalFunction(x,par) + pureReflectionsFunction(x,&par[5]);
}
//-----------------------------------------------------------------------------Main workflow functions---------------------------------------------------------------------------------------------------------
// Module to create TH2D histograms including interest variable
HistogramGroup createHistograms(const std::vector<double>& ptDBinEdges, int xbins, double xmin, double xmax, const std::vector<double>& deltaRBinEdges) {

    // create histograms containers for each case (each container has a number of histograms corresponding to the invariant mass intervals)
    HistogramGroup histograms;

    for (size_t i = 0; i < ptDBinEdges.size() - 1; ++i) {
        // So that the title adapts to fractional binning title
        if (std::fmod(ptDBinEdges[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(ptDBinEdges[i+1], 1.0) != 0) {
                // pure reflections histograms
                histograms.reflections.push_back(new TH2D(Form("histMass_reflections%zu",i+1), Form("%.1f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // pure signals histograms
                histograms.signals.push_back(new TH2D(Form("histMass_signals%zu",i+1), Form("%.1f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // signals+reflections histograms: calculate ratios
                histograms.signals_and_reflections.push_back(new TH2D(Form("histMass_signals_and_reflections%zu",i+1), Form("%.1f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
        
            } else {
                // pure reflections histograms
                histograms.reflections.push_back(new TH2D(Form("histMass_reflections%zu",i+1), Form("%.1f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // pure signals histograms
                histograms.signals.push_back(new TH2D(Form("histMass_signals%zu",i+1), Form("%.1f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // signals+reflections histograms: calculate ratios
                histograms.signals_and_reflections.push_back(new TH2D(Form("histMass_signals_and_reflections%zu",i+1), Form("%.1f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
        
            }
        } else {
            if (std::fmod(ptDBinEdges[i+1], 1.0) != 0) {
                // pure reflections histograms
                histograms.reflections.push_back(new TH2D(Form("histMass_reflections%zu",i+1), Form("%.0f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // pure signals histograms
                histograms.signals.push_back(new TH2D(Form("histMass_signals%zu",i+1), Form("%.0f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // signals+reflections histograms: calculate ratios
                histograms.signals_and_reflections.push_back(new TH2D(Form("histMass_signals_and_reflections%zu",i+1), Form("%.0f < p_{T,D} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
        
            } else {
                // pure reflections histograms
                histograms.reflections.push_back(new TH2D(Form("histMass_reflections%zu",i+1), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // pure signals histograms
                histograms.signals.push_back(new TH2D(Form("histMass_signals%zu",i+1), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // signals+reflections histograms: calculate ratios
                histograms.signals_and_reflections.push_back(new TH2D(Form("histMass_signals_and_reflections%zu",i+1), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
        
            }
        }
        
        // pure reflections histograms
        histograms.reflections[i]->Sumw2();

        // pure signals histograms
        histograms.signals[i]->Sumw2();
        
        // signals+reflections histograms: calculate ratios
        histograms.signals_and_reflections[i]->Sumw2();
        
    }

    std::cout << "Histograms objects created." << std::endl;
    return histograms;
}

// Find the appropriate cut value based on pT
double GetBkgProbabilityCut(double pT, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    for (size_t i = 0; i < bdtPtCuts.size() - 1; ++i) {
        if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) {
            return bdtPtCuts[i].second;
        }
    }
    return 1.0; // Default: accept all if out of range
}

// Module to fill 2D histograms from TTree data
void fillHistograms(TFile* fInputMC, HistogramGroup& histograms, const double& jetptMin, const double& jetptMax, std::vector<double>& ptDBinEdges, std::vector<double>& deltaRBinEdges, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = deltaRBinEdges[deltaRBinEdges.size() - 1]; // on particle level delta R

    // Accessing detector level data TTree
    TTree* tree = (TTree*)fInputMC->Get("DF_merged/O2mcdjetdisttable");

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
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= deltaRBinEdges[0]) && (deltaR < DeltaRcut))) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < ptDBinEdges.size() - 1 && !filled; iEdge++) {
                if ((hfPt >= ptDBinEdges[iEdge]) && (hfPt < ptDBinEdges[iEdge + 1])) {

                    // D0 = +1, D0bar = -1, neither = 0
                    if ((hfMatchedFrom != 0) && (hfSelectedAs != 0)) {
                        if (hfMatchedFrom == hfSelectedAs) {
                            
                            // Get the threshold for this pT range
                            double maxBkgProb = GetBkgProbabilityCut(hfPt, bdtPtCuts);
                            // Fill histogram only if the cut is passed
                            if (hfMlScore0 < maxBkgProb) {
                                // pure signals
                                histograms.signals[iEdge]->Fill(hfMass, deltaR); // 2D case: (hfMass, deltaR);
                            }

                        } else {

                            // Get the threshold for this pT range
                            double maxBkgProb = GetBkgProbabilityCut(hfPt, bdtPtCuts);
                            // Fill histogram only if the cut is passed
                            if (hfMlScore0 < maxBkgProb) {
                                // pure reflections
                                histograms.reflections[iEdge]->Fill(hfMass, deltaR); // 2D case: (hfMass, deltaR);
                            }

                        }
                        
                        // Get the threshold for this pT range
                        double maxBkgProb = GetBkgProbabilityCut(hfPt, bdtPtCuts);
                        // Fill histogram only if the cut is passed
                        if (hfMlScore0 < maxBkgProb) {
                            // signals and reflections altogether, withOUT "neither = 0" entries
                            histograms.signals_and_reflections[iEdge]->Fill(hfMass, deltaR); // 2D case: (hfMass, deltaR);
                        }
                    }
                    // signals and reflections altogether, wiTH "neither = 0" entries (should I include "neither" entries too? Altough they should not appear)
                    //histograms.signals_and_reflections[iEdge]->Fill(hfMass, deltaR);
                    filled = true; // Exit the loop once the correct histogram is found
                }
                
            } // pT,D chosen
            
        } // end of TTree entries loop

        
        
        
    }
    std::cout << "Histograms filled." << std::endl;
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
               "reflections", new TF1("reflectionsFit", pureReflectionsFunction, xmin, xmax, 6));

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
        signalFit->SetParLimits(0, 0., TMath::Infinity());
        signalFit->SetParLimits(1, 0., TMath::Infinity());
        //signalFit->SetParLimits(3, 0.5 * sigma1, 2.0 * sigma1); // Example constraint for sigma1

        histograms.signals_1d[iInterval]->Fit(signalFit, "RQ"); // "Q" option performs quiet fit without drawing the fit function
        fits.signals.push_back(signalFit);
        //fits.signals[iInterval]->Print("V");

        // Also create and store individual components
        TF1* fSignalPrimary = new TF1(Form("signalPrimaryFit_%zu", iInterval), singleGaussianFunction, xmin, xmax, 3);
        fSignalPrimary->SetParameters(fits.signals[iInterval]->GetParameter(0), fits.signals[iInterval]->GetParameter(2), fits.signals[iInterval]->GetParameter(3)); // Amplitude, mean, sigma
        TF1* fSignalSecondary = new TF1(Form("signalSecondaryFit_%zu", iInterval), singleGaussianFunction, xmin, xmax, 3);
        fSignalSecondary->SetParameters(fits.signals[iInterval]->GetParameter(1), fits.signals[iInterval]->GetParameter(2), fits.signals[iInterval]->GetParameter(4)); // Amplitude, mean, sigma
        fits.individualSignals.push_back(std::make_pair(fSignalPrimary,fSignalSecondary));
    }

    // Reflections fits
    for (size_t iInterval = 0; iInterval < ptDBinEdges.size() - 1; iInterval++) {
        TF1* reflectionsFit = new TF1(Form("reflectionsFit_%zu", iInterval), pureReflectionsFunction, xmin, xmax, 6);

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
        reflectionsFit->SetParName(0, "C1");
        reflectionsFit->SetParName(1, "C2");
        reflectionsFit->SetParName(2, "m0_1");
        reflectionsFit->SetParName(3, "m0_2");
        reflectionsFit->SetParName(4, "sigma1");
        reflectionsFit->SetParName(5, "sigma2");
        reflectionsFit->SetParLimits(0, 0., TMath::Infinity());
        reflectionsFit->SetParLimits(1, 0., TMath::Infinity());

        histograms.reflections_1d[iInterval]->Fit(reflectionsFit, "RQ");
        fits.reflections.push_back(reflectionsFit);
        fits.reflections[iInterval]->Print("V");
        std::cout << "A1/A2 = " << fits.reflections[iInterval]->GetParameter(0) / fits.reflections[iInterval]->GetParameter(1) << std::endl;

        // Also create and store individual components
        TF1* fReflectionPrimary = new TF1(Form("reflectionPrimaryFit_%zu", iInterval), singleGaussianFunction, xmin, xmax, 3);
        fReflectionPrimary->SetParameters(fits.reflections[iInterval]->GetParameter(0), fits.reflections[iInterval]->GetParameter(2), fits.reflections[iInterval]->GetParameter(4)); // Amplitude, mean, sigma
        TF1* fReflectionSecondary = new TF1(Form("reflectionSecondaryFit_%zu", iInterval), singleGaussianFunction, xmin, xmax, 3);
        fReflectionSecondary->SetParameters(fits.reflections[iInterval]->GetParameter(1), fits.reflections[iInterval]->GetParameter(3), fits.reflections[iInterval]->GetParameter(5)); // Amplitude, mean, sigma
        fits.individualReflections.push_back(std::make_pair(fReflectionPrimary,fReflectionSecondary));
    }

    // Combined fits: true signal + reflections
    for (size_t iInterval = 0; iInterval < ptDBinEdges.size() - 1; iInterval++) {
        TF1* combinedFit = new TF1(Form("combinedFit_%zu", iInterval), signalAndReflectionsFunction, xmin, xmax, 11);

        // Obtain 1D mass histogram projection
        histograms.signals_and_reflections_1d.push_back( histograms.signals_and_reflections[iInterval]->ProjectionX(Form("histMass_signals_and_reflections_1d_%zu", iInterval+1), 1, -1) );

        // Initialize parameters from previous fits
        combinedFit->SetParameters(fits.signals[iInterval]->GetParameter(0), fits.signals[iInterval]->GetParameter(1), fits.signals[iInterval]->GetParameter(2), fits.signals[iInterval]->GetParameter(3), fits.signals[iInterval]->GetParameter(4), 
                                      fits.reflections[iInterval]->GetParameter(0), fits.reflections[iInterval]->GetParameter(1), fits.reflections[iInterval]->GetParameter(2), fits.reflections[iInterval]->GetParameter(3), fits.reflections[iInterval]->GetParameter(4), fits.reflections[iInterval]->GetParameter(5));
        combinedFit->FixParameter(0,fits.signals[iInterval]->GetParameter(0));
        combinedFit->FixParameter(1,fits.signals[iInterval]->GetParameter(1));
        combinedFit->FixParameter(2,fits.signals[iInterval]->GetParameter(2));
        combinedFit->FixParameter(3,fits.signals[iInterval]->GetParameter(3));
        combinedFit->FixParameter(4,fits.signals[iInterval]->GetParameter(4));
        combinedFit->FixParameter(5,fits.reflections[iInterval]->GetParameter(0));
        combinedFit->FixParameter(6,fits.reflections[iInterval]->GetParameter(1));
        combinedFit->FixParameter(7,fits.reflections[iInterval]->GetParameter(2));
        combinedFit->FixParameter(8,fits.reflections[iInterval]->GetParameter(3));
        combinedFit->FixParameter(9,fits.reflections[iInterval]->GetParameter(4));
        combinedFit->FixParameter(10,fits.reflections[iInterval]->GetParameter(5));
        histograms.signals_and_reflections_1d[iInterval]->Fit(combinedFit, "RQ");
        fits.signals_and_reflections.push_back(combinedFit);
        //fits.signals_and_reflections[iInterval]->Print("V");
    }
    
    std::cout << "Fits performed." << std::endl;

    return fits;
};

void PlotHistograms(const HistogramGroup& histograms, FitsGroup& fits, const double jetptMin, const double jetptMax) {
    
    // creating 1D mass projection histograms
    TH1D* tempHist;
    //std::vector<TH1D*> histograms_1d;

    // Defining canvases before plotting
    int nHistos = histograms.signals.size();
    // Start with a square layout (or close to it)
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));
    TCanvas* cSignals = new TCanvas("cSignals","Pure signal histograms");
    //cSignals->SetCanvasSize(1800,1000);
    //cSignals->Divide(3,static_cast<int>(histograms.signals.size() / 3)); // columns, lines
    cSignals->Divide(nCols,nRows);
    TCanvas* cReflections = new TCanvas("cReflections","Pure reflections histograms");
    //cReflections->SetCanvasSize(1800,1000);
    //cReflections->Divide(3,static_cast<int>(histograms.signals.size() / 3)); // columns, lines
    cReflections->Divide(nCols,nRows);
    TCanvas* cSignalsAndReflections = new TCanvas("cSignalsAndReflections","Signal and reflections histograms");
    //cSignalsAndReflections->SetCanvasSize(1800,1000);
    //cSignalsAndReflections->Divide(3,static_cast<int>(histograms.signals.size() / 3)); // columns, lines
    cSignalsAndReflections->Divide(nCols,nRows);

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
        tempHist->SetMinimum(0);
        tempHist->Draw();
        fits.signals[iHist]->SetLineColor(8); // pastel green
        fits.signals[iHist]->Draw("same");
        fits.individualSignals[iHist].first->SetLineColor(8);
        fits.individualSignals[iHist].second->SetLineColor(8);
        fits.individualSignals[iHist].first->SetLineStyle(kDashed);
        fits.individualSignals[iHist].second->SetLineStyle(kDashed);
        fits.individualSignals[iHist].first->Draw("same");
        fits.individualSignals[iHist].second->Draw("same");
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
        latex->DrawLatex(statBoxPos-0.87, 0.65,"JE_HF_LHC24g5_All_D0"); // Derived data from LHC24d3a (MC)
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
        fits.individualReflections[iHist].first->SetLineColor(46);
        fits.individualReflections[iHist].first->SetLineStyle(kDashed);
        fits.individualReflections[iHist].first->Draw("same");
        fits.individualReflections[iHist].second->SetLineColor(46);
        fits.individualReflections[iHist].second->SetLineStyle(kDashed);
        fits.individualReflections[iHist].second->Draw("same");
        latex->DrawLatex(statBoxPos-0.87, 0.40,"MC derived data:"); // Derived data from LHC24d3a (MC)
        latex->DrawLatex(statBoxPos-0.87, 0.35,"JE_HF_LHC24g5_All_D0"); // Derived data from LHC24d3a (MC)
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
        latex->DrawLatex(statBoxPos-0.87, 0.65,"JE_HF_LHC24g5_All_D0"); // Derived data from LHC24d3a (MC)
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

    // Draw individual gaussian contributions with total fit for each plot
    TCanvas* cIndividualGaussians = new TCanvas("cIndividualGaussians","Individual gaussians with total fit");
    cIndividualGaussians->Divide(nCols,nRows);
    // Loop over pT,D intervals/histograms
    for (size_t iHist = 0; iHist < histograms.signals_and_reflections.size(); iHist++) {
        tempHist = histograms.signals_and_reflections[iHist]->ProjectionX(Form("h_mass_signals_and_reflections_proj_%zu", iHist));
        gStyle->SetOptStat(0); // Turn off the default stats box
        tempHist->SetMarkerStyle(kDot);
        tempHist->SetMarkerColor(kBlack);
        tempHist->SetLineColor(kBlack); // kGray
        tempHist->GetYaxis()->SetTitle("counts");
        cIndividualGaussians->cd(iHist+1);
        tempHist->Draw();

        fits.individualSignals[iHist].first->Draw("same");
        fits.individualSignals[iHist].second->Draw("same");
        fits.individualReflections[iHist].first->Draw("same");
        fits.individualReflections[iHist].second->Draw("same");
        fits.signals_and_reflections[iHist]->Draw("same");
    }

    //
    // Storing images
    //
    TString imagePath = "../../Images/1-SignalTreatment/Reflections/";
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
    cSignals->Print(Form(imagePath + "reflections_fits_%s.pdf(",jetPtRange.Data()));
    cReflections->Print(Form(imagePath + "reflections_fits_%s.pdf",jetPtRange.Data()));
    cSignalsAndReflections->Print(Form(imagePath + "reflections_fits_%s.pdf)",jetPtRange.Data()));
    
    std::cout << "Histograms plotted." << std::endl;
}

// Module to save histograms and fits to output file
void SaveData(TFile* fOutput, const HistogramGroup& histograms, FitsGroup& fits, 
              const std::vector<double>& deltaRBinEdges, const std::vector<double>& ptjetBinEdges, const std::vector<double>& ptDBinEdges) {
    // Save histograms to output file
    fOutput->cd();
    for (size_t i = 0; i < histograms.signals.size(); i++) {

        // Storing reflection histograms
        histograms.signals[i]->Write();
        histograms.reflections[i]->Write();
        histograms.signals_and_reflections[i]->Write();

        // Storing reflection fits
        fits.signals[i]->Write();
        fits.reflections[i]->Write();
        fits.signals_and_reflections[i]->Write();
    }
    std::cout << "Reflection histograms and fits stored in file " << fOutput->GetName() << "." << std::endl;

    // Also store the axes used for the histograms
    // Create a directory for axes
    fOutput->mkdir("axes");
    fOutput->cd("axes");
    // Create TVectorD with same content
    TVectorD vecDeltaR(deltaRBinEdges.size());
    for (size_t i = 0; i < deltaRBinEdges.size(); ++i) {
        vecDeltaR[i] = deltaRBinEdges[i];
    }
    vecDeltaR.Write("deltaRBinEdges");
    TVectorD vecPtJet(ptjetBinEdges.size());
    for (size_t i = 0; i < ptjetBinEdges.size(); ++i) {
        vecPtJet[i] = ptjetBinEdges[i];
    }
    vecPtJet.Write("ptjetBinEdges");
    TVectorD vecPtD(ptDBinEdges.size());
    for (size_t i = 0; i < ptDBinEdges.size(); ++i) {
        vecPtD[i] = ptDBinEdges[i];
    }
    vecPtD.Write("ptDBinEdges");

    // Return to root directory (optional)
    fOutput->cd();
}

void JetPtIterator(const double jetptMin, const double jetptMax, const std::vector<double>& ptjetBinEdges, const bool useEmmaYeatsBins) {

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;

    // pT,jet cuts
    //std::vector<double> ptjetBinEdges = {5., 7., 15., 30., 50.};
    //const double jetptMin = ptjetBinEdges[0]; // GeV
    //const double jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 1]; // GeV
    //const double jetptMin = 5.; // GeV
    //const double jetptMax = 7.; // GeV
    // DeltaR bins: default = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.5}
    std::vector<double> deltaRBinEdges = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5};
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
    // // pT,D bins
    std::vector<double> ptDBinEdges = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.}; // 1., 2., 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 8., 10., 12., 15., 30.}
    // Emma Yeats reported analysis bins:
    if (useEmmaYeatsBins) {
        deltaRBinEdges = {0., 0.01, 0.03, 0.05, 0.12, 0.2}; // originally {0., 0.01, 0.03, 0.05, 0.12}, need to add higher bin for unfolding step
        minDeltaR = deltaRBinEdges[0];
        maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
        ptDBinEdges = {5., 6., 7., 8., 9., 10., 12., 20.};
    }
    
    
    
    // default = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.}
    // fractional = {3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 15., 30.}
    // low bins = {1., 2., 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 15., 30.}
    // higher bin size of last one: {1., 2., 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 8., 10., 12., 15., 40.}
    // bins with plenty of statistics: 3-4, 4-5, 5-6, 6-7, 8-10
    // wide peak bins: 8-10, 10-12, 12-15, 15-30 (very wide)
    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.03 (from first of pair)
    //std::vector<std::pair<double, double>> bdtPtCuts = {
    //    {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}
    //};
    // Dataset: JE_HF_LHC24g5_All_D0
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.4}, {8, 0.6}, {10, 0.8}, {12, 0.8}, {16, 1.0}, {50, 1.0}
    };

    // mass histogram
    int massBins = 50; // default=100 
    double minMass = 1.72; // use from 1.72, used to use 1.67
    double maxMass = 2.06; // old = 2.1, new = 2.05 + 0.01

    // Opening data file
    //TFile* fInputMC = new TFile("../../SimulatedData/Hyperloop_output/McEfficiency/New_with_reflections/Merged_AO2D_HF_LHC24d3a_All.root","read"); // previous file used
    TFile* fInputMC = new TFile("../../SimulatedData/Hyperloop_output/Train_runs/410602_Eff/AO2D_mergedDFs.root","read");
    if (!fInputMC || fInputMC->IsZombie()) {
        std::cerr << "Error: Unable to open the input ROOT file." << std::endl;
    }

    TFile* fOutput = new TFile(Form("reflections_%.0f_to_%.0fGeV.root", jetptMin, jetptMax),"recreate");
    if (!fOutput || fOutput->IsZombie()) {
        std::cerr << "Error: Unable to open the output ROOT file." << std::endl;
    }

    // Create multiple histograms
    HistogramGroup histograms = createHistograms(ptDBinEdges,                                 // the pT,D edges will determine the number of mass histograms
                                                     massBins, minMass, maxMass,              // mass histograms binning
                                                     deltaRBinEdges);                         // deltaR histograms with asymmetrical bin widths

    // Fill histograms
    fillHistograms(fInputMC, histograms, jetptMin, jetptMax, ptDBinEdges, deltaRBinEdges, bdtPtCuts);

    // Perform fits
    //std::vector<TF1*> fittings = performFit(histograms, parametersVectors, minMass, maxMass);
    FitsGroup fits = performFits(histograms, ptDBinEdges, minMass, maxMass);

    // signal/side-band region parameters
    int signalSigmas = 2; // 2
    int startingBackSigma = 4; // 4
    int backgroundSigmas = 4; // 4
    
    // Plot histograms
    PlotHistograms(histograms, fits, jetptMin, jetptMax);
    
    std::cout << "Sanity check on loaded bin edges:" << std::endl;
    std::cout << "pT,jet bin edges: ";
    for (const auto& edge : ptjetBinEdges) {
        std::cout << edge << " ";
    }
    std::cout << std::endl;
    std::cout << "DeltaR bin edges: ";
    for (const auto& edge : deltaRBinEdges) {
        std::cout << edge << " ";
    }
    std::cout << std::endl;
    std::cout << "pT,D bin edges: ";
    for (const auto& edge : ptDBinEdges) {
        std::cout << edge << " ";
    }
    std::cout << std::endl;

    // Storing final histograms to output file
    SaveData(fOutput, histograms, fits, deltaRBinEdges, ptjetBinEdges, ptDBinEdges);
    
    // Clean current iteration data in order to avoid memory leaks
    //ClearData(outputStruct, massHistograms);

    
}

void ReflectionsTreatment() {
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // pT,jet cuts
    std::vector<double> ptjetBinEdges = {5., 7., 10., 15., 30., 50.}; // default = {5., 7., 15., 30., 50.}
    //std::vector<double> ptjetBinEdges = {30., 50.}; // default = {5., 7., 15., 30., 50.}
    // Emma Yeats reported analysis bins:
    bool useEmmaYeatsBins = false;
    if (useEmmaYeatsBins) {
        ptjetBinEdges = {5., 10., 15., 20., 30.}; // originally {10., 15., 20.}, need to add lower and higher bins for unfolding step
    }

    for (size_t iJetPt = 0; iJetPt < ptjetBinEdges.size() - 1; iJetPt++) {
        
        // Apply side-band method to each pT,jet bin
        double jetptMin = ptjetBinEdges[iJetPt];
        double jetptMax = ptjetBinEdges[iJetPt+1];

        JetPtIterator(jetptMin, jetptMax, ptjetBinEdges, useEmmaYeatsBins);
    }

    // Compute the entire range too
    double jetptMin = ptjetBinEdges[0];
    double jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 1];
    JetPtIterator(jetptMin, jetptMax, ptjetBinEdges, useEmmaYeatsBins);

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



