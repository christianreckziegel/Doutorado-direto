/**
 * @file BackgroundSubtraction.C
 * @author Christian Reckziegel
 * @brief Macro for executing the background subtraction procedure of D0 candidates using the side-band technique
**/

using namespace std;

double DeltaPhi(double phi1, double phi2) {
    // Compute the absolute difference between phi1 and phi2
    double dphi = std::abs(phi1 - phi2); 
    if (dphi > M_PI) {
        // subtract 2pi if the difference if bigger than pi
        dphi = dphi - 2*M_PI;
    }

    return dphi;
}

//__________________________________________________________________________________________________________________________
// Fit functions
// Custom fit functions for full fit, background, signal, and reflection components
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

//__________________________________________________________________________________________________________________________
// Module to create TH2D histograms including interest variable: VARIABLE bin sizes
std::vector<std::vector<TH2D*>> createMassHistograms(const std::vector<double>& ptDBinEdges, int xbins, double xmin, double xmax, const std::vector<double>& yBinEdges, const int BDTsteps) {
    std::vector<std::vector<TH2D*>> histograms;

    for (size_t iptDBin = 0; iptDBin < ptDBinEdges.size() - 1; ++iptDBin) {
        // Create vector of histograms for each pT,D bin
        std::vector<TH2D*> histograms_pTD;
        for (size_t iBDTstep = 0; iBDTstep < BDTsteps; iBDTstep++) {
            histograms_pTD.push_back(new TH2D(Form("histMass_%zupTD_%zuBDTstep",iptDBin, iBDTstep), Form("%.0f < p_{T,D} < %.0f GeV/c, BDT score %.2f;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[iptDBin],ptDBinEdges[iptDBin+1], static_cast<double>(iBDTstep+1) / BDTsteps), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
            
            
        }
        histograms.push_back(histograms_pTD);
    }
    
    

    return histograms;
}

// Module to create TH1D histograms for significance values storage
std::vector<TH1D*> createSignificanceHistograms(const std::vector<double>& ptDBinEdges, const int BDTsteps) {
    int binNumber = BDTsteps;
    double binWidth = 1. / BDTsteps;
    double minBin = 0.; // 0. + 1. * binWidth
    double maxBin = 1.; // 1. + 1. * binWidth
    // Create histograms for significance values
    std::vector<TH1D*> histograms;
    for (size_t iptDBin = 0; iptDBin < ptDBinEdges.size() - 1; ++iptDBin) {
        // Create vector of histograms for each pT,D bin
        histograms.push_back(new TH1D(Form("histSignificance_%zupTD",iptDBin), Form("%.0f < p_{T,D} < %.0f GeV/c;BDT cut_{bkg}^{maximum}; Significance = #frac{S}{#sqrt{S+B}}",ptDBinEdges[iptDBin],ptDBinEdges[iptDBin+1]), BDTsteps, minBin, maxBin));
        histograms.back()->Sumw2();
    }
    
    return histograms;
}

// Module to create BDT background histograms
std::vector<TH1D*> createBkgHistograms(const std::vector<double>& ptDBinEdges, const int BDTsteps) {
    // Create histograms for BDT background values
    std::vector<TH1D*> histograms;
    for (size_t iptDBin = 0; iptDBin < ptDBinEdges.size() - 1; ++iptDBin) {
        // Create vector of histograms for each pT,D bin
        histograms.push_back(new TH1D(Form("histBkg_%zupTD",iptDBin), Form("%.0f < p_{T,D} < %.0f GeV/c;BDT cut_{bkg}^{maximum}; BDT_{bkg}",ptDBinEdges[iptDBin],ptDBinEdges[iptDBin+1]), BDTsteps * 100, 0., 1.));
        //histograms.back()->Sumw2();
    }
    return histograms;
}

void fillBDTHistograms(TFile* fDist, TH1D* histogramBkg, const double& jetptMin, const double& jetptMax, const double& ptDMin, const double& ptDMax, size_t BDTsteps) {
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = 0.4; // on particle level delta R

    // Accessing TTree
    TTree* tree = (TTree*)fDist->Get("DF_merged/O2jetdisttable");

    // Assuming histograms and tree data correspond in some way
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing the TTree
    float axisDistance, jetPt, jetEta, jetPhi;
    float hfPt, hfEta, hfPhi, hfMass, hfY, hfMlScore0, hfMlScore1, hfMlScore2;
    int jetNConst;

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
    tree->SetBranchAddress("fHfMlScore0",&hfMlScore0); // background ML score
    tree->SetBranchAddress("fHfMlScore1",&hfMlScore1); // prompt D0 ML score
    tree->SetBranchAddress("fHfMlScore2",&hfMlScore2); // non-prompt D0 ML score

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));
        
        // Fill each histogram with their respective pT intervals
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= 0.) && (deltaR < DeltaRcut))) {
            
            // Apply pT,D range cut
            if ((hfPt >= ptDMin) && (hfPt < ptDMax)) {

                histogramBkg->Fill(hfMlScore0);
            }
            
        }

        
        
        
    }
}
//__________________________________________________________________________________________________________________________
// Module to fill 2D histograms from TTree data
void fillMassHistogram(TFile* fDist, TH2D* histogram2d, const double& jetptMin, const double& jetptMax, const double& ptDMin, const double& ptDMax, const double& maxBkgProb) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = 0.4; // on particle level delta R

    // Accessing TTree
    TTree* tree = (TTree*)fDist->Get("DF_merged/O2jetdisttable");

    // Assuming histograms and tree data correspond in some way
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing the TTree
    float axisDistance, jetPt, jetEta, jetPhi;
    float hfPt, hfEta, hfPhi, hfMass, hfY, hfMlScore0, hfMlScore1, hfMlScore2;
    int jetNConst;

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
    tree->SetBranchAddress("fHfMlScore0",&hfMlScore0); // background ML score
    tree->SetBranchAddress("fHfMlScore1",&hfMlScore1); // prompt D0 ML score
    tree->SetBranchAddress("fHfMlScore2",&hfMlScore2); // non-prompt D0 ML score

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));
        
        // Fill each histogram with their respective pT intervals
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= 0.) && (deltaR < DeltaRcut))) {
            
            // Apply pT,D range cut
            if ((hfPt >= ptDMin) && (hfPt < ptDMax)) {

                // Fill histogram only if the cut is passed
                if (hfMlScore0 < maxBkgProb) {
                    histogram2d->Fill(hfMass, deltaR);
                }
            }
            
        }

        
        
        
    }
    //std::cout << "Histogram filled.\n";
}

// Module to perform fits to histograms
struct FitContainer {
    std::vector<std::vector<TF1*>> fitBackgroundOnly;        // background only fits
    std::vector<std::vector<TF1*>> fitSignalOnly;            // signal only fits
    std::vector<std::vector<TF1*>> fitReflectionsOnly;       // reflections only fits
    std::vector<std::vector<TF1*>> fitTotal;                 // total fits

    // Constructor to initialize with given sizes
    FitContainer(size_t N, size_t M): // N = number of pT bins, M = number of BDT steps
        fitBackgroundOnly(N, std::vector<TF1*>(M, nullptr)),
        fitSignalOnly(N, std::vector<TF1*>(M, nullptr)),
        fitReflectionsOnly(N, std::vector<TF1*>(M, nullptr)),
        fitTotal(N, std::vector<TF1*>(M, nullptr)) {}
};
void performFit(TFile* fReflectionsMC, FitContainer& fittings, TH2D* histogram2d, const double minMass, const double maxMass, size_t iPtDBin, size_t iBDTstep) {
    
    double m_0_reference = 1.86484; // D0 mass in GeV/c^2
    double sigma_reference = 0.012; // D0 width in GeV/c^2
    
    // Creating 1D mass projection histograms
    TH1D* histogram1d;
    histogram1d = histogram2d->ProjectionX(Form("h_mass_proj_%zupTD_%zuBDTstep", iPtDBin, iBDTstep));

    // ----Obtain MC fit parameters----

    // Get TF1 objects from MC file
    TF1* fSignal = (TF1*)fReflectionsMC->Get(Form("signalFit_%zu", iPtDBin));
    //fSignal->Print("V");
    TF1* fReflections = (TF1*)fReflectionsMC->Get(Form("reflectionsFit_%zu", iPtDBin));

    double A1Signal = fSignal->GetParameter(0); // A1Signal
    double A2Signal = fSignal->GetParameter(1); // A2Signal
    double m0Signal = fSignal->GetParameter(2); // m0
    double sigma1Signal = fSignal->GetParameter(3); // sigma1
    double sigma2Signal = fSignal->GetParameter(4); // sigma2
    //std::cout << "A1Signal = " << A1Signal << ", A2Signal = " << A2Signal << ", m0Signal = " << m0Signal << ", sigma1Signal = " << sigma1Signal << ", sigma2Signal = " << sigma2Signal << std::endl;

    double A1Reflection = fReflections->GetParameter(0); // A1Reflection
    double A2Reflection = fReflections->GetParameter(1); // A2Reflection
    double m0_1Reflections = fReflections->GetParameter(2); // m0_1
    double m0_2Reflections = fReflections->GetParameter(3); // m0_2
    double sigma1Reflections = fReflections->GetParameter(4); // sigma1
    double sigma2Reflections = fReflections->GetParameter(5); // sigma2
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
    fittings.fitTotal[iPtDBin][iBDTstep] = new TF1(Form("totalFit_histMass_%zupTD_%zuBDTstep", iPtDBin, iBDTstep), customFitFunction, minMass, maxMass, 13); // 12 parameters = 8 fixed + 5 free

    // Set initial values and fix parameters
    Double_t params[13] = {200000, -10.0, 1000, 1.5, m_0_reference, 0.02, 1.2, 2.0, 1.3, 1.83, 1.85, 0.02, 0.03};
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParameters(params);

    // Set parameter names
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(0, "Background A");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(1, "Background B");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(2, "A1 Signal");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(3, "A1/A2 Signal Ratio");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(4, "Signal Mean m0");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(5, "Signal Width Sigma1");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(6, "Sigma Ratio (Sig1/Sig2)");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(7, "A1Signal/A1Reflection Ratio");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(8, "A1/A2 Reflection Ratio");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(9, "Reflection Mean m0_1");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(10, "Reflection Mean m0_2");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(11, "Reflection Width Sigma_1");
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParName(12, "Reflection Width Sigma_2");

    // Fix the parameters from MC fits
    fittings.fitTotal[iPtDBin][iBDTstep]->FixParameter(3, A1toA2MCSignalRatio);                   // A1toA2MCSignalRatio
    fittings.fitTotal[iPtDBin][iBDTstep]->FixParameter(6, Sigma1toSigma2MCSignalRatio);           // Sigma1toSigma2MCSignalRatio
    fittings.fitTotal[iPtDBin][iBDTstep]->FixParameter(7, A1SignalToA1ReflectionMCRatio);         // A1SignalToA1ReflectionMCRatio
    fittings.fitTotal[iPtDBin][iBDTstep]->FixParameter(8, A1toA2MCReflectionsRatio);              // A1toA2MCReflectionsRatio
    fittings.fitTotal[iPtDBin][iBDTstep]->FixParameter(9, m0_1Reflections);                       // m0_1
    fittings.fitTotal[iPtDBin][iBDTstep]->FixParameter(10, m0_2Reflections);                      // m0_2
    fittings.fitTotal[iPtDBin][iBDTstep]->FixParameter(11, sigma1Reflections);                    // sigma_1
    fittings.fitTotal[iPtDBin][iBDTstep]->FixParameter(12, sigma2Reflections);                    // sigma_2

    // Apply range limits to the parameters
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParLimits(4, 0.95 * m_0_reference, 1.05 * m_0_reference);
    fittings.fitTotal[iPtDBin][iBDTstep]->SetParLimits(5, 0.35 * sigma_reference, 3.0 * sigma_reference);

    // Choose line color
    fittings.fitTotal[iPtDBin][iBDTstep]->SetLineColor(kBlack);

    // Perform fit with "Q" (quiet) option: no drawing of the fit function
    histogram1d->Fit(fittings.fitTotal[iPtDBin][iBDTstep], "Q");

    // Print the fit results
    //fittings.fitTotal[iPtDBin][iBDTstep]->Print("V");
    std::cout << "A1 Signal = " << fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(2) << std::endl;
    std::cout << "Signal Mean m0 = " << fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(4) << std::endl;
    std::cout << "Signal Width Sigma1 = " << fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(5) << std::endl;

    // Perform background only fit to each histogram
    // Getting total parameter values
    double a_par = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(0); // Get the value of parameter 'a'
    double b_par = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(1); // Get the value of parameter 'b'
    fittings.fitBackgroundOnly[iPtDBin][iBDTstep] = new TF1(Form("backgroundOnlyFit_%zupTD_%zuBDTstep", iPtDBin, iBDTstep), backgroundFunction, minMass, maxMass, 2);
    fittings.fitBackgroundOnly[iPtDBin][iBDTstep]->SetParameters(a_par,b_par); // fittings.size()-1 = the latest added to the vector
    fittings.fitBackgroundOnly[iPtDBin][iBDTstep]->SetLineStyle(kDashed);
    fittings.fitBackgroundOnly[iPtDBin][iBDTstep]->SetLineColor(kRed);

    
    //std::cout << "Background only fit performed.\n";

    // Perform signal only fit to each histogram
    // Extract signal-related parameters from the total fit
    double a1Signal = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(2);                      // A1 signal
    A1toA2MCSignalRatio = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(3);           // A2 signal
    m0Signal = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(4);                      // Signal mean m0
    double sigmaSignal = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(5);                   // Signal width sigma
    Sigma1toSigma2MCSignalRatio = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(6);   // Sigma1/Sigma2 signal
    
    // Perfoming fit with acquired parameters from total fit
    fittings.fitSignalOnly[iPtDBin][iBDTstep] = new TF1(Form("signalOnlyFit_%zupTD_%zuBDTstep", iPtDBin, iBDTstep), signalOnlyFunction, minMass, maxMass, 5);
    fittings.fitSignalOnly[iPtDBin][iBDTstep]->SetParameters(a1Signal, A1toA2MCSignalRatio, m0Signal, sigmaSignal, Sigma1toSigma2MCSignalRatio); // fittings.size()-1 = the latest added to the vector
    fittings.fitSignalOnly[iPtDBin][iBDTstep]->SetLineStyle(kDashed);
    fittings.fitSignalOnly[iPtDBin][iBDTstep]->SetLineColor(kBlue);

    //fittings.fitSignalOnly[iPtDBin][iBDTstep]->Print("V");
    //std::cout << "Signal only fit performed.\n";

    // Perform reflections only fit to each histogram
    // Extract reflection-related parameters from the total fit
    double A1SignalToA1ReflectionMCRatios = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(7);    // A1 signal / A1 reflection
    A1Signal = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(2);                          // A1 signal
    A1toA2MCReflectionsRatio = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(8);          // A1 reflection / A2 reflection
    double m0_1Reflection = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(9);                    // Reflection mean m0_1
    double m0_2Reflection = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(10);                   // Reflection mean m0_2
    double sigma1Reflection = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(11);                 // Reflection sigma_1
    double sigma2Reflection = fittings.fitTotal[iPtDBin][iBDTstep]->GetParameter(12);                 // Reflection sigma_2
    
    // Perfoming fit with acquired parameters from total fit
    fittings.fitReflectionsOnly[iPtDBin][iBDTstep] = new TF1(Form("reflectionOnlyFit_%zupTD_%zuBDTstep", iPtDBin, iBDTstep), reflectionOnlyFunction, minMass, maxMass, 7);
    fittings.fitReflectionsOnly[iPtDBin][iBDTstep]->SetParameters(A1SignalToA1ReflectionMCRatios, A1Signal, A1toA2MCReflectionsRatio, m0_1Reflection, m0_2Reflection, sigma1Reflection, sigma2Reflection); // fittings.size()-1 = the latest added to the vector
    fittings.fitReflectionsOnly[iPtDBin][iBDTstep]->SetLineStyle(kDashed);
    fittings.fitReflectionsOnly[iPtDBin][iBDTstep]->SetLineColor(kGreen);
    //fittings.fitReflectionsOnly[iPtDBin][iBDTstep]->Print("V");

   //std::cout << "Reflections only fits performed.\n";
    
    //std::cout << "Histogram fit performed.\n";
}

// Calculate significance using fit areas
double calculateSignificance(FitContainer& fittings, const size_t iPtDBin, const size_t iBDTstep, const size_t signalSigmas) {
    // Get signal yield by integrating Gaussian in a 2σ window
    double mean_signalOnly = fittings.fitSignalOnly[iPtDBin][iBDTstep]->GetParameter(2);
    double sigma1_signalOnly = fittings.fitSignalOnly[iPtDBin][iBDTstep]->GetParameter(3);
    //fittings.fitSignalOnly[iPtDBin][iBDTstep]->Print("V");
    double S = fittings.fitSignalOnly[iPtDBin][iBDTstep]->Integral(mean_signalOnly - signalSigmas*sigma1_signalOnly, mean_signalOnly + signalSigmas*sigma1_signalOnly);
    double B = fittings.fitBackgroundOnly[iPtDBin][iBDTstep]->Integral(mean_signalOnly - signalSigmas*sigma1_signalOnly, mean_signalOnly + signalSigmas*sigma1_signalOnly);
    
    // Binomial-like Approximation: when S is comparable to B, consider both signal and background fluctuations
    double significance = S / sqrt(S+B);
    
    return significance;
}

TH1D* optimalBDTCut(const std::vector<TH1D*>& significanceHistograms, const std::vector<double>& ptDBinEdges, const double& jetptMin, const double& jetptMax) {
    // Create a new histogram for the optimal BDT cut
    TH1D* optimalBDTHistogram = new TH1D("optimalBDTHistogram", "Optimal BDT cut^{maximum}_{bkg};p_{T,D} (GeV/c);Optimal BDT cut^{maximum}_{bkg}", ptDBinEdges.size() - 1, ptDBinEdges.data());
    
    // Loop over the bins of the input histogram
    for (int ipTD = 0; ipTD < significanceHistograms.size(); ++ipTD) {
        // Find BDT cut with maximum significance: get x-axis value with highest number of entries (x value searched from low to high values)
        double maxSignificance = significanceHistograms[ipTD]->GetMaximum();
        int maxBin = significanceHistograms[ipTD]->GetMaximumBin();
        double optimalBDT = significanceHistograms[ipTD]->GetBinLowEdge(maxBin);
        
        // Fill the optimal BDT histogram with the significance values
        int binNumber = optimalBDTHistogram->FindBin(ptDBinEdges[ipTD]);
        optimalBDTHistogram->SetBinContent(binNumber, optimalBDT);
        optimalBDTHistogram->SetBinError(binNumber, 0);
    }
    
    return optimalBDTHistogram;
}
void PlotHistograms(const std::vector<TH1D*>& significanceHistograms, const std::vector<TH1D*>& histogramsBkg, const std::vector<std::vector<TH2D*>>& histograms2d, TH1D* optimalBDTHistogram, const FitContainer& fittings, const double jetptMin, const double jetptMax) {
    
    // creating 1D mass projection canvases and histograms
    std::vector<TCanvas*> cMass;
    std::vector<std::vector<TH1D*>> hMass;

    // Initialize cMass with the same dimensions as histograms2d
    cMass.resize(histograms2d.size(), nullptr);
    hMass.resize(histograms2d.size());
    for (size_t iHisto = 0; iHisto < histograms2d.size(); iHisto++) {
        // Initialize with nullptr
        hMass[iHisto].resize(histograms2d[iHisto].size(), nullptr);

        // Divide each canvas in the number of BDT steps (one canvas per pT,D bin)
        cMass[iHisto] = new TCanvas(Form("cMass_%zupTD", iHisto), Form("cMass_%zupTD", iHisto), 800, 600);
        cMass[iHisto]->Divide(4, 3);
    }

    // obtaining 1D invariant mass histograms from the projection
    for (size_t iPtDBin = 0; iPtDBin < histograms2d.size(); iPtDBin++) {
        for (size_t iBDTstep = 0; iBDTstep < histograms2d[0].size(); iBDTstep++) {
            //hMass[iPtDBin][iBDTstep] = histograms2d[iPtDBin][iBDTstep]->ProjectionX(Form("h_mass_proj_1d_%zupTD_%zuBDTstep", iPtDBin, iBDTstep));
            hMass[iPtDBin][iBDTstep] = (TH1D*)histograms2d[iPtDBin][iBDTstep]->ProjectionX(Form("h_mass_proj_1d_%zupTD_%zuBDTstep", iPtDBin, iBDTstep))->Clone();
            // Calculate how many steps to skip (if BDTsteps > 10)
            size_t drawInterval = (histograms2d[iPtDBin].size() > 10) ? histograms2d[iPtDBin].size() / 10 : 1;
            if (iBDTstep % drawInterval == 0) {
                cMass[iPtDBin]->cd(iBDTstep / drawInterval + 1);
                hMass[iPtDBin][iBDTstep]->Draw();
                fittings.fitTotal[iPtDBin][iBDTstep]->Draw("same");
            }
            
            // if it's the last BDT step, draw the BDT histogram (corresponding to the pT,D bin) too
            if (iBDTstep == histograms2d[iPtDBin].size() - 1) {
                // Draw the last histogram with the fit
                cMass[iPtDBin]->cd(iBDTstep / drawInterval + 2);
                histogramsBkg[iPtDBin]->Draw();
            }
        }
    }

    TCanvas* cBDTBkg = new TCanvas("cBDTBkg", "cBDTBkg", 800, 600);
    cBDTBkg->Divide(4, 3);
    for (size_t iHisto = 0; iHisto < histogramsBkg.size(); iHisto++) {
        cBDTBkg->cd(iHisto + 1);
        histogramsBkg[iHisto]->Draw();
    }
    TCanvas* cBDTBkgSame = new TCanvas("cBDTBkgSame", "cBDTBkgSame", 800, 600);
    for (size_t iHisto = 0; iHisto < histogramsBkg.size(); iHisto++) {
        cBDTBkgSame->cd();
        
        if (iHisto == 0) {
            histogramsBkg[iHisto]->Draw();
        } else {
            histogramsBkg[iHisto]->Draw("same");
        }
        
    }

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.04); // default = 0.05

    // Create a canvas for plotting
    TCanvas* cSignificance = new TCanvas("cSignificance", "Significance", 800, 600);
    cSignificance->SetCanvasSize(1800,1000);
    cSignificance->Divide(4, 3);
    for (size_t iHisto = 0; iHisto < significanceHistograms.size(); iHisto++) {
        cSignificance->cd(iHisto+1);
        significanceHistograms[iHisto]->Draw();

    }
    
    TCanvas* cOptimalBDT = new TCanvas("cOptimalBDT", "Optimal BDT cut", 800, 600);
    cOptimalBDT->SetCanvasSize(1800,1000);
    cOptimalBDT->cd();
    optimalBDTHistogram->Draw();

    //
    // Storing images
    //
    TString imagePath = "../../Images/1-SignalTreatment/SideBand/";
    cSignificance->Update();
    cSignificance->SaveAs(imagePath + "Significance_per_pTD_bin.png");

    //
    // Storing in a single pdf file
    //
    cSignificance->Print(Form("significance_per_pTD_bin_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cBDTBkg->Print(Form("significance_per_pTD_bin_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cOptimalBDT->Print(Form("significance_per_pTD_bin_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    for (size_t iPtDBin = 0; iPtDBin < cMass.size(); iPtDBin++) {
        cMass[iPtDBin]->Update();
        if (iPtDBin == cMass.size() - 1) {
            cMass[iPtDBin]->Print(Form("significance_per_pTD_bin_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
        } else {
            cMass[iPtDBin]->Print(Form("significance_per_pTD_bin_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
        }
    }
    std::cout << "Histograms plotted and pdf stored." << std::endl;
}

void SaveData(const std::vector<std::vector<TH2D*>>& histograms2d, TH1D* optimalBDTHistogram, const std::vector<TH1D*>& significanceHistograms, const double& jetptMin, const double& jetptMax){
    // Open output file
    TFile* outFile = new TFile(Form("significance_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Unable to open the output ROOT file." << std::endl;
    }

    // creating 1D mass projection canvases and histograms
    std::vector<std::vector<TH1D*>> hMass;

    // Initialize cMass with the same dimensions as histograms2d
    hMass.resize(histograms2d.size());
    for (size_t iHisto = 0; iHisto < histograms2d.size(); iHisto++) {
        // Initialize with nullptr
        hMass[iHisto].resize(histograms2d[iHisto].size(), nullptr);
    }
    // obtaining 1D invariant mass histograms from the projection
    for (size_t iPtDBin = 0; iPtDBin < histograms2d.size(); iPtDBin++) {
        outFile->mkdir(Form("hMass_%zu", iPtDBin));
        outFile->cd(Form("hMass_%zu", iPtDBin)); // Set the current directory to the ROOT file 
        for (size_t iBDTstep = 0; iBDTstep < histograms2d[iPtDBin].size(); iBDTstep++) {
            hMass[iPtDBin][iBDTstep] = (TH1D*)histograms2d[iPtDBin][iBDTstep]->ProjectionX(Form("h_mass_1d_%zupTD_%zuBDTstep", iPtDBin, iBDTstep))->Clone();
            hMass[iPtDBin][iBDTstep]->Write();
        }
    }
     

    // Create a new directory for each histogram
    outFile->mkdir("Significance");
    outFile->cd("Significance"); // Set the current directory to the ROOT file
    // Loop over background subtracted histograms
    for (size_t iHisto = 0; iHisto < significanceHistograms.size(); iHisto++) {
        
        
        // store each histogram in file
        significanceHistograms[iHisto]->Write();
    }

    // Store the optimal BDT histogram
    outFile->cd(); // Move back to the root directory of the TFile
    optimalBDTHistogram->Write();

    // Close output file
    outFile->Close();
    delete outFile;
    
    std::cout << "Data saved to ROOT file." << std::endl;
}

void OptimalSignificanceBDTScore(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // jet pT cuts
    std::vector<double> ptjetBinEdges = {5., 7., 15., 30.};
    double jetptMin = ptjetBinEdges[0]; // GeV
    double jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 1]; // GeV
    // deltaR histogram
    int deltaRbins = 10000; // deltaRbins = numberOfPoints, default=10 bins for [0. 0.4]
    std::vector<double> deltaRBinEdges = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4}; // chosen by Nima
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
    // pT,D bins
    std::vector<double> ptDBinEdges = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};

    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.03 (from first of pair)
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}
    };

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;

    // mass histogram
    int massBins = 50; // default=100 
    double minMass = 1.72; // use from 1.72, used to use 1.67
    double maxMass = 2.1;

    // Opening data file
    TFile* fDist = new TFile("../../ExperimentalData/Hyperloop_output/HF_LHC23_pass4_Thin_small_2P3PDstar_DATA/AO2D.root","read");
    if (!fDist || fDist->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT data file." << std::endl;
    }
    // Opening
    TFile* fReflectionsMC = new TFile(Form("../Reflections/reflections_%.0f_to_%.0fGeV.root",jetptMin,jetptMax),"read");
    if (!fReflectionsMC || fReflectionsMC->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT reflections file." << std::endl;
    }

    // Number of steps for BDT score cut
    size_t BDTsteps = 100;
    
    // Create multiple histograms
    std::vector<std::vector<TH2D*>> histograms2d = createMassHistograms(ptDBinEdges,                            // the pT,D edges will determine the number of mass histograms
                                                                        massBins, minMass, maxMass,             // mass histograms binning
                                                                        deltaRBinEdges,                         // deltaR histograms with asymmetrical bin widths
                                                                        BDTsteps);                              // number of BDT score cut steps
    

    // Create histograms for BDT background
    std::vector<TH1D*> histogramsBkg = createBkgHistograms(ptDBinEdges, BDTsteps); // BDT background histograms for each pT,D bin

    // Create data container for fit functions
    FitContainer fittings(ptDBinEdges.size() -1, BDTsteps);
    std::vector<TH1D*> significanceHistograms = createSignificanceHistograms(ptDBinEdges, BDTsteps); // significance histograms for each pT,D bin

    // signal/side-band region parameters
    size_t signalSigmas = 2;
    size_t startingBackSigma = 4;
    size_t backgroundSigmas = 4;
    
    // once per pT,D bin
    // Create a histogram for each pT,D bin
    for (size_t iPtDBin = 0; iPtDBin < ptDBinEdges.size() -1; iPtDBin++) {

        double ptDMin = ptDBinEdges[iPtDBin];
        double ptDMax = ptDBinEdges[iPtDBin + 1];
        
        fillBDTHistograms(fDist, histogramsBkg[iPtDBin], jetptMin, jetptMax, ptDMin, ptDMax, BDTsteps);

        // once per BDT score cut value chosen
        for (size_t iBDTstep = 0; iBDTstep < BDTsteps; iBDTstep++) {
            
            // BDT score cut value
            double maxBkgProb = static_cast<double>(iBDTstep+1) / BDTsteps;

            // fill histogram with chosen BDT score cut value
            fillMassHistogram(fDist, histograms2d[iPtDBin][iBDTstep], jetptMin, jetptMax, ptDMin, ptDMax, maxBkgProb);

            // perform fit
            performFit(fReflectionsMC, fittings, histograms2d[iPtDBin][iBDTstep], minMass, maxMass, iPtDBin, iBDTstep);

            // calculate significance
            double significance = calculateSignificance(fittings, iPtDBin, iBDTstep, signalSigmas);
            
            std::cout << "Significance: " << significance << std::endl;
            // store significance in the histogram
            int binNumber = significanceHistograms[iPtDBin]->FindBin(maxBkgProb);
            significanceHistograms[iPtDBin]->SetBinContent(binNumber, significance);
            significanceHistograms[iPtDBin]->SetBinError(binNumber, 0);
            std::cout << "Bin number = " << binNumber << std::endl;

            std::cout << "Filled histogram for BDT step " << iBDTstep << ", BDT bkg maximum score cut = " << maxBkgProb << std::endl;
        }
        std::cout << "Filled histograms for pT,D range: " << ptDBinEdges[iPtDBin] << " < pT,D < " << ptDBinEdges[iPtDBin+1] << " GeV/c" << std::endl;
    }
    
    // Find optimal BDT cut
    TH1D* optimalBDTHistogram = optimalBDTCut(significanceHistograms, ptDBinEdges, jetptMin, jetptMax);

    // Plot histograms
    PlotHistograms(significanceHistograms, histogramsBkg, histograms2d, optimalBDTHistogram, fittings, jetptMin, jetptMax);
    
    // Storing final histograms to output file
    SaveData(histograms2d, optimalBDTHistogram, significanceHistograms, jetptMin, jetptMax);

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    OptimalSignificanceBDTScore();
    return 0;
}