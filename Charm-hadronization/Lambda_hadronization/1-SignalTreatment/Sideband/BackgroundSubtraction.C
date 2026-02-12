/**
 * Lambda_c hadron analysis
 * @file BackgroundSubtraction.C
 * @brief Background subtraction using sideband method
 * Input: Reflections.root: contain histograms and fit templates of Lambda_c signal and reflections
 * Outputs: BackgroundSubtraction.root
 * 
 * @author: Christian Reckziegel
 * Date: February 2026
 */

 #include "../../commonUtilities.h"
using namespace std;

//_________________________________Template for on-the-run chosen model (policy-based strategy)________________________________________________
struct PowerLawBackgroundPolicy {
    static double eval(double* x, double* par) {
        Double_t m = x[0];
        Double_t a = par[0];
        Double_t b = par[1];

        // Defining the custom function
        Double_t result = a * TMath::Power(m, b);
        return result;
    }
};
struct Poly2BackgroundPolicy {
    static double eval(double* x, double* par) {
        Double_t m = x[0];
        Double_t p0 = par[0];   // constant term
        Double_t p1 = par[1];   // linear term
        Double_t p2 = par[13];  // quadratic term (use a new index!)

        Double_t result = p2 + p1 * m + p0 * m * m;
        return result;
    }
};
struct SignalPolicy {
    static double eval(double* x, double* par) {
        Double_t m = x[0];
        Double_t A1Signal = par[2];                                         // Free parameter
        Double_t A1toA2MCSignalRatio = par[3];                              // Fixed parameter
        Double_t A2Signal = A1Signal / A1toA2MCSignalRatio;                 // Constrained to primary/secondary integral obtained from MC data fit
        Double_t m0 = par[4];                                               // Free parameter
        Double_t sigma1 = par[5];                                           // Free parameter
        Double_t Sigma1toSigma2MCSignalRatio = par[6];                      // Fixed parameter
        Double_t sigma2 = sigma1 / Sigma1toSigma2MCSignalRatio;             // Constrained to primary/secondary width obtained from MC data fit

        // Defining the custom function
        Double_t result = A1Signal * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma1), 2)) + A2Signal * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma2), 2));
        return result;
    }
};
struct SignalSinglePolicy { // need to implement in Reflections.C macro as well
    static double eval(double* x, double* par) {
        Double_t m = x[0];
        Double_t A1Signal = par[2];                                         // Free parameter
        Double_t m0 = par[4];                                               // Free parameter
        Double_t sigma1 = par[5];                                           // Free parameter
        
        // Defining the custom function
        Double_t result = A1Signal * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma1), 2));
        return result;
    }
};
struct ReflectionPolicy {
    static double eval(double* x, double* par) {
        // m0 and sigma values fixed from fits obtained from MC data fits
        Double_t m = x[0];
        Double_t A1SignalToA1ReflectionMCRatio = par[7];                    // Fixed parameter
        // Extract A1Signal from the parameter array (same index as in `signalFunction`)
        Double_t A1Signal = par[2];
        Double_t A1Reflection = A1Signal / A1SignalToA1ReflectionMCRatio;   // Constrained to signal/reflection integral obtained from MC data fits
        Double_t A1toA2MCReflectionsRatio = par[8];                         // Fixed parameter
        Double_t A2Reflection = A1Reflection / A1toA2MCReflectionsRatio;    // Constrained to primary/secondary integral obtained from MC data fit
        Double_t m0_1 = par[9];                                             // Fixed parameter
        Double_t m0_2 = par[10];                                            // Fixed parameter
        Double_t sigma_1 = par[11];                                         // Fixed parameter
        Double_t sigma_2 = par[12];                                         // Fixed parameter

        // Defining the custom function
        Double_t result = A1Reflection * TMath::Exp(-TMath::Power((m - m0_1) / (2 * sigma_1), 2)) + A2Reflection * TMath::Exp(-TMath::Power((m - m0_2) / (2 * sigma_2), 2));
        return result;
    }
};
// Compile time polymorphism via policy-based design for fit models
using FullModelPowerLaw = FitModel<SignalPolicy, ReflectionPolicy, PowerLawBackgroundPolicy>;               // Signal + Reflection + power law Background
using FullModelPoly2 = FitModel<SignalPolicy, ReflectionPolicy, Poly2BackgroundPolicy>;                     // Signal + Reflection + 2nd order polynomial Background
using SigRefModel = FitModel<SignalPolicy, ReflectionPolicy>;                                               // Signal + Reflection only
using FullSingleGaussianModel = FitModel<SignalSinglePolicy, ReflectionPolicy, PowerLawBackgroundPolicy>;   // Single Gaussian Signal + Reflection + power law Background
// Custom fit functions for background component
Double_t backgroundFunctionPoly2(Double_t* x, Double_t* par) {
    Double_t m = x[0];
    Double_t p0 = par[2];   // constant term
    Double_t p1 = par[1];   // linear term
    Double_t p2 = par[0];  // quadratic term (use a new index!)

    // Defining the custom function
    Double_t result = p0 + p1 * m + p2 * m * m;
    return result;
}
Double_t backgroundFunction(Double_t* x, Double_t* par) {
    Double_t m = x[0];
    Double_t a = par[0];
    Double_t b = par[1];

    // Defining the custom function
    Double_t result = a * TMath::Power(m, b);
    return result;
}
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

struct SubtractionResult {
    std::vector<TH1D*> sidebandHist;                                    // 1D sideband background histograms vector
    std::vector<TH1D*> signalHist;                                      // 1D signal histograms vector
    std::vector<TH1D*> subtractedHist;                                  // 1D sideband subtracted histograms vector
    TH1D* hSubtractedFullJetPt;                                         // final deltaR distribution for all pT,HF summed
    TH1D* hSignificance;                                                // final significance distribution for all pT,HF summed
    std::vector<std::pair<double, double>> signal_background_values;    // vector of signal and background values for each pT,HF bin
    std::vector<std::array<double, 2>> scalingFactorsArrays;           // alpha and beta scaling factors for each pT,HF bin
};
struct SidebandData {
    std::vector<TH2D*> histograms2d;
    std::vector<TH1D*> histograms1d;
    FitContainer fittings;
    SubtractionResult subtractionResults;
};

// Module to create TH2D histograms including interest variable: VARIABLE bin sizes
SidebandData createHistograms(const std::vector<double>& ptHFBinEdges, int xbins, double xmin, double xmax, const std::vector<double>& yBinEdges, const double& jetptMin) {
    // Change binning to low statistics range [30,50] GeV/c
    if (jetptMin >= 30.) {
        xbins = xbins / 2;
    }
    
    SidebandData dataContainer;
    for (size_t i = 0; i < ptHFBinEdges.size() - 1; ++i) {
        // So that the title adapts to fractional binning title
        if (std::fmod(ptHFBinEdges[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(ptHFBinEdges[i+1], 1.0) != 0) {
                dataContainer.histograms2d.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.1f < #it{p}_{T, #Lambda_{c}^{+}} < %.1f GeV/#it{c};#it{M}(#piKPr) (GeV/#it{c}^{2});#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            } else {
                dataContainer.histograms2d.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.1f < #it{p}_{T, #Lambda_{c}^{+}} < %.0f GeV/#it{c};#it{M}(#piKPr) (GeV/#it{c}^{2});#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            }
        } else {
            if (std::fmod(ptHFBinEdges[i+1], 1.0) != 0) {
                dataContainer.histograms2d.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < #it{p}_{T, #Lambda_{c}^{+}} < %.1f GeV/#it{c};#it{M}(#piKPr) (GeV/#it{c}^{2});#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            } else {
                dataContainer.histograms2d.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < #it{p}_{T, #Lambda_{c}^{+}} < %.0f GeV/#it{c};#it{M}(#piKPr) (GeV/#it{c}^{2});#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            }
        }
        dataContainer.histograms2d[i]->Sumw2();
    }

    return dataContainer;
}
//__________________________________________________________________________________________________________________________

void fillHistograms(TFile* fDist, SidebandData& dataContainer, double jetptMin, double jetptMax, const BinningStruct& binning) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on particle level delta R

    // Accessing TTree
    TTree* tree = (TTree*)fDist->Get("DF_2261906155510144/O2jetdisttable");

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
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= binning.deltaRBinEdges_detector[0]) && (deltaR < DeltaRcut))) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < binning.ptHFBinEdges_detector.size() - 1 && !filled; iEdge++) {
                if ((hfPt >= binning.ptHFBinEdges_detector[iEdge]) && (hfPt < binning.ptHFBinEdges_detector[iEdge + 1])) {
                    // Get the threshold for this pT range
                    double maxBkgProb = GetBkgProbabilityCut(hfPt, binning.bdtPtCuts);

                    // Fill histogram only if the cut is passed
                    if (hfMlScore0 < maxBkgProb) {
                        dataContainer.histograms2d[iEdge]->Fill(hfMass, deltaR);
                    }
                    filled = true; // Exit the loop once the correct histogram is found (alternative: break)
                }
                
            }
            
        } // kinematic range choice

    } // end of TTree entries loop
    cout << "Histograms filled.\n";

    // Creating 1D mass projections
    TH1D* tempHist;
    // Obtaining 1D invariant mass histograms from the projection
    for (size_t iHist = 0; iHist < dataContainer.histograms2d.size(); iHist++) {
        //
        tempHist = dataContainer.histograms2d[iHist]->ProjectionX(Form("h_mass_proj_%zu", iHist));
        dataContainer.histograms1d.push_back(tempHist);
    }
}

FitContainer performFit(TFile* fReflectionsMC, SidebandData& dataContainer, double& minMass, double& maxMass, const FitModelType& modelToUse, const double& jetptMin, const double& jetptMax) {
    
    double m_0_reference = 2.28646; // Lc mass in GeV/c^2
    double sigma_reference = 0.012;

    // --- Total fits: loop through pT,HF intervals/histograms and perform fits
    for (size_t iHisto = 0; iHisto < dataContainer.histograms2d.size(); ++iHisto) {
        
        // Get TF1 objects from MC file
        TF1* fSignal = (TF1*)fReflectionsMC->Get(Form("signalFit_%zu", iHisto));
        TF1* fReflections = (TF1*)fReflectionsMC->Get(Form("reflectionsFit_%zu", iHisto));
        std::cout << " --- Printing signal fit parameters. ---" << std::endl;
        fSignal->Print("V");
        std::cout << " --- Signal fit parameters printed. ---" << std::endl;
        std::cout << " --- Printing reflections fit parameters. ---" << std::endl;
        fReflections->Print("V");
        std::cout << " --- Reflections fit parameters printed. ---" << std::endl;

        // Extract parameters from MC pure signal fits to set constraints for data fits
        // Define signal fit variables
        double A1Signal; // A1Signal
        double A2Signal; // A2Signal
        double m0Signal; // m0
        double sigma1Signal; // sigma1
        double sigma2Signal; // sigma2
        // Fetch that Gaussian "1" will always be the primary (with bigger amplitude)
        if (fSignal->GetParameter(0) > fSignal->GetParameter(1)) {
            A1Signal = fSignal->GetParameter(0); // A1Signal
            A2Signal = fSignal->GetParameter(1); // A2Signal
            sigma1Signal = fSignal->GetParameter(3); // sigma1
            sigma2Signal = fSignal->GetParameter(4); // sigma2
        } else {
            A1Signal = fSignal->GetParameter(1); // A1Signal
            A2Signal = fSignal->GetParameter(0); // A2Signal
            sigma1Signal = fSignal->GetParameter(4); // sigma1
            sigma2Signal = fSignal->GetParameter(3); // sigma2
        }
        m0Signal = fSignal->GetParameter(2); // m0
        //std::cout << "A1Signal = " << A1Signal << ", A2Signal = " << A2Signal << ", m0Signal = " << m0Signal << ", sigma1Signal = " << sigma1Signal << ", sigma2Signal = " << sigma2Signal << std::endl;

        // Extract parameters from MC pure reflection fits to set constraints for data fits
        // Define signal fit variables
        double A1Reflection; // A1Reflection
        double A2Reflection; // A2Reflection
        double m0_1Reflections; // m0_1
        double m0_2Reflections; // m0_2
        double sigma1Reflections; // sigma1
        double sigma2Reflections; // sigma2
        // Fetch that Gaussian "1" will always be the primary (with bigger amplitude)
        if (fSignal->GetParameter(0) > fSignal->GetParameter(1)) {
            A1Reflection = fReflections->GetParameter(5); // A1Reflection
            A2Reflection = fReflections->GetParameter(6); // A2Reflection
            m0_1Reflections = fReflections->GetParameter(7); // m0_1
            m0_2Reflections = fReflections->GetParameter(8); // m0_2
            sigma1Reflections = fReflections->GetParameter(9); // sigma1
            sigma2Reflections = fReflections->GetParameter(10); // sigma2
        } else {
            A1Reflection = fReflections->GetParameter(6); // A1Reflection
            A2Reflection = fReflections->GetParameter(5); // A2Reflection
            m0_1Reflections = fReflections->GetParameter(8); // m0_1
            m0_2Reflections = fReflections->GetParameter(7); // m0_2
            sigma1Reflections = fReflections->GetParameter(10); // sigma1
            sigma2Reflections = fReflections->GetParameter(9); // sigma2
        }
        //std::cout << "A1Reflection = " << A1Reflection << ", A2Reflection = " << A2Reflection << ", m0_1Reflections = " << m0_1Reflections << ", m0_2Reflections = " << m0_2Reflections << ", sigma1Reflections = " << sigma1Reflections << ", sigma2Reflections = " << sigma2Reflections << std::endl;

        // -> signal constrains
        double A1toA2MCSignalRatio = A1Signal / A2Signal;
        double Sigma1toSigma2MCSignalRatio = sigma1Signal / sigma2Signal;

        // -> reflection constrains
        double A1toA2MCReflectionsRatio = A1Reflection / A2Reflection;
        double Sigma1toSigma2MCReflectionsRatio = sigma1Reflections / sigma2Reflections;
        double A1SignalToA1ReflectionMCRatio = A1Signal / A1Reflection;

        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        if (modelToUse == FitModelType::FullPowerLaw) { // 12 parameters = 8 fixed + 5 free
            dataContainer.fittings.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<FullModelPowerLaw>, minMass, maxMass, 13));
            // Apply positive amplitude constraint to power law background only
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(0, 0., TMath::Infinity()); // only accepts non-negative background fits
            // dataContainer.fittings.fitTotal[iHisto]->SetParLimits(1, -1e10, 0.); // only descending inclination background fits
        } else if (modelToUse == FitModelType::FullPoly2) { // 13 parameters = 8 fixed + 6 free
            dataContainer.fittings.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<FullModelPoly2>, minMass, maxMass, 14));
            dataContainer.fittings.fitTotal[iHisto]->SetParName(13, "Background C (m^2 term)");  // Free
        } else if (modelToUse == FitModelType::SignalReflectionsOnly) {
            dataContainer.fittings.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<SigRefModel>, minMass, maxMass, 13));
        }
        // Set initial values and fix parameters (original background amplitude guess = 200000)
        double guessBackgroundA = dataContainer.histograms1d[iHisto]->GetBinContent(dataContainer.histograms1d[iHisto]->FindBin(m_0_reference)) / 5; // initial guess for background amplitude at the signal peak
        double guessBackgroundB = -10.0; // initial guess for background slope (power law exponent or linear term)
        double guessA1Signal = dataContainer.histograms1d[iHisto]->GetBinContent(dataContainer.histograms1d[iHisto]->FindBin(m_0_reference)) / 2; // initial guess for A1Signal amplitude from MC fit
        Double_t params[14] = {guessBackgroundA, guessBackgroundB, guessA1Signal, A1toA2MCSignalRatio, m_0_reference, 0.02, 1.2, 2.0, 1.3, 1.83, 1.85, 0.02, 0.03, 2}; // last parameter only for poly2 background
        dataContainer.fittings.fitTotal[iHisto]->SetParameters(params);

        // Set parameter names
        dataContainer.fittings.fitTotal[iHisto]->SetParName(0, "Background A"); // Free
        dataContainer.fittings.fitTotal[iHisto]->SetParName(1, "Background B"); // Free
        dataContainer.fittings.fitTotal[iHisto]->SetParName(2, "A1 Signal"); // Free
        dataContainer.fittings.fitTotal[iHisto]->SetParName(3, "A1/A2 Signal Ratio");
        dataContainer.fittings.fitTotal[iHisto]->SetParName(4, "Signal Mean m0"); // Free
        dataContainer.fittings.fitTotal[iHisto]->SetParName(5, "Signal Width Sigma1"); // Free
        dataContainer.fittings.fitTotal[iHisto]->SetParName(6, "Sigma Ratio (Sig1/Sig2)");
        dataContainer.fittings.fitTotal[iHisto]->SetParName(7, "A1Signal/A1Reflection Ratio");
        dataContainer.fittings.fitTotal[iHisto]->SetParName(8, "A1/A2 Reflection Ratio");
        dataContainer.fittings.fitTotal[iHisto]->SetParName(9, "Reflection Mean m0_1");
        dataContainer.fittings.fitTotal[iHisto]->SetParName(10, "Reflection Mean m0_2");
        dataContainer.fittings.fitTotal[iHisto]->SetParName(11, "Reflection Width Sigma_1");
        dataContainer.fittings.fitTotal[iHisto]->SetParName(12, "Reflection Width Sigma_2");

        // Apply range limits to the parameters
        dataContainer.fittings.fitTotal[iHisto]->SetParLimits(0, 0.001, TMath::Infinity()); // only accepts positive background fits
        dataContainer.fittings.fitTotal[iHisto]->SetParLimits(1, -TMath::Infinity(), -1.01); // only descending inclination background fits
        dataContainer.fittings.fitTotal[iHisto]->SetParLimits(2, 3.0, TMath::Infinity()); // only accepts non-negative primary gaussian amplitude
        dataContainer.fittings.fitTotal[iHisto]->SetParLimits(4, m_0_reference, m_0_reference); // mean invariant mass should be around the one from literature, between [0.995 * m_0_reference, 1.005 * m_0_reference]
        // if (jetptMin == 30. && jetptMax == 50.) {
        //     dataContainer.fittings.fitTotal[iHisto]->SetParLimits(4, m_0_reference, m_0_reference);
        // } else {
        //     dataContainer.fittings.fitTotal[iHisto]->SetParLimits(4, 0.995 * m_0_reference, 1.005 * m_0_reference);
        // }
        dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.35 * sigma_reference, 4.0 * sigma_reference);

        // Fix the parameters from MC fits
        dataContainer.fittings.fitTotal[iHisto]->FixParameter(3, A1toA2MCSignalRatio);   // A1toA2MCSignalRatio
        dataContainer.fittings.fitTotal[iHisto]->FixParameter(6, Sigma1toSigma2MCSignalRatio);   // Sigma1toSigma2MCSignalRatio
        dataContainer.fittings.fitTotal[iHisto]->FixParameter(7, A1SignalToA1ReflectionMCRatio);   // A1SignalToA1ReflectionMCRatio
        dataContainer.fittings.fitTotal[iHisto]->FixParameter(8, A1toA2MCReflectionsRatio);   // A1toA2MCReflectionsRatio
        dataContainer.fittings.fitTotal[iHisto]->FixParameter(9, m0_1Reflections);  // m0_1
        dataContainer.fittings.fitTotal[iHisto]->FixParameter(10, m0_2Reflections); // m0_2
        dataContainer.fittings.fitTotal[iHisto]->FixParameter(11, sigma1Reflections); // sigma_1
        dataContainer.fittings.fitTotal[iHisto]->FixParameter(12, sigma2Reflections); // sigma_2

        // Perform fit with "Q" (quiet) option: no drawing of the fit function
        dataContainer.histograms1d[iHisto]->Fit(dataContainer.fittings.fitTotal[iHisto], "RQN");
        dataContainer.fittings.fitTotal[iHisto]->Print("V");

        // Check if signal gaussian acomodates the literature D0 rest mass
        double primaryMean = dataContainer.fittings.fitTotal[iHisto]->GetParameter(4);
        double primarySigma = dataContainer.fittings.fitTotal[iHisto]->GetParameter(5);
        if ((m_0_reference < (primaryMean - 2*primarySigma)) || (m_0_reference > (primaryMean + 2*primarySigma))) { // fit failed: PDG mass outside of signal region
            // if literature mass is outside, clear histogram as it should be excluded
            //histograms[iHisto]->Reset();
            //histograms[iHisto]->SetTitle(TString(histograms[iHisto]->GetTitle()) + " (m_{D^{0}^{literature}} #notin [#bar{m}_{primary}-2#sigma_{primary};#bar{m}_{primary}+2#sigma_{primary}])");
            dataContainer.histograms1d[iHisto]->SetTitle(TString(dataContainer.histograms1d[iHisto]->GetTitle()) + " (Fit failed)");
            //histograms2d[iHisto]->Reset("ICES");
            //histograms2d[iHisto]->Reset();
            //histograms2d[iHisto]->SetTitle(TString(histograms[iHisto]->GetTitle()) + " (m_{D^{0}^{literature}} #notin [#bar{m}_{primary}-2#sigma_{primary};#bar{m}_{primary}+2#sigma_{primary}])");
            dataContainer.histograms2d[iHisto]->SetTitle(TString(dataContainer.histograms1d[iHisto]->GetTitle()) + " (Fit failed)");
            dataContainer.fittings.workingFits.push_back(false);
        } else { // fit worked: PDG mass inside of signal region
            dataContainer.fittings.workingFits.push_back(true);
        }
    }
    std::cout << "Total fits performed.\n";

    // --- Background only fits: loop through pT,HF intervals/histograms and perform fits
    // Perform background only fit to each histogram
    int numberOfFits = dataContainer.fittings.fitTotal.size();
    for (size_t iHisto = 0; iHisto < numberOfFits; iHisto++) {
        if (modelToUse == FitModelType::FullPowerLaw) { // backgroundFunctionPowerLaw
            // Getting total parameter values (f(x) = a * x^b)
            double a_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(0); // Get the value of parameter 'a'
            double b_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(1); // Get the value of parameter 'b'
            dataContainer.fittings.fitBackgroundOnly.push_back(new TF1(Form("backgroundOnlyFit_%zu", iHisto), backgroundFunction, minMass, maxMass, 2));
            dataContainer.fittings.fitBackgroundOnly[iHisto]->FixParameter(0, a_par);
            dataContainer.fittings.fitBackgroundOnly[iHisto]->FixParameter(1, b_par);
        } else if (modelToUse == FitModelType::FullPoly2) { // backgroundFunctionPoly2
            // Getting total parameter values (f(x) = a*x² + b*x + c)
            double a_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(0); // Get the value of parameter 'a'
            double b_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(1); // Get the value of parameter 'b'
            double c_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(13); // Get the value of parameter 'c'
            dataContainer.fittings.fitBackgroundOnly.push_back(new TF1(Form("backgroundOnlyFit_%zu", iHisto), backgroundFunctionPoly2, minMass, maxMass, 3));
            dataContainer.fittings.fitBackgroundOnly[iHisto]->FixParameter(0, a_par);
            dataContainer.fittings.fitBackgroundOnly[iHisto]->FixParameter(1, b_par);
            dataContainer.fittings.fitBackgroundOnly[iHisto]->FixParameter(2, c_par);
        }
        
        // Getting total parameter values
        //double a_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(0); // Get the value of parameter 'a'
        //double b_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(1); // Get the value of parameter 'b'
        //dataContainer.fittings.fitBackgroundOnly.push_back(new TF1(Form("backgroundOnlyFit_%zu", iHisto), backgroundFunction, minMass, maxMass, 2));
        //dataContainer.fittings.fitBackgroundOnly[iHisto]->SetParameters(a_par,b_par); // dataContainer.fittings.size()-1 = the latest added to the vector
        dataContainer.fittings.fitBackgroundOnly[iHisto]->SetLineStyle(kDashed);
        dataContainer.fittings.fitBackgroundOnly[iHisto]->SetLineColor(kRed+1);

    }
    std::cout << "Background only fits performed.\n";

    // Perform signal only fit to each histogram
    for (size_t iHisto = 0; iHisto < numberOfFits; iHisto++) {
        // Extract signal-related parameters from the total fit
        double a1Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(2); // A1 signal
        double A1toA2MCSignalRatio = dataContainer.fittings.fitTotal[iHisto]->GetParameter(3); // A2 signal
        double m0Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(4); // Signal mean m0
        double sigmaSignal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(5); // Signal width sigma
        double Sigma1toSigma2MCSignalRatio = dataContainer.fittings.fitTotal[iHisto]->GetParameter(6); // Sigma1/Sigma2 signal
        //std::cout << "A1Signal = " << a1Signal << ", A1toA2MCSignalRatio = " << A1toA2MCSignalRatio << ", m0Signal = " << m0Signal << ", sigmaSignal = " << sigmaSignal << ", Sigma1toSigma2MCSignalRatio = " << Sigma1toSigma2MCSignalRatio << std::endl;

        // Perfoming fit with acquired parameters from total fit
        dataContainer.fittings.fitSignalOnly.push_back(new TF1(Form("signalOnlyFit_%zu", iHisto), signalOnlyFunction, minMass, maxMass, 5));
        dataContainer.fittings.fitSignalOnly[iHisto]->FixParameter(0, a1Signal);
        dataContainer.fittings.fitSignalOnly[iHisto]->FixParameter(1, A1toA2MCSignalRatio);
        dataContainer.fittings.fitSignalOnly[iHisto]->FixParameter(2, m0Signal);
        dataContainer.fittings.fitSignalOnly[iHisto]->FixParameter(3, sigmaSignal);
        dataContainer.fittings.fitSignalOnly[iHisto]->FixParameter(4, Sigma1toSigma2MCSignalRatio);
        dataContainer.fittings.fitSignalOnly[iHisto]->SetLineStyle(kDashed);
        dataContainer.fittings.fitSignalOnly[iHisto]->SetLineColor(kBlue+1);
        //dataContainer.fittings.fitSignalOnly[iHisto]->Print("V");
    }
    std::cout << "Signal only fits performed.\n";

    // Perform reflections only fit to each histogram
    for (size_t iHisto = 0; iHisto < numberOfFits; iHisto++) {
        // Extract reflection-related parameters from the total fit
        double A1SignalToA1ReflectionMCRatios = dataContainer.fittings.fitTotal[iHisto]->GetParameter(7); // A1 signal / A1 reflection
        double A1Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(2); // A1 signal
        double A1toA2MCReflectionsRatio = dataContainer.fittings.fitTotal[iHisto]->GetParameter(8); // A1 reflection / A2 reflection
        double m0_1Reflection = dataContainer.fittings.fitTotal[iHisto]->GetParameter(9); // Reflection mean m0_1
        double m0_2Reflection = dataContainer.fittings.fitTotal[iHisto]->GetParameter(10); // Reflection mean m0_2
        double sigma1Reflection = dataContainer.fittings.fitTotal[iHisto]->GetParameter(11); // Reflection sigma_1
        double sigma2Reflection = dataContainer.fittings.fitTotal[iHisto]->GetParameter(12); // Reflection sigma_2
        //std:cout << "A1SignalToA1ReflectionMCRatios = " << A1SignalToA1ReflectionMCRatios << ", A1Signal = " << A1Signal << ", A1toA2MCReflectionsRatio = " << A1toA2MCReflectionsRatio << ", m0_1Reflection = " << m0_1Reflection << ", m0_2Reflection = " << m0_2Reflection << ", sigma1Reflection = " << sigma1Reflection << ", sigma2Reflection = " << sigma2Reflection << std::endl;

        // Perfoming fit with acquired parameters from total fit
        dataContainer.fittings.fitReflectionsOnly.push_back(new TF1(Form("reflectionOnlyFit_%zu", iHisto), reflectionOnlyFunction, minMass, maxMass, 7));
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(0, A1SignalToA1ReflectionMCRatios);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(1, A1Signal);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(2, A1toA2MCReflectionsRatio);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(3, m0_1Reflection);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(4, m0_2Reflection);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(5, sigma1Reflection);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(6, sigma2Reflection);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->SetLineStyle(kDashed);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->SetLineColor(kGreen+3);
        //dataContainer.fittings.fitReflectionsOnly[iHisto]->Print("V");
        //std::cout << "A1 = " << A1Signal / A1SignalToA1ReflectionMCRatios << ", A2 = " << (A1Signal / A1SignalToA1ReflectionMCRatios) / A1toA2MCReflectionsRatio << std::endl;
        //std::cout << "A1/A2 = " << A1toA2MCReflectionsRatio << std::endl;

    }
    std::cout << "Reflections only fits performed.\n";
    
    return dataContainer.fittings;
}

SubtractionResult SideBand(SidebandData& dataContainer, const BinningStruct& binning, double& signalSigmas, int& startingBackSigma, int& backgroundSigmas) {
    
    // Creating histograms for collecting data
    TH1D* h_sideBand;
    TH1D* h_signal;
    TH1D* h_back_subtracted;

    // Loop through pT,HF intervals/histograms and perform side band subtraction
    for (size_t iHisto = 0; iHisto < dataContainer.histograms1d.size(); ++iHisto) {
        std::cout << "\nPerforming side-band subtraction for fit number " << iHisto << std::endl;
        // Get total fit parameters
        double m_0 = dataContainer.fittings.fitTotal[iHisto]->GetParameter(4); // Get the value of parameter 'm_0'
        double signalSigma1 = dataContainer.fittings.fitTotal[iHisto]->GetParameter(5); // Get the value of parameter 'sigma1'
        double signalSigma2 = dataContainer.fittings.fitTotal[iHisto]->GetParameter(5) / dataContainer.fittings.fitTotal[iHisto]->GetParameter(6); // signalSigma2 = sigma1 / sigmaRatio12
        double sigma = signalSigma1;

        // Calculate sideband regions
        std::pair<std::array<double, 2>, std::array<double, 2>> sidebandRanges = calculateSidebandRegions(iHisto, dataContainer.fittings, dataContainer.histograms1d[iHisto], startingBackSigma, backgroundSigmas);

        // Calculate scaling factor to apply to the sideband subtraction original method
        std::array<double, 2> scallingFactors = calculateScalingFactor(iHisto, dataContainer.fittings, sidebandRanges.first, sidebandRanges.second, signalSigmas);
        dataContainer.subtractionResults.scalingFactorsArrays.push_back(scallingFactors);

        // Create signal histogram
        int lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(m_0 - signalSigmas * sigma);
        int highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(m_0 + signalSigmas * sigma);
        h_signal = dataContainer.histograms2d[iHisto]->ProjectionY(Form("h_signal_proj_%zu",iHisto), lowBin, highBin);
        dataContainer.subtractionResults.signalHist.push_back(h_signal);

        // Create side-band histogram
        if (sidebandRanges.first[0] == 0. && sidebandRanges.first[1] == 0.) {
            // Use only the right sideband
            lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[0]);
            highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[1]);
            h_sideBand = dataContainer.histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_%zu",iHisto), lowBin, highBin); // sum the right sideband
            double rightSBHistogram = dataContainer.histograms1d[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram

            double rightSB = dataContainer.fittings.fitTotal[iHisto]->Integral(sidebandRanges.second[0],sidebandRanges.second[1]);
            double fitSigIntegral = dataContainer.fittings.fitTotal[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

            lowBin = dataContainer.histograms1d[iHisto]->FindBin(m_0 - signalSigmas * sigma);
            highBin = dataContainer.histograms1d[iHisto]->FindBin(m_0 + signalSigmas * sigma);
            double signalHistogram = dataContainer.histograms1d[iHisto]->Integral(lowBin, highBin);
        } else {
            // Use both sidebands
            lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.first[0]);
            highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.first[1]);
            h_sideBand = dataContainer.histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_left_%zu",iHisto), lowBin, highBin); // sum the left sideband
            double leftSBHistogram = dataContainer.histograms1d[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram
            std::cout << "Sideband integral (only left) = " << h_sideBand->Integral() << std::endl;
            std::cout << "leftSBHistogram = " << leftSBHistogram << std::endl;

            lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[0]);
            highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[1]);
            h_sideBand->Add(dataContainer.histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_%zu",iHisto), lowBin, highBin)); // sum the right sideband
            double rightSBHistogram = dataContainer.histograms1d[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram

            double leftSB = dataContainer.fittings.fitTotal[iHisto]->Integral(sidebandRanges.first[0],sidebandRanges.first[1]); // in the opposite order because the histogram is decreasing
            double rightSB = dataContainer.fittings.fitTotal[iHisto]->Integral(sidebandRanges.second[0],sidebandRanges.second[1]);
            double fitSigIntegral = dataContainer.fittings.fitTotal[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

            lowBin = dataContainer.histograms1d[iHisto]->FindBin(m_0 - signalSigmas * sigma);
            highBin = dataContainer.histograms1d[iHisto]->FindBin(m_0 + signalSigmas * sigma);
            double signalHistogram = dataContainer.histograms1d[iHisto]->Integral(lowBin, highBin);
        }

        // Scale the sideband histogram by ratio of it in the signal region Bs/(B1+B2)
        std::cout << "Sideband integral (before beta scaling) = " << h_sideBand->Integral() << std::endl;
        h_sideBand->Scale(scallingFactors[1]);
        std::cout << "Sideband integral (after beta scaling) = " << h_sideBand->Integral() << std::endl;
        dataContainer.subtractionResults.sidebandHist.push_back(h_sideBand);

        // Subtract background histogram from signal histogram
        h_back_subtracted = (TH1D*)h_signal->Clone(Form("h_back_subtracted_%zu",iHisto));
        std::cout << "Before subtraction: Signal integral = " << h_signal->Integral() << ", Sideband integral = " << h_sideBand->Integral() << std::endl;
        h_back_subtracted->Add(h_sideBand,-1.0);

        // Scale by reflections correction
        h_back_subtracted->Scale(scallingFactors[0]);

        // Account for two sigma only area used for signal region
        double coverage = TMath::Erf(signalSigmas / sqrt(2)); // coverage of a Gaussian in a +/- signalSigmas window
        h_back_subtracted->Scale(1 / coverage);
        // If histogram isn't empty, store it in the container
        bool isHistoEmpty = (h_back_subtracted->GetEntries() != 0) ? true : false;
        dataContainer.subtractionResults.subtractedHist.push_back(h_back_subtracted);

        // Calculating significance
        // Get signal yield by integrating Gaussian in a 2σ window
        double mean_signalOnly = dataContainer.fittings.fitSignalOnly[iHisto]->GetParameter(2);
        double sigma1_signalOnly = dataContainer.fittings.fitSignalOnly[iHisto]->GetParameter(3);
        //fittings.fitSignalOnly[iHisto]->Print("V");
        double S = dataContainer.fittings.fitSignalOnly[iHisto]->Integral(mean_signalOnly - 2*sigma1_signalOnly, mean_signalOnly + 2*sigma1_signalOnly);
        double B = dataContainer.fittings.fitBackgroundOnly[iHisto]->Integral(mean_signalOnly - 2*sigma1_signalOnly, mean_signalOnly + 2*sigma1_signalOnly);
        double significance = 0.;
        if (!std::isnan(S)) {
            significance = S / sqrt(S+B); // Binomial-like Approximation: when S is comparable to B, consider both signal and background fluctuations
        }
        
        std::cout << "For histogram index " << iHisto << ": S = " << S << ", B = " << B << "\t Significance = " << significance << std::endl;
        // Store signal and background values for each pT,D bin
        dataContainer.subtractionResults.signal_background_values.push_back(std::make_pair(S, B));
        // Storing significance in the output struct object
        if (iHisto == 0) {
            // Set a significance value for each corresponding pT,D bin
            dataContainer.subtractionResults.hSignificance = new TH1D("hSignificance", "Estimated significance for each m_{invariant} distribution bin;p_{T, D^{0}} (GeV/#it{c});Significance = #frac{S}{#sqrt{S+B}}", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data());
            dataContainer.subtractionResults.hSignificance->SetBinContent(iHisto+1, significance);
            dataContainer.subtractionResults.hSignificance->SetBinError(iHisto+1, 0);
        } else {
            dataContainer.subtractionResults.hSignificance->SetBinContent(iHisto+1, significance);
            dataContainer.subtractionResults.hSignificance->SetBinError(iHisto+1, 0);
        }
    } // end of pT,HF loop for sideband subtraction

    // Final adjustments to subtracted histograms: cleaning and storage
    for (size_t iHisto = 0; iHisto < dataContainer.subtractionResults.subtractedHist.size(); iHisto++) {
        
        // Check if histogram should be erased before storing based on fit performed
        bool eraseHist = eraseHistogram(dataContainer.fittings.workingFits, iHisto);
        if (eraseHist) {
            dataContainer.subtractionResults.subtractedHist[iHisto]->Reset("ICES");
        } else {
            // Set negative count bin entries to 0
            for (int iBin = 1; iBin <= dataContainer.subtractionResults.subtractedHist[iHisto]->GetNbinsX(); iBin++) {
                if (dataContainer.subtractionResults.subtractedHist[iHisto]->GetBinContent(iBin) < 0) {
                    dataContainer.subtractionResults.subtractedHist[iHisto]->SetBinContent(iBin,0);
                    dataContainer.subtractionResults.subtractedHist[iHisto]->SetBinError(iBin,0);
                }
            }
        }
        
        // Obtain summed histogram for all pT,D bins after sideband subtraction
        if (iHisto == 0) {
            dataContainer.subtractionResults.hSubtractedFullJetPt = (TH1D*)dataContainer.subtractionResults.subtractedHist[iHisto]->Clone("hSubtractedFullJetPt");
        } else {
            dataContainer.subtractionResults.hSubtractedFullJetPt->Add(dataContainer.subtractionResults.subtractedHist[iHisto]);
        }
    }

    return dataContainer.subtractionResults;
}

void PlotHistograms(const SidebandData& dataContainer, const BinningStruct& binning, double jetptMin, double jetptMax, const std::vector<double>& ptHFBinEdges) {
    std::cout << "Plotting histograms..." << std::endl;


    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.04); // default = 0.05

    // Create a canvas for plotting
    int nHistos = dataContainer.histograms1d.size();
    // Start with a square layout (or close to it)
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));
    TCanvas* c1d_fit = new TCanvas("c1d_fit", "1D histograms with Fit", 800, 600);
    c1d_fit->SetCanvasSize(1800,1000);
    c1d_fit->Divide(nCols,nRows); // columns, lines
    TCanvas* c_2d = new TCanvas("c_2d", "2D histograms", 800, 600);
    c_2d->SetCanvasSize(1800,1000);
    c_2d->Divide(nCols,nRows); // columns, lines

    // Loop through all histograms and fitting functions
    for(size_t iHisto = 0; iHisto < dataContainer.histograms1d.size(); ++iHisto) {
        //
        c1d_fit->cd(iHisto+1);
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        dataContainer.histograms1d[iHisto]->SetMarkerStyle(kDot); //kFullDotMedium
        dataContainer.histograms1d[iHisto]->SetMarkerColor(kBlack);
        dataContainer.histograms1d[iHisto]->SetLineColor(kBlack);
        dataContainer.histograms1d[iHisto]->GetYaxis()->SetTitle("counts");
        dataContainer.histograms1d[iHisto]->SetMinimum(0);
        
        dataContainer.histograms1d[iHisto]->Draw();
        dataContainer.fittings.fitTotal[iHisto]->SetLineColor(kBlack);
        dataContainer.fittings.fitTotal[iHisto]->SetLineStyle(kSolid);
        dataContainer.fittings.fitTotal[iHisto]->SetLineWidth(1);
        dataContainer.fittings.fitTotal[iHisto]->Draw("same");
        dataContainer.fittings.fitSignalOnly[iHisto]->Draw("same");
        dataContainer.fittings.fitBackgroundOnly[iHisto]->Draw("same");
        dataContainer.fittings.fitReflectionsOnly[iHisto]->Draw("same");

        // double A1Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(2); // Get the value of parameter 'A' from primary signal gaussian
        // double m0_1Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(4); // Get the value of parameter 'm0' from primary signal gaussian
        // double sigma1Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(5); // Get the value of parameter 'sigma' from primary signal gaussian
        // double A2Signal = A1Signal / dataContainer.fittings.fitTotal[iHisto]->GetParameter(3); // Get the value of parameter 'A' from secondary signal gaussian
        // double sigma2Signal = sigma1Signal / dataContainer.fittings.fitTotal[iHisto]->GetParameter(6); // Get the value of parameter 'sigma' from secondary signal gaussian
        double chi2 = dataContainer.fittings.fitTotal[iHisto]->GetChisquare();
        double degOfFreedom = dataContainer.fittings.fitTotal[iHisto]->GetNDF();
        
        // latex->DrawLatex(statBoxPos-0.52, 0.80, Form("A_{1}^{signal} = %.2f, #bar{m_{1}} = %.2f, #sigma_{1} = %.2f GeV/#it{c}^{2}", A1Signal, m0_1Signal,sigma1Signal));
        // latex->DrawLatex(statBoxPos-0.52, 0.75, Form("A_{2}^{signal} = %.2f, #bar{m_{2}} = %.2f, #sigma_{2} = %.2f GeV/#it{c}^{2}", A2Signal, m0_1Signal,sigma2Signal));
        latex->DrawLatex(statBoxPos-0.28, 0.85, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));
        latex->DrawLatex(statBoxPos-0.87, 0.34, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax)); // Display jet pT cut applied
        latex->DrawLatex(statBoxPos-0.87, 0.28, "MC: LHC24g5_All"); // HF_LHC24g5_All
        latex->DrawLatex(statBoxPos-0.87, 0.23, "Data: LHC22o_pass7"); // HF_LHC22o_pass7_minBias_small_2P3PDstar

        // Drawing 2D histograms
        c_2d->cd(iHisto+1);
        dataContainer.histograms2d[iHisto]->Draw("colz");

    }

    // Plotting extracted raw yields DeltaR observable
    TCanvas* cSideBand = new TCanvas("cSideBand", "delta R for side-band", 800, 600);
    cSideBand->SetCanvasSize(1800,1000);
    cSideBand->Divide(nCols,nRows); // columns, lines
    TCanvas* cSignal = new TCanvas("cSignal", "delta R for signal", 800, 600);
    cSignal->SetCanvasSize(1800,1000);
    cSignal->Divide(nCols,nRows); // columns, lines
    TCanvas* cSubtracted = new TCanvas("cSubtracted", "delta R for side-band subtracted signal", 800, 600);
    cSubtracted->SetCanvasSize(1800,1000);
    cSubtracted->Divide(nCols,nRows); // columns, lines
    TCanvas* cSigPlusBack = new TCanvas("cSigPlusBack", "delta R for side-band and signal in the same plot", 800, 600);
    cSigPlusBack->SetCanvasSize(1800,1000);
    cSigPlusBack->Divide(nCols,nRows); // columns, lines

    TLegend* legend = new TLegend(0.6,0.57,0.9,0.77);
    legend->AddEntry(dataContainer.subtractionResults.sidebandHist[0],"Sideband", "lpe");
    legend->AddEntry(dataContainer.subtractionResults.signalHist[0],"Signal", "lpe");
    legend->AddEntry(dataContainer.subtractionResults.subtractedHist[0],"Signal (minus background)", "lpe");
    std::cout << "Starting subtracted histograms plotting..." << std::endl;
    for (size_t iHisto = 0; iHisto < dataContainer.subtractionResults.subtractedHist.size(); iHisto++) {
        
        cSideBand->cd(iHisto+1);
        //dataContainer.subtractionResults.sidebandHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        dataContainer.subtractionResults.sidebandHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        dataContainer.subtractionResults.sidebandHist[iHisto]->GetYaxis()->SetTitle("yields");
        dataContainer.subtractionResults.sidebandHist[iHisto]->SetMarkerStyle(kFullTriangleUp);
        dataContainer.subtractionResults.sidebandHist[iHisto]->SetMarkerColor(kAzure);
        dataContainer.subtractionResults.sidebandHist[iHisto]->SetLineColor(kAzure);
        double statBoxPos = gPad->GetUxmax();
        dataContainer.subtractionResults.sidebandHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));

        cSignal->cd(iHisto+1);
        //dataContainer.subtractionResults.signalHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        dataContainer.subtractionResults.signalHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        dataContainer.subtractionResults.signalHist[iHisto]->GetYaxis()->SetTitle("yields");
        dataContainer.subtractionResults.signalHist[iHisto]->SetMarkerStyle(kFullSquare);
        dataContainer.subtractionResults.signalHist[iHisto]->SetMarkerColor(kViolet);
        dataContainer.subtractionResults.signalHist[iHisto]->SetLineColor(kViolet);
        dataContainer.subtractionResults.signalHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));
        std::cout << "Plotting subtracted histogram number " << iHisto << std::endl;
        cSubtracted->cd(iHisto+1);
        //dataContainer.subtractionResults.subtractedHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        dataContainer.subtractionResults.subtractedHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        dataContainer.subtractionResults.subtractedHist[iHisto]->GetYaxis()->SetTitle("yields");
        dataContainer.subtractionResults.subtractedHist[iHisto]->SetMarkerStyle(kFullCircle);
        dataContainer.subtractionResults.subtractedHist[iHisto]->SetMarkerColor(kPink);
        dataContainer.subtractionResults.subtractedHist[iHisto]->SetLineColor(kPink);
        dataContainer.subtractionResults.subtractedHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));

        cSigPlusBack->cd(iHisto+1);
        //dataContainer.subtractionResults.signalHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        dataContainer.subtractionResults.signalHist[iHisto]->Draw();
        dataContainer.subtractionResults.sidebandHist[iHisto]->Draw("same");
        dataContainer.subtractionResults.subtractedHist[iHisto]->Draw("same");
        legend->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.5, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));
    }

    // Plot scaling factors
    TCanvas* cscallingFactorsArrays = new TCanvas("cscallingFactorsArrays", "Scaling Factors Arrays", 1800, 1000);
    cscallingFactorsArrays->cd();
    // Fill two histograms with the scaling factors
    int nBins = binning.ptHFBinEdges_detector.size();
    TH1D* hAlpha = new TH1D("hAlpha", "Scaling factors;p_{T,D^{0}} (GeV/c); s.f.", binning.ptHFBinEdges_detector.size()-1, binning.ptHFBinEdges_detector.data());
    TH1D* hBeta = new TH1D("hBeta", "Scaling factor;p_{T,D^{0}} (GeV/c); s.f.", binning.ptHFBinEdges_detector.size()-1, binning.ptHFBinEdges_detector.data());
    for (int iBin = 0; iBin < nBins; ++iBin) {
        hAlpha->SetBinContent(iBin + 1, dataContainer.subtractionResults.scalingFactorsArrays[iBin][0]); // alpha
        hBeta->SetBinContent(iBin + 1, dataContainer.subtractionResults.scalingFactorsArrays[iBin][1]); // beta
    }
    hBeta->SetMarkerStyle(kFullCircle);
    hBeta->SetMarkerColor(kBlue+1);
    hBeta->SetLineColor(kBlue+1);
    hBeta->SetLineWidth(2);
    hAlpha->SetMarkerStyle(kFullCircle);
    hAlpha->SetMarkerColor(kRed+1);
    hAlpha->SetLineColor(kRed+1);
    hAlpha->SetLineWidth(2);
    if (hAlpha->GetBinContent(hAlpha->GetMaximumBin()) > hBeta->GetBinContent(hBeta->GetMaximumBin())) {
        hAlpha->SetMinimum(0);
        hAlpha->Draw();
        hBeta->Draw("same");
    } else {
        hBeta->SetMinimum(0);
        hBeta->Draw();
        hAlpha->Draw("same");
    }
    TLegend* lScalingFactors = new TLegend(0.1,0.3,0.18,0.4);
    lScalingFactors->AddEntry(hAlpha,"#alpha","l");
    lScalingFactors->AddEntry(hBeta,"#beta","l");
    lScalingFactors->Draw();



    //
    // Storing images
    //
    TString imagePath = "../../Images/1-SignalTreatment/SideBand/";

    //
    // Storing in a single pdf file
    //
    c1d_fit->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf(",jetptMin,jetptMax));
    c_2d->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf",jetptMin,jetptMax));
    cSigPlusBack->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf",jetptMin,jetptMax));
    cscallingFactorsArrays->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf)",jetptMin,jetptMax));
}

void SaveData(SidebandData& dataContainer, const BinningStruct& binning, double jetptMin, double jetptMax) {
    // Open output file
    TFile* fOutput = new TFile(Form("backSub_%d_to_%d_jetpt_with_reflections.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    if (!fOutput || fOutput->IsZombie()) {
        std::cerr << "Error: Unable to open the output ROOT file." << std::endl;
    }

    // Loop over background subtracted histograms
    for (size_t iHisto = 0; iHisto < dataContainer.subtractionResults.subtractedHist.size(); iHisto++) {
        // store each histogram in file
        dataContainer.subtractionResults.subtractedHist[iHisto]->Write();
    }
    // Significance as a function of pT,D mass histogram
    dataContainer.subtractionResults.hSignificance->Write();
    // Final summed pT,D histogram
    dataContainer.subtractionResults.hSubtractedFullJetPt->Write();

    // Store binning information in the output file for later use
    storeBinningInFile(fOutput, binning);
}

void JetPtIterator(const double jetptMin, const double jetptMax, const BinningStruct& binning) {
    std::cout << "============================================= " << jetptMin << " GeV/c < pT,jet < " << jetptMax << " GeV/c =============================================" << std::endl;

    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.05 (from first of pair)
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {1, 0.05}, {2, 0.08}, {3, 0.08}, {4, 0.1}, {5, 0.2}, {6, 0.25}, {7, 0.3}, {8, 0.4}, {12, 0.5}, {24, 0.7}, {100, 0.7}
    }; // on dataset HF_LHC22o_pass7_minBias_small_2P3PDstar

    // Lc mass in GeV/c^2
    double m_0_parameter = 2.28646;
    double sigmaInitial = 0.012;

    // mass histogram
    int massBins = 50; // default=50 
    double minMass = 2.1; // default = 2.1
    double maxMass = 2.4822; // default = 2.49, due to drop of counts on last bin (of various pT,HF ranges) changed to 2.49-0.0078=2.4822

    // Opening data file
    TFile* fDist = new TFile("../../Data/Experimental/Train_606159/AO2D.root","read");
    if (!fDist || fDist->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT data file." << std::endl;
    }
    // Opening
    TFile* fReflectionsMC = new TFile(Form("../Reflections/reflections_%.0f_to_%.0fGeV.root",jetptMin,jetptMax),"read");
    if (!fReflectionsMC || fReflectionsMC->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT reflections file." << std::endl;
    }
    // Load ΔR bin edges
    double minDeltaR = binning.deltaRBinEdges_detector[0];
    double maxDeltaR = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1];
    // Load pTHF bin edges
    double minPtHF = binning.ptHFBinEdges_detector[0];
    double maxPtHF = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1];

    // Create multiple histograms
    SidebandData dataContainer = createHistograms(binning.ptHFBinEdges_detector,                    // the pT,HF edges will determine the number of mass histograms
                                                     massBins, minMass, maxMass,                    // mass histograms binning
                                                     binning.deltaRBinEdges_detector, jetptMin);    // deltaR histograms with asymmetrical bin widths
    //
    // Fill histograms
    fillHistograms(fDist, dataContainer, jetptMin, jetptMax, binning);

    // Perform fits
    FitModelType modelToUse = FitModelType::FullPoly2; // options: FullPowerLaw, FullPoly2, SignalReflectionsOnly, SignalOnly, ReflectionsOnly
    FitContainer fittings = performFit(fReflectionsMC, dataContainer, minMass, maxMass, modelToUse, jetptMin, jetptMax);

    // signal/side-band region parameters
    double signalSigmas = 2; // default delta = 2
    int startingBackSigma = 4; // default position = 4
    int backgroundSigmas = 4; // default delta = 4

    // Subtract side-band from signal
    SubtractionResult subtractionResult = SideBand(dataContainer, binning, signalSigmas, startingBackSigma, backgroundSigmas);

    // Plot histograms
    PlotHistograms(dataContainer, binning, jetptMin, jetptMax, binning.ptHFBinEdges_detector);

    // Storing final histograms to output file
    SaveData(dataContainer, binning, jetptMin, jetptMax);
}

void create3DBackgroundSubtracted(const BinningStruct& binning) {

    // jet pT cuts
    double jetptMin = binning.ptjetBinEdges_detector[0]; // GeV
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1]; // GeV
    // deltaR histogram
    double minDeltaR = binning.deltaRBinEdges_detector[0];
    double maxDeltaR = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1];

    TH3D* h3D = new TH3D("h3DBackgroundSubtracted", "Background subtracted; p_{T,jet} (GeV/c); #DeltaR; p_{T,D^{0}} (GeV/c)",
                     binning.ptjetBinEdges_detector.size()-1, binning.ptjetBinEdges_detector.data(),
                     binning.deltaRBinEdges_detector.size()-1, binning.deltaRBinEdges_detector.data(),
                     binning.ptHFBinEdges_detector.size()-1, binning.ptHFBinEdges_detector.data());
    h3D->Sumw2();

    for (size_t iJetBin = 0; iJetBin < binning.ptjetBinEdges_detector.size() - 1; iJetBin++) {
        TFile* fJetRange = new TFile(Form("backSub_%0.f_to_%0.f_jetpt_with_reflections.root",binning.ptjetBinEdges_detector[iJetBin],binning.ptjetBinEdges_detector[iJetBin+1]),"read");
        if (!fJetRange || fJetRange->IsZombie()) {
            std::cerr << "Error opening file " << Form("backSub_%0.f_to_%0.f_jetpt_with_reflections.root",binning.ptjetBinEdges_detector[iJetBin],binning.ptjetBinEdges_detector[iJetBin+1]) << std::endl;
            continue;
        }

        for (int iHist = 0; iHist < 100; ++iHist) {
            TString histName = Form("h_back_subtracted_%d", iHist);
            TH1D* hDeltaR = (TH1D*)fJetRange->Get(histName);
            if (!hDeltaR) break; // No more histograms
            //std::cout << "Histogram x axis range: " << hDeltaR->GetXaxis()->GetXmin() << " to " << hDeltaR->GetXaxis()->GetXmax() << std::endl;

            // Get bin centers
            double ptJetCenter = 0.5 * (binning.ptjetBinEdges_detector[iJetBin] + binning.ptjetBinEdges_detector[iJetBin+1]);
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

            double ptHFlow = -1, ptHFhigh = -1;
            if (sscanf(title.Data(), "%lf < #it{p}_{T, D^{0}} < %lf GeV/#it{c}", &ptHFlow, &ptHFhigh) != 2) {
                std::cerr << "Could not parse pT,D range from histogram title: " << title << std::endl;
                std::cout << "Trying to parse title: '" << title << "'" << std::endl;
                std::cout << "hDeltaR histogram -> " << hDeltaR->GetName() << std::endl;
                std::cout << "In file " << fJetRange->GetName() << std::endl;
                continue;
            }
            double ptHFCenter = 0.5 * (ptHFlow + ptHFhigh);


            for (int iBin = 1; iBin <= hDeltaR->GetNbinsX(); ++iBin) {
                double deltaRcenter = hDeltaR->GetBinCenter(iBin);
                double content = hDeltaR->GetBinContent(iBin);
                double error = hDeltaR->GetBinError(iBin);

                // Fill the TH3D using centers
                h3D->Fill(ptJetCenter, deltaRcenter, ptHFCenter, content);
                // Optionally, if you want to preserve error propagation:
                int binX = h3D->GetXaxis()->FindBin(ptJetCenter);
                int binY = h3D->GetYaxis()->FindBin(deltaRcenter);
                int binZ = h3D->GetZaxis()->FindBin(ptHFCenter);
                h3D->SetBinError(binX, binY, binZ, error); // only if needed
            }
        }
    }

    h3D->Draw("colz");

    TFile* fOutput = new TFile(Form("full_merged_ranges_back_sub.root"), "RECREATE");
    h3D->Write();

    // Store binning information in the output file for later use
    storeBinningInFile(fOutput, binning);

    std::cout << "3D histogram created and saved to " << fOutput->GetName() << " with 3D histogram." << std::endl;
}

void BackgroundSubtraction() {
    
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Opening
    double jetptMin = 5.; // default = 5 GeV, for Emma's use 5 GeV
    bool useEmmaYeatsBins = false;
    double jetptMax = useEmmaYeatsBins ? 10. : 7.; // default = 7GeV, for Emma's use 10 GeV
    TFile* fBinning = new TFile(Form("../Reflections/binningInfo.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    // Load binning from reflections file
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    // binning.ptjetBinEdges_detector = {7., 10.}; // for quick tests

    for (size_t iJetPt = 0; iJetPt < binning.ptjetBinEdges_detector.size() - 1; iJetPt++) {
        // Apply side-band method to each pT,jet bin
        std::cout << "Processing pT,jet bin: " << binning.ptjetBinEdges_detector[iJetPt] << " to " << binning.ptjetBinEdges_detector[iJetPt+1] << " GeV/c" << std::endl;
        jetptMin = binning.ptjetBinEdges_detector[iJetPt];
        jetptMax = binning.ptjetBinEdges_detector[iJetPt+1];
        JetPtIterator(jetptMin, jetptMax, binning);
    }

    // Compute the entire range too
    jetptMin = binning.ptjetBinEdges_detector[0];
    jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    JetPtIterator(jetptMin, jetptMax, binning);

    // Create 3D final histogram with pT,jet vs DeltaR vs pT,D
    create3DBackgroundSubtracted(binning);

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl;
}

int main() {
    BackgroundSubtraction();
    return 0;
}