/**
 * D0 meson analysis
 * @file BackgroundSubtraction.C
 * @brief Background subtraction using sideband method
 * Input: Reflections.root - contain histograms and fit templates of Lambda_c signal and reflections
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
using StandardSideBandSubtraction = FitModel<SignalSinglePolicy,PowerLawBackgroundPolicy>;                  // Standard side-band subtraction procedure
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
    std::vector<TH1D*> subtractedHistCopy;                              // 1D sideband subtracted histograms vector copy for scaling factor variations
    TH1D* hSubtractedFullJetPt;                                         // final deltaR distribution for all pT,HF summed
    TH1D* hSignificance;                                                // final significance distribution for all pT,HF summed
    std::vector<std::pair<double, double>> signal_background_values;    // vector of signal and background values for each pT,HF bin
    std::vector<std::array<double, 2>> scalingFactorsArrays;            // alpha and beta scaling factors for each pT,HF bin
    std::vector<TH1D*> sidebandHistDrawing;                             // histogram used for drawing the side-band shaded invariant mass distribution -in red
    std::vector<TH1D*> signalHistDrawing;                               // histogram used for drawing the signal shaded invariant mass distribution - in blue
};
struct SidebandData {
    double D0_reference_mass = 1.86484; // D0 mass in GeV/c^2
    std::vector<TH2D*> histograms2d;
    std::vector<TH1D*> histograms1d;
    FitContainer fittings;
    SubtractionResult subtractionResults;
};

// Module to create TH2D histograms including interest variable: VARIABLE bin sizes
SidebandData createHistograms(const BinningStruct& binning, const double& massBinDensity, 
                               const double& minMass, const std::vector<double>& maxMass) {

    // Validate that maxMass vector has the correct size
    if (maxMass.size() != binning.ptHFBinEdges_detector.size() - 1) {
        std::cerr << "Error: maxMass vector size (" << maxMass.size() 
                  << ") does not match number of pT,D bins (" 
                  << binning.ptHFBinEdges_detector.size() - 1 << ").\n";
        return SidebandData();
    }

    SidebandData dataContainer;

    for (size_t i = 0; i < binning.ptHFBinEdges_detector.size() - 1; ++i) {

        double ptDlow  = binning.ptHFBinEdges_detector[i];
        double ptDhigh = binning.ptHFBinEdges_detector[i + 1];
        double massMax = maxMass[i];

        // Compute number of bins from density and mass range
        // massBinDensity = nBins / (maxMass - minMass) → nBins = density * range
        int nMassBins = std::max(1, static_cast<int>(std::round(massBinDensity * (massMax - minMass))));

        std::cout << "pT,D bin [" << ptDlow << ", " << ptDhigh << "] GeV/c: "
                  << "mass range [" << minMass << ", " << massMax << "] GeV/c^2, "
                  << "nMassBins = " << nMassBins 
                  << " (bin width = " << (massMax - minMass) / nMassBins << " GeV/c^2)\n";

        // Build title string handling integer/fractional bin edges
        auto formatEdge = [](double val) -> std::string {
            if (std::fmod(val, 1.0) == 0.0) {
                return Form("%.0f", val);
            } else {
                return Form("%.1f", val);
            }
        };

        TString title = Form("%s < #it{p}_{T, D^{0}} < %s GeV/#it{c};"
                             "#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",
                             formatEdge(ptDlow).c_str(),
                             formatEdge(ptDhigh).c_str());

        dataContainer.histograms2d.push_back(
            new TH2D(Form("histMass%zu", i + 1), title, nMassBins, minMass, massMax, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));

        dataContainer.histograms2d[i]->Sumw2();
    }

    return dataContainer;
}
// Overloaded version with histogram name as argument for more flexibility
SidebandData createHistograms(const BinningStruct& binning, TString histName = "") {

    // Validate that maxMass vector has the correct size
    if (binning.maxMass.size() != binning.ptHFBinEdges_detector.size() - 1) {
        std::cerr << "Error: maxMass vector size (" << binning.maxMass.size() 
                  << ") does not match number of pT,D bins (" 
                  << binning.ptHFBinEdges_detector.size() - 1 << ").\n";
        return SidebandData();
    }

    SidebandData dataContainer;

    for (size_t i = 0; i < binning.ptHFBinEdges_detector.size() - 1; ++i) {

        double ptDlow  = binning.ptHFBinEdges_detector[i];
        double ptDhigh = binning.ptHFBinEdges_detector[i + 1];
        double massMin = binning.minMass[i];
        double massMax = binning.maxMass[i];

        // Compute number of bins from density and mass range
        // massBinDensity = nBins / (maxMass - minMass) → nBins = density * range
        int nMassBins = std::max(1, static_cast<int>(std::round(binning.massBinDensity * (massMax - massMin))));

        std::cout << "pT,D bin [" << ptDlow << ", " << ptDhigh << "] GeV/c: "
                  << "mass range [" << massMin << ", " << massMax << "] GeV/c^2, "
                  << "nMassBins = " << nMassBins 
                  << " (bin width = " << (massMax - massMin) / nMassBins << " GeV/c^2)\n";

        // Build title string handling integer/fractional bin edges
        auto formatEdge = [](double val) -> std::string {
            if (std::fmod(val, 1.0) == 0.0) {
                return Form("%.0f", val);
            } else {
                return Form("%.1f", val);
            }
        };

        TString title = Form("%s < #it{p}_{T, D^{0}} < %s GeV/#it{c};"
                             "#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",
                             formatEdge(ptDlow).c_str(),
                             formatEdge(ptDhigh).c_str());

        dataContainer.histograms2d.push_back(
            new TH2D(Form("%s_histMass%zu", histName.Data(), i + 1),
                     title,
                     nMassBins, massMin, massMax,
                     binning.deltaRBinEdges_detector.size() - 1,
                     binning.deltaRBinEdges_detector.data()));

        dataContainer.histograms2d[i]->Sumw2();
    }

    return dataContainer;
}


//__________________________________________________________________________________________________________________________
void fillHistograms(TFile* fDist, SidebandData& dataContainer, const double& jetptMin, const double& jetptMax, const BinningStruct& binning,
                    const std::vector<std::pair<double,double>>& bdtPtCustomCuts = {}) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on particle level delta R

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

    int nEntries = tree->GetEntries(); //  * 2 / 10
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double deltaR = axisDistance;
        
        // Fill each histogram with their respective pT intervals
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= binning.deltaRBinEdges_detector[0]) && (deltaR < DeltaRcut))) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < binning.ptHFBinEdges_detector.size() - 1 && !filled; iEdge++) {
                if ((hfPt >= binning.ptHFBinEdges_detector[iEdge]) && (hfPt < binning.ptHFBinEdges_detector[iEdge + 1])) {
                    // Get the threshold for this pT range
                    double maxBkgProb;
                    if (!bdtPtCustomCuts.empty()) {
                        maxBkgProb = GetBkgProbabilityCut(hfPt, bdtPtCustomCuts); // use custom cut for this pT,D bin
                    } else {
                        maxBkgProb = GetBkgProbabilityCut(hfPt, binning.bdtPtCuts); // normal behaviour
                    }

                    // Fill histogram only if the cut is passed
                    if (hfMlScore0 < maxBkgProb) {
                        // use Emma's run 2 cuts in case useEmmaYeatsBins is true
                        if (!binning.useEmmaYeatsBins || passEmmaCut(jetPt, hfPt)) {
                            dataContainer.histograms2d[iEdge]->Fill(hfMass, deltaR);
                        }
                        
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
        double binWidth = tempHist->GetXaxis()->GetBinWidth(1);
        tempHist->GetYaxis()->SetTitle(Form("Entries / %.0f MeV/#it{c}^{2}", binWidth*1000)); // Set y-axis title with bin width in MeV/c^2
        dataContainer.histograms1d.push_back(tempHist);
    }
}

FitContainer performFit(TFile* fReflectionsMC, SidebandData& dataContainer, const BinningStruct& binning, const FitModelType& modelToUse, const double& jetptMin, const double& jetptMax, const double& signalSigmas, const int& startingBackSigma, const int& backgroundSigmas) {
    
    double m_0_reference = dataContainer.D0_reference_mass; // D0 mass in GeV/c^2 (hard-coded)
    double sigma_reference = 0.012;

    // --- Total fits: loop through pT,HF intervals/histograms and perform fits
    for (size_t iHisto = 0; iHisto < dataContainer.histograms2d.size(); ++iHisto) {
        double minMass = dataContainer.histograms1d[iHisto]->GetXaxis()->GetXmin();
        double maxMass = dataContainer.histograms1d[iHisto]->GetXaxis()->GetXmax();

        // Get TF1 objects from MC file
        TF1* fSignal = (TF1*)fReflectionsMC->Get(Form("signalFit_%zu", iHisto));
        TF1* fReflections = (TF1*)fReflectionsMC->Get(Form("reflectionsFit_%zu", iHisto));
        // std::cout << " --- Printing signal fit parameters. ---" << std::endl;
        // fSignal->Print("V");
        // std::cout << " --- Signal fit parameters printed. ---" << std::endl;
        // std::cout << " --- Printing reflections fit parameters. ---" << std::endl;
        // fReflections->Print("V");
        // std::cout << " --- Reflections fit parameters printed. ---" << std::endl;

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
        } else if (modelToUse == FitModelType::StandardSideBand) {
            dataContainer.fittings.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<StandardSideBandSubtraction>, minMass, maxMass, 13));
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
        dataContainer.fittings.fitTotal[iHisto]->SetParLimits(0, 0.01, TMath::Infinity()); // only accepts positive background fits
        dataContainer.fittings.fitTotal[iHisto]->SetParLimits(1, -TMath::Infinity(), -1.01); // only descending inclination background fits
        dataContainer.fittings.fitTotal[iHisto]->SetParLimits(2, 3.0, TMath::Infinity()); // only accepts non-negative primary gaussian amplitude
        dataContainer.fittings.fitTotal[iHisto]->SetParLimits(4, 0.99 * m_0_reference, 1.01 * m_0_reference); // mean invariant mass should be around the one from literature, between [0.995 * m_0_reference, 1.005 * m_0_reference]
        dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.35 * sigma_reference, 4.0 * sigma_reference);
        
        if (iHisto == 0) {          // 1 < pT,D < 2 GeV/c
            sigma_reference = 0.0105;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0105
        } else if (iHisto == 1) {   // 2 < pT,D < 3 GeV/c
            sigma_reference = 0.0123;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0123
        } else if (iHisto == 2) {   // 3 < pT,D < 4 GeV/c
            sigma_reference = 0.0139;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0139
        } else if (iHisto == 3) {   // 4 < pT,D < 5 GeV/c
            sigma_reference = 0.0157;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0157
        } else if (iHisto == 4) {   // 5 < pT,D < 6 GeV/c
            sigma_reference = 0.0172;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0172
        } else if (iHisto == 5) {   // 6 < pT,D < 7 GeV/c
            sigma_reference = 0.0185;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0185
        } else if (iHisto == 6) {   // 7 < pT,D < 8 GeV/c
            sigma_reference = 0.0199;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0199
        } else if (iHisto == 7) {   // 8 < pT,D < 12 GeV/c
            sigma_reference = 0.0223;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0223
        } else if (iHisto == 8) {   // 12 < pT,D < 16 GeV/c
            sigma_reference = 0.0240;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.3 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0240
        } else if (iHisto == 9) {   // 16 < pT,D < 24 GeV/c
            sigma_reference = 0.0240;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.25 * sigma_reference); // sigma_1 = 0.0240
        } else if (iHisto == 10) {  // 24 < pT,D < 36 GeV/c
            sigma_reference = 0.0373;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.2 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0373
        }
        // sigma_reference = 0.012;
        // dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.35 * sigma_reference, 2.0 * sigma_reference);
        
        
        if (binning.useEmmaYeatsBins) {
            if (jetptMin == 5. && jetptMax == 7.) {
                if (iHisto == 4) { // 1 < pT,D < 2 GeV/c
                    dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.35 * sigma_reference, 2.5 * sigma_reference);
                }
            }
            
        }
        
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
        // dataContainer.fittings.fitTotal[iHisto]->Print("V");

        // Check if signal gaussian acomodates the literature D0 rest mass
        double primaryMean = dataContainer.fittings.fitTotal[iHisto]->GetParameter(4);
        double primarySigma = dataContainer.fittings.fitTotal[iHisto]->GetParameter(5);
        if (didFitFailed(dataContainer.histograms1d[iHisto], dataContainer.fittings.fitTotal[iHisto], m_0_reference, signalSigmas, startingBackSigma)) {
            // if literature mass is outside, clear histogram as it should be excluded
            //histograms[iHisto]->Reset();
            //dataContainer.histograms1d[iHisto]->SetTitle(TString(dataContainer.histograms1d[iHisto]->GetTitle()) + " (Fit failed)");
            
            dataContainer.fittings.workingFits.push_back(false);
        } else {
            dataContainer.fittings.workingFits.push_back(true);
            std::cout << "Fit parameters for histogram " << iHisto << ": Background A = " << dataContainer.fittings.fitTotal[iHisto]->GetParameter(0) 
                      << ", Background B = " << dataContainer.fittings.fitTotal[iHisto]->GetParameter(1) 
                      << ", A1 Signal = " << dataContainer.fittings.fitTotal[iHisto]->GetParameter(2) 
                      << ", m0 Signal = " << dataContainer.fittings.fitTotal[iHisto]->GetParameter(4) 
                      << ", Sigma1 Signal = " << dataContainer.fittings.fitTotal[iHisto]->GetParameter(5) 
                      << std::endl;
        }
    }
    std::cout << "Total fits performed.\n";

    // --- Background only fits: loop through pT,HF intervals/histograms and perform fits
    // Perform background only fit to each histogram
    int numberOfFits = dataContainer.fittings.fitTotal.size();
    for (size_t iHisto = 0; iHisto < numberOfFits; iHisto++) {
        double minMass = dataContainer.histograms1d[iHisto]->GetXaxis()->GetXmin();
        double maxMass = dataContainer.histograms1d[iHisto]->GetXaxis()->GetXmax();
        if (modelToUse == FitModelType::FullPowerLaw || modelToUse == FitModelType::StandardSideBand) { // backgroundFunctionPowerLaw
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
        double minMass = dataContainer.histograms1d[iHisto]->GetXaxis()->GetXmin();
        double maxMass = dataContainer.histograms1d[iHisto]->GetXaxis()->GetXmax();
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

    // Significance calculation from fit integrals
    for (size_t iHisto = 0; iHisto < dataContainer.fittings.fitSignalOnly.size(); ++iHisto) {
        double sig = 0., S = 0., B = 0.;

        // Only compute for working fits that have both total and background fits
        bool hasSignalFit = (iHisto < dataContainer.fittings.fitSignalOnly.size()) && dataContainer.fittings.fitSignalOnly[iHisto];
        bool hasBkgFit = (iHisto < dataContainer.fittings.fitBackgroundOnly.size()) && dataContainer.fittings.fitBackgroundOnly[iHisto];
        bool isWorking = (iHisto < dataContainer.fittings.workingFits.size()) && dataContainer.fittings.workingFits[iHisto];

        if (isWorking && hasSignalFit && hasBkgFit) {
            TF1* fSignal = dataContainer.fittings.fitSignalOnly[iHisto];
            TF1* fBkg = dataContainer.fittings.fitBackgroundOnly[iHisto];
            double m0 = fSignal->GetParameter(2); // total->4, pure signal->2
            double sigma = fSignal->GetParameter(3); // total->5, pure signal->3

            // Calculate areas
            double totalInSR = fSignal->Integral(m0 - signalSigmas*sigma, m0 + signalSigmas*sigma);
            double bkgInSR = fBkg->Integral(  m0 - signalSigmas*sigma, m0 + signalSigmas*sigma);
            S = totalInSR;
            B = bkgInSR;
            if (S > 0. && (S + B) > 0.) {
                sig = S / std::sqrt(S + B);
            }
            std::cout << "  [Significance] pT,D bin " << iHisto << ": S=" << S << "  B=" << B << "  sig=" << sig << "\n";
        } else {
            std::cout << "  [Significance] pT,D bin " << iHisto << ": skipped (fit failed or missing)\n";
        }

        dataContainer.fittings.significance.push_back(sig);
        dataContainer.fittings.signalYield.push_back(S);
        dataContainer.fittings.backgroundYield.push_back(B);
    }
    std::cout << "Significance computed for all pT,D bins.\n";

    // Perform reflections only fit to each histogram
    for (size_t iHisto = 0; iHisto < numberOfFits; iHisto++) {
        double minMass = dataContainer.histograms1d[iHisto]->GetXaxis()->GetXmin();
        double maxMass = dataContainer.histograms1d[iHisto]->GetXaxis()->GetXmax();
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

//__________________________________________________________________________________________________________________________
struct BDTScanResult {
    int nPtDBins;
    std::vector<std::vector<double>> bdtThresholds;   // [iScan][iPtD] ← change to 2D
    std::vector<std::vector<double>> significance;    // [iScan][iPtD]
    std::vector<double>              optimalCuts;     // [iPtD]
    std::vector<double>              maxSignificance; // [iPtD]
    std::vector<std::vector<double>> bdtThresholdsFine;  // [iScan][iPtD]
    std::vector<std::vector<double>> significanceFine;   // [iScan][iPtD]
};
void plotBDTScanFits(const SidebandData& dataContainer, const FitContainer& fittings, const std::vector<std::pair<double,double>>& bdtCuts, const BinningStruct& binning,
                     int scanStep, bool isFine, const TString& pdfDir, const FitModelType& modelToUse, double signalSigmas, int startingBackSigma) {

    const int nPtDBins = dataContainer.histograms1d.size();
    int nCols = std::min(4, nPtDBins);
    int nRows = std::ceil((double)nPtDBins / nCols);

    TString stepLabel = isFine ? Form("Fine step %d", scanStep+1)
                               : Form("Coarse step %d", scanStep+1);

    TCanvas* c = new TCanvas(Form("cBDTFits_%s_%d", isFine?"fine":"coarse", scanStep),
                              stepLabel, 400*nCols, 350*nRows);
    c->Divide(nCols, nRows);

    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.038);

    for (int iPtD = 0; iPtD < nPtDBins; ++iPtD) {
        c->cd(iPtD + 1);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.15);
        gStyle->SetOptStat(0);

        TH1D* h = dataContainer.histograms1d[iPtD];
        if (!h) continue;

        // Draw histogram
        h->SetMarkerStyle(kDot);
        h->SetMarkerColor(kBlack);
        h->SetLineColor(kBlack);
        h->SetMinimum(0);
        h->Draw("E");

        // Draw fits if working
        bool isWorking = (iPtD < (int)fittings.workingFits.size())
                       && fittings.workingFits[iPtD];

        if (isWorking) {
            TF1* fTotal = (iPtD < (int)fittings.fitTotal.size())
                        ? fittings.fitTotal[iPtD] : nullptr;
            TF1* fBkg   = (iPtD < (int)fittings.fitBackgroundOnly.size())
                        ? fittings.fitBackgroundOnly[iPtD] : nullptr;
            TF1* fSig   = (iPtD < (int)fittings.fitSignalOnly.size())
                        ? fittings.fitSignalOnly[iPtD] : nullptr;
            TF1* fRefl  = (iPtD < (int)fittings.fitReflectionsOnly.size())
                        ? fittings.fitReflectionsOnly[iPtD] : nullptr;

            // Compute shaded sideband and signal regions from fit parameters
            if (fTotal) {
                double m0    = fTotal->GetParameter(4);
                double sigma = fTotal->GetParameter(5);
                double binW  = h->GetBinWidth(1);
                double yMax  = h->GetMaximum();

                // Signal region shading (blue)
                double srLow  = m0 - signalSigmas * sigma;
                double srHigh = m0 + signalSigmas * sigma;
                int srLowBin  = safeLowBin(h, srLow);
                int srHighBin = safeHighBin(h, srHigh);
                TH1D* hSigShade = (TH1D*)h->Clone(
                    Form("hSigShade_%d_%d", iPtD, scanStep));
                hSigShade->Reset();
                for (int iBin = srLowBin; iBin <= srHighBin; ++iBin) {
                    hSigShade->SetBinContent(iBin, h->GetBinContent(iBin));
                }
                hSigShade->SetFillColorAlpha(kBlue, 0.2);
                hSigShade->SetLineColor(kBlue);
                hSigShade->Draw("hist same");

                // Sideband region shading (red)
                // Left sideband: [m0 - sbEnd*sigma, m0 - startingBackSigma*sigma]
                // Right sideband: [m0 + startingBackSigma*sigma, m0 + sbEnd*sigma]
                // Clamp to histogram range
                double sbLow1  = std::max(m0 - 9.*sigma, h->GetXaxis()->GetXmin());
                double sbHigh1 = m0 - startingBackSigma * sigma;
                double sbLow2  = m0 + startingBackSigma * sigma;
                double sbHigh2 = std::min(m0 + 9.*sigma, h->GetXaxis()->GetXmax());

                TH1D* hSBShade = (TH1D*)h->Clone(
                    Form("hSBShade_%d_%d", iPtD, scanStep));
                hSBShade->Reset();
                int sb1Low  = safeLowBin(h, sbLow1);
                int sb1High = safeHighBin(h, sbHigh1);
                int sb2Low  = safeLowBin(h, sbLow2);
                int sb2High = safeHighBin(h, sbHigh2);
                for (int iBin = sb1Low; iBin <= sb1High; ++iBin) {
                    hSBShade->SetBinContent(iBin, h->GetBinContent(iBin));
                }
                for (int iBin = sb2Low; iBin <= sb2High; ++iBin) {
                    hSBShade->SetBinContent(iBin, h->GetBinContent(iBin));
                }
                hSBShade->SetFillColorAlpha(kRed, 0.2);
                hSBShade->SetLineColor(kRed);
                hSBShade->Draw("hist same");

                // Draw fit components — same style as actual analysis
                fTotal->SetLineColor(kBlack);
                fTotal->SetLineStyle(kSolid);
                fTotal->SetLineWidth(1);
                fTotal->Draw("same");

                if (fSig) {
                    fSig->SetLineColor(kBlue);
                    fSig->SetLineStyle(kDashed);
                    fSig->SetLineWidth(1);
                    fSig->Draw("same");
                }
                if (fBkg) {
                    fBkg->SetLineColor(kRed);
                    fBkg->SetLineStyle(kDashed);
                    fBkg->SetLineWidth(1);
                    fBkg->Draw("same");
                }
                if (fRefl && modelToUse != FitModelType::StandardSideBand) {
                    fRefl->SetLineColor(kGreen+2);
                    fRefl->SetLineStyle(kDashed);
                    fRefl->SetLineWidth(1);
                    fRefl->Draw("same");
                }

                // Annotations
                double sigma2 = sigma / fTotal->GetParameter(6);
                double chi2   = fTotal->GetChisquare();
                int    ndf    = fTotal->GetNDF();
                latex->SetTextColor(kBlack);
                latex->DrawLatex(0.14, 0.89,
                    Form("#Chi^{2}_{red} = %.3f", ndf > 0 ? chi2/ndf : -1.));
                latex->DrawLatex(0.14, 0.84,
                    Form("#sigma_{primary} = %.4f", sigma));
                latex->DrawLatex(0.14, 0.79,
                    Form("#sigma_{secondary} = %.4f", sigma2));
            }

            // Significance annotation
            double cut = (iPtD < (int)bdtCuts.size()) ? bdtCuts[iPtD].second : 0.;
            double sig = (iPtD < (int)fittings.significance.size())
                       ? fittings.significance[iPtD] : 0.;
            double S   = (iPtD < (int)fittings.signalYield.size())
                       ? fittings.signalYield[iPtD] : 0.;
            double B   = (iPtD < (int)fittings.backgroundYield.size())
                       ? fittings.backgroundYield[iPtD] : 0.;
            latex->SetTextColor(kBlue);
            latex->DrawLatex(0.14, 0.73, Form("BDT cut = %.4f", cut));
            latex->DrawLatex(0.14, 0.68, Form("S = %.0f, B = %.0f", S, B));
            latex->DrawLatex(0.14, 0.63, Form("Sig = %.2f", sig));

        } else {
            // Fit failed annotation
            double cut = (iPtD < (int)bdtCuts.size()) ? bdtCuts[iPtD].second : 0.;
            latex->SetTextColor(kRed);
            latex->DrawLatex(0.14, 0.73, Form("BDT cut = %.4f", cut));
            latex->DrawLatex(0.14, 0.68, "Fit failed");
        }

        // pT,D label
        latex->SetTextColor(kBlack);
        latex->DrawLatex(0.14, 0.57,
            Form("%.0f < p_{T,D} < %.0f GeV/c",
                 binning.ptHFBinEdges_detector[iPtD],
                 binning.ptHFBinEdges_detector[iPtD+1]));

        // Redraw histogram on top of shading
        h->Draw("E same");
    }

    // Save
    TString scanType = isFine ? "fine" : "coarse";
    TString fileName = pdfDir + Form("BDT_scan_%s_step%02d.pdf",
                                      scanType.Data(), scanStep + 1);
    for (int i = 1; i <= nCols * nRows; ++i) {
        TVirtualPad* pad = c->GetPad(i);
        if (pad) { pad->Modified(); pad->Update(); }
    }
    c->cd(0);
    c->Modified();
    c->Update();
    c->SaveAs(fileName);
    delete c;
    delete latex;

    std::cout << "[BDT Optimization] Saved: " << fileName << "\n";
}
void plotSignificanceScan(const BDTScanResult& scanResult, const BinningStruct& binning, const TString& pdfDir) { // ← directory, not file

    const int nPtDBins = scanResult.nPtDBins;
    const int nScans   = scanResult.bdtThresholds.size();

    // One pad per pT,D bin
    int nCols = std::min(4, nPtDBins);
    int nRows = std::ceil((double)nPtDBins / nCols);
    TCanvas* cOptimalBdtCuts = new TCanvas("cBDTScan", "BDT Significance Scan", 400*nCols, 350*nRows);
    cOptimalBdtCuts->Divide(nCols, nRows);

    for (int iPtD = 0; iPtD < nPtDBins; ++iPtD) {
        cOptimalBdtCuts->cd(iPtD + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

        // Build TGraph: x = BDT threshold, y = significance
        TGraph* graph = new TGraph(nScans);
        graph->SetName(Form("gBDTScan_ptD%d", iPtD));
        graph->SetTitle(Form("%.0f < p_{T,D^{0}} < %.0f GeV/c;" "BDT background score threshold;" "S/#sqrt{S+B}", binning.ptHFBinEdges_detector[iPtD], binning.ptHFBinEdges_detector[iPtD+1]));

        for (int iScan = 0; iScan < nScans; ++iScan) {
            double threshold = scanResult.bdtThresholds[iScan][iPtD];
            double sig = (iPtD < (int)scanResult.significance[iScan].size()) ? scanResult.significance[iScan][iPtD] : 0.;
            graph->SetPoint(iScan, threshold, sig);
        }

        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.2);
        graph->SetMarkerColor(kBlue+1);
        graph->SetLineColor(kBlue+1);
        graph->SetLineWidth(1);
        graph->Draw("APL");

        // Mark the optimal cut with a vertical line
        double optCut = (iPtD < (int)scanResult.optimalCuts.size())
                      ? scanResult.optimalCuts[iPtD] : 0.;
        double maxSig = (iPtD < (int)scanResult.maxSignificance.size())
                      ? scanResult.maxSignificance[iPtD] : 0.;

        if (optCut > 0. && graph->GetHistogram()) {
            double yMin = graph->GetHistogram()->GetMinimum();
            double yMax = graph->GetHistogram()->GetMaximum() * 1.15;
            graph->GetHistogram()->GetYaxis()->SetRangeUser(yMin, yMax);

            TLine* lOpt = new TLine(optCut, yMin, optCut, yMax);
            lOpt->SetLineColor(kRed);
            lOpt->SetLineWidth(1);
            lOpt->SetLineStyle(2);
            lOpt->Draw("same");

            // Label the optimal cut value
            TLatex* tex = new TLatex();
            tex->SetNDC();
            tex->SetTextSize(0.04);
            tex->SetTextColor(kRed);
            tex->DrawLatex(0.18, 0.82, Form("Opt. cut = %.4f", optCut));
            tex->DrawLatex(0.18, 0.74, Form("Max sig = %.2f", maxSig));
        }

        // After drawing coarse scan graph gCoarse:
        int nFine = scanResult.bdtThresholdsFine.size();
        if (nFine > 0) {
            TGraph* gFine = new TGraph(nFine);
            for (int j = 0; j < nFine; ++j) {
                gFine->SetPoint(j,
                    scanResult.bdtThresholdsFine[j][iPtD],
                    scanResult.significanceFine[j][iPtD]);
            }
            gFine->SetMarkerStyle(21);
            gFine->SetMarkerColor(kRed);
            gFine->SetLineColor(kRed);
            gFine->Draw("P same"); // points only, no connecting line
        }
    }

    
    //
    // Storing in a single pdf file
    //
    TString imagePath = "../../Images/1-SignalTreatment/SideBand/";
    double jetptMin = binning.ptjetBinEdges_detector.front();
    double jetptMax = binning.ptjetBinEdges_detector.back();
    TString fileName = pdfDir + "BDT_significance_scan.pdf";

    for (int i = 1; i <= nCols * nRows; ++i) {
        TVirtualPad* pad = cOptimalBdtCuts->GetPad(i);
        if (pad) { pad->Modified(); pad->Update(); }
    }
    cOptimalBdtCuts->cd(0);
    cOptimalBdtCuts->Modified();
    cOptimalBdtCuts->Update();
    cOptimalBdtCuts->SaveAs(fileName);
    delete cOptimalBdtCuts;

    std::cout << "[BDT Optimization] Saved: " << fileName << "\n";
}

BDTScanResult optimizeBDTCuts(TFile* fDist, TFile* fReflectionsMCFullRange, BinningStruct& binning, const FitModelType& modelToUse, double& signalSigmas, int& startingBackSigma, int& backgroundSigmas) {
    
    BDTScanResult scanResult;
    
    // 0 - Check if file with optimal BDT cuts already exists. If yes, read the cuts from the file and skip the optimization procedure. If not, proceed with the optimization and save the results in a file for future use.
    // --- 0: Check if file already exists ---
    const std::string resultFile = "bdtOptimalCuts.root"; // outputPath + "/bdtOptimalCuts.root"
    if (!gSystem->AccessPathName(resultFile.c_str())) {
        std::cout << "[BDT Optimization] Found existing file: " << resultFile << ". Loading cuts into binning struct...\n";
        TFile* fCuts = TFile::Open(resultFile.c_str(), "READ");
        TVectorD* vEdges = (TVectorD*)fCuts->Get("bdt/bdtPtEdges");
        TVectorD* vCuts  = (TVectorD*)fCuts->Get("bdt/bdtCutValues");
        if (!vCuts || !vEdges) {
            std::cerr << "[BDT Optimization] Error: could not read optimal cuts from file.\n";
            fCuts->Close();
            return scanResult;
        }
        // Direct reconstruction — no special casing needed
        int nFullEntries = vCuts->GetNoElements(); // = 13
        binning.bdtPtCuts.clear();
        for (int i = 0; i < nFullEntries; ++i) {
            binning.bdtPtCuts.emplace_back((*vEdges)[i], (*vCuts)[i]);
        }
        std::cout << "[BDT Optimization] Loaded bdtPtCuts with "
                << binning.bdtPtCuts.size() << " entries.\n";
        fCuts->Close();
        return scanResult;
    }
    std::cout << "[BDT Optimization] No existing file found. Running optimization...\n";
    
    // --- Configuration ---
    const int nPtHFBins = binning.ptHFBinEdges_detector.size() - 1;
    const int nCoarseScanSteps = 20;  // coarse scan first: 10 steps
    const int nFineScan = 5;  // fine scan: 20 steps around the best region
    // Use integrated pT,jet range for optimization
    const double jetptMin = binning.ptjetBinEdges_detector.front();
    const double jetptMax = binning.ptjetBinEdges_detector.back();

    TString imagePath = "../../Images/1-SignalTreatment/SideBand/";
    // Create directory if it doesn't exist
    TString pdfDir = imagePath + Form("BDT_scan_%.0f_to_%.0fGeV/", jetptMin, jetptMax);
    gSystem->Exec(Form("mkdir -p %s", pdfDir.Data()));

    // 1 - Loop through course scan of BDT scores range
    for (size_t iCoarseScore = 0; iCoarseScore < nCoarseScanSteps; iCoarseScore++) {
        // 2 - Create histograms for this score cut
        SidebandData dataContainer = createHistograms(binning, Form("bdt_optimal_%zu", iCoarseScore));

        // 3 - Build custom BDT cuts
        std::vector<std::pair<double,double>> bdtPtCustomCuts;
        // loops over analysis pT,D bins (11 entries, starting from pT=1):
        for (int iPtHFBin = 0; iPtHFBin < nPtHFBins; ++iPtHFBin) {
            double ptDcenter = 0.5 * (binning.ptHFBinEdges_detector[iPtHFBin] + binning.ptHFBinEdges_detector[iPtHFBin+1]);
            double maxCut = GetBkgProbabilityCut(ptDcenter, binning.bdtPtCuts);
            // double cut = (iCoarseScore + 1) * (binning.bdtPtCuts[iPtHFBin].second / nCoarseScanSteps); // fraction varies with scan step
            // fraction = (iCoarseScore + 1) / nCoarseScanSteps
            double cut = (iCoarseScore + 1) * maxCut / nCoarseScanSteps;
            bdtPtCustomCuts.push_back({binning.ptHFBinEdges_detector[iPtHFBin], cut});
            //bdtPtCustomCuts.push_back({binning.bdtPtCuts[iPtHFBin].first, (iCoarseScore + 1) * (binning.bdtPtCuts[iPtHFBin].second / nCoarseScanSteps)});
        }
        bdtPtCustomCuts.push_back({50., -1});

        // 4 - Fill histograms with custom BDT cuts
        fillHistograms(fDist, dataContainer, jetptMin, jetptMax, binning, bdtPtCustomCuts);

        // 5 - Perform fits and calculate significance for this set of score cuts
        FitContainer fittings = performFit(fReflectionsMCFullRange, dataContainer, binning, modelToUse, jetptMin, jetptMax, signalSigmas, startingBackSigma, backgroundSigmas);

        // Open PDF on first page, middle on others: invariant mass plots for this scan step
        plotBDTScanFits(dataContainer, fittings, bdtPtCustomCuts, binning, iCoarseScore, false, pdfDir, modelToUse, signalSigmas, startingBackSigma);

        // 6 - Store significance and threshold for this scan step
        std::vector<double> thresholdsThisStep;
        for (int iPtD = 0; iPtD < nPtHFBins; ++iPtD) {
            thresholdsThisStep.push_back(bdtPtCustomCuts[iPtD].second);
        }
        scanResult.bdtThresholds.push_back(thresholdsThisStep);
        scanResult.significance.push_back(fittings.significance);

        // Print results for this step
        std::cout << "[BDT Optimization] Coarse step " << iCoarseScore+1
                  << "/" << nCoarseScanSteps << "\n";
        for (int iPtD = 0; iPtD < nPtHFBins; ++iPtD) {
            std::cout << "  pT,D bin " << iPtD
                      << " [" << binning.ptHFBinEdges_detector[iPtD]
                      << "-"  << binning.ptHFBinEdges_detector[iPtD+1] << "]"
                      << "  cut=" << bdtPtCustomCuts[iPtD].second
                      << "  sig=" << fittings.significance[iPtD] << "\n";
        }

        // Cleanup histograms and fits for this step
        for (auto* h : dataContainer.histograms2d) { if (h) delete h; }
        for (auto* h : dataContainer.histograms1d) { if (h) delete h; }
        for (auto* f : fittings.fitTotal)          { if (f) delete f; }
        for (auto* f : fittings.fitBackgroundOnly) { if (f) delete f; }
        for (auto* f : fittings.fitSignalOnly)      { if (f) delete f; }
        for (auto* f : fittings.fitReflectionsOnly) { if (f) delete f; }
    } // end coarse scan loop

    // 7 - Find the best coarse step per pT,HF bin
    // Initialize optimal cuts and significances
    scanResult.nPtDBins       = nPtHFBins;
    scanResult.optimalCuts    = std::vector<double>(nPtHFBins, 0.);
    scanResult.maxSignificance = std::vector<double>(nPtHFBins, 0.);

    std::vector<int> bestCoarseStep(nPtHFBins, 0);
    for (int iPtD = 0; iPtD < nPtHFBins; ++iPtD) {
        for (int iScan = 0; iScan < nCoarseScanSteps; ++iScan) {
            if (scanResult.significance[iScan][iPtD] >
                scanResult.significance[bestCoarseStep[iPtD]][iPtD]) {
                bestCoarseStep[iPtD] = iScan;
            }
        }
        // Store best coarse result as starting point
        scanResult.optimalCuts[iPtD] = scanResult.bdtThresholds[bestCoarseStep[iPtD]][iPtD];
        scanResult.maxSignificance[iPtD] = scanResult.significance[bestCoarseStep[iPtD]][iPtD];

        std::cout << "[BDT Optimization] pT,D bin " << iPtD
                  << ": best coarse step = " << bestCoarseStep[iPtD]
                  << ", cut = " << scanResult.optimalCuts[iPtD]
                  << ", sig = " << scanResult.maxSignificance[iPtD] << "\n";
    }

    // 8 - Fine scan: each pT,D bin scans independently around its best coarse step
    std::cout << "[BDT Optimization] === Fine scan ===\n";
    for (int iFine = 0; iFine <= nFineScan; ++iFine) {

        // Build fine cuts — each pT,D bin uses its own fine range
        std::vector<std::pair<double,double>> bdtPtCustomCuts;
        for (int iPtD = 0; iPtD < nPtHFBins; ++iPtD) {
            // Fine range spans one coarse step below and above the best
            double fineMin = (bestCoarseStep[iPtD] > 0) ? scanResult.bdtThresholds[bestCoarseStep[iPtD] - 1][iPtD] : 0.;
            double fineMax = (bestCoarseStep[iPtD] < nCoarseScanSteps - 1) ? scanResult.bdtThresholds[bestCoarseStep[iPtD] + 1][iPtD] : binning.bdtPtCuts[iPtD].second;
            double cut = fineMin + (iFine + 1) * (fineMax - fineMin) / (nFineScan + 2);
            bdtPtCustomCuts.push_back({binning.bdtPtCuts[iPtD].first, cut});
        }
        bdtPtCustomCuts.push_back({50., -1});

        SidebandData dataContainer = createHistograms(binning, Form("bdt_fine_%d", iFine));
        fillHistograms(fDist, dataContainer, jetptMin, jetptMax, binning, bdtPtCustomCuts);
        FitContainer fittings = performFit(fReflectionsMCFullRange, dataContainer, binning, modelToUse, jetptMin, jetptMax, signalSigmas, startingBackSigma, backgroundSigmas);

        std::cout << "[BDT Optimization] Fine step " << iFine+1 << "/" << nFineScan+1 << "\n";

        plotBDTScanFits(dataContainer, fittings, bdtPtCustomCuts, binning, iFine, true, pdfDir, modelToUse, signalSigmas, startingBackSigma);

        // Update optimal cut per pT,D bin if this step gives better significance
        for (int iPtD = 0; iPtD < nPtHFBins; ++iPtD) {
            double sig = (iPtD < (int)fittings.significance.size()) ? fittings.significance[iPtD] : 0.;
            double cut = bdtPtCustomCuts[iPtD].second;
            std::cout << "  pT,D bin " << iPtD << "  cut=" << cut << "  sig=" << sig << "\n";
            if (sig > scanResult.maxSignificance[iPtD]) {
                scanResult.maxSignificance[iPtD] = sig;
                scanResult.optimalCuts[iPtD]     = cut;
            }
        }

        for (auto* h : dataContainer.histograms2d) { if (h) delete h; }
        for (auto* h : dataContainer.histograms1d) { if (h) delete h; }
        for (auto* f : fittings.fitTotal)          { if (f) delete f; }
        for (auto* f : fittings.fitBackgroundOnly) { if (f) delete f; }
        for (auto* f : fittings.fitSignalOnly)      { if (f) delete f; }
        for (auto* f : fittings.fitReflectionsOnly) { if (f) delete f; }
    } // end fine scan loop

    // 9 - Add significance scan plots as final pages
    plotSignificanceScan(scanResult, binning, pdfDir);

    // 10 - Build full bdtPtCuts FIRST
    binning.bdtPtCuts.clear();
    binning.bdtPtCuts.emplace_back(0., -1.);                                         // pre-range
    for (int i = 0; i < nPtHFBins + 1; ++i) {
        binning.bdtPtCuts.emplace_back(binning.ptHFBinEdges_detector[i], scanResult.optimalCuts[i]);
    }
    binning.bdtPtCuts.emplace_back(50., -1.0); // sentinel

    // 11 - Save optimal cuts to file
    TFile* fOut = TFile::Open(resultFile.c_str(), "RECREATE");
    TDirectory* dBDT = fOut->mkdir("bdt");
    dBDT->cd();

    // Store the FULL bdtPtCuts structure including pre-range and sentinel
    // At this point binning.bdtPtCuts already has the correct structure from step 11
    int nFullEntries = binning.bdtPtCuts.size(); // = nPtHFBins + 2 = 13
    TVectorD vEdges(nFullEntries);
    TVectorD vCuts(nFullEntries);
    for (int i = 0; i < nFullEntries; ++i) {
        vEdges[i] = binning.bdtPtCuts[i].first;
        vCuts[i]  = binning.bdtPtCuts[i].second;
    }
    vEdges.Write("bdtPtEdges");
    vCuts.Write("bdtCutValues");

    // Also save max significance for the analysis bins only
    TVectorD vSig(nPtHFBins);
    for (int i = 0; i < nPtHFBins; ++i) {
        vSig[i] = scanResult.maxSignificance[i];
    }
    vSig.Write("maxSignificance");
    fOut->Close();
    std::cout << "[BDT Optimization] Optimal cuts saved to: " << resultFile << "\n";

    // Print summary table
    std::cout << "\n[BDT Optimization] === Final summary ===\n";
    std::cout << std::setw(20) << "pT,D bin" << std::setw(15) << "Optimal cut" << std::setw(15) << "Max sig\n";
    for (int i = 0; i < nPtHFBins; ++i) {
        std::cout << std::setw(6)  << binning.ptHFBinEdges_detector[i]
                  << " - " << std::setw(6) << binning.ptHFBinEdges_detector[i+1]
                  << std::setw(15) << std::fixed << std::setprecision(4)
                  << scanResult.optimalCuts[i]
                  << std::setw(15) << std::fixed << std::setprecision(2)
                  << scanResult.maxSignificance[i] << "\n";
    }
  
    std::cout << "[BDT Optimization] binning.bdtPtCuts updated with optimal cuts.\n";

    return scanResult;
}

void createShadedRegions(SidebandData& dataContainer, const std::pair<std::array<double, 2>, std::array<double, 2>>& sidebandRanges, const int& iHisto, const double& signalSigmas, const double& sigma) {
    
    // Create shaded regions for the sidebands
    TH1D* hInvMassLeftSidebandDrawing;
    TH1D* hInvMassSidebandDrawing;

    // Left sideband
    if (!(sidebandRanges.first[0] == 0.) && !(sidebandRanges.first[1] == 0.)) {
        double leftSidebandLowBin = dataContainer.histograms1d[iHisto]->GetXaxis()->FindBin(sidebandRanges.first[0]);
        double leftSidebandHighBin = dataContainer.histograms1d[iHisto]->GetXaxis()->FindBin(sidebandRanges.first[1]);
        hInvMassLeftSidebandDrawing = (TH1D*)dataContainer.histograms1d[iHisto]->Clone(Form("hSidebandHistDrawingLeft_%d", iHisto));
        hInvMassLeftSidebandDrawing->Reset();
        for (int iBin = leftSidebandLowBin; iBin <= leftSidebandHighBin; ++iBin) {
            hInvMassLeftSidebandDrawing->SetBinContent(iBin, dataContainer.histograms1d[iHisto]->GetBinContent(iBin));
            hInvMassLeftSidebandDrawing->SetBinError(iBin, dataContainer.histograms1d[iHisto]->GetBinError(iBin));
        }
        hInvMassLeftSidebandDrawing->SetFillColorAlpha(kRed, 0.2);
        hInvMassLeftSidebandDrawing->SetLineColor(kRed);

        // Right sideband
        if (!(sidebandRanges.second[0] == 0.) && !(sidebandRanges.second[1] == 0.)) {
            double rightSidebandLowBin = dataContainer.histograms1d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[0]);
            double rightSidebandHighBin = dataContainer.histograms1d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[1]);
            hInvMassSidebandDrawing = (TH1D*)dataContainer.histograms1d[iHisto]->Clone(Form("hSidebandHistDrawingRight_%d", iHisto));
            hInvMassSidebandDrawing->Reset();
            for (int iBin = rightSidebandLowBin; iBin <= rightSidebandHighBin; ++iBin) {
                hInvMassSidebandDrawing->SetBinContent(iBin, dataContainer.histograms1d[iHisto]->GetBinContent(iBin));
                hInvMassSidebandDrawing->SetBinError(iBin, dataContainer.histograms1d[iHisto]->GetBinError(iBin));
            }
            hInvMassSidebandDrawing->SetFillColorAlpha(kRed, 0.2);
            hInvMassSidebandDrawing->SetLineColor(kRed);

            // Add both sidebands to the container for drawing
            hInvMassSidebandDrawing->Add(hInvMassLeftSidebandDrawing);
            dataContainer.subtractionResults.sidebandHistDrawing.push_back(hInvMassSidebandDrawing);
        } else {
            dataContainer.subtractionResults.sidebandHistDrawing.push_back(hInvMassLeftSidebandDrawing);
        }

    } else {
        // Right sideband
        if (!(sidebandRanges.second[0] == 0.) && !(sidebandRanges.second[1] == 0.)) {
            double rightSidebandLowBin = dataContainer.histograms1d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[0]);
            double rightSidebandHighBin = dataContainer.histograms1d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[1]);
            hInvMassSidebandDrawing = (TH1D*)dataContainer.histograms1d[iHisto]->Clone(Form("hSidebandHistDrawingRight_%d", iHisto));
            hInvMassSidebandDrawing->Reset();
            for (int iBin = rightSidebandLowBin; iBin <= rightSidebandHighBin; ++iBin) {
                hInvMassSidebandDrawing->SetBinContent(iBin, dataContainer.histograms1d[iHisto]->GetBinContent(iBin));
                hInvMassSidebandDrawing->SetBinError(iBin, dataContainer.histograms1d[iHisto]->GetBinError(iBin));
            }
            hInvMassSidebandDrawing->SetFillColorAlpha(kRed, 0.2);
            hInvMassSidebandDrawing->SetLineColor(kRed);
            dataContainer.subtractionResults.sidebandHistDrawing.push_back(hInvMassSidebandDrawing);
        } else {
            std::cout << "Error: No valid sideband ranges provided for histogram index " << iHisto << ".\n";
        }
    }
    
    // Create shaded region for the signal
    double m_0 = dataContainer.fittings.fitTotal[iHisto]->GetParameter(4); // Get the value of parameter 'm_0' from the total fit
    double signalLowBin = dataContainer.histograms1d[iHisto]->GetXaxis()->FindBin(m_0 - signalSigmas * sigma);
    double signalHighBin = dataContainer.histograms1d[iHisto]->GetXaxis()->FindBin(m_0 + signalSigmas * sigma);
    TH1D* hInvMassSignalDrawing = (TH1D*)dataContainer.histograms1d[iHisto]->Clone(Form("hSignalHistDrawing_%d", iHisto));
    hInvMassSignalDrawing->Reset();
    for (int iBin = signalLowBin; iBin <= signalHighBin; ++iBin) {
        hInvMassSignalDrawing->SetBinContent(iBin, dataContainer.histograms1d[iHisto]->GetBinContent(iBin));
        hInvMassSignalDrawing->SetBinError(iBin, dataContainer.histograms1d[iHisto]->GetBinError(iBin));
    }
    hInvMassSignalDrawing->SetFillColorAlpha(kBlue, 0.2);
    hInvMassSignalDrawing->SetLineColor(kBlue);
    dataContainer.subtractionResults.signalHistDrawing.push_back(hInvMassSignalDrawing);
    
}
// To store information per jet pT interval
struct JetPtContainerSideband {
    std::vector<SidebandData> jetPtTotalSidebandData;
};
SubtractionResult SideBand(SidebandData& dataContainer, const FitModelType& modelToUse, const BinningStruct& binning, double& signalSigmas, int& startingBackSigma, int& backgroundSigmas, JetPtContainerSideband& jetPtContainer) {
    
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
        if (dataContainer.histograms1d[iHisto]->GetEntries() == 0) {
            scallingFactors = {0., 0.};
        }
        
        createShadedRegions(dataContainer,sidebandRanges, iHisto, signalSigmas, sigma);

        dataContainer.subtractionResults.scalingFactorsArrays.push_back(scallingFactors);

        // Create DeltaR signal histogram
        // int lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(m_0 - signalSigmas * sigma);
        // int highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(m_0 + signalSigmas * sigma);
        int lowBin = safeLowBin(dataContainer.histograms2d[iHisto], m_0 - signalSigmas * sigma);
        int highBin = safeHighBin(dataContainer.histograms2d[iHisto], m_0 + signalSigmas * sigma);
        h_signal = dataContainer.histograms2d[iHisto]->ProjectionY(Form("h_signal_proj_%zu",iHisto), lowBin, highBin);
        dataContainer.subtractionResults.signalHist.push_back(h_signal);

        // Create DeltaR side-band histogram
        if (sidebandRanges.first[0] == 0. && sidebandRanges.first[1] == 0.) {
            // Use only the right sideband
            std::cout << "Using only the right sideband..." << std::endl;
            // lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[0]);
            // highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[1]);
            lowBin = safeLowBin(dataContainer.histograms2d[iHisto], sidebandRanges.second[0]);
            highBin = safeHighBin(dataContainer.histograms2d[iHisto], sidebandRanges.second[1]);
            h_sideBand = dataContainer.histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_%zu",iHisto), lowBin, highBin); // sum the right sideband
            double rightSBHistogram = dataContainer.histograms1d[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram

            double rightSB = dataContainer.fittings.fitTotal[iHisto]->Integral(sidebandRanges.second[0],sidebandRanges.second[1]);
            double fitSigIntegral = dataContainer.fittings.fitTotal[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

            // lowBin = dataContainer.histograms1d[iHisto]->FindBin(m_0 - signalSigmas * sigma);
            // highBin = dataContainer.histograms1d[iHisto]->FindBin(m_0 + signalSigmas * sigma);
            lowBin = safeLowBin(dataContainer.histograms1d[iHisto], m_0 - signalSigmas * sigma);
            highBin = safeHighBin(dataContainer.histograms1d[iHisto], m_0 + signalSigmas * sigma);
            double signalHistogram = dataContainer.histograms1d[iHisto]->Integral(lowBin, highBin);
        } else {
            // Use both sidebands
            std::cout << "Using both sidebands..." << std::endl;
            // lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.first[0]);
            // highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.first[1]);
            lowBin = safeLowBin(dataContainer.histograms2d[iHisto], sidebandRanges.first[0]);
            highBin = safeHighBin(dataContainer.histograms2d[iHisto], sidebandRanges.first[1]);
            h_sideBand = dataContainer.histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_left_%zu",iHisto), lowBin, highBin); // sum the left sideband
            double leftSBHistogram = dataContainer.histograms1d[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram
            // std::cout << "Sideband integral (only left) = " << h_sideBand->Integral() << std::endl;
            // std::cout << "leftSBHistogram = " << leftSBHistogram << std::endl;

            // lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[0]);
            // highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[1]);
            lowBin = safeLowBin(dataContainer.histograms2d[iHisto], sidebandRanges.second[0]);
            highBin = safeHighBin(dataContainer.histograms2d[iHisto], sidebandRanges.second[1]);
            h_sideBand->Add(dataContainer.histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_%zu",iHisto), lowBin, highBin)); // sum the right sideband
            double rightSBHistogram = dataContainer.histograms1d[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram

            double leftSB = dataContainer.fittings.fitTotal[iHisto]->Integral(sidebandRanges.first[0],sidebandRanges.first[1]); // in the opposite order because the histogram is decreasing
            double rightSB = dataContainer.fittings.fitTotal[iHisto]->Integral(sidebandRanges.second[0],sidebandRanges.second[1]);
            double fitSigIntegral = dataContainer.fittings.fitTotal[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

            // lowBin = dataContainer.histograms1d[iHisto]->FindBin(m_0 - signalSigmas * sigma);
            // highBin = dataContainer.histograms1d[iHisto]->FindBin(m_0 + signalSigmas * sigma);
            lowBin = safeLowBin(dataContainer.histograms1d[iHisto], m_0 - signalSigmas * sigma);
            highBin = safeHighBin(dataContainer.histograms1d[iHisto], m_0 + signalSigmas * sigma);
            double signalHistogram = dataContainer.histograms1d[iHisto]->Integral(lowBin, highBin);
        }

        // Scale the sideband histogram by ratio of it in the signal region Bs/(B1+B2)
        // std::cout << "Sideband integral (before beta scaling) = " << h_sideBand->Integral() << std::endl;
        // std::cout << "beta = " << scallingFactors[1] << std::endl;
        h_sideBand->Scale(scallingFactors[1]);
        // std::cout << "Sideband integral (after beta scaling) = " << h_sideBand->Integral() << std::endl;
        dataContainer.subtractionResults.sidebandHist.push_back(h_sideBand);

        // Subtract background histogram from signal histogram
        h_back_subtracted = (TH1D*)h_signal->Clone(Form("h_back_subtracted_%zu",iHisto));
        h_back_subtracted->Add(h_sideBand,-1.0);
        
        // There is no alpha scaling parameter when there is no reflections
        if (modelToUse != FitModelType::StandardSideBand) {
            // Scale by reflections correction
            h_back_subtracted->Scale(scallingFactors[0]);
        }
        
        

        // Account for two sigma only area used for signal region
        double coverage = TMath::Erf(signalSigmas / sqrt(2)); // coverage of a Gaussian in a +/- signalSigmas window
        h_back_subtracted->Scale(1 / coverage);
        // If histogram isn't empty, store it in the container
        bool isHistoEmpty = (h_back_subtracted->GetEntries() != 0) ? true : false;
        dataContainer.subtractionResults.subtractedHist.push_back(h_back_subtracted);

        // Calculate areas for the invariant mass distribution
        // lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.first[0]);
        // highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.first[1]);
        lowBin = safeLowBin(dataContainer.histograms2d[iHisto], sidebandRanges.first[0]);
        highBin = safeHighBin(dataContainer.histograms2d[iHisto], sidebandRanges.first[1]);
        double leftBandHistoArea = dataContainer.histograms1d[iHisto]->Integral(lowBin,highBin, "width");
        double leftBandHistoArea_no_width = dataContainer.histograms1d[iHisto]->Integral(lowBin,highBin);

        TH1D* h_sideBand_left_test = dataContainer.histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_left_test_%zu",iHisto), lowBin, highBin);
        std::cout << "Left areas: \t DeltaR_area = " << h_sideBand_left_test->Integral() << ";\t m_inv_area = " << leftBandHistoArea_no_width << std::endl;

        // lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[0]);
        // highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[1]);// right sideband limit
        lowBin = safeLowBin(dataContainer.histograms2d[iHisto], sidebandRanges.second[0]);
        highBin = safeHighBin(dataContainer.histograms2d[iHisto], sidebandRanges.second[1]);
        // if (highBin > dataContainer.histograms2d[iHisto]->GetXaxis()->GetNbins()) {
        //     highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->GetNbins();
        // }
        
        double rightBandHistoArea = dataContainer.histograms1d[iHisto]->Integral(lowBin,highBin, "width");
        // lowBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(m_0 - signalSigmas * sigma);
        // highBin = dataContainer.histograms2d[iHisto]->GetXaxis()->FindBin(m_0 + signalSigmas * sigma);
        lowBin = safeLowBin(dataContainer.histograms2d[iHisto], m_0 - signalSigmas * sigma);
        highBin = safeHighBin(dataContainer.histograms2d[iHisto], m_0 + signalSigmas * sigma);
        double signalregionHistoArea = dataContainer.histograms1d[iHisto]->Integral(lowBin,highBin, "width");
        std::cout << "leftBandHistoArea=" << leftBandHistoArea << "\t rightBandHistoArea=" << rightBandHistoArea << "\t left+right=" << leftBandHistoArea+rightBandHistoArea << "\t signalregionHistoArea=" << signalregionHistoArea << std::endl;
        std::cout << "Ys / B (calculated from histograms) = " << signalregionHistoArea / (leftBandHistoArea+rightBandHistoArea) << std::endl;
        std::cout << "Y(x=2.28646) evaluated with the fit function = " << dataContainer.fittings.fitTotal[iHisto]->Eval(2.28646) << std::endl;
        std::cout << "2nd attempt of area calculation:" << std::endl;
        double totalFitIntegralLeft = dataContainer.fittings.fitTotal[iHisto]->Integral(sidebandRanges.first[0],sidebandRanges.first[1]);
        double totalFitIntegralRight = dataContainer.fittings.fitTotal[iHisto]->Integral(sidebandRanges.second[0],sidebandRanges.second[1]);
        double totalFitIntegralSignal = dataContainer.fittings.fitTotal[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);
        // now individually each component of the fit
        double signalOnlyIntegral = dataContainer.fittings.fitSignalOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);
        double backgroundOnlyIntegral = dataContainer.fittings.fitBackgroundOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);
        double reflectionsOnlyIntegral = dataContainer.fittings.fitReflectionsOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);
        std::cout << "totalFitIntegralLeft = " << totalFitIntegralLeft << "\t totalFitIntegralRight = " << totalFitIntegralRight << "\t totalFitIntegralSignal = " << totalFitIntegralSignal << std::endl;
        std::cout << "signalOnlyIntegral = " << signalOnlyIntegral << "\t backgroundOnlyIntegral = " << backgroundOnlyIntegral << "\t reflectionsOnlyIntegral = " << reflectionsOnlyIntegral << std::endl;
        double backgroundOnlyLeft = dataContainer.fittings.fitBackgroundOnly[iHisto]->Integral(sidebandRanges.first[0],sidebandRanges.first[1]);
        double backgroundOnlyRight = dataContainer.fittings.fitBackgroundOnly[iHisto]->Integral(sidebandRanges.second[0],sidebandRanges.second[1]);
        std::cout << "backgroundOnlyLeft = " << backgroundOnlyLeft << "\t backgroundOnlyRight = " << backgroundOnlyRight << std::endl;
        double reflectionsOnlyLeft = dataContainer.fittings.fitReflectionsOnly[iHisto]->Integral(sidebandRanges.first[0],sidebandRanges.first[1]);
        double reflectionsOnlyRight = dataContainer.fittings.fitReflectionsOnly[iHisto]->Integral(sidebandRanges.second[0],sidebandRanges.second[1]);
        std::cout << "reflectionsOnlyLeft = " << reflectionsOnlyLeft << "\t reflectionsOnlyRight = " << reflectionsOnlyRight << std::endl;
        double signalOnlyMiddle = dataContainer.fittings.fitSignalOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);
        double backgroundOnlyMiddle = dataContainer.fittings.fitBackgroundOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);
        double reflectionsOnlyMiddle = dataContainer.fittings.fitReflectionsOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);
        std::cout << "signalOnlyMiddle = " << signalOnlyMiddle << "\t backgroundOnlyMiddle = " << backgroundOnlyMiddle << "\t reflectionsOnlyMiddle = " << reflectionsOnlyMiddle << std::endl;

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
        
        // Copy for later comparison
        dataContainer.subtractionResults.subtractedHistCopy.push_back((TH1D*)dataContainer.subtractionResults.subtractedHist[iHisto]->Clone(Form("hSubtractedFullJetPtCopy_%zu", iHisto)));

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
    // Store fit results for this jet pT interval
    jetPtContainer.jetPtTotalSidebandData.push_back(dataContainer);
    return dataContainer.subtractionResults;
}

void PlotHistograms(const SidebandData& dataContainer, const FitModelType& modelToUse, const BinningStruct& binning, double jetptMin, double jetptMax, const std::vector<double>& ptHFBinEdges) {
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
    TCanvas* c1d_fit = new TCanvas("c1d_fit", "1D histograms with Fit", 1800,1000);
    c1d_fit->Divide(nCols,nRows); // columns, lines
    TCanvas* c_2d = new TCanvas("c_2d", "2D histograms", 1800,1000);
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
        // dataContainer.histograms1d[iHisto]->GetYaxis()->SetTitle("counts");
        dataContainer.histograms1d[iHisto]->SetMinimum(0);
        
        dataContainer.histograms1d[iHisto]->Draw();
        dataContainer.subtractionResults.sidebandHistDrawing[iHisto]->Draw("hist same");
        dataContainer.subtractionResults.signalHistDrawing[iHisto]->Draw("hist same");
        dataContainer.fittings.fitTotal[iHisto]->SetLineColor(kBlack);
        dataContainer.fittings.fitTotal[iHisto]->SetLineStyle(kSolid);
        dataContainer.fittings.fitTotal[iHisto]->SetLineWidth(1);
        dataContainer.fittings.fitTotal[iHisto]->Draw("same");
        dataContainer.fittings.fitSignalOnly[iHisto]->Draw("same");
        dataContainer.fittings.fitBackgroundOnly[iHisto]->Draw("same");
        if (modelToUse != FitModelType::StandardSideBand) {
            dataContainer.fittings.fitReflectionsOnly[iHisto]->Draw("same");
        }

        // double A1Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(2); // Get the value of parameter 'A' from primary signal gaussian
        // double m0_1Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(4); // Get the value of parameter 'm0' from primary signal gaussian
        double sigma1Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(5); // Get the value of parameter 'sigma' from primary signal gaussian
        // double A2Signal = A1Signal / dataContainer.fittings.fitTotal[iHisto]->GetParameter(3); // Get the value of parameter 'A' from secondary signal gaussian
        double sigma2Signal = sigma1Signal / dataContainer.fittings.fitTotal[iHisto]->GetParameter(6); // Get the value of parameter 'sigma' from secondary signal gaussian
        double chi2 = dataContainer.fittings.fitTotal[iHisto]->GetChisquare();
        double degOfFreedom = dataContainer.fittings.fitTotal[iHisto]->GetNDF();
        
        // latex->DrawLatex(statBoxPos-0.52, 0.80, Form("A_{1}^{signal} = %.2f, #bar{m_{1}} = %.2f, #sigma_{1} = %.2f GeV/#it{c}^{2}", A1Signal, m0_1Signal,sigma1Signal));
        // latex->DrawLatex(statBoxPos-0.52, 0.75, Form("A_{2}^{signal} = %.2f, #bar{m_{2}} = %.2f, #sigma_{2} = %.2f GeV/#it{c}^{2}", A2Signal, m0_1Signal,sigma2Signal));
        latex->DrawLatex(statBoxPos-0.87, 0.23, Form("Data: %s", binning.inputDATA.first.Data()));
        latex->DrawLatex(statBoxPos-0.87, 0.28, Form("MC: %s",binning.inputMC.first.Data()));
        latex->DrawLatex(statBoxPos-0.87, 0.34, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax)); // Display jet pT cut applied
        latex->DrawLatex(statBoxPos-0.30, 0.80, Form("#sigma_{primary}^{signal} = %.3f", sigma1Signal));
        latex->DrawLatex(statBoxPos-0.32, 0.75, Form("#sigma_{secondary}^{signal} = %.3f", sigma2Signal));
        latex->DrawLatex(statBoxPos-0.28, 0.85, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));
        
        

        // Drawing 2D histograms
        c_2d->cd(iHisto+1);
        dataContainer.histograms2d[iHisto]->Draw("colz");

    }

    // Plotting extracted raw yields DeltaR observable
    TCanvas* cSideBand = new TCanvas("cSideBand", "delta R for side-band", 1800,1000);
    cSideBand->Divide(nCols,nRows); // columns, lines
    TCanvas* cSignal = new TCanvas("cSignal", "delta R for signal", 1800,1000);
    cSignal->Divide(nCols,nRows); // columns, lines
    TCanvas* cSubtracted = new TCanvas("cSubtracted", "delta R for side-band subtracted signal", 1800,1000);
    cSubtracted->Divide(nCols,nRows); // columns, lines
    TCanvas* cSigPlusBack = new TCanvas("cSigPlusBack", "delta R for side-band and signal in the same plot", 1800,1000);
    cSigPlusBack->Divide(nCols,nRows); // columns, lines

    TLegend* legend = new TLegend(0.6,0.57,0.9,0.77);
    legend->AddEntry(dataContainer.subtractionResults.sidebandHist[0],"Scaled sideband", "lpe");
    legend->AddEntry(dataContainer.subtractionResults.signalHist[0],"Raw yield", "lpe");
    legend->AddEntry(dataContainer.subtractionResults.subtractedHist[0],"Real signal (minus background)", "lpe");
    std::cout << "Starting subtracted histograms plotting..." << std::endl;
    for (size_t iHisto = 0; iHisto < dataContainer.subtractionResults.subtractedHist.size(); iHisto++) {
        
        cSideBand->cd(iHisto+1);
        //dataContainer.subtractionResults.sidebandHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        dataContainer.subtractionResults.sidebandHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        // dataContainer.subtractionResults.sidebandHist[iHisto]->GetYaxis()->SetTitle("yields");
        dataContainer.subtractionResults.sidebandHist[iHisto]->SetMarkerStyle(kFullTriangleUp);
        dataContainer.subtractionResults.sidebandHist[iHisto]->SetMarkerColor(kAzure);
        dataContainer.subtractionResults.sidebandHist[iHisto]->SetLineColor(kAzure);
        double statBoxPos = gPad->GetUxmax();
        dataContainer.subtractionResults.sidebandHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));

        cSignal->cd(iHisto+1);
        //dataContainer.subtractionResults.signalHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        dataContainer.subtractionResults.signalHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        // dataContainer.subtractionResults.signalHist[iHisto]->GetYaxis()->SetTitle("yields");
        dataContainer.subtractionResults.signalHist[iHisto]->SetMarkerStyle(kFullSquare);
        dataContainer.subtractionResults.signalHist[iHisto]->SetMarkerColor(kViolet);
        dataContainer.subtractionResults.signalHist[iHisto]->SetLineColor(kViolet);
        dataContainer.subtractionResults.signalHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));
        std::cout << "Plotting subtracted histogram number " << iHisto << std::endl;
        cSubtracted->cd(iHisto+1);
        //dataContainer.subtractionResults.subtractedHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        dataContainer.subtractionResults.subtractedHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        // dataContainer.subtractionResults.subtractedHist[iHisto]->GetYaxis()->SetTitle("yields");
        dataContainer.subtractionResults.subtractedHist[iHisto]->SetMarkerStyle(kFullCircle);
        dataContainer.subtractionResults.subtractedHist[iHisto]->SetMarkerColor(kPink);
        dataContainer.subtractionResults.subtractedHist[iHisto]->SetLineColor(kPink);
        dataContainer.subtractionResults.subtractedHist[iHisto]->Draw();
        dataContainer.subtractionResults.subtractedHistCopy[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        // dataContainer.subtractionResults.subtractedHistCopy[iHisto]->GetYaxis()->SetTitle("yields");
        dataContainer.subtractionResults.subtractedHistCopy[iHisto]->SetMarkerStyle(kFullCircle);
        dataContainer.subtractionResults.subtractedHistCopy[iHisto]->SetMarkerColor(kPink);
        dataContainer.subtractionResults.subtractedHistCopy[iHisto]->SetLineColor(kPink);
        dataContainer.subtractionResults.subtractedHistCopy[iHisto]->Draw("same");
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));

        cSigPlusBack->cd(iHisto+1);
        //dataContainer.subtractionResults.signalHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        double maxY = std::max(dataContainer.subtractionResults.signalHist[iHisto]->GetMaximum(), dataContainer.subtractionResults.sidebandHist[iHisto]->GetMaximum());
        maxY = std::max(maxY, dataContainer.subtractionResults.subtractedHistCopy[iHisto]->GetMaximum());
        double minY = std::min(dataContainer.subtractionResults.signalHist[iHisto]->GetMinimum(), dataContainer.subtractionResults.sidebandHist[iHisto]->GetMinimum());
        minY = std::min(minY, dataContainer.subtractionResults.subtractedHistCopy[iHisto]->GetMinimum());
        dataContainer.subtractionResults.signalHist[iHisto]->GetYaxis()->SetRangeUser(minY*0.9, maxY*1.1);
        dataContainer.subtractionResults.signalHist[iHisto]->Draw();
        dataContainer.subtractionResults.sidebandHist[iHisto]->Draw("same");
        dataContainer.subtractionResults.subtractedHist[iHisto]->Draw("same"); // negative bins set to 0 version
        //dataContainer.subtractionResults.subtractedHistCopy[iHisto]->Draw("same"); // not reset version
        legend->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.5, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));
    }
    
    // Plot scaling factors
    TCanvas* cscallingFactorsArrays = new TCanvas("cscallingFactorsArrays", "Scaling Factors Arrays", 1800, 1000);
    cscallingFactorsArrays->cd();
    // Fill two histograms with the scaling factors
    int nBins = binning.ptHFBinEdges_detector.size() - 1;
    TH1D* hAlpha = new TH1D("hAlpha", "Scaling factors;p_{T,D^{0}} (GeV/c); s.f.", binning.ptHFBinEdges_detector.size()-1, binning.ptHFBinEdges_detector.data());
    TH1D* hBeta = new TH1D("hBeta", "Scaling factor;p_{T,D^{0}} (GeV/c); s.f.", binning.ptHFBinEdges_detector.size()-1, binning.ptHFBinEdges_detector.data());
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (modelToUse != FitModelType::StandardSideBand) {
            hAlpha->SetBinContent(iBin + 1, dataContainer.subtractionResults.scalingFactorsArrays[iBin][0]); // alpha
        }
        hBeta->SetBinContent(iBin + 1, dataContainer.subtractionResults.scalingFactorsArrays[iBin][1]); // beta
        std::cout << "Beta value for " << iBin+1 << "th bin: " << dataContainer.subtractionResults.scalingFactorsArrays[iBin][1] << std::endl;
    }
    hBeta->SetMarkerStyle(kFullCircle);
    hBeta->SetMarkerColor(kBlue+1);
    hBeta->SetLineColor(kBlue+1);
    hBeta->SetLineWidth(2);
    hAlpha->SetMarkerStyle(kFullCircle);
    hAlpha->SetMarkerColor(kRed+1);
    hAlpha->SetLineColor(kRed+1);
    hAlpha->SetLineWidth(2);
    if (modelToUse != FitModelType::StandardSideBand) {
        if (hAlpha->GetBinContent(hAlpha->GetMaximumBin()) > hBeta->GetBinContent(hBeta->GetMaximumBin())) {
            hAlpha->SetMinimum(0);
            hAlpha->Draw();
            hBeta->Draw("same");
        } else {
            hBeta->SetMinimum(0);
            hBeta->Draw();
            hAlpha->Draw("same");
        }
    } else {
        hBeta->SetMinimum(0);
        hBeta->Draw();
    }
    TLegend* lScalingFactors = new TLegend(0.1,0.3,0.18,0.4);
    if (modelToUse != FitModelType::StandardSideBand) {
        lScalingFactors->AddEntry(hAlpha,"#alpha","l");
    }
    lScalingFactors->AddEntry(hBeta,"#beta","l");
    lScalingFactors->Draw();



    //
    // Storing images
    //
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../../Images/1-SignalTreatment/SideBand/" + sEmmaBins + "/" + binning.dataPeriod + "/";   
    std::cout << "Image path is: " << imagePath << std::endl;
    //
    // Storing in a single pdf file
    //
    c1d_fit->Print(imagePath + Form("sb_subtraction_deltaR_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    c_2d->Print(imagePath + Form("sb_subtraction_deltaR_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cSigPlusBack->Print(imagePath + Form("sb_subtraction_deltaR_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cscallingFactorsArrays->Print(imagePath + Form("sb_subtraction_deltaR_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
}

void SaveData(SidebandData& dataContainer, const BinningStruct& binning, double jetptMin, double jetptMax) {
    // Open output file
    TFile* fOutput = new TFile(Form("backSub_%d_to_%d_jetpt_" + binning.dataPeriod + ".root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
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

void widthsPerJetPtBin(JetPtContainerSideband& jetPtContainer, const BinningStruct& binning) {
    
    // Create a canvas for plotting
    int nHistos = binning.ptHFBinEdges_detector.size() - 1;
    // Start with a square layout (or close to it)
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));
    TCanvas* cWidths = new TCanvas("cWidths", "Widths of signal Gaussian for each jet pT bin", 1800, 1000);
    cWidths->Divide(nCols,nRows); // columns, lines

    // Create histograms to store the widths for each jet pT bin
    std::vector<TH1D*> hWidthsPerJetPtBin(binning.ptHFBinEdges_detector.size() - 1);

    for (size_t iHFPtBin = 0; iHFPtBin < binning.ptHFBinEdges_detector.size() - 1; iHFPtBin++) {
        // Create the histogram for the current HF pT bin
        hWidthsPerJetPtBin[iHFPtBin] = new TH1D(Form("hWidthsPerJet_HFPtBin_%zu", iHFPtBin),Form("%.0f < p_{T,D^{0}} < %.0f GeV/c;p_{T,jet} (GeV/c);#sigma (GeV/c^{2})", binning.ptHFBinEdges_detector[iHFPtBin], binning.ptHFBinEdges_detector[iHFPtBin + 1]), binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data());
        
        // Fill histogram with the widths for each jet pT bin
        for (size_t iJetPtBin = 0; iJetPtBin < jetPtContainer.jetPtTotalSidebandData.size(); iJetPtBin++) {

            // Do not plot the points corresponding to failed fits
            if (!jetPtContainer.jetPtTotalSidebandData[iJetPtBin].fittings.workingFits[iHFPtBin]) {
                continue;
            }

            double jetPtBinCenter = (binning.ptjetBinEdges_detector[iJetPtBin] + binning.ptjetBinEdges_detector[iJetPtBin + 1]) / 2.0;
            double width1 = jetPtContainer.jetPtTotalSidebandData[iJetPtBin].fittings.fitTotal[iHFPtBin]->GetParameter(5);
            double width1Error = jetPtContainer.jetPtTotalSidebandData[iJetPtBin].fittings.fitTotal[iHFPtBin]->GetParError(5);
            double widthRatio = jetPtContainer.jetPtTotalSidebandData[iJetPtBin].fittings.fitTotal[iHFPtBin]->GetParameter(6);
            if (widthRatio == 0) {
                continue; // skip this point to avoid division by zero
            }
            
            double width2 = width1 / widthRatio;
            double widthRatioError = jetPtContainer.jetPtTotalSidebandData[iJetPtBin].fittings.fitTotal[iHFPtBin]->GetParError(6);
            double width2Error = width2 * sqrt(pow(width1Error/width1, 2) + pow(widthRatioError/widthRatio, 2)); // propagating the errors sigma2 = sigma1 / ratio
            int bin = hWidthsPerJetPtBin[iHFPtBin]->FindBin(jetPtBinCenter);
            hWidthsPerJetPtBin[iHFPtBin]->SetBinContent(bin, width1);
            hWidthsPerJetPtBin[iHFPtBin]->SetBinError(bin, width1Error);
        }

        // Plot them on their respective pad of the canvas
        cWidths->cd(iHFPtBin+1);
        hWidthsPerJetPtBin[iHFPtBin]->SetMarkerStyle(kFullCircle);
        hWidthsPerJetPtBin[iHFPtBin]->SetMarkerColor(kYellow+2);
        hWidthsPerJetPtBin[iHFPtBin]->SetLineColor(kYellow+2);
        hWidthsPerJetPtBin[iHFPtBin]->Draw();

    }

    // Store canvas on PDF file
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../../Images/1-SignalTreatment/SideBand/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    cWidths->Print(imagePath + Form("signalGaussianWidths_perJetPtBin_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",binning.ptjetBinEdges_detector[0], binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1]));
}

void JetPtIterator(TFile* fDist, JetPtContainerSideband& jetPtContainer, const double jetptMin, const double jetptMax, const BinningStruct& binning,
                   const FitModelType& modelToUse, double& signalSigmas, int& startingBackSigma, int& backgroundSigmas) {
    std::cout << "============================================= " << jetptMin << " GeV/c < pT,jet < " << jetptMax << " GeV/c =============================================" << std::endl;

    // mass histogram
    int massBins = 200; // default=50
    double massBinDensity = 50 / (2.06 - 1.72);
    double minMass = 1.72; // default = 1.72, in order to exclude D⁰->KPiPi decay entries on the left side
    //double maxMass = 2.5; // default = 2.06, due to drop of counts on last bin (of various pT,HF ranges) changed
    std::vector<double> maxMass = {2.02, 2.045, 2.06, 2.08, 2.1, 2.14, 2.20, 2.28, 2.47, 2.47, 2.47};
    // Opening template file
    TFile* fReflectionsMC = new TFile(Form("../Reflections/reflections_%.0f_to_%.0fGeV_" + binning.dataPeriod + ".root",jetptMin,jetptMax),"read");
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
    SidebandData dataContainer = createHistograms(binning);
    
    // Fill histograms
    fillHistograms(fDist, dataContainer, jetptMin, jetptMax, binning);

    // Perform fits
    FitContainer fittings = performFit(fReflectionsMC, dataContainer, binning, modelToUse, jetptMin, jetptMax, signalSigmas, startingBackSigma, backgroundSigmas);

    // Subtract side-band from signal
    SubtractionResult subtractionResult = SideBand(dataContainer, modelToUse, binning, signalSigmas, startingBackSigma, backgroundSigmas, jetPtContainer);

    // Plot histograms
    PlotHistograms(dataContainer, modelToUse, binning, jetptMin, jetptMax, binning.ptHFBinEdges_detector);

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
        TFile* fJetRange = new TFile(Form("backSub_%0.f_to_%0.f_jetpt_" + binning.dataPeriod + ".root",binning.ptjetBinEdges_detector[iJetBin],binning.ptjetBinEdges_detector[iJetBin+1]),"read");
        if (!fJetRange || fJetRange->IsZombie()) {
            std::cerr << "Error opening file " << Form("backSub_%0.f_to_%0.f_jetpt.root",binning.ptjetBinEdges_detector[iJetBin],binning.ptjetBinEdges_detector[iJetBin+1]) << std::endl;
            continue;
        }

        for (int iHist = 0; iHist < binning.ptHFBinEdges_detector.size() - 1; ++iHist) {
            TString histName = Form("h_back_subtracted_%d", iHist);
            TH1D* hDeltaR = (TH1D*)fJetRange->Get(histName);
            if (!hDeltaR) {
                std::cout << "Warning: histogram " << histName << " not found!" << std::endl;
                break; // No more histograms
            }
            
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
                std::cerr << "Could not parse pT,HF range from histogram title: " << title << std::endl;
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

                if (std::isnan(content)) {
                    content = 0.;
                    error = 0.;
                }

                // Fill the TH3D using centers
                // h3D->Fill(ptJetCenter, deltaRcenter, ptHFCenter, content);
                h3D->SetBinContent(h3D->GetXaxis()->FindBin(ptJetCenter), h3D->GetYaxis()->FindBin(deltaRcenter), h3D->GetZaxis()->FindBin(ptHFCenter), content);
                // Optionally, if you want to preserve error propagation:
                int binX = h3D->GetXaxis()->FindBin(ptJetCenter);
                int binY = h3D->GetYaxis()->FindBin(deltaRcenter);
                int binZ = h3D->GetZaxis()->FindBin(ptHFCenter);
                h3D->SetBinError(binX, binY, binZ, error); // only if needed
            }
        }
    }

    h3D->Draw("colz");

    TFile* fOutput = new TFile(Form("full_merged_ranges_back_sub_%s.root", binning.dataPeriod.Data()), "RECREATE");
    h3D->Write();

    // Storing images
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../../Images/1-SignalTreatment/SideBand/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    // Storing in a single pdf file
    h3D->Print(imagePath + Form("3d_deltaR_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));

    // Store binning information in the output file for later use
    storeBinningInFile(fOutput, binning);

    std::cout << "3D histogram created and saved to " << fOutput->GetName() << " with 3D histogram." << std::endl;
}

void BackgroundSubtraction() {
    
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Opening binning file
    double jetptMin, jetptMax;
    TFile* fBinning = new TFile(Form("../Reflections/binningInfo.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    // Load binning from reflections file
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    // Opening data file
    TFile* fDist = new TFile("../../" + binning.inputDATA.second + "/AO2D_mergedDFs.root","read");
    if (!fDist || fDist->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT data file." << std::endl;
    }
    // Opening template file
    jetptMin = binning.ptjetBinEdges_detector[0];
    jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    TFile* fReflectionsMCFullRange = new TFile(Form("../Reflections/reflections_%.0f_to_%.0fGeV_" + binning.dataPeriod + ".root",jetptMin,jetptMax),"read");
    if (!fReflectionsMCFullRange || fReflectionsMCFullRange->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT reflections file." << std::endl;
    }
    FitModelType modelToUse = FitModelType::FullPowerLaw; // options: StandardSideBand, FullPowerLaw, FullPoly2, SignalReflectionsOnly, SignalOnly, ReflectionsOnly
    // signal/side-band region parameters
    double signalSigmas = 2; // default delta = 2
    int startingBackSigma = 4; // default position = 4
    int backgroundSigmas = 4; // default delta = 4
    //BDTScanResult scanResult = optimizeBDTCuts(fDist, fReflectionsMCFullRange, binning, modelToUse, signalSigmas, startingBackSigma, backgroundSigmas);
    // binning.ptjetBinEdges_detector = {7., 10.}; // for quick tests

    // Container to store results for each jet pT bin
    JetPtContainerSideband jetPtContainer;

    for (size_t iJetPt = 0; iJetPt < binning.ptjetBinEdges_detector.size() - 1; iJetPt++) {
        // Apply side-band method to each pT,jet bin
        std::cout << "Processing pT,jet bin: " << binning.ptjetBinEdges_detector[iJetPt] << " to " << binning.ptjetBinEdges_detector[iJetPt+1] << " GeV/c" << std::endl;
        jetptMin = binning.ptjetBinEdges_detector[iJetPt];
        jetptMax = binning.ptjetBinEdges_detector[iJetPt+1];
        JetPtIterator(fDist, jetPtContainer, jetptMin, jetptMax, binning, modelToUse, signalSigmas, startingBackSigma, backgroundSigmas);
    }

    // Fit widths as a function of jet pT for each pT,HF bin and store results in the container
    if (!binning.useEmmaYeatsBins) {
        widthsPerJetPtBin(jetPtContainer, binning);
    }

    // Compute the entire range too
    jetptMin = binning.ptjetBinEdges_detector[0];
    jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    JetPtIterator(fDist, jetPtContainer, jetptMin, jetptMax, binning, modelToUse, signalSigmas, startingBackSigma, backgroundSigmas);

    // Create 3D final histogram with pT,jet vs DeltaR vs pT,HF
    create3DBackgroundSubtracted(binning);

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl;

    std::time_t now = std::time(nullptr);
    std::cout << "Finished at: " << std::ctime(&now);
}

int main() {
    BackgroundSubtraction();
    return 0;
}