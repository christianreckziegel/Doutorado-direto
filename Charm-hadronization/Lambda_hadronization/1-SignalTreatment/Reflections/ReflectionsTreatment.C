/**
 * @file ReflectionsTreatment.C
 * @author Christian Reckziegel
 * @brief Macro for obtaining the fit parameters, integrals and reflection scaling parameter for each invariant mass histogram fit
**/

#include "../../commonUtilities.h"

using namespace std;

//-----------------------------------------------------------------------------Fits---------------------------------------------------------------------------------------------------------
// Defining fit functions using policy-based design for flexibility and modularity in fit model construction
// Define the enum class with descriptive names
enum class FitModelType {
    FullPowerLaw,                   // Signal + Reflections + power law Background,
    FullPoly2,                      // Signal + Reflections + 2nd order polynominal Background,
    FullSingleGaussian,             // Signal as single Gaussian + Reflections + power law background,
    SignalReflectionsOnly,          // Signal + Reflections only
    SignalOnly,                     // Only Signal
    ReflectionsOnly                 // Only Reflections
};
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
        Double_t C_1 = par[0];
        Double_t C_2 = par[1];
        Double_t m0 = par[2]; // constrain: same mean value parameter for both single gaussians
        Double_t sigma_1 = par[3];
        Double_t sigma_2 = par[4];

        // Defining the custom function
        Double_t result = C_1 * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma_1), 2)) + C_2 * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma_2), 2));
        return result;
    }
};
struct SignalSinglePolicy { // need to implement in Reflections.C macro as well
    static double eval(double* x, double* par) {
        Double_t m = x[0];
        Double_t C_1 = par[0];                                         // Free parameter
        Double_t m0 = par[2];                                               // Free parameter
        Double_t sigma1 = par[3];                                           // Free parameter
        
        // Defining the custom function
        Double_t result = C_1 * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma1), 2));
        return result;
    }
};
struct ReflectionPolicy {
    static double eval(double* x, double* par) {
        // m0 and sigma values fixed from fits obtained from MC data fits
        Double_t m = x[0];
        Double_t C_1 = par[0];      // total parameter index: 5
    Double_t C_2 = par[5];      // total parameter index: 6
    Double_t m0_1 = par[6];     // total parameter index: 7
    Double_t m0_2 = par[7];     // total parameter index: 8
    Double_t sigma_1 = par[8];  // total parameter index: 9
    Double_t sigma_2 = par[9];  // total parameter index: 10

    // Defining the custom function
    Double_t result = C_1 * TMath::Exp(-TMath::Power((m - m0_1) / (2 * sigma_1), 2)) + C_2 * TMath::Exp(-TMath::Power((m - m0_2) / (2 * sigma_2), 2));
    return result;
    }
};
// Parameters indices description:
// 0-4: signal parameters (C1, C2, m0, sigma1, sigma2)
// 5-10: reflection parameters (C1, C2, m0_1, m0_2, sigma1, sigma2)
// template for combinations of policies
template <typename... Policies>
struct FitModel {
    static double eval(double* x, double* par) {
        return (Policies::eval(x, par) + ...); // fold expression in C++17
    }
};
// ROOT TF1 expects a regular function pointer
template <typename Model>
double fitWrapper(double* x, double* par) {
    return Model::eval(x, par);
}
// Compile time polymorphism via policy-based design for fit models
using SigRefModel = FitModel<SignalPolicy, ReflectionPolicy>;                                               // Signal + Reflection only
using SingleGausRefModel = FitModel<SignalSinglePolicy, ReflectionPolicy>;                                  // Single Gaussian Signal + Reflection
using PureSignalModel = FitModel<SignalPolicy>;                                                             // Single Gaussian Signal only
using PureReflectionsModel = FitModel<ReflectionPolicy>;                                                     // Reflection only
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

    // Individual components (primary and secondary Gaussians)
    std::vector<std::pair<TF1*,TF1*>> individualSignals;
    std::vector<std::pair<TF1*,TF1*>> individualReflections;
    std::vector<std::pair<TF1*,TF1*>> individualSignals_and_reflections;
};
struct ReflectionsData {

    // Histograms
    HistogramGroup histograms;

    // Fits
    FitsGroup fits;
};
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
Double_t pureReclectionsFunction(Double_t* x, Double_t* par) {

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
    
    return pureSignalFunction(x,par) + pureReclectionsFunction(x,&par[5]);
}
//-----------------------------------------------------------------------------Main workflow functions---------------------------------------------------------------------------------------------------------
// Module to create TH2D histograms including interest variable
HistogramGroup createHistograms(const std::vector<double>& ptHFBinEdges, int massBins, double minMass, double maxMass, const std::vector<double>& deltaRBinEdges) {

    // create histograms containers for each case (each container has a number of histograms corresponding to the invariant mass intervals)
    HistogramGroup histograms;

    for (size_t i = 0; i < ptHFBinEdges.size() - 1; ++i) {
        // So that the title adapts to fractional binning title
        if (std::fmod(ptHFBinEdges[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(ptHFBinEdges[i+1], 1.0) != 0) {
                // pure reflections histograms
                histograms.reflections.push_back(new TH2D(Form("histMass_reflections%zu",i+1), Form("%.1f < p_{T,#Lambda_{c}} < %.1f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // pure signals histograms
                histograms.signals.push_back(new TH2D(Form("histMass_signals%zu",i+1), Form("%.1f < p_{T,#Lambda_{c}} < %.1f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // signals+reflections histograms: calculate ratios
                histograms.signals_and_reflections.push_back(new TH2D(Form("histMass_signals_and_reflections%zu",i+1), Form("%.1f < p_{T,#Lambda_{c}} < %.1f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
        
            } else {
                // pure reflections histograms
                histograms.reflections.push_back(new TH2D(Form("histMass_reflections%zu",i+1), Form("%.1f < p_{T,#Lambda_{c}} < %.0f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // pure signals histograms
                histograms.signals.push_back(new TH2D(Form("histMass_signals%zu",i+1), Form("%.1f < p_{T,#Lambda_{c}} < %.0f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // signals+reflections histograms: calculate ratios
                histograms.signals_and_reflections.push_back(new TH2D(Form("histMass_signals_and_reflections%zu",i+1), Form("%.1f < p_{T,#Lambda_{c}} < %.0f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
        
            }
        } else {
            if (std::fmod(ptHFBinEdges[i+1], 1.0) != 0) {
                // pure reflections histograms
                histograms.reflections.push_back(new TH2D(Form("histMass_reflections%zu",i+1), Form("%.0f < p_{T,#Lambda_{c}} < %.1f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // pure signals histograms
                histograms.signals.push_back(new TH2D(Form("histMass_signals%zu",i+1), Form("%.0f < p_{T,#Lambda_{c}} < %.1f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // signals+reflections histograms: calculate ratios
                histograms.signals_and_reflections.push_back(new TH2D(Form("histMass_signals_and_reflections%zu",i+1), Form("%.0f < p_{T,#Lambda_{c}} < %.1f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
        
            } else {
                // pure reflections histograms
                histograms.reflections.push_back(new TH2D(Form("histMass_reflections%zu",i+1), Form("%.0f < p_{T,#Lambda_{c}} < %.0f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // pure signals histograms
                histograms.signals.push_back(new TH2D(Form("histMass_signals%zu",i+1), Form("%.0f < p_{T,#Lambda_{c}} < %.0f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
                // signals+reflections histograms: calculate ratios
                histograms.signals_and_reflections.push_back(new TH2D(Form("histMass_signals_and_reflections%zu",i+1), Form("%.0f < p_{T,#Lambda_{c}} < %.0f GeV/c;#it{M}(#piKPr) GeV/c^{2};#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()));
        
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

// Module to fill 2D histograms from TTree data
void fillHistograms(TFile* fInputMC, HistogramGroup& histograms, const double& jetptMin, const double& jetptMax, std::vector<double>& ptHFBinEdges, std::vector<double>& deltaRBinEdges, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level Lc
    const double DeltaRcut = deltaRBinEdges[deltaRBinEdges.size() - 1]; // on particle level delta R
    double hfptMin = ptHFBinEdges[0]; //ptHFBinEdges[0] - should start from 0 or from the lowest pT,Lc value?
    double hfptMax = ptHFBinEdges[ptHFBinEdges.size() - 1];
    const double hfPtMincut = hfptMin; // on particle level HF
    const double hfPtMaxcut = hfptMax; // on particle level HF

    // Accessing detector level data TTree
    TTree* tree = (TTree*)fInputMC->Get("DF_merged/O2matchtable");

    // Assuming histograms and tree data correspond in some way
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    int MCDhfMatchedFrom, MCDhfSelectedAs;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;

    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    MCPhfMass = 2.28646; // Lc rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&MCPhfY);
    tree->SetBranchAddress("fMcHfPrompt",&MCPhfprompt);
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

        // correct selection BIT shift info: Checks whether BIT(i) is set, regardless of other bits
        /*if (hfSelectedAs & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as Lc
            hfSelectedAs = 1;
        } else if (hfSelectedAs & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as Lcbar
            hfSelectedAs = -1;
        }*/

        // calculating delta R
        double deltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));
        
        // Fill each histogram with their respective pT intervals
        if ((abs(MCDhfEta) < etaCut) && (abs(MCDhfY) < yCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((deltaR >= deltaRBinEdges[0]) && (deltaR < DeltaRcut))) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < ptHFBinEdges.size() - 1 && !filled; iEdge++) {
                if ((MCDhfPt >= ptHFBinEdges[iEdge]) && (MCDhfPt < ptHFBinEdges[iEdge + 1])) {

                    // Lc = +1, Lcbar = -1, neither = 0
                    if ((MCDhfMatchedFrom != 0) && (MCDhfSelectedAs != 0)) {
                        if (MCDhfMatchedFrom == MCDhfSelectedAs) {
                            
                            // Get the threshold for this pT range
                            double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
                            // Fill histogram only if the cut is passed
                            if (MCDhfMlScore0 < maxBkgProb) {
                                // pure signals
                                histograms.signals[iEdge]->Fill(MCDhfMass, deltaR); // 2D case: (hfMass, deltaR);
                            }

                        } else {

                            // Get the threshold for this pT range
                            double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
                            // Fill histogram only if the cut is passed
                            if (MCDhfMlScore0 < maxBkgProb) {
                                // pure reflections
                                histograms.reflections[iEdge]->Fill(MCDhfMass, deltaR); // 2D case: (hfMass, deltaR);
                            }

                        }
                        
                        // Get the threshold for this pT range
                        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);
                        // Fill histogram only if the cut is passed
                        if (MCDhfMlScore0 < maxBkgProb) {
                            // signals and reflections altogether, withOUT "neither = 0" entries
                            histograms.signals_and_reflections[iEdge]->Fill(MCDhfMass, deltaR); // 2D case: (hfMass, deltaR);
                        }
                    }
                    // signals and reflections altogether, wiTH "neither = 0" entries (should I include "neither" entries too? Altough they should not appear)
                    //histograms.signals_and_reflections[iEdge]->Fill(hfMass, deltaR);
                    filled = true; // Exit the loop once the correct histogram is found
                }
            } // pT,Lc chosen
        } // end of TTree entries loop
    }
    std::cout << "Histograms filled." << std::endl;

    // Obtain 1D histograms from 2D histograms by projecting onto the mass axis (X axis)
    for (size_t iHist = 0; iHist < histograms.signals.size(); iHist++) {
        // pure reflections 1D histograms
        histograms.reflections_1d.push_back(histograms.reflections[iHist]->ProjectionX(Form("h_mass_reflections_proj_%zu", iHist)));
        // histograms.reflections_1d[iHist]->Sumw2();
        // pure signals 1D histograms
        histograms.signals_1d.push_back(histograms.signals[iHist]->ProjectionX(Form("h_mass_signals_proj_%zu", iHist)));
        // histograms.signals_1d[iHist]->Sumw2();
        // signals+reflections 1D histograms
        histograms.signals_and_reflections_1d.push_back(histograms.signals_and_reflections[iHist]->ProjectionX(Form("h_mass_signals_and_reflections_proj_%zu", iHist)));
        // histograms.signals_and_reflections_1d[iHist]->Sumw2();
    }

    std::cout << "1D Histograms projected from 2D histograms." << std::endl;
}

FitsGroup performFits(HistogramGroup& histograms, std::vector<double>& ptHFBinEdges, const double& jetptMin, const double& jetptMax, const double& minMass, const double& maxMass) {

    FitsGroup fittings;

    if ((histograms.signals.size() != histograms.reflections.size()) || (histograms.signals.size() != histograms.signals_and_reflections.size()) || (histograms.reflections.size() != histograms.signals_and_reflections.size())) {
        std::cout << "Warning: different size of histogram vectors (signals, reflections, signal_and_reflections)" << std::endl;
    }
    
    // --- Signal fits
    for (size_t iInterval = 0; iInterval < ptHFBinEdges.size() - 1; iInterval++) {
        
        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        fittings.signals.push_back(new TF1(Form("signalFit_%zu", iInterval), fitWrapper<PureSignalModel>, minMass, maxMass, 5));

        // Initialize parameters from the histogram distribution characteristics
        double C1, C2, m0, sigma1, sigma2;
        C1 = 0.7 * histograms.signals_1d[iInterval]->GetMaximum(); // Assume first Gaussian dominates
        C2 = 0.3 * histograms.signals_1d[iInterval]->GetMaximum(); // Second Gaussian smaller
        m0 = histograms.signals_1d[iInterval]->GetMean();       // Mean of the distribution
        sigma1 = histograms.signals_1d[iInterval]->GetRMS() / 2; // First Gaussian narrower
        sigma2 = histograms.signals_1d[iInterval]->GetRMS();    // Second Gaussian wider
        fittings.signals[iInterval]->SetParameters(C1, C2, m0, sigma1, sigma2);
        fittings.signals[iInterval]->SetParName(0, "C1");
        fittings.signals[iInterval]->SetParName(1, "C2");
        fittings.signals[iInterval]->SetParName(2, "m0");
        fittings.signals[iInterval]->SetParName(3, "sigma1");
        fittings.signals[iInterval]->SetParName(4, "sigma2");
        fittings.signals[iInterval]->SetParLimits(0, 0., TMath::Infinity());
        fittings.signals[iInterval]->SetParLimits(1, 0., TMath::Infinity());
        //fittings.signals[iInterval]->SetParLimits(3, 0.5 * sigma1, 2.0 * sigma1); // Example constraint for sigma1

        histograms.signals_1d[iInterval]->Fit(fittings.signals[iInterval], "RQ"); // "Q" option performs quiet fit without drawing the fit function
        fittings.signals.push_back(fittings.signals[iInterval]);
        //fits.signals[iInterval]->Print("V");

        // Also create and store individual components
        TF1* fSignalPrimary = new TF1(Form("signalPrimaryFit_%zu", iInterval), singleGaussianFunction, minMass, maxMass, 3);
        fSignalPrimary->SetParameters(fittings.signals[iInterval]->GetParameter(0), fittings.signals[iInterval]->GetParameter(2), fittings.signals[iInterval]->GetParameter(3)); // Amplitude, mean, sigma
        TF1* fSignalSecondary = new TF1(Form("signalSecondaryFit_%zu", iInterval), singleGaussianFunction, minMass, maxMass, 3);
        fSignalSecondary->SetParameters(fittings.signals[iInterval]->GetParameter(1), fittings.signals[iInterval]->GetParameter(2), fittings.signals[iInterval]->GetParameter(4)); // Amplitude, mean, sigma
        fittings.individualSignals.push_back(std::make_pair(fSignalPrimary,fSignalSecondary));
    }

    // --- Reflections fits
    for (size_t iInterval = 0; iInterval < ptHFBinEdges.size() - 1; iInterval++) {
        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        fittings.reflections.push_back(new TF1(Form("reflectionsFit_%zu", iInterval), fitWrapper<PureReflectionsModel>, minMass, maxMass, 6));

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
        fittings.reflections[iInterval]->SetParameters(C1, C2, m0_1, m0_2, sigma1, sigma2);
        fittings.reflections[iInterval]->SetParName(0, "C1");
        fittings.reflections[iInterval]->SetParName(1, "C2");
        fittings.reflections[iInterval]->SetParName(2, "m0_1");
        fittings.reflections[iInterval]->SetParName(3, "m0_2");
        fittings.reflections[iInterval]->SetParName(4, "sigma1");
        fittings.reflections[iInterval]->SetParName(5, "sigma2");
        fittings.reflections[iInterval]->SetParLimits(0, 0., TMath::Infinity());
        fittings.reflections[iInterval]->SetParLimits(1, 0., TMath::Infinity());

        histograms.reflections_1d[iInterval]->Fit(fittings.reflections[iInterval], "RQ");
        fittings.reflections.push_back(fittings.reflections[iInterval]);
        //fittings.reflections[iInterval]->Print("V");
        std::cout << "A1/A2 = " << fittings.reflections[iInterval]->GetParameter(0) / fittings.reflections[iInterval]->GetParameter(1) << std::endl;

        // Also create and store individual components
        TF1* fReflectionPrimary = new TF1(Form("reflectionPrimaryFit_%zu", iInterval), singleGaussianFunction, minMass, maxMass, 3);
        fReflectionPrimary->SetParameters(fittings.reflections[iInterval]->GetParameter(0), fittings.reflections[iInterval]->GetParameter(2), fittings.reflections[iInterval]->GetParameter(4)); // Amplitude, mean, sigma
        TF1* fReflectionSecondary = new TF1(Form("reflectionSecondaryFit_%zu", iInterval), singleGaussianFunction, minMass, maxMass, 3);
        fReflectionSecondary->SetParameters(fittings.reflections[iInterval]->GetParameter(1), fittings.reflections[iInterval]->GetParameter(3), fittings.reflections[iInterval]->GetParameter(5)); // Amplitude, mean, sigma
        fittings.individualReflections.push_back(std::make_pair(fReflectionPrimary,fReflectionSecondary));
    }

    // --- Combined fits: true signal + reflections
    for (size_t iInterval = 0; iInterval < ptHFBinEdges.size() - 1; iInterval++) {
        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        fittings.signals_and_reflections.push_back(new TF1(Form("combinedFit_%zu", iInterval), fitWrapper<SigRefModel>, minMass, maxMass, 11));

        // Obtain 1D mass histogram projection
        histograms.signals_and_reflections_1d.push_back( histograms.signals_and_reflections[iInterval]->ProjectionX(Form("histMass_signals_and_reflections_1d_%zu", iInterval+1), 1, -1) );

        // Initialize parameters from previous fits
        fittings.signals_and_reflections[iInterval]->SetParameters(fittings.signals[iInterval]->GetParameter(0), fittings.signals[iInterval]->GetParameter(1), fittings.signals[iInterval]->GetParameter(2), fittings.signals[iInterval]->GetParameter(3), fittings.signals[iInterval]->GetParameter(4), 
                                      fittings.reflections[iInterval]->GetParameter(0), fittings.reflections[iInterval]->GetParameter(1), fittings.reflections[iInterval]->GetParameter(2), fittings.reflections[iInterval]->GetParameter(3), fittings.reflections[iInterval]->GetParameter(4), fittings.reflections[iInterval]->GetParameter(5));
        fittings.signals_and_reflections[iInterval]->FixParameter(0,fittings.signals[iInterval]->GetParameter(0));
        fittings.signals_and_reflections[iInterval]->FixParameter(1,fittings.signals[iInterval]->GetParameter(1));
        fittings.signals_and_reflections[iInterval]->FixParameter(2,fittings.signals[iInterval]->GetParameter(2));
        fittings.signals_and_reflections[iInterval]->FixParameter(3,fittings.signals[iInterval]->GetParameter(3));
        fittings.signals_and_reflections[iInterval]->FixParameter(4,fittings.signals[iInterval]->GetParameter(4));
        fittings.signals_and_reflections[iInterval]->FixParameter(5,fittings.reflections[iInterval]->GetParameter(0));
        fittings.signals_and_reflections[iInterval]->FixParameter(6,fittings.reflections[iInterval]->GetParameter(1));
        fittings.signals_and_reflections[iInterval]->FixParameter(7,fittings.reflections[iInterval]->GetParameter(2));
        fittings.signals_and_reflections[iInterval]->FixParameter(8,fittings.reflections[iInterval]->GetParameter(3));
        fittings.signals_and_reflections[iInterval]->FixParameter(9,fittings.reflections[iInterval]->GetParameter(4));
        fittings.signals_and_reflections[iInterval]->FixParameter(10,fittings.reflections[iInterval]->GetParameter(5));
        histograms.signals_and_reflections_1d[iInterval]->Fit(fittings.signals_and_reflections[iInterval], "RQ");
        fittings.signals_and_reflections.push_back(fittings.signals_and_reflections[iInterval]);
        //fittings.signals_and_reflections[iInterval]->Print("V");
    }
    
    std::cout << "Fits performed." << std::endl;

    return fittings;
}

void PlotHistograms(const HistogramGroup& histograms, FitsGroup& fittings, const double jetptMin, const double jetptMax) {
    //
    // Defining canvases before plotting
    int nHistos = histograms.signals.size();
    // Start with a square layout (or close to it)
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));
    TCanvas* cSignals = new TCanvas("cSignals","Pure signal histograms");
    //cSignals->SetCanvasSize(1800,1000);
    cSignals->Divide(nCols,nRows);
    TCanvas* cReflections = new TCanvas("cReflections","Pure reflections histograms");
    //cReflections->SetCanvasSize(1800,1000);
    cReflections->Divide(nCols,nRows);
    TCanvas* cSignalsAndReflections = new TCanvas("cSignalsAndReflections","Signal and reflections histograms");
    //cSignalsAndReflections->SetCanvasSize(1800,1000);
    cSignalsAndReflections->Divide(nCols,nRows);

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.05);

    // Loop over pT,D intervals/histograms
    for (size_t iHist = 0; iHist < histograms.signals.size(); iHist++) {

        //
        // Pure signals
        //
        cSignals->cd(iHist+1);
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        histograms.signals_1d[iHist]->SetMarkerStyle(kDot); // kFullDotMedium
        histograms.signals_1d[iHist]->SetMarkerColor(kBlack); // kBlack
        histograms.signals_1d[iHist]->SetLineColor(kBlack); // kGray
        histograms.signals_1d[iHist]->GetYaxis()->SetTitle("counts");
        histograms.signals_1d[iHist]->SetMinimum(0);
        histograms.signals_1d[iHist]->Draw();
        fittings.signals[iHist]->SetLineColor(8); // pastel green
        fittings.signals[iHist]->Draw("same");
        fittings.individualSignals[iHist].first->SetLineColor(8);
        fittings.individualSignals[iHist].second->SetLineColor(8);
        fittings.individualSignals[iHist].first->SetLineStyle(kDashed);
        fittings.individualSignals[iHist].second->SetLineStyle(kDashed);
        fittings.individualSignals[iHist].first->Draw("same");
        fittings.individualSignals[iHist].second->Draw("same");
        //fittings[iHist]->Draw("same");
        //int backgroundHist = histograms.size() + iHist;
        //fittings[backgroundHist]->Draw("same");
        //double m_0 = fittings[iHist]->GetParameter(3); // Get the value of parameter 'm_0'
        //double sigma = fittings[iHist]->GetParameter(4); // Get the value of parameter 'sigma'
        // double chi2 = fittings.signals[iHist]->GetChisquare();
        // double degOfFreedom = fittings.signals[iHist]->GetNDF();
        // original position at statBoxPos-0.35, 0.70 with 0.03 of size
        // latex->DrawLatex(statBoxPos-0.4, 0.65, Form("m_{0} = %.3f GeV/c^{2}", fittings.signals[iHist]->GetParameter(2))); // Display parameter 'm_0' value
        latex->DrawLatex(statBoxPos-0.87, 0.70,"MC derived data:"); // Derived data from LHC24d3a (MC)
        latex->DrawLatex(statBoxPos-0.87, 0.65,"HF_LHC24g5_All"); // Derived data from LHC24d3a (MC)
        //latex->DrawLatex(statBoxPos-0.3, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax)); // Display jet pT cut applied
        // latex->DrawLatex(statBoxPos-0.4, 0.58, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));

        // Drawing 2D histograms
        //cSignals_2d->cd(iHist+1);
        //histograms2d[iHist]->GetYaxis()->SetTitle("#DeltaR");
        //histograms2d[iHist]->Draw("colz");

        //
        // Pure reflections
        //
        cReflections->cd(iHist+1);
        statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        histograms.reflections_1d[iHist]->SetMarkerStyle(kDot);
        histograms.reflections_1d[iHist]->SetMarkerColor(kBlack);
        histograms.reflections_1d[iHist]->SetLineColor(kBlack);
        histograms.reflections_1d[iHist]->GetYaxis()->SetTitle("counts");
        histograms.reflections_1d[iHist]->Draw();
        // chi2 = fittings.reflections[iHist]->GetChisquare();
        // degOfFreedom = fittings.reflections[iHist]->GetNDF();
        fittings.reflections[iHist]->SetLineColor(46); // pastel red
        fittings.reflections[iHist]->Draw("same");
        fittings.individualReflections[iHist].first->SetLineColor(46);
        fittings.individualReflections[iHist].first->SetLineStyle(kDashed);
        fittings.individualReflections[iHist].first->Draw("same");
        fittings.individualReflections[iHist].second->SetLineColor(46);
        fittings.individualReflections[iHist].second->SetLineStyle(kDashed);
        fittings.individualReflections[iHist].second->Draw("same");
        latex->DrawLatex(statBoxPos-0.87, 0.40,"MC derived data:"); // Derived data from LHC24d3a (MC)
        latex->DrawLatex(statBoxPos-0.87, 0.35,"HF_LHC24g5_All"); // Derived data from LHC24d3a (MC)
        //latex->DrawLatex(statBoxPos-0.87, 0.30, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));

        //
        // Signal and reflections
        //
        cSignalsAndReflections->cd(iHist+1);
        //statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        histograms.signals_and_reflections_1d[iHist]->SetMarkerStyle(kDot);
        histograms.signals_and_reflections_1d[iHist]->SetMarkerColor(kBlack);
        histograms.signals_and_reflections_1d[iHist]->SetLineColor(kBlack); // kGray
        histograms.signals_and_reflections_1d[iHist]->GetYaxis()->SetTitle("counts");
        histograms.signals_and_reflections_1d[iHist]->Draw();
        // chi2 = fittings.signals_and_reflections[iHist]->GetChisquare();
        // degOfFreedom = fittings.signals_and_reflections[iHist]->GetNDF();
        fittings.signals_and_reflections[iHist]->SetLineColor(38); // pastel blue
        latex->DrawLatex(statBoxPos-0.87, 0.70,"MC derived data:"); // Derived data from LHC24d3a (MC)
        latex->DrawLatex(statBoxPos-0.87, 0.65,"HF_LHC24g5_All"); // Derived data from LHC24d3a (MC)
        // latex->DrawLatex(statBoxPos-0.87, 0.60, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));

        // Creating legend for the three components in the final fit
        TLegend* legend = new TLegend(0.6,0.57,0.85,0.77);
        legend->AddEntry(fittings.signals[iHist],"Pure signal", "l");
        legend->AddEntry(fittings.reflections[iHist],"Reflections", "l");
        legend->AddEntry(fittings.signals_and_reflections[iHist],"Total fit", "l");

        fittings.signals[iHist]->Draw("same");
        fittings.reflections[iHist]->Draw("same");
        fittings.signals_and_reflections[iHist]->Draw("same");
        legend->Draw();

    }


    //
    // Storing images
    //
    TString imagePath = "../../Images/1-SignalTreatment/Reflections/";
    TString jetPtRange = Form("%.0f_to_%.0fGeV", jetptMin, jetptMax);

    //
    // Storing in a single pdf file
    //
    cSignals->Print(Form(imagePath + "reflections_fits_%s.pdf(",jetPtRange.Data()));
    cReflections->Print(Form(imagePath + "reflections_fits_%s.pdf",jetPtRange.Data()));
    cSignalsAndReflections->Print(Form(imagePath + "reflections_fits_%s.pdf)",jetPtRange.Data()));
    
    std::cout << "Histograms plotted." << std::endl;
}

// Module to save histograms and fits to output file
void SaveData(TFile* fOutput, const HistogramGroup& histograms, FitsGroup& fits, const BinningStruct& binning) {
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
    storeBinningInFile(fOutput, binning);
}

void JetPtIterator(const double jetptMin, const double jetptMax, const std::vector<double>& ptjetBinEdges, const bool useEmmaYeatsBins) {
    std::cout << "============================================= " << jetptMin << " GeV/c < pT,jet < " << jetptMax << " GeV/c =============================================" << std::endl;

    // Lc mass in GeV/c^2
    double m_0_parameter = 2.28646;
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
    std::vector<double> ptHFBinEdges = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.}; // 1., 2., 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 8., 10., 12., 15., 30.}
    // Emma Yeats reported analysis bins:
    if (useEmmaYeatsBins) {
        deltaRBinEdges = {0., 0.01, 0.03, 0.05, 0.12, 0.2}; // originally {0., 0.01, 0.03, 0.05, 0.12}, need to add higher bin for unfolding step
        minDeltaR = deltaRBinEdges[0];
        maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
        ptHFBinEdges = {5., 6., 7., 8., 9., 10., 12., 20.};
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
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {1, 0.05}, {2, 0.08}, {3, 0.08}, {4, 0.1}, {5, 0.2}, {6, 0.25}, {7, 0.3}, {8, 0.4}, {12, 0.5}, {24, 0.7}, {100, 0.7}
    }; // on dataset HF_LHC22o_pass7_minBias_small_2P3PDstar
    // --- Binning objects
    BinningStruct binning;
    binning.ptjetBinEdges_detector = ptjetBinEdges;
    binning.deltaRBinEdges_detector = deltaRBinEdges;
    binning.ptHFBinEdges_detector = ptHFBinEdges;
    binning.bdtPtCuts = bdtPtCuts;
    binning.ptjetBinEdges_particle = ptjetBinEdges;
    binning.deltaRBinEdges_particle = deltaRBinEdges;
    binning.ptHFBinEdges_particle = ptHFBinEdges;

    // mass histogram
    int massBins = 50; // default=50 
    double minMass = 2.1; // default = 2.1
    double maxMass = 2.49; // default = 2.49

    // Opening data file
    //TFile* fInputMC = new TFile("../../SimulatedData/Hyperloop_output/McEfficiency/New_with_reflections/Merged_AO2D_HF_LHC24d3a_All.root","read"); // previous file used
    TFile* fInputMC = new TFile("../../Data/MonteCarlo/Train_607573/AO2D_mergedDFs.root","read");
    if (!fInputMC || fInputMC->IsZombie()) {
        std::cerr << "Error: Unable to open the input ROOT file." << std::endl;
    }

    TFile* fOutput = new TFile(Form("reflections_%.0f_to_%.0fGeV.root", jetptMin, jetptMax),"recreate");
    if (!fOutput || fOutput->IsZombie()) {
        std::cerr << "Error: Unable to open the output ROOT file." << std::endl;
    }

    // Create multiple histograms
    HistogramGroup histograms = createHistograms(ptHFBinEdges,                                 // the pT,D edges will determine the number of mass histograms
                                                     massBins, minMass, maxMass,              // mass histograms binning
                                                     deltaRBinEdges);                         // deltaR histograms with asymmetrical bin widths

    // Fill histograms
    fillHistograms(fInputMC, histograms, jetptMin, jetptMax, ptHFBinEdges, deltaRBinEdges, bdtPtCuts);

    // Perform fits
    FitsGroup fittings = performFits(histograms, ptHFBinEdges, jetptMin, jetptMax, minMass, maxMass);

    // signal/side-band region parameters
    int signalSigmas = 2; // 2
    int startingBackSigma = 4; // 4
    int backgroundSigmas = 4; // 4
    
    // Plot histograms
    PlotHistograms(histograms, fittings, jetptMin, jetptMax);

    // Storing final histograms to output file
    SaveData(fOutput, histograms, fittings, binning);

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



