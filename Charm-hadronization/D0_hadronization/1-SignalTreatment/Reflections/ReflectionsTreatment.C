/**
 * D0 meson analysis
 * @file ReflectionsTreatment.C
 * @author Christian Reckziegel
 * @brief Macro for obtaining the fit parameters, integrals and reflection scaling parameter for each invariant mass histogram fit
**/

#include "../../commonUtilities.h"

using namespace std;

//-----------------------------------------------------------------------------Fits---------------------------------------------------------------------------------------------------------
// Defining fit functions using policy-based design for flexibility and modularity in fit model construction
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
        Double_t C_1 = par[5];      // total parameter index: 5
    Double_t C_2 = par[6];          // total parameter index: 6
    Double_t m0_1 = par[7];         // total parameter index: 7
    Double_t m0_2 = par[8];         // total parameter index: 8
    Double_t sigma_1 = par[9];      // total parameter index: 9
    Double_t sigma_2 = par[10];     // total parameter index: 10
    // Defining the custom function
    Double_t result = C_1 * TMath::Exp(-TMath::Power((m - m0_1) / (2 * sigma_1), 2)) + C_2 * TMath::Exp(-TMath::Power((m - m0_2) / (2 * sigma_2), 2));
    return result;
    }
};
// Parameters indices description:
// 0-4: signal parameters (C1, C2, m0, sigma1, sigma2)
// 5-10: reflection parameters (C1, C2, m0_1, m0_2, sigma1, sigma2)
// Compile time polymorphism via policy-based design for fit models
using SigRefModel = FitModel<SignalPolicy, ReflectionPolicy>;                                               // Signal + Reflection only
using SingleGausRefModel = FitModel<SignalSinglePolicy, ReflectionPolicy>;                                  // Single Gaussian Signal + Reflection
using PureSignalModel = FitModel<SignalPolicy>;                                                             // Single Gaussian Signal only
using PureReflectionsModel = FitModel<ReflectionPolicy>;                                                    // Reflection only
//-----------------------------------------------------------------------------Structs---------------------------------------------------------------------------------------------------------
// Histograms containers for each case (each container has a number of histograms corresponding to the invariant mass intervals)
struct HistogramGroup {
    std::vector<TH2D*> signals;
    std::vector<TH2D*> reflections;
    std::vector<TH2D*> contamination_d0Tokpipi0;
    std::vector<TH2D*> signals_and_reflections;

    std::vector<TH1D*> signals_1d;
    std::vector<TH1D*> reflections_1d;
    std::vector<TH1D*> contamination_d0Tokpipi0_1d;
    std::vector<TH1D*> signals_and_reflections_1d;
};
struct FitsGroup {
    std::vector<TF1*> signals;
    std::vector<TF1*> reflections;
    std::vector<TF1*> contamination_d0Tokpipi0;
    std::vector<TF1*> signals_and_reflections;

    // Individual components (primary and secondary Gaussians)
    std::vector<std::pair<TF1*,TF1*>> individualSignals;
    std::vector<std::pair<TF1*,TF1*>> individualReflections;
    std::vector<std::pair<TF1*,TF1*>> individualContamination_d0Tokpipi0;
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
HistogramGroup createHistograms(const BinningStruct& binning) {

    HistogramGroup histograms;

    // Helper to format a bin edge as integer or 1-decimal depending on value
    auto formatEdge = [](double val) -> std::string {
        return (std::fmod(val, 1.0) == 0.0) ? Form("%.0f", val) : Form("%.1f", val);
    };

    for (size_t i = 0; i < binning.ptHFBinEdges_detector.size() - 1; ++i) {
        double massMin = binning.minMass[i];
        double massMax = binning.maxMass[i];
        int nMassBins = std::max(1, static_cast<int>( std::round(binning.massBinDensity * (massMax - massMin))));

        // Build title once — works for any combination of integer/fractional edges
        TString title = Form("%s < p_{T,D^{0}} < %s GeV/c;#it{M}(K#pi) GeV/c^{2};#DeltaR", formatEdge(binning.ptHFBinEdges_detector[i]).c_str(), formatEdge(binning.ptHFBinEdges_detector[i+1]).c_str());

        // Common axis parameters
        int nDRBins = binning.deltaRBinEdges_detector.size() - 1;
        //double* dREdges = binning.deltaRBinEdges_detector.data();

        // Create all three histogram types with the same title and axes
        auto makeHisto = [&](const char* prefix) -> TH2D* {
            
            TH2D* h = new TH2D(Form("%s%zu", prefix, i+1), title, nMassBins, massMin, massMax, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());
            h->Sumw2();
            return h;
        };

        //histograms.contamination_d0Tokpipi0.push_back(makeHisto("histMass_contamination_d0Tokpipi0"));
        histograms.reflections.push_back(makeHisto("histMass_reflections"));
        histograms.signals.push_back(makeHisto("histMass_signals"));
        histograms.signals_and_reflections.push_back(makeHisto("histMass_signals_and_reflections"));
    }

    std::cout << "Histograms objects created." << std::endl;
    return histograms;
}

// Module to fill 2D histograms from TTree data
void fillHistograms(TFile* fInputMC, HistogramGroup& histograms, const double& jetptMin, const double& jetptMax, const BinningStruct& binning) {
    std::cout << "File name is " << fInputMC->GetName() << std::endl;
    for (size_t iPtHFBin = 0; iPtHFBin < binning.ptHFBinEdges_detector.size(); iPtHFBin++) {
        std::cout << binning.bdtPtCuts[iPtHFBin].first << "-" << binning.bdtPtCuts[iPtHFBin+1].first << "\t" << binning.bdtPtCuts[iPtHFBin].second << std::endl;
    }
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCDetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDyCut = 0.8; // on detector level Lc
    const double DeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on particle level delta R
    double hfptMin = binning.ptHFBinEdges_detector[0]; //ptHFBinEdges[0] - should start from 0 or from the lowest pT,Lc value?
    double hfptMax = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1];
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
    int MCPjetNConst;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    int MCDhfMatchedFrom, MCDhfSelectedAs;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;
    int MCDjetNConst;
    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetNConst); // float
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
    tree->SetBranchAddress("fJetNConst",&MCDjetNConst); // int
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

        // calculating delta R
        double MCDDeltaR = MCDaxisDistance;
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb);
        bool recoAngularRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut);
        bool recoJetPtRange = ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax));
        bool recoHfPtRange = ((MCDhfPt >= binning.ptHFBinEdges_detector[0]) && (MCDhfPt < binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1]));
        bool recoDeltaRRange = ((MCDDeltaR >= binning.deltaRBinEdges_detector[0]) && (MCDDeltaR < DeltaRcut));
        bool recoLevelRange = passBDTcut && recoAngularRange && recoJetPtRange && recoHfPtRange && recoDeltaRRange;

        // only matched detector level entries
        if (MCPjetNConst < 0) {
            continue;
        }

        // Fill each histogram with their respective pT intervals
        if (recoLevelRange) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < binning.ptHFBinEdges_detector.size() - 1 && !filled; iEdge++) {
                if ((MCDhfPt >= binning.ptHFBinEdges_detector[iEdge]) && (MCDhfPt < binning.ptHFBinEdges_detector[iEdge + 1])) {

                    if (isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs)) {
                        // use Emma's run 2 cuts in case useEmmaYeatsBins is true
                        if (!binning.useEmmaYeatsBins || passEmmaCut(MCDjetPt, MCDhfPt)) {
                            histograms.signals[iEdge]->Fill(MCDhfMass, MCDDeltaR);
                        }
                    } else if (isReflection(MCDhfMatchedFrom, MCDhfSelectedAs)) {
                        // use Emma's run 2 cuts in case useEmmaYeatsBins is true
                        if (!binning.useEmmaYeatsBins || passEmmaCut(MCDjetPt, MCDhfPt)) {
                            histograms.reflections[iEdge]->Fill(MCDhfMass, MCDDeltaR);
                        }
                    }
                    if (isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs) || isReflection(MCDhfMatchedFrom, MCDhfSelectedAs)) {
                        // use Emma's run 2 cuts in case useEmmaYeatsBins is true
                        if (!binning.useEmmaYeatsBins || passEmmaCut(MCDjetPt, MCDhfPt)) {
                            histograms.signals_and_reflections[iEdge]->Fill(MCDhfMass, MCDDeltaR);
                        }
                    }
                    
                    filled = true; // Exit the loop once the correct histogram is found
                }
            } // pT,HF chosen
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

    // Check for discontinuities in the histograms, e.g., empty bins, unexpected spikes and print warnings if found (currently only empty bins surrounded by filled bins)
    for (size_t iHist = 0; iHist < histograms.signals_1d.size(); iHist++) {
        if (checkHistoDiscontinuity(histograms.signals_1d[iHist])) {
            // std::cerr << "Warning: Discontinuity found in signal histogram " << iHist << std::endl;
        }
        if (checkHistoDiscontinuity(histograms.reflections_1d[iHist])) {
            // std::cerr << "Warning: Discontinuity found in reflection histogram " << iHist << std::endl;
        }
        if (checkHistoDiscontinuity(histograms.signals_and_reflections_1d[iHist])) {
            // std::cerr << "Warning: Discontinuity found in signals+reflections histogram " << iHist << std::endl;
        }
    }
}
// To store information per jet pT interval
struct JetPtContainerReflections {
    std::vector<HistogramGroup> jetPtTotalHistograms;
    std::vector<FitsGroup> jetPtTotalFits;
};
FitsGroup performFits(HistogramGroup& histograms, const BinningStruct& binning, const double& jetptMin, const double& jetptMax, JetPtContainerReflections& jetPtContainer) {

    FitsGroup fittings;

    if ((histograms.signals.size() != histograms.reflections.size()) || (histograms.signals.size() != histograms.signals_and_reflections.size()) || (histograms.reflections.size() != histograms.signals_and_reflections.size())) {
        std::cout << "Warning: different size of histogram vectors (signals, reflections, signal_and_reflections)" << std::endl;
    }
    
    // --- Signal fits
    for (size_t iInterval = 0; iInterval < binning.ptHFBinEdges_detector.size() - 1; iInterval++) {
        double minMass = binning.minMass[iInterval];
        double maxMass = binning.maxMass[iInterval];

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
        fittings.signals[iInterval]->SetParLimits(2, 1.72, 2.1);
        fittings.signals[iInterval]->SetParLimits(3, 0.006, TMath::Infinity());     // sigma1
        fittings.signals[iInterval]->SetParLimits(4, 0.006, TMath::Infinity());     // sigma2

        histograms.signals_1d[iInterval]->Fit(fittings.signals[iInterval], "RQN"); // "Q" option performs quiet fit without drawing the fit function
        //fittings.signals.push_back(fittings.signals[iInterval]);
        //fits.signals[iInterval]->Print("V");

        // Also create and store individual components
        TF1* fSignalPrimary = new TF1(Form("signalPrimaryFit_%zu", iInterval), singleGaussianFunction, minMass, maxMass, 3);
        // fSignalPrimary->SetParameters(fittings.signals[iInterval]->GetParameter(0), fittings.signals[iInterval]->GetParameter(2), fittings.signals[iInterval]->GetParameter(3)); // Amplitude, mean, sigma
        fSignalPrimary->FixParameter(0, fittings.signals[iInterval]->GetParameter(0)); // Fix amplitude to the value from the combined fit
        fSignalPrimary->FixParameter(1, fittings.signals[iInterval]->GetParameter(2)); // Fix mean to the value from the combined fit
        fSignalPrimary->FixParameter(2, fittings.signals[iInterval]->GetParameter(3)); // Fix sigma to the value from the combined fit
        fSignalPrimary->SetParName(0, "C1_primary");
        fSignalPrimary->SetParName(1, "m0_primary");
        fSignalPrimary->SetParName(2, "sigma_primary");
        TF1* fSignalSecondary = new TF1(Form("signalSecondaryFit_%zu", iInterval), singleGaussianFunction, minMass, maxMass, 3);
        // fSignalSecondary->SetParameters(fittings.signals[iInterval]->GetParameter(1), fittings.signals[iInterval]->GetParameter(2), fittings.signals[iInterval]->GetParameter(4)); // Amplitude, mean, sigma
        fSignalSecondary->FixParameter(0, fittings.signals[iInterval]->GetParameter(1)); // Fix amplitude to the value from the combined fit
        fSignalSecondary->FixParameter(1, fittings.signals[iInterval]->GetParameter(2)); // Fix mean to the value from the combined fit
        fSignalSecondary->FixParameter(2, fittings.signals[iInterval]->GetParameter(4)); // Fix sigma to the value from the combined fit
        fSignalSecondary->SetParName(0, "C2_secondary");
        fSignalSecondary->SetParName(1, "m0_secondary");
        fSignalSecondary->SetParName(2, "sigma_secondary");
        fittings.individualSignals.push_back(std::make_pair(fSignalPrimary,fSignalSecondary));
    }

    // --- Reflections fits
    for (size_t iInterval = 0; iInterval < binning.ptHFBinEdges_detector.size() - 1; iInterval++) {
        double minMass = binning.minMass[iInterval];
        double maxMass = binning.maxMass[iInterval];

        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        fittings.reflections.push_back(new TF1(Form("reflectionsFit_%zu", iInterval), fitWrapper<PureReflectionsModel>, minMass, maxMass, 11));

        // Obtain 1D mass histogram projection
        //histograms.reflections_1d.push_back( histograms.reflections[iInterval]->ProjectionX(Form("histMass_reflections_1d_%zu", iInterval+1), 1, -1) );

        // Initialize parameters from the histogram distribution characteristics
        double C1, C2, m0_1, m0_2, sigma1, sigma2;
        C1 = 0.6 * histograms.reflections_1d[iInterval]->GetMaximum();   // Dominant reflection, default to 60% of the histogram maximum
        C2 = 0.4 * histograms.reflections_1d[iInterval]->GetMaximum();   // Secondary reflection, default to 40% of the histogram maximum
        m0_1 = histograms.reflections_1d[iInterval]->GetMean() - histograms.reflections_1d[iInterval]->GetRMS() / 2; // Offset from mean
        m0_2 = histograms.reflections_1d[iInterval]->GetMean() + histograms.reflections_1d[iInterval]->GetRMS() / 2; // Offset from mean
        sigma1 = histograms.reflections_1d[iInterval]->GetRMS() / 2;  // Narrower
        sigma2 = histograms.reflections_1d[iInterval]->GetRMS();      // Wider
        fittings.reflections[iInterval]->SetParameters(0., 0., 0., 0., 0., C1, C2, m0_1, m0_2, sigma1, sigma2);
        fittings.reflections[iInterval]->SetParName(5, "C1");
        fittings.reflections[iInterval]->SetParName(6, "C2");
        fittings.reflections[iInterval]->SetParName(7, "m0_1");
        fittings.reflections[iInterval]->SetParName(8, "m0_2");
        fittings.reflections[iInterval]->SetParName(9, "sigma1");
        fittings.reflections[iInterval]->SetParName(10, "sigma2");
        fittings.reflections[iInterval]->SetParLimits(5, 0., TMath::Infinity());
        fittings.reflections[iInterval]->SetParLimits(6, 0., TMath::Infinity()); // amplitudes
        fittings.reflections[iInterval]->SetParLimits(7, 1.72, 2.1);
        fittings.reflections[iInterval]->SetParLimits(8, 1.72, 2.1);
        fittings.reflections[iInterval]->SetParLimits(9, 0.0041, TMath::Infinity());
        fittings.reflections[iInterval]->SetParLimits(10, 0.0041, TMath::Infinity()); // sigmas

        histograms.reflections_1d[iInterval]->Fit(fittings.reflections[iInterval], "RQN"); // "Q" option sets minimum printing, "R" option fit in the specified range, "N" option does not store the fit result in the histogram (to avoid overwriting the previous fit results)

        // Reflections parameters:
        // std::cout << "Reflection parameter 0 = " << fittings.reflections[iInterval]->GetParameter(0) << std::endl;
        // std::cout << "Reflection parameter 1 = " << fittings.reflections[iInterval]->GetParameter(1) << std::endl;
        // std::cout << "Reflection parameter 2 = " << fittings.reflections[iInterval]->GetParameter(2) << std::endl;
        // std::cout << "Reflection parameter 3 = " << fittings.reflections[iInterval]->GetParameter(3) << std::endl;
        // std::cout << "Reflection parameter 4 = " << fittings.reflections[iInterval]->GetParameter(4) << std::endl;
        // std::cout << "Reflection parameter 5 = " << fittings.reflections[iInterval]->GetParameter(5) << std::endl;
        // std::cout << "Reflection parameter 6 = " << fittings.reflections[iInterval]->GetParameter(6) << std::endl;
        // std::cout << "Reflection parameter 7 = " << fittings.reflections[iInterval]->GetParameter(7) << std::endl;
        // std::cout << "Reflection parameter 8 = " << fittings.reflections[iInterval]->GetParameter(8) << std::endl;
        // std::cout << "Reflection parameter 9 = " << fittings.reflections[iInterval]->GetParameter(9) << std::endl;
        // std::cout << "Reflection parameter 10 = " << fittings.reflections[iInterval]->GetParameter(10) << std::endl;
        //fittings.reflections.push_back(fittings.reflections[iInterval]);
        fittings.reflections[iInterval]->Print("V");
        std::cout << "A1/A2 = " << fittings.reflections[iInterval]->GetParameter(0) / fittings.reflections[iInterval]->GetParameter(1) << std::endl;

        // Also create and store individual components
        TF1* fReflectionPrimary = new TF1(Form("reflectionPrimaryFit_%zu", iInterval), singleGaussianFunction, minMass, maxMass, 3);
        // fReflectionPrimary->SetParameters(fittings.reflections[iInterval]->GetParameter(0), fittings.reflections[iInterval]->GetParameter(2), fittings.reflections[iInterval]->GetParameter(4)); // Amplitude, mean, sigma
        fReflectionPrimary->FixParameter(0, fittings.reflections[iInterval]->GetParameter(5)); // Fix amplitude to the value from the combined fit
        fReflectionPrimary->FixParameter(1, fittings.reflections[iInterval]->GetParameter(7)); // Fix mean to the value from the combined fit
        fReflectionPrimary->FixParameter(2, fittings.reflections[iInterval]->GetParameter(9)); // Fix sigma to the value from the combined fit
        fReflectionPrimary->SetParName(0, "C1_primary");
        fReflectionPrimary->SetParName(1, "m0_primary");
        fReflectionPrimary->SetParName(2, "sigma_primary");
        fReflectionPrimary->Print("V");
        TF1* fReflectionSecondary = new TF1(Form("reflectionSecondaryFit_%zu", iInterval), singleGaussianFunction, minMass, maxMass, 3);
        // fReflectionSecondary->SetParameters(fittings.reflections[iInterval]->GetParameter(1), fittings.reflections[iInterval]->GetParameter(3), fittings.reflections[iInterval]->GetParameter(5)); // Amplitude, mean, sigma
        fReflectionSecondary->FixParameter(0, fittings.reflections[iInterval]->GetParameter(6)); // Fix amplitude to the value from the combined fit
        fReflectionSecondary->FixParameter(1, fittings.reflections[iInterval]->GetParameter(8)); // Fix mean to the value from the combined fit
        fReflectionSecondary->FixParameter(2, fittings.reflections[iInterval]->GetParameter(10)); // Fix sigma to the value from the combined fit
        fReflectionSecondary->SetParName(0, "C2_secondary");
        fReflectionSecondary->SetParName(1, "m0_secondary");
        fReflectionSecondary->SetParName(2, "sigma_secondary");
        fReflectionSecondary->Print("V");
        fittings.individualReflections.push_back(std::make_pair(fReflectionPrimary,fReflectionSecondary));
    }

    // --- Combined fits: true signal + reflections
    for (size_t iInterval = 0; iInterval < binning.ptHFBinEdges_detector.size() - 1; iInterval++) {
        double minMass = binning.minMass[iInterval];
        double maxMass = binning.maxMass[iInterval];

        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        fittings.signals_and_reflections.push_back(new TF1(Form("combinedFit_%zu", iInterval), fitWrapper<SigRefModel>, minMass, maxMass, 11));

        // Obtain 1D mass histogram projection
        histograms.signals_and_reflections_1d.push_back( histograms.signals_and_reflections[iInterval]->ProjectionX(Form("histMass_signals_and_reflections_1d_%zu", iInterval+1), 1, -1) );

        // Initialize parameters from previous fits
        fittings.signals_and_reflections[iInterval]->SetParameters(fittings.signals[iInterval]->GetParameter(0), fittings.signals[iInterval]->GetParameter(1), fittings.signals[iInterval]->GetParameter(2), fittings.signals[iInterval]->GetParameter(3), fittings.signals[iInterval]->GetParameter(4), 
                                      fittings.reflections[iInterval]->GetParameter(5), fittings.reflections[iInterval]->GetParameter(6), fittings.reflections[iInterval]->GetParameter(7), fittings.reflections[iInterval]->GetParameter(8), fittings.reflections[iInterval]->GetParameter(9), fittings.reflections[iInterval]->GetParameter(10));
        fittings.signals_and_reflections[iInterval]->FixParameter(0,fittings.signals[iInterval]->GetParameter(0));
        fittings.signals_and_reflections[iInterval]->FixParameter(1,fittings.signals[iInterval]->GetParameter(1));
        fittings.signals_and_reflections[iInterval]->FixParameter(2,fittings.signals[iInterval]->GetParameter(2));
        fittings.signals_and_reflections[iInterval]->FixParameter(3,fittings.signals[iInterval]->GetParameter(3));
        fittings.signals_and_reflections[iInterval]->FixParameter(4,fittings.signals[iInterval]->GetParameter(4));
        fittings.signals_and_reflections[iInterval]->FixParameter(5,fittings.reflections[iInterval]->GetParameter(5));
        fittings.signals_and_reflections[iInterval]->FixParameter(6,fittings.reflections[iInterval]->GetParameter(6));
        fittings.signals_and_reflections[iInterval]->FixParameter(7,fittings.reflections[iInterval]->GetParameter(7));
        fittings.signals_and_reflections[iInterval]->FixParameter(8,fittings.reflections[iInterval]->GetParameter(8));
        fittings.signals_and_reflections[iInterval]->FixParameter(9,fittings.reflections[iInterval]->GetParameter(9));
        fittings.signals_and_reflections[iInterval]->FixParameter(10,fittings.reflections[iInterval]->GetParameter(10));
        histograms.signals_and_reflections_1d[iInterval]->Fit(fittings.signals_and_reflections[iInterval], "RQN");
        //fittings.signals_and_reflections.push_back(fittings.signals_and_reflections[iInterval]);
        //fittings.signals_and_reflections[iInterval]->Print("V");
    }
    
    std::cout << "Fits performed." << std::endl;

    // Store all fits for this jet pT interval in the jetPtContainer
    jetPtContainer.jetPtTotalFits.push_back(fittings);

    return fittings;
}

void PlotHistograms(const HistogramGroup& histograms, FitsGroup& fittings, const BinningStruct& binning, const double jetptMin, const double jetptMax) {
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
    TCanvas* cContamination_d0Tokpipi0 = new TCanvas("cContamination_d0Tokpipi0","Contamination histograms");
    //cContamination_d0Tokpipi0->SetCanvasSize(1800,1000);
    cContamination_d0Tokpipi0->Divide(nCols,nRows);

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.04);

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
        //double m_0 = fittings[iHist]->GetParameter(3); // Get the value of parameter 'm_0'
        //double sigma = fittings[iHist]->GetParameter(4); // Get the value of parameter 'sigma'
        double chi2 = fittings.signals[iHist]->GetChisquare();
        double degOfFreedom = fittings.signals[iHist]->GetNDF();
        double sigmaSig1 = fittings.signals[iHist]->GetParameter(3);
        double sigmaSig2 = fittings.signals[iHist]->GetParameter(4);
        // original position at statBoxPos-0.35, 0.70 with 0.03 of size
        latex->DrawLatex(statBoxPos-0.4, 0.65, Form("m_{0} = %.3f GeV/c^{2}", fittings.signals[iHist]->GetParameter(2))); // Display parameter 'm_0' value
        latex->DrawLatex(statBoxPos-0.87, 0.70,"MC derived data:");
        latex->DrawLatex(statBoxPos-0.87, 0.65, binning.inputMC.first); // Derived data from anchored MC
        latex->DrawLatex(statBoxPos-0.4, 0.72, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax)); // Display jet pT cut applied
        latex->DrawLatex(statBoxPos-0.4, 0.58, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));
        latex->DrawLatex(statBoxPos-0.35, 0.85, Form("#sigma_{1} = %.3f",sigmaSig1));
        latex->DrawLatex(statBoxPos-0.35, 0.80, Form("#sigma_{2} = %.3f",sigmaSig2));

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
        histograms.reflections_1d[iHist]->GetYaxis()->SetRangeUser(0., histograms.reflections_1d[iHist]->GetMaximum() * 1.3); // Set x-axis range to start from 0
        histograms.reflections_1d[iHist]->Draw();
        chi2 = fittings.reflections[iHist]->GetChisquare();
        degOfFreedom = fittings.reflections[iHist]->GetNDF();
        double sigmaRef1 = fittings.reflections[iHist]->GetParameter(9);
        double sigmaRef2 = fittings.reflections[iHist]->GetParameter(10);
        fittings.reflections[iHist]->SetLineColor(46); // pastel red
        fittings.reflections[iHist]->Draw("same");
        fittings.individualReflections[iHist].first->SetLineColor(46);
        fittings.individualReflections[iHist].first->SetLineStyle(kDashed);
        fittings.individualReflections[iHist].first->Draw("same");
        fittings.individualReflections[iHist].second->SetLineColor(46);
        fittings.individualReflections[iHist].second->SetLineStyle(kDashed);
        fittings.individualReflections[iHist].second->Draw("same");
        latex->DrawLatex(statBoxPos-0.87, 0.40,"MC derived data:");
        latex->DrawLatex(statBoxPos-0.87, 0.35, binning.inputMC.first); // Derived data from anchored MC
        latex->DrawLatex(statBoxPos-0.35, 0.80, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));
        latex->DrawLatex(statBoxPos-0.35, 0.75, Form("#sigma_{1} = %.3f",sigmaRef1));
        latex->DrawLatex(statBoxPos-0.35, 0.70, Form("#sigma_{2} = %.3f",sigmaRef2));

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
        chi2 = fittings.signals_and_reflections[iHist]->GetChisquare();
        degOfFreedom = fittings.signals_and_reflections[iHist]->GetNDF();
        fittings.signals_and_reflections[iHist]->SetLineColor(38); // pastel blue
        latex->DrawLatex(statBoxPos-0.87, 0.70,"MC derived data:");
        latex->DrawLatex(statBoxPos-0.87, 0.65, binning.inputMC.first); // Derived data from anchored MC
        latex->DrawLatex(statBoxPos-0.87, 0.60, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));

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
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../../Images/1-SignalTreatment/Reflections/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    TString jetPtRange = Form("%.0f_to_%.0fGeV", jetptMin, jetptMax);
    

    //
    // Storing in a single pdf file
    //
    cSignals->Print(Form(imagePath + "reflections_fits_" + sEmmaBins + "_%s.pdf(",jetPtRange.Data()));
    cReflections->Print(Form(imagePath + "reflections_fits_" + sEmmaBins + "_%s.pdf",jetPtRange.Data()));
    cSignalsAndReflections->Print(Form(imagePath + "reflections_fits_" + sEmmaBins + "_%s.pdf)",jetPtRange.Data()));
    
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

void JetPtIterator(JetPtContainerReflections& jetPtContainer, const double jetptMin, const double jetptMax, const BinningStruct& binning) {
    std::cout << "============================================= " << jetptMin << " GeV/c < pT,jet < " << jetptMax << " GeV/c =============================================" << std::endl;

    // D0 mass in GeV/c^2 (hard coded parameter)
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;

    // mass histogram
    int massBins = 50; // default=50 
    double minMass = 1.72; // default = 1.72
    double maxMass = 2.06; // default = 2.06

    // Opening data file
    //TFile* fInputMC = new TFile("../../SimulatedData/Hyperloop_output/McEfficiency/New_with_reflections/Merged_AO2D_HF_LHC24d3a_All.root","read"); // previous file used
    TFile* fInputMC = new TFile("../../" + binning.inputMC.second + "/AO2D_mergedDFs.root","read");
    if (!fInputMC || fInputMC->IsZombie()) {
        std::cerr << "Error: Unable to open the input ROOT file." << std::endl;
    }

    TFile* fOutput = new TFile(Form("reflections_%.0f_to_%.0fGeV_" + binning.dataPeriod + ".root", jetptMin, jetptMax),"recreate");
    if (!fOutput || fOutput->IsZombie()) {
        std::cerr << "Error: Unable to open the output ROOT file." << std::endl;
    }

    // Create multiple histograms
    HistogramGroup histograms = createHistograms(binning);

    // Fill histograms
    fillHistograms(fInputMC, histograms, jetptMin, jetptMax, binning);

    // Perform fits
    FitsGroup fittings = performFits(histograms, binning, jetptMin, jetptMax, jetPtContainer);

    // signal/side-band region parameters
    int signalSigmas = 2; // 2
    int startingBackSigma = 4; // 4
    int backgroundSigmas = 4; // 4
    
    // Plot histograms
    PlotHistograms(histograms, fittings, binning, jetptMin, jetptMax);

    // Storing final histograms to output file
    SaveData(fOutput, histograms, fittings, binning);

}

void ReflectionsTreatment() {
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // // --- Binning objects
    // BinningStruct binning;

    // // Emma Yeats reported analysis bins:
    // bool useEmmaYeatsBins = false;
    // if (useEmmaYeatsBins) {
    //     binning.useEmmaYeatsBins = useEmmaYeatsBins;
    //     // pT,jet cuts
    //     binning.ptjetBinEdges_detector = {5., 7., 10., 20., 50.};
    //     // DeltaR bins
    //     binning.deltaRBinEdges_detector = {0., 0.01, 0.03, 0.05, 0.12, 0.2};
    //     // pT,HF bins, must be the same from the BDT models
    //     // binning.ptHFBinEdges_detector = {2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36.};
    // } else {
    //     // pT,jet cuts
    //     binning.ptjetBinEdges_detector = {5., 7., 10., 16., 36., 50.}; // default = {5., 7., 10., 15., 30., 50.}
    //     // DeltaR bins
    //     binning.deltaRBinEdges_detector = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5};
    //     // pT,HF bins, must be the same from the BDT models
    //     // binning.ptHFBinEdges_detector = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.};
    // }
    // double minDeltaR = binning.deltaRBinEdges_detector[0];
    // double maxDeltaR = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1];    
    // // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.05 (from first of pair)
    // std::vector<std::pair<double, double>> bdtPtCuts = {
    //     {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.6}, {8, 0.8}, {10, 0.8}, {12, 1}, {16, 1}, {50, -1}
    // }; // on dataset JE_HF_LHC23_pass4_Thin_2P3PDstar_D0CJ_4_D0_1
    // // std::vector<std::pair<double, double>> bdtPtCuts = {
    // //     {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}, {48, -1}
    // // }; // from Nima's AN
    // binning.bdtPtCuts = bdtPtCuts;
    // binning.ptHFEfficiencyBinEdges = efficiencyBinEdges(binning.bdtPtCuts);
    // binning.ptjetBinEdges_particle = binning.ptjetBinEdges_detector;
    // binning.deltaRBinEdges_particle = binning.deltaRBinEdges_detector;
    // // the HF binning depends on the ones used in the BDT models, copy ptHFEfficiencyBinEdges minus first and last edge
    // binning.ptHFBinEdges_particle = std::vector<double>(binning.ptHFEfficiencyBinEdges.begin() + 1, binning.ptHFEfficiencyBinEdges.end() - 1);
    // binning.ptHFBinEdges_particle.push_back(36.); // add last edge to cover the entire range, as done on Nima's AN
    // //binning.ptHFBinEdges_particle = {1., 2., 3., 4., 5., 6., 7., 8., 12., 16., 24., 36.}; // override with custom binning for HF pT, same as used for BDT models
    // binning.ptHFBinEdges_detector = binning.ptHFBinEdges_particle;
    // // Define mass range here — this is the ONLY place you need to change it
    // binning.massBinDensity = 50.0 / (2.06 - 1.72); // = 147.06 bins/GeV/c²
    // // Per pT,D upper mass edge — adjust per bin if needed
    // // Default: all bins use 2.06 GeV/c²
    // binning.minMass.assign(binning.ptHFBinEdges_detector.size() - 1, 1.72);
    // binning.maxMass.assign(binning.ptHFBinEdges_detector.size() - 1, 2.06);
    // binning.maxMass = {2.02, 2.045, 2.06, 2.08, 2.1, 2.14, 2.20, 2.28, 2.47, 2.47, 2.47};
    // // Override specific bins if needed, e.g.:
    // // binning.maxMass[0] = 1.98; // first pT,D bin has less room on the right

    // // Choose file period
    // binning.dataPeriod = "2023";
    // if (binning.dataPeriod == "2023") {
    //     // DATA
    //     binning.inputDATA.first = "JE_HF_LHC23_pass4_Thin_2P3PDstar_D0CJ_4_D0_1";
    //     binning.inputDATA.second = "Data/Experimental/Train_643652";
    //     // Anchored MC
    //     binning.inputMC.first = "HF_LHC24h1c_All_D0";
    //     binning.inputMC.second = "Data/MonteCarlo/Train_645447";
    // } else if (binning.dataPeriod == "2022") {
    //     // DATA
    //     binning.inputDATA.first = "JE_HF_LHC22o_pass7_minBias_2P3PDstar_D0CJ_4_D0_1";
    //     binning.inputDATA.second = "Data/Experimental/Train_659513";
    //     // Anchored MC
    //     binning.inputMC.first = "HF_LHC24g5_All_D0";
    //     binning.inputMC.second = "Data/MonteCarlo/Train_659747";
    // }
    TString dataPeriod = "2023";
    // Open binning information file with optimal BDT score thresholds
    TFile* fBinning = new TFile("../BDTOptimization/binningInfo_" + dataPeriod + ".root", "read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the binning info ROOT file." << std::endl;
    }
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    std::cout << "ptHFBinEdges = [";
    for (size_t ipTHF = 0; ipTHF < binning.ptHFBinEdges_particle.size(); ipTHF++) {
        std::cout << binning.ptHFBinEdges_particle[ipTHF] << ",\t";
    }
    std::cout << "]" << std::endl;    
    
    // // Emma Yeats reported analysis bins:
    binning.useEmmaYeatsBins = true;
    if (binning.useEmmaYeatsBins) {
        // pT,jet cuts
        binning.ptjetBinEdges_detector = {5., 7., 10., 20., 50.};
        binning.ptjetBinEdges_particle = binning.ptjetBinEdges_detector;
        // DeltaR bins
        binning.deltaRBinEdges_detector = {0., 0.01, 0.03, 0.05, 0.12, 0.2};
        binning.deltaRBinEdges_particle = binning.deltaRBinEdges_detector;
        // pT,HF bins, must be the same from the BDT models
        // binning.ptHFBinEdges_detector = {2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36.};
    }

    // Container to store the histograms and fits for each pT,jet bin
    JetPtContainerReflections jetPtContainer;    

    for (size_t iJetPt = 0; iJetPt < binning.ptjetBinEdges_detector.size() - 1; iJetPt++) {
        
        // Apply side-band method to each pT,jet bin
        double jetptMin = binning.ptjetBinEdges_detector[iJetPt];
        double jetptMax = binning.ptjetBinEdges_detector[iJetPt+1];

        JetPtIterator(jetPtContainer, jetptMin, jetptMax, binning);
    }

    // Compute the entire range too
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    JetPtIterator(jetPtContainer, jetptMin, jetptMax, binning);

    // Also store separate binning file without mention to data period in order to automatize workflow
    TFile* fBinningNoPeriod = new TFile("binningInfo.root", "recreate");
    storeBinningInFile(fBinningNoPeriod, binning);

    fBinning->Close();

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl;

    std::time_t now = std::time(nullptr);
    std::cout << "Finished at: " << std::ctime(&now);
}

int main(){
    ReflectionsTreatment();
    return 0;
}



