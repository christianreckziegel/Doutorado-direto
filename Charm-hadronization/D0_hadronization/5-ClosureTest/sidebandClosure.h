/*
 * Macro for performing side-band subtraction procedure to second closure test
 * 
 * 
 * 
 * Author: Christian Reckziegel
**/


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
struct SignalMCTemplatePolicy {
    static double eval(double* x, double* par) {
        Double_t m = x[0];
        Double_t A1Signal = par[0];                                         // Free parameter
        Double_t A2Signal = par[1];                                         // Fixed parameter
        Double_t m0 = par[2];                                               // Free parameter
        Double_t sigma1 = par[3];                                           // Free parameter
        Double_t sigma2 = par[4];                                           // Fixed parameter

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
struct ReflectionMCTemplatePolicy {
    static double eval(double* x, double* par) {
        // m0 and sigma values fixed from fits obtained from MC data fits
        Double_t m = x[0];
        Double_t A1Reflection = par[5];
        Double_t A2Reflection = par[6];
        Double_t m0_1 = par[7];
        Double_t m0_2 = par[8];
        Double_t sigma_1 = par[9];
        Double_t sigma_2 = par[10];

        // Defining the custom function
        Double_t result = A1Reflection * TMath::Exp(-TMath::Power((m - m0_1) / (2 * sigma_1), 2)) + A2Reflection * TMath::Exp(-TMath::Power((m - m0_2) / (2 * sigma_2), 2));
        return result;
    }
};

// Compile time polymorphism via policy-based design for fit models
using FullModelPowerLaw = FitModel<SignalPolicy, ReflectionPolicy, PowerLawBackgroundPolicy>;               // Signal + Reflection + power law Background
using FullModelPoly2 = FitModel<SignalPolicy, ReflectionPolicy, Poly2BackgroundPolicy>;                     // Signal + Reflection + 2nd order polynomial Background
using SigRefModel = FitModel<SignalPolicy, ReflectionPolicy>;                                               // Signal + Reflection only for "DATA"
using PureSignalModel = FitModel<SignalPolicy>;                                                             // Signal only for "DATA"
using PureReflectionsModel = FitModel<ReflectionPolicy>;                                                    // Reflection only for "DATA" template

using SigRefTemplateModel = FitModel<SignalMCTemplatePolicy, ReflectionMCTemplatePolicy>;                   // Signal + Reflection only for MC template
using PureSignalTemplateModel = FitModel<SignalMCTemplatePolicy>;                                           // Signal only for MC template
using PureReflectionsTemplateModel = FitModel<ReflectionMCTemplatePolicy>;                                  // Reflection only for MC template

using PowerLawBackgroundModel = FitModel<PowerLawBackgroundPolicy>;                                         // Power law background only
using Poly2BackgroundModel = FitModel<Poly2BackgroundPolicy>;                                               // 2nd order polynomial background only
using FullSingleGaussianModel = FitModel<SignalSinglePolicy, ReflectionPolicy, PowerLawBackgroundPolicy>;   // Single Gaussian Signal + Reflection + power law Background
using StandardSideBandSubtraction = FitModel<SignalSinglePolicy,PowerLawBackgroundPolicy>;                  // Standard side-band subtraction procedure
//_________________________________________________________________________________
// Histograms containers for each case
struct HistogramGroup {
    std::vector<TH2D*> signals;
    std::vector<TH2D*> reflections;
    std::vector<TH2D*> signals_and_reflections;

    std::vector<TH1D*> signals_1d;
    std::vector<TH1D*> reflections_1d;
    std::vector<TH1D*> signals_and_reflections_1d;
};
struct SidebandClosureResult {
    TH3D* hBackgroundSubtracted;                            // final 3D distribution
    // For each jet pT bin, the set of pT,D bin indices that survived
    std::vector<std::vector<bool>> workingFitsPerJetPt;     // [iJetPt][iPtD]
};
std::pair<FitContainer, HistogramGroup> calculateFitTemplates(TFile* fClosureInput, const double& jetptMin, const double& jetptMax, const BinningStruct& binning, const int& massBins, const double& minMass, const double& maxMass) {
    // Create template fits container
    FitContainer fTemplateFits;
    std::cout << "MC templates: creating template fits for side-band subtraction procedure..." << std::endl;
    // Create template histograms for template fits from detector level correction data sample
    HistogramGroup histogramTemplates;
    for (size_t i = 0; i < binning.ptHFBinEdges_detector.size() - 1; ++i) {
        // So that the title adapts to fractional binning title
        if (std::fmod(binning.ptHFBinEdges_detector[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(binning.ptHFBinEdges_detector[i+1], 1.0) != 0) {
                // pure reflections histograms
                histogramTemplates.reflections.push_back(new TH2D(Form("histMassTemplate_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D^{0}} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // pure signals histograms
                histogramTemplates.signals.push_back(new TH2D(Form("histMassTemplate_signals%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D^{0}} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // signals+reflections histograms: calculate ratios
                histogramTemplates.signals_and_reflections.push_back(new TH2D(Form("histMassTemplate_signals_and_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D^{0}} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
            } else {
                // pure reflections histograms
                histogramTemplates.reflections.push_back(new TH2D(Form("histMassTemplate_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D^{0}} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // pure signals histograms
                histogramTemplates.signals.push_back(new TH2D(Form("histMassTemplate_signals%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D^{0}} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // signals+reflections histograms: calculate ratios
                histogramTemplates.signals_and_reflections.push_back(new TH2D(Form("histMassTemplate_signals_and_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < p_{T,D^{0}} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
            }
        } else {
            if (std::fmod(binning.ptHFBinEdges_detector[i+1], 1.0) != 0) {
                // pure reflections histograms
                histogramTemplates.reflections.push_back(new TH2D(Form("histMassTemplate_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D^{0}} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // pure signals histograms
                histogramTemplates.signals.push_back(new TH2D(Form("histMassTemplate_signals%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D^{0}} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // signals+reflections histograms: calculate ratios
                histogramTemplates.signals_and_reflections.push_back(new TH2D(Form("histMassTemplate_signals_and_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D^{0}} < %.1f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
            } else {
                // pure reflections histograms
                histogramTemplates.reflections.push_back(new TH2D(Form("histMassTemplate_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D^{0}} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // pure signals histograms
                histogramTemplates.signals.push_back(new TH2D(Form("histMassTemplate_signals%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D^{0}} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // signals+reflections histograms: calculate ratios
                histogramTemplates.signals_and_reflections.push_back(new TH2D(Form("histMassTemplate_signals_and_reflections%zu_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < p_{T,D^{0}} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
            }
        }
        // pure reflections histograms
        histogramTemplates.reflections[i]->Sumw2();
        histogramTemplates.reflections[i]->SetLineColor(kGreen+1);
        histogramTemplates.reflections[i]->SetDirectory(0);
        // pure signals histograms
        histogramTemplates.signals[i]->Sumw2();
        histogramTemplates.signals[i]->SetLineColor(kBlue+1);
        histogramTemplates.signals[i]->SetDirectory(0);
        // signals+reflections histograms: calculate ratios
        histogramTemplates.signals_and_reflections[i]->Sumw2();
        histogramTemplates.signals_and_reflections[i]->SetLineColor(kBlack+1);
    }
    std::cout << "MC templates: histograms created." << std::endl;
    // Fill histograms with appropriate cuts
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on particle level delta R
    
    // Accessing detector level data TTree
    TTree* tree = (TTree*)fClosureInput->Get("CorrectionTree");
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
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetNConst); // float
    MCPhfMass = 1.86484; // D0 rest mass in GeV/c^2
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

        // Convert to enum class for type-safe comparison
        D0Species hfMatchedFrom = intToD0Species(MCDhfMatchedFrom);
        D0Species hfSelectedAs = intToD0Species(MCDhfSelectedAs);

        // only matched detector level entries
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isPrompt = MCDhfprompt && MCPhfprompt;                                 // both need to be prompt in order to accept the candidate
        if ((!isPrompt) || !MCDhfmatch) { // detector level
            continue;
        }
        
        // calculating delta R
        double MCDdeltaR = MCDaxisDistance;

        // Fill each histogram with their respective pT intervals
        if ((abs(MCDjetEta) < etaCut) && (abs(MCDhfY) < yCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDdeltaR >= binning.deltaRBinEdges_detector[0]) && (MCDdeltaR < binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]))) {
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < binning.ptHFBinEdges_detector.size() - 1 && !filled; iEdge++) {
                if ((MCDhfPt >= binning.ptHFBinEdges_detector[iEdge]) && (MCDhfPt < binning.ptHFBinEdges_detector[iEdge + 1])) {
                    // D0 = +1, D0bar = -1, neither = 0
                    if ((hfMatchedFrom != D0Species::NEITHER) && (hfSelectedAs != D0Species::NEITHER)) {
                        if (hfMatchedFrom == hfSelectedAs) {
                            // Get the threshold for this pT range
                            double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
                            // Fill histogram only if the cut is passed
                            if (MCDhfMlScore0 < maxBkgProb) {
                                // use Emma's run 2 cuts in case useEmmaYeatsBins is true
                                if (!binning.useEmmaYeatsBins || passEmmaCut(MCDjetPt, MCDhfPt)) {
                                    // pure signals
                                    histogramTemplates.signals[iEdge]->Fill(MCDhfMass, MCDdeltaR); // 2D case: (MCDhfMass, MCDdeltaR);
                                }
                            }
                        } else {
                            // Get the threshold for this pT range
                            double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
                            // Fill histogram only if the cut is passed
                            if (MCDhfMlScore0 < maxBkgProb) {
                                // use Emma's run 2 cuts in case useEmmaYeatsBins is true
                                if (!binning.useEmmaYeatsBins || passEmmaCut(MCDjetPt, MCDhfPt)) {
                                    // pure reflections
                                    histogramTemplates.reflections[iEdge]->Fill(MCDhfMass, MCDdeltaR); // 2D case: (MCDhfMass, MCDdeltaR);
                                }
                            }
                        }
                        // Get the threshold for this pT range
                        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
                        // Fill histogram only if the cut is passed
                        if (MCDhfMlScore0 < maxBkgProb) {
                            // use Emma's run 2 cuts in case useEmmaYeatsBins is true
                            if (!binning.useEmmaYeatsBins || passEmmaCut(MCDjetPt, MCDhfPt)) {
                                // signals and reflections altogether, withOUT "neither = 0" entries
                                histogramTemplates.signals_and_reflections[iEdge]->Fill(MCDhfMass, MCDdeltaR); // 2D case: (MCDhfMass, MCDdeltaR);
                            }
                        }
                    }
                    filled = true; // Exit the loop once the correct histogram is found
                }
            }
        }
    }
    std::cout << "MC templates: histograms filled, now performing fits..." << std::endl;
    // Fit signal and reflection models to the histograms
    // Signal fits
    for (size_t iInterval = 0; iInterval < binning.ptHFBinEdges_detector.size() - 1; iInterval++) {
        TF1* signalFit = new TF1(Form("signalFit_%zu_%.0f_to_%.0fGeV", iInterval, jetptMin, jetptMax), fitWrapper<PureSignalTemplateModel>, minMass, maxMass, 5);

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
        signalFit->SetParName(0, "C1_signal");
        signalFit->SetParName(1, "C2_signal");
        signalFit->SetParName(2, "m0_signal");
        signalFit->SetParName(3, "sigma1_signal");
        signalFit->SetParName(4, "sigma2_signal");
        signalFit->SetParLimits(0, 0., TMath::Infinity()); // Example constraint for C1_signal: 0, 0., TMath::Infinity()
        signalFit->SetParLimits(1, 0., TMath::Infinity()); // Example constraint for C2_signal: 0, 0., TMath::Infinity()
        if ((jetptMin == 30.) && (jetptMax == 50.)) {
            if (iInterval == 6) { // pT,HF \in [7;8] GeV/c
                signalFit->SetParLimits(0, 0., 10.); // Example constraint for C1_signal: 0, 0., TMath::Infinity()
                // signalFit->SetParLimits(1, 0., 10.);
            }
            
        }
        

        histogramTemplates.signals_1d[iInterval]->Fit(signalFit, "RQ0"); // "Q" option performs quiet fit without drawing the fit function
        fTemplateFits.fitSignalOnly.push_back(signalFit);
        fTemplateFits.fitSignalOnly[iInterval]->SetLineColor(8); // pastel green
        //fTemplateFits.fitSignalOnly[iInterval]->Print("V");
    }

    // Reflections fits
    for (size_t iInterval = 0; iInterval < binning.ptHFBinEdges_detector.size() - 1; iInterval++) {
        TF1* reflectionsFit = new TF1(Form("reflectionsFit_%zu_%.0f_to_%.0fGeV", iInterval, jetptMin, jetptMax), fitWrapper<PureReflectionsTemplateModel>, minMass, maxMass, 11);

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
        reflectionsFit->SetParName(5, "A1_reflections");
        reflectionsFit->SetParName(6, "A2_reflections");
        reflectionsFit->SetParName(7, "m0_1_reflections");
        reflectionsFit->SetParName(8, "m0_2_reflections");
        reflectionsFit->SetParName(9, "Sigma1_reflections");
        reflectionsFit->SetParName(10, "Sigma2_reflections");
        reflectionsFit->SetParameters(0., 0., 0., 0., 0., C1, C2, m0_1, m0_2, sigma1, sigma2);
        reflectionsFit->SetParLimits(5, 0., TMath::Infinity());
        reflectionsFit->SetParLimits(6, 0., TMath::Infinity());

        histogramTemplates.reflections_1d[iInterval]->Fit(reflectionsFit, "RQ0");
        fTemplateFits.fitReflectionsOnly.push_back(reflectionsFit);
        fTemplateFits.fitReflectionsOnly[iInterval]->SetLineColor(46); // pastel red
        fTemplateFits.fitReflectionsOnly[iInterval]->Print("V");
        std::cout << "A1/A2 = " << fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(5) / fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(6) << std::endl;
    }

    // Combined fits: true signal + reflections
    for (size_t iInterval = 0; iInterval < binning.ptHFBinEdges_detector.size() - 1; iInterval++) {
        TF1* combinedFit = new TF1(Form("combinedFit_%zu_%.0f_to_%.0fGeV", iInterval, jetptMin, jetptMax), fitWrapper<SigRefTemplateModel>, minMass, maxMass, 11);
        // Obtain 1D mass histogram projection
        histogramTemplates.signals_and_reflections_1d.push_back( histogramTemplates.signals_and_reflections[iInterval]->ProjectionX(Form("histMass_signals_and_reflections_1d_%zu_%.0f_to_%.0fGeV", iInterval+1, jetptMin, jetptMax), 1, -1) );

        // Set parameter names
        combinedFit->SetParName(0, "A1_signal");
        combinedFit->SetParName(1, "A2_signal");
        combinedFit->SetParName(2, "m0_signal");
        combinedFit->SetParName(3, "Sigma1_signal");
        combinedFit->SetParName(4, "Sigma2_signal");
        combinedFit->SetParName(5, "A1_reflections");
        combinedFit->SetParName(6, "A2_reflections");
        combinedFit->SetParName(7, "m0_1_reflections");
        combinedFit->SetParName(8, "m0_2_reflections");
        combinedFit->SetParName(9, "Sigma1_reflections");
        combinedFit->SetParName(10, "Sigma2_reflections");

        // Initialize parameters from previous fits
        combinedFit->SetParameters(fTemplateFits.fitSignalOnly[iInterval]->GetParameter(0), fTemplateFits.fitSignalOnly[iInterval]->GetParameter(1), fTemplateFits.fitSignalOnly[iInterval]->GetParameter(2), fTemplateFits.fitSignalOnly[iInterval]->GetParameter(3), fTemplateFits.fitSignalOnly[iInterval]->GetParameter(4), 
                                      fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(5), fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(6), fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(7), fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(8), fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(9), fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(10));
        combinedFit->FixParameter(0,fTemplateFits.fitSignalOnly[iInterval]->GetParameter(0));
        combinedFit->FixParameter(1,fTemplateFits.fitSignalOnly[iInterval]->GetParameter(1));
        combinedFit->FixParameter(2,fTemplateFits.fitSignalOnly[iInterval]->GetParameter(2));
        combinedFit->FixParameter(3,fTemplateFits.fitSignalOnly[iInterval]->GetParameter(3));
        combinedFit->FixParameter(4,fTemplateFits.fitSignalOnly[iInterval]->GetParameter(4));
        combinedFit->FixParameter(5,fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(5));
        combinedFit->FixParameter(6,fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(6));
        combinedFit->FixParameter(7,fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(7));
        combinedFit->FixParameter(8,fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(8));
        combinedFit->FixParameter(9,fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(9));
        combinedFit->FixParameter(10,fTemplateFits.fitReflectionsOnly[iInterval]->GetParameter(10));
        histogramTemplates.signals_and_reflections_1d[iInterval]->Fit(combinedFit, "RQ0");
        fTemplateFits.fitSignalsAndReflections.push_back(combinedFit);
        fTemplateFits.fitSignalsAndReflections[iInterval]->SetLineColor(38); // pastel blue
        //fTemplateFits.fitSignalsAndReflections[iInterval]->Print("V");
    }
    std::cout << "MC templates: fits performed." << std::endl;
    std::pair<FitContainer, HistogramGroup> fTemplateFitsAndCanvas = std::make_pair(fTemplateFits, histogramTemplates);
    return fTemplateFitsAndCanvas;
};

void PlotHistograms(HistogramGroup histogramTemplates, FitContainer fTemplateFits, const std::vector<TH2D*>& hInvariantMass2D, const std::vector<TH1D*>& hInvariantMass1D, const std::vector<TH2D*>& hInvariantMass2DOracle, const std::vector<TH1D*>& hInvariantMass1DOracle,
     FitContainer& fDataFits, const std::vector<TH1D*>& hsideBandSubtracted, TH2D* hDeltaR_vs_PtD, const double& jetptMin, const double& jetptMax, const FitModelType& modelToUse, const BinningStruct& binning) {
    //
    int nHistos = histogramTemplates.signals_and_reflections_1d.size();
    // Start with a square layout (or close to it)
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));

    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->SetTextAlign(31); // align right
    
    //
    // MC template fits and histograms
    //
    TCanvas* cTemplateFits = new TCanvas(Form("cTemplateFits_%.0f_%.0f", jetptMin, jetptMax),Form("Template histograms and fits %.0f < pT,jet < %0.f GeV/c", jetptMin, jetptMax),1920,1080);
    cTemplateFits->Divide(nCols,nRows); // columns, lines
    TCanvas* cTemplateFits2D = new TCanvas(Form("cTemplateFits2D_%.0f_%.0f", jetptMin, jetptMax),Form("2D Template histograms %.0f < pT,jet < %0.f GeV/c", jetptMin, jetptMax),1920,1080);
    cTemplateFits2D->Divide(nCols,nRows); // columns, lines
    // Loop through all histograms (and fitting functions in the future)
    for(size_t iHisto = 0; iHisto < nHistos; ++iHisto) {
        
        // Drawing 2D histograms
        cTemplateFits2D->cd(iHisto+1);
        histogramTemplates.signals_and_reflections[iHisto]->Draw("colz");

        // Drawing 1D histograms with fits
        cTemplateFits->cd(iHisto+1);
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        histogramTemplates.signals_and_reflections_1d[iHisto]->Draw();
        fTemplateFits.fitSignalOnly[iHisto]->Draw("same");
        fTemplateFits.fitReflectionsOnly[iHisto]->Draw("same");
        fTemplateFits.fitSignalsAndReflections[iHisto]->Draw("same");
        latex->DrawLatex(statBoxPos-0.1,0.42,Form("%.1f < #it{p}_{T,jet} < %.1f GeV/c", jetptMin, jetptMax));
        double chi2red = fTemplateFits.fitSignalsAndReflections[iHisto]->GetChisquare();
        int ndf = fTemplateFits.fitSignalsAndReflections[iHisto]->GetNDF();
        latex->DrawLatex(statBoxPos-0.1, 0.36,Form("#Chi^{2}_{red} = %.2f",chi2red/ndf));
        latex->DrawLatex(statBoxPos-0.1, 0.26,Form("#sigma_{primary}^{signal} = %.2f",fTemplateFits.fitSignalOnly[iHisto]->GetParameter(3)));
        latex->DrawLatex(statBoxPos-0.1, 0.16,Form("#sigma_{secondary}^{signal} = %.2f",fTemplateFits.fitSignalOnly[iHisto]->GetParameter(4)));
        // Creating legend for the three components in the final fit
        TLegend* legendMCTemplate = new TLegend(0.6,0.57,0.85,0.77);
        legendMCTemplate->AddEntry(fTemplateFits.fitSignalOnly[iHisto],"Pure signal", "l");
        legendMCTemplate->AddEntry(fTemplateFits.fitReflectionsOnly[iHisto],"Reflections", "l");
        legendMCTemplate->AddEntry(fTemplateFits.fitSignalsAndReflections[iHisto],"Total fit", "l");
        legendMCTemplate->Draw();
    }
    
    //
    // "Data" fits and histograms
    //
    TCanvas* cDataFits = new TCanvas(Form("cDataFits_%.0f_%.0f", jetptMin, jetptMax),Form("Data histograms and fits %.0f < pT,jet < %0.f GeV/c", jetptMin, jetptMax),1920,1080);
    cDataFits->Divide(nCols,nRows); // columns, lines
    TCanvas* cDataFits2D = new TCanvas(Form("cDataFits2D_%.0f_%.0f", jetptMin, jetptMax),Form("2D Data histograms %.0f < pT,jet < %0.f GeV/c", jetptMin, jetptMax),1920,1080);
    cDataFits2D->Divide(nCols,nRows); // columns, lines
    // Loop through all histograms (and fitting functions in the future)
    for(size_t iHisto = 0; iHisto < nHistos; ++iHisto) {
        
        // 2D histograms
        cDataFits2D->cd(iHisto+1);
        hInvariantMass2D[iHisto]->Draw("colz");

        // 1D histograms with fits
        cDataFits->cd(iHisto+1);
        if (hInvariantMass1D[iHisto]->GetEntries() > 0) {
            hInvariantMass1D[iHisto]->Draw();
            if (modelToUse != FitModelType::SignalReflectionsOnly)
                fDataFits.fitBackgroundOnly[iHisto]->Draw("same");
            fDataFits.fitSignalOnly[iHisto]->Draw("same");
            fDataFits.fitReflectionsOnly[iHisto]->Draw("same");
            fDataFits.fitTotal[iHisto]->Draw("same");
            double chi2red = fDataFits.fitTotal[iHisto]->GetChisquare();
            int ndf = fDataFits.fitTotal[iHisto]->GetNDF();
            double statBoxPos = gPad->GetUxmax();
            if (ndf > 0) {
                latex->DrawLatex(statBoxPos-0.2, 0.38, Form("#Chi^{2}_{red} = %.2f", chi2red/ndf));
            }
            latex->DrawLatex(statBoxPos-0.2, 0.48,Form("#sigma_{primary}^{signal} = %.2f",fDataFits.fitSignalOnly[iHisto]->GetParameter(5)));
            latex->DrawLatex(statBoxPos-0.2, 0.58,Form("#sigma_{secondary}^{signal} = %.2f",fDataFits.fitSignalOnly[iHisto]->GetParameter(5)/fDataFits.fitSignalOnly[iHisto]->GetParameter(6)));
        } else {
            // Use DrawFrame instead of Draw() — avoids degenerate [0,0] y-axis
            // gPad->DrawFrame(minMass, 0, maxMass, 1);
            TLatex* tl = new TLatex();
            tl->SetNDC(); tl->SetTextSize(0.05);
            tl->DrawLatex(0.15, 0.5, "(No data / Fit failed)");
        }
    }
    
    // Side-band subtracted histograms
    TCanvas* c1d_sideBandSubtracted = new TCanvas(Form("c1d_sideBandSubtracted_%.0f_%.0f", jetptMin, jetptMax),Form("Side-band subtracted histograms %.0f < pT,jet < %0.f GeV/c", jetptMin, jetptMax),1920,1080);
    c1d_sideBandSubtracted->Divide(nCols,nRows); // columns, lines
    for(size_t iHisto = 0; iHisto < nHistos; ++iHisto) {
        // c1d_sideBandSubtracted->cd(iHisto+1);
        
        // // MC truth histogram
        // TH1D* hsideBandSubtractedMCTruth = (TH1D*) hInvariantMass2DOracle[iHisto]->ProjectionY(Form("hsideBandSubtractedMCTruthProjY_%zu_%.0f_to_%.0fGeV", iHisto, jetptMin, jetptMax));
        // hsideBandSubtractedMCTruth->SetLineColor(kBlue+2);
        // hsideBandSubtractedMCTruth->SetLineWidth(1);
        // hsideBandSubtractedMCTruth->SetLineStyle(kDashed);

        // // Background subtracted histogram from data
        // hsideBandSubtracted[iHisto]->SetLineColor(kRed+2);
        // hsideBandSubtracted[iHisto]->SetLineWidth(1);
        // hsideBandSubtracted[iHisto]->SetTitle(Form("Side-band subtracted histogram for %.1f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c}", binning.ptHFBinEdges_detector[iHisto], binning.ptHFBinEdges_detector[iHisto + 1]));
        // hsideBandSubtracted[iHisto]->GetXaxis()->SetTitle("#Delta R");
        // hsideBandSubtracted[iHisto]->GetYaxis()->SetTitle("Counts");
        // hsideBandSubtracted[iHisto]->GetYaxis()->SetRangeUser(std::min(hsideBandSubtracted[iHisto]->GetMinimum(), hsideBandSubtractedMCTruth->GetMinimum()) * 0.8, std::max(hsideBandSubtracted[iHisto]->GetMaximum(), hsideBandSubtractedMCTruth->GetMaximum()) * 1.2); // Set y-axis range to fit both histograms
        // hsideBandSubtracted[iHisto]->Draw();
        // hsideBandSubtractedMCTruth->Draw("same");

        // =====================================================
        // Select cell in divided canvas
        // =====================================================

        c1d_sideBandSubtracted->cd(iHisto+1);

        // =====================================================
        // Create TOP and BOTTOM pads
        // =====================================================

        TPad* padTop = new TPad(
            Form("padTop_%zu", iHisto),"",0.0, 0.30,1.0, 1.0);

        TPad* padBottom = new TPad(Form("padBottom_%zu", iHisto),"",0.0, 0.0,1.0, 0.30);

        padTop->SetBottomMargin(0.02);

        padBottom->SetTopMargin(0.02);
        padBottom->SetBottomMargin(0.30);

        padTop->Draw();
        padBottom->Draw();

        // =====================================================
        // TOP PAD
        // =====================================================

        padTop->cd();

        // MC truth histogram
        TH1D* hMC = (TH1D*)hInvariantMass2DOracle[iHisto]->ProjectionY(Form("hsideBandSubtractedMCTruthProjY_%zu_%.0f_to_%.0fGeV",iHisto, jetptMin, jetptMax));

        hMC->SetLineColor(kRed+2);
        hMC->SetLineWidth(1);
        hMC->SetLineStyle(kDashed);

        // Data histogram
        TH1D* hData = hsideBandSubtracted[iHisto];

        hData->SetLineColor(kBlue+2);
        hData->SetLineWidth(1);

        hData->SetTitle(Form("%.1f < #it{p}_{T,D^{0}} < %.1f GeV/#it{c}",binning.ptHFBinEdges_detector[iHisto],binning.ptHFBinEdges_detector[iHisto + 1]));

        // Hide x-axis labels in top pad
        hData->GetXaxis()->SetLabelSize(0);
        hData->GetXaxis()->SetTitleSize(0);

        hData->GetYaxis()->SetTitle("Counts");

        double ymin = std::min(hData->GetMinimum(),hMC->GetMinimum()) * 0.8;

        double ymax = std::max(hData->GetMaximum(),hMC->GetMaximum()) * 1.2;

        hData->GetYaxis()->SetRangeUser(ymin, ymax);

        hData->Draw("E");
        hMC->Draw("HIST SAME");
        TLegend* legsubtracted = new TLegend(0.6,0.57,0.85,0.77);
        legsubtracted->AddEntry(hData,"Subtracted", "l");
        legsubtracted->AddEntry(hMC,"Truth", "l");
        legsubtracted->Draw();

        // =====================================================
        // BOTTOM PAD (RATIO)
        // =====================================================

        padBottom->cd();

        TH1D* hRatio = (TH1D*) hData->Clone(Form("hRatio_%zu", iHisto));

        hRatio->Divide(hMC);

        hRatio->SetTitle("");

        hRatio->GetYaxis()->SetTitle("Data/MC");
        hRatio->GetXaxis()->SetTitle("#Delta R");

        hRatio->GetYaxis()->SetNdivisions(505);

        // Axis formatting for small ratio pad
        hRatio->GetXaxis()->SetTitleSize(0.10);
        hRatio->GetXaxis()->SetLabelSize(0.08);

        hRatio->GetYaxis()->SetTitleSize(0.08);
        hRatio->GetYaxis()->SetTitleOffset(0.5);
        hRatio->GetYaxis()->SetLabelSize(0.08);

        hRatio->SetMinimum(0.5);
        hRatio->SetMaximum(1.5);

        hRatio->SetMarkerStyle(20);
        hRatio->SetMarkerSize(0.7);

        hRatio->Draw("EP");

        // Ratio = 1 line
        TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(),1.0,hRatio->GetXaxis()->GetXmax(),1.0);

        line->SetLineStyle(kDashed);
        line->Draw("same");
    }
    
    // Delta R vs PtD histogram
    TCanvas* cDeltaR_vs_PtD = new TCanvas(Form("cDeltaR_vs_PtD_%.0f_%.0f", jetptMin, jetptMax),Form("#Delta R vs PtD histogram %.0f < pT,jet < %0.f GeV/c", jetptMin, jetptMax),1920,1080);
    hDeltaR_vs_PtD->Draw("colz");
    
    //
    // Storing images in a single pdf file
    //
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../Images/5-ClosureTest/Second/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    cTemplateFits2D->Print(imagePath + Form("sb_subtraction_deltaR_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cTemplateFits->Print(imagePath + Form("sb_subtraction_deltaR_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cDataFits2D->Print(imagePath + Form("sb_subtraction_deltaR_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cDataFits->Print(imagePath + Form("sb_subtraction_deltaR_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    c1d_sideBandSubtracted->Print(imagePath + Form("sb_subtraction_deltaR_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cDeltaR_vs_PtD->Print(imagePath + Form("sb_subtraction_deltaR_" + sEmmaBins + "_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
}

std::pair<TH2D*, std::vector<bool>> AnalyzeJetPtRange(TFile* fClosureInput, const double& jetptMin, const double& jetptMax, 
                                                      const BinningStruct& binning,
                                                      const int& signalSigmas, const int& startingBackSigma, const int& backgroundSigmas,
                                                      const FitModelType& modelToUse) {
    //
    TH2D* hDeltaR_vs_PtD;
    // mass histogram
    int massBins = 50; // default=100 
    double minMass = 1.72; // use from 1.72, used to use 1.67
    double maxMass = 2.06; // old = 2.1, new = 2.05 + 0.01
    // Mass and std reference values for fits
    double m_0_reference = 1.86484; // D0 mass in GeV/c^2
    double sigma_reference = 0.012;
    std::cout << "Starting fit templates calculation for jet pT range: " << jetptMin << " - " << jetptMax << " GeV/c.\n";
    // ----- Obtain template fits from MC particle level data
    std::pair<FitContainer, HistogramGroup> fTemplateFitsAndCanvas = calculateFitTemplates(fClosureInput, jetptMin, jetptMax, binning, massBins, minMass, maxMass);
    FitContainer fTemplateFits = fTemplateFitsAndCanvas.first;
    HistogramGroup histogramTemplates = fTemplateFitsAndCanvas.second;
    std::cout << "Finished calculating template fits.\n";

    //
    // Perform side-band subtraction procedure on input data (not template)
    //
    // ----- Create and fill invariant mass distributions with detector level data
    std::vector<TH2D*> hInvariantMass2D;
    std::vector<TH2D*> hInvariantMass2DOracle;
    for (size_t i = 0; i < binning.ptHFBinEdges_detector.size() - 1; ++i) {
        // So that the title adapts to fractional binning title
        if (std::fmod(binning.ptHFBinEdges_detector[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(binning.ptHFBinEdges_detector[i+1], 1.0) != 0) {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // MC truth version, only real D⁰s (no reflections, no background)
                hInvariantMass2DOracle.push_back(new TH2D(Form("histMassOracle%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
            } else {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // MC truth version, only real D⁰s (no reflections, no background)
                hInvariantMass2DOracle.push_back(new TH2D(Form("histMassOracle%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
            }
        } else {
            if (std::fmod(binning.ptHFBinEdges_detector[i+1], 1.0) != 0) {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // MC truth version, only real D⁰s (no reflections, no background)
                hInvariantMass2DOracle.push_back(new TH2D(Form("histMassOracle%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
            } else {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
                // MC truth version, only real D⁰s (no reflections, no background)
                hInvariantMass2DOracle.push_back(new TH2D(Form("histMassOracle%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",binning.ptHFBinEdges_detector[i],binning.ptHFBinEdges_detector[i+1]), massBins, minMass, maxMass, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data()));
            }
        }
        hInvariantMass2D[i]->Sumw2();
        hInvariantMass2DOracle[i]->Sumw2();
    }
    std::cout << "Finished creating 2D invariant mass histograms.\n";
    // ----- Fill data
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCDyCut = 0.8; // on detector level D0
    const double MCDDeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double MCDHfPtMincut = binning.ptHFBinEdges_detector[0];;
    const double MCDHfPtMaxcut = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1];;
    const double MCDjetptMin = jetptMin;
    const double MCDjetptMax = jetptMax;
    
    // Accessing TTree
    TTree* tree = (TTree*)fClosureInput->Get("InputTree");
    // Check for correct access
    if (!tree) {
        std::cout << "Error opening correction data tree.\n";
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
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetNConst); // float
    MCPhfMass = 1.86484; // D0 rest mass in GeV/c^2
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
    std::cout << "Finished setting branch addresses for the TTree.\n";
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // Convert to enum class for type-safe comparison
        D0Species hfMatchedFrom = intToD0Species(MCDhfMatchedFrom);
        D0Species hfSelectedAs = intToD0Species(MCDhfSelectedAs);
        bool isReflection = (hfMatchedFrom != hfSelectedAs) ? true : false;

        // Apply prompt (real+reflections, only detector level) selection (no b feed-down correction is performed on the 2nd closure test)
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isRealD0 = isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs);
        if ((!MCDhfprompt) || !MCDhfmatch) {
            continue;
        }
        
        // calculating delta R
        double MCDdeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb) ? true : false;
        bool recoJetPtRange = (MCDjetPt >= MCDjetptMin) && (MCDjetPt < MCDjetptMax);
        bool recoHfPtRange = (MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut);
        bool recoDeltaRRange = (MCDdeltaR >= binning.deltaRBinEdges_detector[0]) && (MCDdeltaR < MCDDeltaRcut);
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && recoJetPtRange && recoHfPtRange && recoDeltaRRange; // currently used

        // Fill each histogram with their respective pT intervals
        if (recoLevelRange && passBDTcut) {
        //if ((abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= MCDjetptMin) && (MCDjetPt < MCDjetptMax)) && ((MCDdeltaR >= binning.deltaRBinEdges_detector[0]) && (MCDdeltaR < MCDDeltaRcut))) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < binning.ptHFBinEdges_detector.size() - 1 && !filled; iEdge++) {
                if ((MCDhfPt >= binning.ptHFBinEdges_detector[iEdge]) && (MCDhfPt < binning.ptHFBinEdges_detector[iEdge + 1])) {

                    // use Emma's run 2 cuts in case useEmmaYeatsBins is true
                    if (!binning.useEmmaYeatsBins || passEmmaCut(MCDjetPt, MCDhfPt)) {
                        hInvariantMass2D[iEdge]->Fill(MCDhfMass, MCDdeltaR);
                        // Fill the "oracle" histogram with true D0s (no reflections, no background)
                        if (isRealD0) {
                            hInvariantMass2DOracle[iEdge]->Fill(MCDhfMass, MCDdeltaR);
                        }
                        
                    }
                    filled = true; // Exit the loop once the correct histogram is found (alternative: break)
                }
            }
        }
        
    }
    std::cout << "Finished filling 2D invariant mass histograms with data from the TTree.\n";
    // ----- Obtain 1D projections: invariant mass distributions
    std::vector<TH1D*> hInvariantMass1D;
    std::vector<TH1D*> hInvariantMass1DOracle; // only real D⁰s (no reflections, no background)
    for (size_t iHisto = 0; iHisto < hInvariantMass2D.size(); iHisto++) {

        hInvariantMass1D.push_back(hInvariantMass2D[iHisto]->ProjectionX(Form("h_mass_proj_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax)));
        hInvariantMass1DOracle.push_back(hInvariantMass2DOracle[iHisto]->ProjectionX(Form("h_mass_proj_oracle_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax)));
    }
    std::cout << "Finished filling histograms with data and obtaining 1D projections.\n";
    // ----- Fit input detector level data using template constraints
    FitContainer fDataFits;
    for (size_t iHisto = 0; iHisto < hInvariantMass2D.size(); iHisto++) {

        double A1Signal = fTemplateFits.fitSignalOnly[iHisto]->GetParameter(0); // A1Signal
        double A2Signal = fTemplateFits.fitSignalOnly[iHisto]->GetParameter(1); // A2Signal
        double m0Signal = fTemplateFits.fitSignalOnly[iHisto]->GetParameter(2); // m0
        double sigma1Signal = fTemplateFits.fitSignalOnly[iHisto]->GetParameter(3); // sigma1
        double sigma2Signal = fTemplateFits.fitSignalOnly[iHisto]->GetParameter(4); // sigma2
        //std::cout << "A1Signal = " << A1Signal << ", A2Signal = " << A2Signal << ", m0Signal = " << m0Signal << ", sigma1Signal = " << sigma1Signal << ", sigma2Signal = " << sigma2Signal << std::endl;

        double A1Reflection = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(5); // A1Reflection
        double A2Reflection = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(6); // A2Reflection
        double m0_1Reflections = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(7); // m0_1
        double m0_2Reflections = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(8); // m0_2
        double sigma1Reflections = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(9); // sigma1
        double sigma2Reflections = fTemplateFits.fitReflectionsOnly[iHisto]->GetParameter(10); // sigma2
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
        if (modelToUse == FitModelType::FullPowerLaw) {
            // Signal + Reflections + Background
            fDataFits.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<FullModelPowerLaw>, minMass, maxMass, 13));
        } else if (modelToUse == FitModelType::SignalReflectionsOnly) {
            // Signal + Reflections only (no background)
            fDataFits.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<SigRefModel>, minMass, maxMass, 13));
        } else if (modelToUse == FitModelType::SignalOnly) {
            // Signal only (no reflections, no background)
            fDataFits.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<PureSignalModel>, minMass, maxMass, 13));
        }

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
        fDataFits.fitTotal[iHisto]->SetParLimits(0, 0., TMath::Infinity()); // only accepts non-negative background fits
        fDataFits.fitTotal[iHisto]->SetParLimits(4, 0.99 * m_0_reference, 1.01 * m_0_reference); // D0 invariant mass
        fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 3.0 * sigma_reference); // primary gaussian standard deviation
        fDataFits.fitTotal[iHisto]->SetParLimits(7, 1.0, TMath::Infinity()); // signal amplitude always larger than reflections'
        // Different sigma for each pT,HF bin
        if (iHisto == 0) {          // 1 < pT,D < 2 GeV/c
            sigma_reference = 0.0105;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0105
        } else if (iHisto == 1) {   // 2 < pT,D < 3 GeV/c
            sigma_reference = 0.0123;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0123
        } else if (iHisto == 2) {   // 3 < pT,D < 4 GeV/c
            sigma_reference = 0.0139;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 2.0 * sigma_reference); // sigma_1 = 0.0139
        } else if (iHisto == 3) {   // 4 < pT,D < 5 GeV/c
            sigma_reference = 0.0157;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 2.0 * sigma_reference); // sigma_1 = 0.0157
        } else if (iHisto == 4) {   // 5 < pT,D < 6 GeV/c
            sigma_reference = 0.0172;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0172
        } else if (iHisto == 5) {   // 6 < pT,D < 7 GeV/c
            sigma_reference = 0.0185;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0185
        } else if (iHisto == 6) {   // 7 < pT,D < 8 GeV/c
            sigma_reference = 0.0199;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0199
        } else if (iHisto == 7) {   // 8 < pT,D < 12 GeV/c
            sigma_reference = 0.0223;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0223
        } else if (iHisto == 8) {   // 12 < pT,D < 16 GeV/c
            sigma_reference = 0.0240;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.3 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0240
        } else if (iHisto == 9) {   // 16 < pT,D < 24 GeV/c
            sigma_reference = 0.0240;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 1.0 * sigma_reference, 2.0 * sigma_reference); // sigma_1 = 0.0240
        } else if (iHisto == 10) {  // 24 < pT,D < 36 GeV/c
            sigma_reference = 0.0373;
            fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.2 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0373
        }
        // sigma_reference = 0.012;
        // fDataFits.fitTotal[iHisto]->SetParLimits(5, 0.35 * sigma_reference, 2.0 * sigma_reference);
        // Choose line color
        fDataFits.fitTotal[iHisto]->SetLineColor(kBlack);
        // Perform fit with "Q" (quiet) option: no drawing of the fit function
        hInvariantMass1D[iHisto]->Fit(fDataFits.fitTotal[iHisto], "Q0");

        // Check if signal gaussian acomodates the literature D0 rest mass
        double primaryMean = fDataFits.fitTotal[iHisto]->GetParameter(4);
        double primarySigma = fDataFits.fitTotal[iHisto]->GetParameter(5);
        double D0LiteratureMass = 1.86484; // GeV/c²
        if (didFitFailed(hInvariantMass1D[iHisto], fDataFits.fitTotal[iHisto], D0LiteratureMass, signalSigmas, startingBackSigma)) { // fit failed: PDG mass outside of signal region
            // if literature mass is outside, clear histogram as it should be excluded
            //histograms[iHisto]->Reset();
            // hInvariantMass1D[iHisto]->SetTitle(TString(hInvariantMass1D[iHisto]->GetTitle()) + " (Fit failed)");
            //histograms2d[iHisto]->Reset("ICES");
            //histograms2d[iHisto]->Reset();
            hInvariantMass2D[iHisto]->SetTitle(TString(hInvariantMass2D[iHisto]->GetTitle()) + " (Fit failed)");

            fDataFits.workingFits.push_back(false);
        } else { // fit worked: PDG mass inside of signal region
            fDataFits.workingFits.push_back(true);
        }
        
        // Getting total parameter values
        double a_par = fDataFits.fitTotal[iHisto]->GetParameter(0); // Get the value of parameter 'a'
        double b_par = fDataFits.fitTotal[iHisto]->GetParameter(1); // Get the value of parameter 'b'
        fDataFits.fitBackgroundOnly.push_back(new TF1(Form("backgroundOnlyFit_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax), fitWrapper<PowerLawBackgroundModel>, minMass, maxMass, 2));
        fDataFits.fitBackgroundOnly[iHisto]->SetParameters(a_par,b_par); // fDataFits.size()-1 = the latest added to the vector
        fDataFits.fitBackgroundOnly[iHisto]->SetLineStyle(kDashed);
        fDataFits.fitBackgroundOnly[iHisto]->SetLineColor(kRed+1);

        // Extract signal-related parameters from the total fit
        double a1Signal = fDataFits.fitTotal[iHisto]->GetParameter(2);              // A1 signal
        A1toA2MCSignalRatio = fDataFits.fitTotal[iHisto]->GetParameter(3);          // A2 signal
        m0Signal = fDataFits.fitTotal[iHisto]->GetParameter(4);                     // Signal mean m0
        double sigmaSignal = fDataFits.fitTotal[iHisto]->GetParameter(5);           // Signal width sigma
        Sigma1toSigma2MCSignalRatio = fDataFits.fitTotal[iHisto]->GetParameter(6);  // Sigma1/Sigma2 signal
        // Perfoming fit with acquired parameters from total fit
        fDataFits.fitSignalOnly.push_back(new TF1(Form("signalOnlyFit_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax), fitWrapper<PureSignalModel>, minMass, maxMass, 7));
        fDataFits.fitSignalOnly[iHisto]->SetParameters(a1Signal, A1toA2MCSignalRatio, m0Signal, sigmaSignal, Sigma1toSigma2MCSignalRatio); // fDataFits.size()-1 = the latest added to the vector
        fDataFits.fitSignalOnly[iHisto]->FixParameter(2, a1Signal);
        fDataFits.fitSignalOnly[iHisto]->FixParameter(3, A1toA2MCSignalRatio);
        fDataFits.fitSignalOnly[iHisto]->FixParameter(4, m0Signal);
        fDataFits.fitSignalOnly[iHisto]->FixParameter(5, sigmaSignal);
        fDataFits.fitSignalOnly[iHisto]->FixParameter(6, Sigma1toSigma2MCSignalRatio);
        fDataFits.fitSignalOnly[iHisto]->SetLineStyle(kDashed);
        fDataFits.fitSignalOnly[iHisto]->SetLineColor(kBlue+1);

        // Extract reflection-related parameters from the total fit
        double A1SignalToA1ReflectionMCRatios = fDataFits.fitTotal[iHisto]->GetParameter(7);    // A1 signal / A1 reflection
        A1Signal = fDataFits.fitTotal[iHisto]->GetParameter(2);                                 // A1 signal
        A1toA2MCReflectionsRatio = fDataFits.fitTotal[iHisto]->GetParameter(8);                 // A1 reflection / A2 reflection
        double m0_1Reflection = fDataFits.fitTotal[iHisto]->GetParameter(9);                    // Reflection mean m0_1
        double m0_2Reflection = fDataFits.fitTotal[iHisto]->GetParameter(10);                   // Reflection mean m0_2
        double sigma1Reflection = fDataFits.fitTotal[iHisto]->GetParameter(11);                 // Reflection sigma_1
        double sigma2Reflection = fDataFits.fitTotal[iHisto]->GetParameter(12);                 // Reflection sigma_2
        //std:cout << "A1SignalToA1ReflectionMCRatios = " << A1SignalToA1ReflectionMCRatios << ", A1Signal = " << A1Signal << ", A1toA2MCReflectionsRatio = " << A1toA2MCReflectionsRatio << ", m0_1Reflection = " << m0_1Reflection << ", m0_2Reflection = " << m0_2Reflection << ", sigma1Reflection = " << sigma1Reflection << ", sigma2Reflection = " << sigma2Reflection << std::endl;

        // Perfoming fit with acquired parameters from total fit
        fDataFits.fitReflectionsOnly.push_back(new TF1(Form("reflectionOnlyFit_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax), fitWrapper<PureReflectionsModel>, minMass, maxMass, 13));
        fDataFits.fitReflectionsOnly[iHisto]->SetParameters(A1SignalToA1ReflectionMCRatios, A1Signal, A1toA2MCReflectionsRatio, m0_1Reflection, m0_2Reflection, sigma1Reflection, sigma2Reflection); // fDataFits.size()-1 = the latest added to the vector
        fDataFits.fitReflectionsOnly[iHisto]->FixParameter(7, A1SignalToA1ReflectionMCRatios);
        fDataFits.fitReflectionsOnly[iHisto]->FixParameter(2, A1Signal);
        fDataFits.fitReflectionsOnly[iHisto]->FixParameter(8, A1toA2MCReflectionsRatio);
        fDataFits.fitReflectionsOnly[iHisto]->FixParameter(9, m0_1Reflection);
        fDataFits.fitReflectionsOnly[iHisto]->FixParameter(10, m0_2Reflection);
        fDataFits.fitReflectionsOnly[iHisto]->FixParameter(11, sigma1Reflection);
        fDataFits.fitReflectionsOnly[iHisto]->FixParameter(12, sigma2Reflection);
        fDataFits.fitReflectionsOnly[iHisto]->SetLineStyle(kDashed);
        fDataFits.fitReflectionsOnly[iHisto]->SetLineColor(kGreen+3);
    }
    std::cout << "Finished performing fits on DATA histograms.\n";
    // Side-band method
    std::vector<TH1D*> hsideBandSubtracted;
    // Creating histograms for collecting data
    TH1D* tempHist; // temporary histogram for collecting data
    TH1D* h_sideBand;
    TH1D* h_signal;
    TH1D* h_back_subtracted;
    for (size_t iHisto = 0; iHisto < hInvariantMass2D.size(); iHisto++) {
        
        // Get total fit parameters
        double m_0 = fDataFits.fitTotal[iHisto]->GetParameter(4); // Get the value of parameter 'm_0'
        double sigma = fDataFits.fitTotal[iHisto]->GetParameter(5); // Get the value of parameter 'sigma1'

        // Obtain signal histogram
        // int lowBin = hInvariantMass2D[iHisto]->GetXaxis()->FindBin(m_0 - signalSigmas * sigma);
        // int highBin = hInvariantMass2D[iHisto]->GetXaxis()->FindBin(m_0 + signalSigmas * sigma);
        int lowBin = safeLowBin(hInvariantMass2D[iHisto], m_0 - signalSigmas * sigma);
        int highBin = safeHighBin(hInvariantMass2D[iHisto], m_0 + signalSigmas * sigma);
        h_signal = hInvariantMass2D[iHisto]->ProjectionY(Form("h_signal_proj_%zu_%0.f_to_%0.fGeV",iHisto, jetptMin, jetptMax), lowBin, highBin);
        if ((modelToUse == FitModelType::SignalReflectionsOnly) || (modelToUse == FitModelType::SignalOnly)) { // no background, just reflections removal
            std::cout << "Background fit empty for histogram " << iHisto << ", skipping side-band subtraction." << std::endl;
            // Subtracted histogram = original signal histogram
            h_back_subtracted = (TH1D*)h_signal->Clone(Form("h_back_subtracted_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax));

            // Calculate fit areas under signal region for scaling factors
            double Rs = fDataFits.fitReflectionsOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // Reflections component in signal region area
            double S = fDataFits.fitSignalOnly[iHisto]->Integral(m_0 - signalSigmas * sigma, m_0 + signalSigmas * sigma); // Signal component in signal region area
            if (std::isnan(S)) {
                std::cout << "Warning: Signal area S is NaN for histogram " << iHisto << ". Setting S to 0 and resetting histogram." << std::endl;
                S = 0.0;
                h_back_subtracted->Reset();
            } else if (std::isnan(Rs)) {
                std::cout << "Warning: Reflections area Rs is NaN for histogram " << iHisto << ". Setting Rs to 0." << std::endl;
                Rs = 0.0;
                // Scale by S/(S+Rs) to remove reflections
                h_back_subtracted->Scale(S / (S + Rs));
                // Account for two sigma_primary only area used for signal region
                double coverage = TMath::Erf(signalSigmas / sqrt(2)); // coverage of a Gaussian in a +/- signalSigmas window
                double twoSigmaArea = fDataFits.fitSignalOnly[iHisto]->Integral(m_0 - signalSigmas * sigma, m_0 + signalSigmas * sigma);
                double totalBoundsArea = fDataFits.fitSignalOnly[iHisto]->Integral(minMass, maxMass);
                coverage = twoSigmaArea / totalBoundsArea;
                std::cout << "2 sigmas for signal region has coverage of " << coverage << std::endl;
                //std::cout << "Total scaling factor of " <<  (S / (S + Rs)) / TMath::Erf(sqrt(2)) << std::endl;
                h_back_subtracted->Scale(1 / coverage);
            } else if (std::isnan(S) && std::isnan(Rs)) {
                std::cout << "Warning: Both Signal area S and Reflections area Rs are NaN for histogram " << iHisto << ". Setting both to 1 and resetting histogram." << std::endl;
                S = 1.0;
                Rs = 1.0;
                h_back_subtracted->Reset();
            } else {
                std::cout << "Signal area S = " << S << ", Reflections area Rs = " << Rs << " for histogram " << iHisto << std::endl;
                // Scale by S/(S+Rs) to remove reflections
                h_back_subtracted->Scale(S / (S + Rs));
                // Account for two sigma only area used for signal region
                double coverage = TMath::Erf(signalSigmas / sqrt(2)); // coverage of a Gaussian in a +/- signalSigmas window
                double twoSigmaArea = fDataFits.fitSignalOnly[iHisto]->Integral(m_0 - signalSigmas * sigma, m_0 + signalSigmas * sigma);
                double totalBoundsArea = fDataFits.fitSignalOnly[iHisto]->Integral(minMass, maxMass);
                coverage = twoSigmaArea / totalBoundsArea;
                std::cout << "2 sigmas for signal region has coverage of " << coverage << std::endl;
                //std::cout << "Total scaling factor of " <<  (S / (S + Rs)) / TMath::Erf(sqrt(2)) << std::endl;
                h_back_subtracted->Scale(1 / coverage);
            }
            
            hsideBandSubtracted.push_back(h_back_subtracted);

        } else {
            std::array<double, 2> leftRange = {0., 0.}; // remains zero if no sigmas fit inside the left range
            std::array<double, 2> rightRange = {0., 0.};
            
            // Count for how many sigmas is there room inside the left side range
            double leftSidebandRange = (m_0 - startingBackSigma * sigma) - hInvariantMass1D[iHisto]->GetBinLowEdge(1);
            leftSidebandRange = std::max(leftSidebandRange, 0.0); // Ensure non-negative range
            double leftSigmas = leftSidebandRange / sigma;// get the fractional number
            std::cout << leftSigmas << " sigmas fit inside the left side range for the background distribution estimation." << std::endl;

            // Count for how many sigmas is there room inside the right side range
            double rightSidebandRange = hInvariantMass1D[iHisto]->GetBinLowEdge(hInvariantMass1D[iHisto]->GetNbinsX()+1) - (m_0 + startingBackSigma*sigma);
            rightSidebandRange = std::max(rightSidebandRange, 0.0); // Ensure non-negative range
            double rightSigmas = rightSidebandRange / sigma;// get the fractional number
            std::cout << rightSigmas << " sigmas fit inside the right side range for the background distribution estimation." << std::endl;

            // Ensure only a maximum of 4 sigmas (backgroundSigmas) are used
            if (leftSigmas > backgroundSigmas) {
                leftSigmas = backgroundSigmas;
            }
            if (rightSigmas > backgroundSigmas) {
                rightSigmas = backgroundSigmas;
            }

            // Calculate left range limits
            leftRange[0] = m_0 - (startingBackSigma + leftSigmas) * sigma;
            leftRange[1] = m_0 - startingBackSigma * sigma;

            // Calculate right range limits
            rightRange[0] = m_0 + startingBackSigma * sigma;
            rightRange[1] = m_0 + (startingBackSigma + rightSigmas) * sigma;

            // If there is no room for the left side range, use only right side-band
            if ((leftRange[0] - leftRange[1]) == 0.) {
                // use only right side-band
                // lowBin = hInvariantMass2D[iHisto]->GetXaxis()->FindBin(rightRange[0]); // ToDo: bug fix, FindBin can be looking for underflow/overflow of side-bands
                // highBin = hInvariantMass2D[iHisto]->GetXaxis()->FindBin(rightRange[1]);
                lowBin = safeLowBin(hInvariantMass2D[iHisto], rightRange[0]);
                highBin = safeHighBin(hInvariantMass2D[iHisto], rightRange[1]);
                h_sideBand = hInvariantMass2D[iHisto]->ProjectionY(Form("h_sideband_proj_temp_%zu_%0.f_to_%0.fGeV",iHisto, jetptMin, jetptMax), lowBin, highBin); // sum the right sideband
                double rightSBHistogram = hInvariantMass1D[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram
            } else {
                // Use both sidebands
                // lowBin = hInvariantMass2D[iHisto]->GetXaxis()->FindBin(leftRange[0]);
                // highBin = hInvariantMass2D[iHisto]->GetXaxis()->FindBin(leftRange[1]);
                lowBin = safeLowBin(hInvariantMass2D[iHisto], leftRange[0]);
                highBin = safeHighBin(hInvariantMass2D[iHisto], leftRange[1]);
                h_sideBand = hInvariantMass2D[iHisto]->ProjectionY(Form("h_sideband_proj_temp_left_%zu_%0.f_to_%0.fGeV",iHisto, jetptMin, jetptMax), lowBin, highBin); // sum the left sideband
                double leftSBHistogram = hInvariantMass1D[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram

                // lowBin = hInvariantMass2D[iHisto]->GetXaxis()->FindBin(rightRange[0]);
                // highBin = hInvariantMass2D[iHisto]->GetXaxis()->FindBin(rightRange[1]);
                lowBin = safeLowBin(hInvariantMass2D[iHisto], rightRange[0]);
                highBin = safeHighBin(hInvariantMass2D[iHisto], rightRange[1]);
                h_sideBand->Add(hInvariantMass2D[iHisto]->ProjectionY(Form("h_sideband_proj_temp_%zu_%0.f_to_%0.fGeV",iHisto, jetptMin, jetptMax), lowBin, highBin)); // sum the right sideband
                double rightSBHistogram = hInvariantMass1D[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram
            }
            
            // Calculate scaling factors for sideband subtraction
            // Calculate scaling factor to apply to the sideband subtraction original method
            std::array<double, 2> scallingFactors = calculateScalingFactor(iHisto, fDataFits, leftRange, rightRange, signalSigmas);
            std::cout << "scallingFactors[0] = " << scallingFactors[0] << ", scallingFactors[1] = " << scallingFactors[1] << std::endl;
            
            // Scale the sideband histogram by ratio of it in the signal region Bs/(B1+B2)
            h_sideBand->Scale(scallingFactors[1]);

            // Subtract the sideband histogram from the signal histogram
            h_back_subtracted = (TH1D*)h_signal->Clone(Form("h_back_subtracted_%zu_%0.f_to_%0.fGeV", iHisto, jetptMin, jetptMax));
            h_back_subtracted->Add(h_sideBand,-1.0);

            // Scale by reflections correction
            h_back_subtracted->Scale(scallingFactors[0]);

            // Account for two sigma only area used for signal region
            double coverage = TMath::Erf(signalSigmas / sqrt(2)); // coverage of a Gaussian in a +/- signalSigmas window
            double twoSigmaArea = fDataFits.fitSignalOnly[iHisto]->Integral(m_0 - signalSigmas * sigma, m_0 + signalSigmas * sigma);
            double totalBoundsArea = fDataFits.fitSignalOnly[iHisto]->Integral(minMass, maxMass);
            coverage = twoSigmaArea / totalBoundsArea;
            std::cout << "2 sigmas for signal region has coverage of " << coverage << std::endl;
            // std::cout << "Total scaling factor of " <<  (S / (S + Rs)) / TMath::Erf(sqrt(2)) << std::endl;
            h_back_subtracted->Scale(1 / coverage);
            
            hsideBandSubtracted.push_back(h_back_subtracted);
        }
    }
    // Final adjustments to subtracted histograms
    for (size_t iHisto = 0; iHisto < hsideBandSubtracted.size(); iHisto++) {
        // Check if histogram should be erased before storing based on fit performed
        bool eraseHist = eraseHistogram(fDataFits.workingFits, iHisto);
        if (eraseHist) {
            hsideBandSubtracted[iHisto]->Reset("ICES");
        } else {
            // Set negative count bin entries to 0
            for (int iBin = 1; iBin <= hsideBandSubtracted[iHisto]->GetNbinsX(); iBin++) {
                if (hsideBandSubtracted[iHisto]->GetBinContent(iBin) < 0) {
                    hsideBandSubtracted[iHisto]->SetBinContent(iBin,0);
                    hsideBandSubtracted[iHisto]->SetBinError(iBin,0);
                }
            }
        }
    }
    std::cout << "[workingFits] pattern: ";
    for (size_t i = 0; i < fDataFits.workingFits.size(); ++i) {
        std::cout << (fDataFits.workingFits[i] ? "OK" : "FAIL");
        if (i < fDataFits.workingFits.size() - 1) std::cout << " | ";
    }
    std::cout << std::endl;
    std::cout << "Finished performing side-band subtraction.\n";
    // ----- Create 2D distribution of DeltaR vs. pT,D0
    int nPtBins = binning.ptHFBinEdges_detector.size() - 1;
    int nDeltaRBins = binning.deltaRBinEdges_detector.size() - 1;
    hDeltaR_vs_PtD = new TH2D(Form("hDeltaR_vs_PtD_%0.f_to_%0.fGeV", jetptMin, jetptMax), "#DeltaR vs #it{p}_{T,D^{0}};#DeltaR;#it{p}_{T,D^{0}} (GeV/#it{c})",
                                  nDeltaRBins, binning.deltaRBinEdges_detector.data(),  // X: ΔR
                                  nPtBins, binning.ptHFBinEdges_detector.data()); // Y: pT,D);
    for (size_t iHisto = 0; iHisto < hsideBandSubtracted.size(); ++iHisto) {
        TH1D* h1D = hsideBandSubtracted[iHisto];
        
        for (int iBin = 1; iBin <= h1D->GetNbinsX(); ++iBin) {
            double deltaR_center = h1D->GetBinCenter(iBin);
            double ptD_center = 0.5*(binning.ptHFBinEdges_detector[iHisto] + binning.ptHFBinEdges_detector[iHisto+1]);
            double content = h1D->GetBinContent(iBin);
            double contentError = h1D->GetBinError(iBin);
            // if (hsideBandSubtracted[iHisto]->GetEntries() == 0) {
            //     std::cout << "Histo with no entries: content = " << content << ", contentError = " << contentError << std::endl;
            //     content = 0;
            //     contentError = 0;
            // } else if (content < 0) {
            //     std::cout << "Histo with negative entries: content = " << content << ", contentError = " << contentError << std::endl;
            //     content = 0;
            //     contentError = 0;
            // } else if (std::isnan(content) || std::isnan(contentError)) {
            //     std::cout << "Histo with nan entries: content = " << content << ", contentError = " << contentError << std::endl;
            // }
            

            // Find the corresponding bin numbers in the 2D histogram
            int xBin = hDeltaR_vs_PtD->GetXaxis()->FindBin(deltaR_center);
            int yBin = hDeltaR_vs_PtD->GetYaxis()->FindBin(ptD_center);

            hDeltaR_vs_PtD->SetBinContent(xBin, yBin, content);
            hDeltaR_vs_PtD->SetBinError(xBin, yBin, contentError);
        }
    }
    std::cout << "Plotting and storing data to PDF file..." << std::endl;
    // Plot histograms and fits for visualization in a file
    PlotHistograms(histogramTemplates, fTemplateFits, hInvariantMass2D, hInvariantMass1D, hInvariantMass2DOracle, hInvariantMass1DOracle, fDataFits, hsideBandSubtracted, hDeltaR_vs_PtD, jetptMin, jetptMax, modelToUse, binning);
    std::cout << "Sideband closure test completed for pT,jet range: " << jetptMin << " - " << jetptMax << " GeV/c" << std::endl;
    return std::make_pair(hDeltaR_vs_PtD, fDataFits.workingFits);
}
SidebandClosureResult create3DBackgroundSubtracted(const std::vector<TH2D*>& outputHistograms, SidebandClosureResult& sidebandDataContainer, const BinningStruct& binning) {

    int nDeltaRBins = binning.deltaRBinEdges_detector.size() - 1;
    int nPtDBins   = binning.ptHFBinEdges_detector.size() - 1;
    int nPtJetBins = binning.ptjetBinEdges_detector.size() - 1;

    sidebandDataContainer.hBackgroundSubtracted = new TH3D("hBackgroundSubtracted",
                        "Background-subtracted 3D histogram;p_{T,jet} (GeV/c);#DeltaR;p_{T,D} (GeV/c)",
                        nPtJetBins, binning.ptjetBinEdges_detector.data(),   // X axis → pT,jet
                        nDeltaRBins, binning.deltaRBinEdges_detector.data(), // Y axis → ΔR
                        nPtDBins, binning.ptHFBinEdges_detector.data());     // Z axis → pT,D
    for (size_t iPtJet = 0; iPtJet < outputHistograms.size(); ++iPtJet) {
        TH2D* h2D = outputHistograms[iPtJet];

        double ptJetCenter = 0.5 * (binning.ptjetBinEdges_detector[iPtJet] + binning.ptjetBinEdges_detector[iPtJet + 1]);

        for (int iX = 1; iX <= h2D->GetNbinsX(); ++iX) {       // ΔR bins
            double deltaR_center = h2D->GetXaxis()->GetBinCenter(iX);

            for (int iY = 1; iY <= h2D->GetNbinsY(); ++iY) {   // pT,D bins
                double ptD_center = h2D->GetYaxis()->GetBinCenter(iY);
                double content    = h2D->GetBinContent(iX, iY);
                double error      = h2D->GetBinError(iX, iY);

                if (std::isnan(content)) {
                    content = 0.;
                    error = 0.;
                }

                // Fill the 3D histogram
                int binX3D = sidebandDataContainer.hBackgroundSubtracted->GetXaxis()->FindBin(ptJetCenter);
                int binY3D = sidebandDataContainer.hBackgroundSubtracted->GetYaxis()->FindBin(deltaR_center);
                int binZ3D = sidebandDataContainer.hBackgroundSubtracted->GetZaxis()->FindBin(ptD_center);
                sidebandDataContainer.hBackgroundSubtracted->SetBinContent(binX3D, binY3D, binZ3D, content);
                sidebandDataContainer.hBackgroundSubtracted->SetBinError(binX3D, binY3D, binZ3D, error);

            }
        }
    }

    // Clean NaN values in the 3D histogram (if any)
    cleanNaNs(sidebandDataContainer.hBackgroundSubtracted);

    return sidebandDataContainer;
}

SidebandClosureResult SidebandClosure(TFile* fClosureInput, const BinningStruct& binning, const FitModelType& modelToUse) {

    // Hardcoded sideband subtraction parameters
    int signalSigmas = 2; // Number of sigmas for signal distribution
    int startingBackSigma = 4; // beginning of background distribution
    int backgroundSigmas = 4; // (maximum) range of background distribution after its starting point
    // One TH2D for each pT,jet range computed
    std::vector<TH2D*> outputHistograms;

    // Testing range of pT,jet bins
    //ptjetBinEdges_detector = {30., 50.};
    
    // Create container to store 3D histogram and working fits vector for each pT,jet range
    SidebandClosureResult sidebandDataContainer;

    for (size_t iPtJetRange = 0; iPtJetRange < binning.ptjetBinEdges_detector.size() - 1; iPtJetRange++) {
        std::cout << std::endl << std::endl << "Processing pT,jet range: " << binning.ptjetBinEdges_detector[iPtJetRange] << " - " << binning.ptjetBinEdges_detector[iPtJetRange + 1] << " GeV/c" << std::endl;
        std::pair<TH2D*, std::vector<bool>> hDeltaR_vs_PtD = AnalyzeJetPtRange(fClosureInput, binning.ptjetBinEdges_detector[iPtJetRange], binning.ptjetBinEdges_detector[iPtJetRange + 1], binning,
                                                 signalSigmas, startingBackSigma, backgroundSigmas,
                                                 modelToUse);
        std::cout << "pT,jet range: " << binning.ptjetBinEdges_detector[iPtJetRange] << " - " << binning.ptjetBinEdges_detector[iPtJetRange + 1] << " GeV/c processed. Storing histogram..." << std::endl;
        outputHistograms.push_back(hDeltaR_vs_PtD.first);
        sidebandDataContainer.workingFitsPerJetPt.push_back(hDeltaR_vs_PtD.second);
    }
    std::cout << "All pT,jet ranges processed. Creating 3D histogram..." << std::endl;
    create3DBackgroundSubtracted(outputHistograms, sidebandDataContainer, binning);
    TCanvas* cBackgroundSubtracted = new TCanvas("cBackgroundSubtracted", "Background-subtracted 3D histogram",1920,1080);
    sidebandDataContainer.hBackgroundSubtracted->Draw("colz");
    return sidebandDataContainer;
}
