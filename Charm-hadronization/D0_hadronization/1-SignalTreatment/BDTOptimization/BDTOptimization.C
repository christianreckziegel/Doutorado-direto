/**
 * D0 meson analysis
 * @file BDTOptimization.C
 * @author Christian Reckziegel
 * @brief Macro for obtaining the optimal BDT score values (highest significance S/sqrt(S+B))
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
// Data models
using FullModelPowerLaw = FitModel<SignalPolicy, ReflectionPolicy, PowerLawBackgroundPolicy>;               // Signal + Reflection + power law Background
using FullModelPoly2 = FitModel<SignalPolicy, ReflectionPolicy, Poly2BackgroundPolicy>;                     // Signal + Reflection + 2nd order polynomial Background
using SigRefModel = FitModel<SignalPolicy, ReflectionPolicy>;                                               // Signal + Reflection only for "DATA"
using PureSignalModel = FitModel<SignalPolicy>;                                                             // Signal only for "DATA"
using PureReflectionsModel = FitModel<ReflectionPolicy>;                                                    // Reflection only for "DATA" template
using PurePowerLawModel = FitModel<PowerLawBackgroundPolicy>;
using PurePoly2Model = FitModel<Poly2BackgroundPolicy>;
// Template models
using SigRefTemplateModel = FitModel<SignalMCTemplatePolicy, ReflectionMCTemplatePolicy>;                   // Signal + Reflection only for MC template
using PureSignalTemplateModel = FitModel<SignalMCTemplatePolicy>;                                           // Signal only for MC template
using PureSingleGaussianTemplateModel = FitModel<SignalSinglePolicy>;                                       // Single Gaussian for individual contributions of double Gaussian
using PureReflectionsTemplateModel = FitModel<ReflectionMCTemplatePolicy>;                                  // Reflection only for MC template
//-----------------------------------------------------------------------------Structs---------------------------------------------------------------------------------------------------------
struct FitsGroup {
    // MC fits
    std::vector<TF1*> fMCPureReflectionsFit;
    std::vector<TF1*> fMCPureSignalsFit;
    std::vector<TF1*> fMCSignalAndReflectionsFit;

    std::vector<std::pair<TF1*,TF1*>> fIndividualMCPureReflectionsFit;          // primary Gaussian, secondary Gaussian
    std::vector<std::pair<TF1*,TF1*>> fIndividualMCPureSignalsFit;              // primary Gaussian, secondary Gaussian
    std::vector<std::pair<TF1*,TF1*>> fIndividualMCSignalAndReflectionsFit;     // primary Gaussian, secondary Gaussian

    // Data fits
    std::vector<TF1*> fitTotal;
    std::vector<TF1*> fitBackgroundOnly;
    std::vector<TF1*> fitSignalOnly;
    std::vector<TF1*> fitReflectionsOnly;

    // List of working fits
    std::vector<bool> workingFits;

    // Optimal BDT evaluation
    std::vector<double> significance;    // S/sqrt(S+B) per pT,HF bin
    std::vector<double> signalYield;     // S per pT,HF bin
    std::vector<double> backgroundYield; // B per pT,HF bin
};
struct BDTScanContainer {
    // --- MC ---
    // Histograms
    std::vector<std::pair<TH2D*,TH1D*>> hMCReflections;
    std::vector<std::pair<TH2D*,TH1D*>> hMCSignals;
    std::vector<std::pair<TH2D*,TH1D*>> hMCSignalAndReflections;

    // DATA
    // Histograms
    std::vector<std::pair<TH2D*,TH1D*>> hDataTotal;
    std::vector<TH1D*> sidebandHistDrawing;                             // histogram used for drawing the side-band shaded invariant mass distribution -in red
    std::vector<TH1D*> signalHistDrawing;                               // histogram used for drawing the signal shaded invariant mass distribution - in blue

    FitsGroup fittings;

    // Clear data for next iteration
    void ClearData() {
    // Histograms
        hMCReflections.clear();
        hMCSignals.clear();
        hMCSignalAndReflections.clear();
        hDataTotal.clear();       
        sidebandHistDrawing.clear();
        signalHistDrawing.clear();

        // Fits
        fittings.fMCPureReflectionsFit.clear();
        fittings.fMCPureSignalsFit.clear();    
        fittings.fMCSignalAndReflectionsFit.clear();
        fittings.fIndividualMCPureReflectionsFit.clear();
        fittings.fIndividualMCPureSignalsFit.clear(); 
        fittings.fIndividualMCSignalAndReflectionsFit.clear();
        fittings.fitTotal.clear();      
        fittings.fitBackgroundOnly.clear();
        fittings.fitSignalOnly.clear(); 
        fittings.fitReflectionsOnly.clear();
        fittings.workingFits.clear();    
        fittings.significance.clear();   
        fittings.signalYield.clear();    
        fittings.backgroundYield.clear();
    }
};
//-----------------------------------------------------------------------------Main workflow functions---------------------------------------------------------------------------------------------------------

std::vector<std::pair<TH2D*,TH1D*>> createHistograms(const BinningStruct& binning, TString histName = "") {

    // Validate that maxMass vector has the correct size
    if (binning.maxMass.size() != binning.ptHFBinEdges_detector.size() - 1) {
        std::cerr << "Error: maxMass vector size (" << binning.maxMass.size() << ") does not match number of pT,D bins (" << binning.ptHFBinEdges_detector.size() - 1 << ").\n";
        return {std::make_pair(nullptr,nullptr)};
    }

    std::vector<std::pair<TH2D*,TH1D*>> histogramContainer;

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

        // Create all three histogram types with the same title and axes
        auto makeHisto = [&](const char* prefix) -> TH2D* {
            
            TH2D* h = new TH2D(Form("%s%zu", prefix, i+1), title, nMassBins, massMin, massMax, binning.deltaRBinEdges_detector.size() - 1, binning.deltaRBinEdges_detector.data());
            h->Sumw2();
            return h;
        };

        TH2D* hTemp2d = makeHisto(Form("%s_histMass%zu", histName.Data(), i + 1));
        histogramContainer.push_back(std::make_pair(hTemp2d, hTemp2d->ProjectionX(Form("%s_projX",hTemp2d->GetName()))));
    }

    return histogramContainer;
}

struct BDTScanResult {
    std::vector<std::vector<double>> bdtThresholds;         // [iScan][iPtD]
    std::vector<std::vector<double>> significance;          // [iScan][iPtD]
    std::vector<double> optimalCuts;                        // [iPtD]
    std::vector<double> maxSignificance;                    // [iPtD]
    std::vector<std::vector<double>> bdtThresholdsFine;     // [iScan][iPtD]
    std::vector<std::vector<double>> significanceFine;      // [iScan][iPtD]

    // Run 2 style efficiency - one for each scan entry
    std::vector<TH1D*> hMcpPt;
    std::vector<TH1D*> hMcdPt;
    std::vector<TH1D*> efficiencyHist;
};
// Create custom BDT cuts object for each iteration
std::vector<std::pair<double, double>> createCustomBdtCuts(BDTScanResult& bdtResults, const BinningStruct& binning, const int& iteration, const int& numIterations) {
    // Start from the maximum allowed loose range of scores
    std::vector<std::pair<double, double>> bdtCustomCuts = binning.bdtPtCuts;

    // Story copy of score iteration for each pT,HF
    std::vector<double> iterScorePerPtHF;

    for (size_t iHFBin = 0; iHFBin < bdtCustomCuts.size(); iHFBin++) {
        iterScorePerPtHF.emplace_back(bdtCustomCuts[iHFBin].second * iteration / numIterations); // needs to be executed first not to double modify the cut
        bdtCustomCuts[iHFBin].second = bdtCustomCuts[iHFBin].second * iteration / numIterations;
        
    }
    // Store for current iteration
    bdtResults.bdtThresholds.emplace_back(iterScorePerPtHF);
    
    bool printCustomCuts = true;
    if (printCustomCuts) {
        std::cout << "Interval (GeV/c)\tscore" << std::endl;
        for (size_t iHFBin = 0; iHFBin < bdtCustomCuts.size()-1; iHFBin++) {
            std::cout << "\t" << bdtCustomCuts[iHFBin].first << "-" << bdtCustomCuts[iHFBin+1].first << "\t\t\t" << bdtCustomCuts[iHFBin].second << std::endl;
        }
        
    }
    
    return bdtCustomCuts;
}

void fillMCHistograms(TFile* fInputMC, const BinningStruct& binning, BDTScanContainer& dataContainer, const std::vector<std::pair<double,double>>& customBdtCuts) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCDetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDyCut = 0.8; // on detector level Lc
    const double DeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    double hfptMin = binning.ptHFBinEdges_detector[0]; //ptHFBinEdges[0] - should start from 0 or from the lowest pT,HF value?
    double hfptMax = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1];
    const double hfPtMincut = hfptMin; // on particle level HF
    const double hfPtMaxcut = hfptMax; // on particle level HF
    const double jetptMin = binning.ptjetBinEdges_detector[0];
    const double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

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
        
        double MCDDeltaR = MCDaxisDistance;
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, customBdtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb);
        bool recoAngularRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut);
        bool recoJetPtRange = ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax));
        bool recoHfPtRange = ((MCDhfPt >= binning.ptHFBinEdges_detector[0]) && (MCDhfPt < binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1]));
        bool recoDeltaRRange = ((MCDDeltaR >= binning.deltaRBinEdges_detector[0]) && (MCDDeltaR < DeltaRcut));
        bool recoLevelRange = passBDTcut && recoAngularRange && recoJetPtRange && recoHfPtRange && recoDeltaRRange;

        // only matched detector level entries
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        if (!MCDhfmatch) {
            continue;
        }

        if (recoLevelRange) {
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < binning.ptHFBinEdges_detector.size() - 1 && !filled; iEdge++) {
                if ((MCDhfPt >= binning.ptHFBinEdges_detector[iEdge]) && (MCDhfPt < binning.ptHFBinEdges_detector[iEdge + 1])) {
                    
                    if (isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs)) {                                                      // It is a signal entry
                        dataContainer.hMCSignals[iEdge].first->Fill(MCDhfMass, MCDDeltaR);
                    } else if (isReflection(MCDhfMatchedFrom, MCDhfSelectedAs)) {                                               // It is a reflection entry
                        dataContainer.hMCReflections[iEdge].first->Fill(MCDhfMass, MCDDeltaR);
                    }
                    if (isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs) || isReflection(MCDhfMatchedFrom, MCDhfSelectedAs)) {   // It is a signal or reflection entry
                        dataContainer.hMCSignalAndReflections[iEdge].first->Fill(MCDhfMass, MCDDeltaR);
                    }
                    
                }
            }
        }
    } // end of entries loop

    // Create 1D projections
    for (size_t iEdge = 0; iEdge < binning.ptHFBinEdges_detector.size() - 1; iEdge++) {
        dataContainer.hMCSignals[iEdge].second = dataContainer.hMCSignals[iEdge].first->ProjectionX(Form("%s_projX",dataContainer.hMCSignals[iEdge].first->GetName()));
        dataContainer.hMCReflections[iEdge].second = dataContainer.hMCReflections[iEdge].first->ProjectionX(Form("%s_projX",dataContainer.hMCReflections[iEdge].first->GetName()));
        dataContainer.hMCSignalAndReflections[iEdge].second = dataContainer.hMCSignalAndReflections[iEdge].first->ProjectionX(Form("%s_projX",dataContainer.hMCSignalAndReflections[iEdge].first->GetName()));
    }

    std::cout << "MC template histograms filled." << std::endl;
}

void fillDataHistograms(TFile* fInputData, const BinningStruct& binning, BDTScanContainer& dataContainer, const std::vector<std::pair<double,double>>& customBdtCuts) {
    //
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCDetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDyCut = 0.8; // on detector level D0
    const double DeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on particle level delta R
    const double jetptMin = binning.ptjetBinEdges_detector[0];
    const double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

    // Accessing TTree
    TTree* tree = (TTree*)fInputData->Get("DF_merged/O2jetdisttable");

    // Assuming histograms and tree data correspond in some way
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing the TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY, MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;
    int MCDjetNConst;

    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fJetNConst",&MCDjetNConst);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfMlScore0",&MCDhfMlScore0); // background ML score
    tree->SetBranchAddress("fHfMlScore1",&MCDhfMlScore1); // prompt D0 ML score
    tree->SetBranchAddress("fHfMlScore2",&MCDhfMlScore2); // non-prompt D0 ML score

    int nEntries = tree->GetEntries(); //  * 2 / 10
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        double MCDDeltaR = MCDaxisDistance;
        // Building selections
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, customBdtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb);
        bool recoAngularRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut);
        bool recoJetPtRange = ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax));
        bool recoHfPtRange = ((MCDhfPt >= binning.ptHFBinEdges_detector[0]) && (MCDhfPt < binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1]));
        bool recoDeltaRRange = ((MCDDeltaR >= binning.deltaRBinEdges_detector[0]) && (MCDDeltaR < DeltaRcut));
        bool recoLevelRange = passBDTcut && recoAngularRange && recoJetPtRange && recoHfPtRange && recoDeltaRRange;

        // Fill each histogram with their respective pT intervals
        if (recoLevelRange) {
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < binning.ptHFBinEdges_detector.size() - 1 && !filled; iEdge++) {
                if ((MCDhfPt >= binning.ptHFBinEdges_detector[iEdge]) && (MCDhfPt < binning.ptHFBinEdges_detector[iEdge + 1])) {

                    // use Emma's run 2 cuts in case useEmmaYeatsBins is true
                    // if (!binning.useEmmaYeatsBins || passEmmaCut(jetPt, hfPt)) {
                    //     dataContainer.histograms2d[iEdge]->Fill(hfMass, deltaR);
                    // }
                    dataContainer.hDataTotal[iEdge].first->Fill(MCDhfMass, MCDDeltaR);
                    filled = true; // Exit the loop once the correct histogram is found (alternative: break)
                }
                
            }
        }
    } // end of entries loop

    // Create 1D projections
    for (size_t iEdge = 0; iEdge < binning.ptHFBinEdges_detector.size() - 1; iEdge++) {
        dataContainer.hDataTotal[iEdge].second = dataContainer.hDataTotal[iEdge].first->ProjectionX(Form("%s_projX",dataContainer.hDataTotal[iEdge].first->GetName()));
    }

    std::cout << "Data histograms filled." << std::endl;
}

void performMCFits(const BinningStruct& binning, BDTScanContainer& dataContainer) {
    
    // --- Perform individual fits ---
    for (size_t iInterval = 0; iInterval < binning.ptHFBinEdges_detector.size() - 1; iInterval++) {
        double minMass = binning.minMass[iInterval];
        double maxMass = binning.maxMass[iInterval];

        // --- Pure signal fits ---
        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        dataContainer.fittings.fMCPureSignalsFit.emplace_back(new TF1(Form("fMCPureSignalsFit_%zu", iInterval), fitWrapper<PureSignalTemplateModel>, minMass, maxMass, 5));
        // Initialize parameters from the histogram distribution characteristics
        double C1, C2, m0, sigma1, sigma2;
        C1 = 0.7 * dataContainer.hMCSignals[iInterval].second->GetMaximum(); // Assume first Gaussian dominates
        C2 = 0.3 * dataContainer.hMCSignals[iInterval].second->GetMaximum(); // Second Gaussian smaller
        m0 = dataContainer.hMCSignals[iInterval].second->GetMean();       // Mean of the distribution
        sigma1 = dataContainer.hMCSignals[iInterval].second->GetRMS() / 2; // First Gaussian narrower
        sigma2 = dataContainer.hMCSignals[iInterval].second->GetRMS();    // Second Gaussian wider
        dataContainer.fittings.fMCPureSignalsFit[iInterval]->SetParameters(C1, C2, m0, sigma1, sigma2);
        dataContainer.fittings.fMCPureSignalsFit[iInterval]->SetParName(0, "C1");
        dataContainer.fittings.fMCPureSignalsFit[iInterval]->SetParName(1, "C2");
        dataContainer.fittings.fMCPureSignalsFit[iInterval]->SetParName(2, "m0");
        dataContainer.fittings.fMCPureSignalsFit[iInterval]->SetParName(3, "sigma1");
        dataContainer.fittings.fMCPureSignalsFit[iInterval]->SetParName(4, "sigma2");
        dataContainer.fittings.fMCPureSignalsFit[iInterval]->SetParLimits(0, 0., TMath::Infinity());
        dataContainer.fittings.fMCPureSignalsFit[iInterval]->SetParLimits(1, 0., TMath::Infinity());
        dataContainer.hMCSignals[iInterval].second->Fit(dataContainer.fittings.fMCPureSignalsFit[iInterval], "RQN"); // "Q" option performs quiet fit without drawing the fit function
        // Also create and store individual components
        TF1* fSignalPrimary = new TF1(Form("signalPrimaryFit_%zu", iInterval), fitWrapper<PureSingleGaussianTemplateModel>, minMass, maxMass, 3);
        fSignalPrimary->FixParameter(0, dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(0)); // Fix amplitude to the value from the combined fit
        fSignalPrimary->FixParameter(1, dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(2)); // Fix mean to the value from the combined fit
        fSignalPrimary->FixParameter(2, dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(3)); // Fix sigma to the value from the combined fit
        fSignalPrimary->SetParName(0, "C1_primary");
        fSignalPrimary->SetParName(1, "m0_primary");
        fSignalPrimary->SetParName(2, "sigma_primary");
        TF1* fSignalSecondary = new TF1(Form("signalSecondaryFit_%zu", iInterval), fitWrapper<PureSingleGaussianTemplateModel>, minMass, maxMass, 3);
        fSignalSecondary->FixParameter(0, dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(1)); // Fix amplitude to the value from the combined fit
        fSignalSecondary->FixParameter(1, dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(2)); // Fix mean to the value from the combined fit
        fSignalSecondary->FixParameter(2, dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(4)); // Fix sigma to the value from the combined fit
        fSignalSecondary->SetParName(0, "C2_secondary");
        fSignalSecondary->SetParName(1, "m0_secondary");
        fSignalSecondary->SetParName(2, "sigma_secondary");
        dataContainer.fittings.fIndividualMCPureSignalsFit.emplace_back(std::make_pair(fSignalPrimary,fSignalSecondary));

        // --- Pure reflection fits ---
        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        dataContainer.fittings.fMCPureReflectionsFit.emplace_back(new TF1(Form("fMCPureReflectionsFit_%zu", iInterval), fitWrapper<PureReflectionsTemplateModel>, minMass, maxMass, 11));
        // Initialize parameters from the histogram distribution characteristics
        double m0_1, m0_2;
        C1 = 0.6 * dataContainer.hMCReflections[iInterval].second->GetMaximum();   // Dominant reflection, default to 60% of the histogram maximum
        C2 = 0.4 * dataContainer.hMCReflections[iInterval].second->GetMaximum();   // Secondary reflection, default to 40% of the histogram maximum
        m0_1 = dataContainer.hMCReflections[iInterval].second->GetMean() - dataContainer.hMCReflections[iInterval].second->GetRMS() / 2; // Offset from mean
        m0_2 = dataContainer.hMCReflections[iInterval].second->GetMean() + dataContainer.hMCReflections[iInterval].second->GetRMS() / 2; // Offset from mean
        sigma1 = dataContainer.hMCReflections[iInterval].second->GetRMS() / 2;  // Narrower
        sigma2 = dataContainer.hMCReflections[iInterval].second->GetRMS();      // Wider
        dataContainer.fittings.fMCPureReflectionsFit[iInterval]->SetParameters(0., 0., 0., 0., 0., C1, C2, m0_1, m0_2, sigma1, sigma2);
        dataContainer.fittings.fMCPureReflectionsFit[iInterval]->SetParName(5, "C1");
        dataContainer.fittings.fMCPureReflectionsFit[iInterval]->SetParName(6, "C2");
        dataContainer.fittings.fMCPureReflectionsFit[iInterval]->SetParName(7, "m0_1");
        dataContainer.fittings.fMCPureReflectionsFit[iInterval]->SetParName(8, "m0_2");
        dataContainer.fittings.fMCPureReflectionsFit[iInterval]->SetParName(9, "sigma1");
        dataContainer.fittings.fMCPureReflectionsFit[iInterval]->SetParName(10, "sigma2");
        dataContainer.fittings.fMCPureReflectionsFit[iInterval]->SetParLimits(5, 0., TMath::Infinity());
        dataContainer.fittings.fMCPureReflectionsFit[iInterval]->SetParLimits(6, 0., TMath::Infinity()); // amplitudes
        dataContainer.hMCReflections[iInterval].second->Fit(dataContainer.fittings.fMCPureReflectionsFit[iInterval], "RQN");
        // Also create and store individual components
        TF1* fReflectionPrimary = new TF1(Form("reflectionPrimaryFit_%zu", iInterval), fitWrapper<PureSingleGaussianTemplateModel>, minMass, maxMass, 3);
        fReflectionPrimary->FixParameter(0, dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(5)); // Fix amplitude to the value from the combined fit
        fReflectionPrimary->FixParameter(1, dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(7)); // Fix mean to the value from the combined fit
        fReflectionPrimary->FixParameter(2, dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(9)); // Fix sigma to the value from the combined fit
        fReflectionPrimary->SetParName(0, "C1_primary");
        fReflectionPrimary->SetParName(1, "m0_primary");
        fReflectionPrimary->SetParName(2, "sigma_primary");
        TF1* fReflectionSecondary = new TF1(Form("reflectionSecondaryFit_%zu", iInterval), fitWrapper<PureSingleGaussianTemplateModel>, minMass, maxMass, 3);
        fReflectionSecondary->FixParameter(0, dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(6)); // Fix amplitude to the value from the combined fit
        fReflectionSecondary->FixParameter(1, dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(8)); // Fix mean to the value from the combined fit
        fReflectionSecondary->FixParameter(2, dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(10)); // Fix sigma to the value from the combined fit
        fReflectionSecondary->SetParName(0, "C2_secondary");
        fReflectionSecondary->SetParName(1, "m0_secondary");
        fReflectionSecondary->SetParName(2, "sigma_secondary");
        dataContainer.fittings.fIndividualMCPureReflectionsFit.emplace_back(std::make_pair(fReflectionPrimary,fReflectionSecondary));

        // ---Combined fits ---
        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        dataContainer.fittings.fMCSignalAndReflectionsFit.emplace_back(new TF1(Form("fMCSignalAndReflectionsFit_%zu", iInterval), fitWrapper<SigRefTemplateModel>, minMass, maxMass, 11));
        // Initialize parameters from previous fits
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(0,dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(0));
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(1,dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(1));
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(2,dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(2));
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(3,dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(3));
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(4,dataContainer.fittings.fMCPureSignalsFit[iInterval]->GetParameter(4));
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(5,dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(5));
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(6,dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(6));
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(7,dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(7));
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(8,dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(8));
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(9,dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(9));
        dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval]->FixParameter(10,dataContainer.fittings.fMCPureReflectionsFit[iInterval]->GetParameter(10));
        dataContainer.hMCSignalAndReflections[iInterval].second->Fit(dataContainer.fittings.fMCSignalAndReflectionsFit[iInterval], "RQN");
    }

    std::cout << "MC template fits performed." << std::endl;
}

void performDataFits(const BinningStruct& binning, BDTScanContainer& dataContainer, const FitModelType& modelToUse) {
    //
    double m_0_reference = 1.86484; // D0 mass in GeV/c^2
    double sigma_reference = 0.012;
    const double jetptMin = binning.ptjetBinEdges_detector[0];
    const double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

    // --- Total fits: loop through pT,HF intervals/histograms and perform fits
    for (size_t iHisto = 0; iHisto < dataContainer.hDataTotal.size(); ++iHisto) {
        double minMass = binning.minMass[iHisto];
        double maxMass = binning.maxMass[iHisto];

        // Get TF1 objects from MC file
        TF1* fSignal = dataContainer.fittings.fMCPureSignalsFit[iHisto];
        TF1* fReflections = dataContainer.fittings.fMCPureReflectionsFit[iHisto];

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

        // -> signal constrains
        double A1toA2MCSignalRatio = A1Signal / A2Signal;
        double Sigma1toSigma2MCSignalRatio = sigma1Signal / sigma2Signal;

        // -> reflection constrains
        double A1toA2MCReflectionsRatio = A1Reflection / A2Reflection;
        double Sigma1toSigma2MCReflectionsRatio = sigma1Reflections / sigma2Reflections;
        double A1SignalToA1ReflectionMCRatio = A1Signal / A1Reflection;

        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        if (modelToUse == FitModelType::FullPowerLaw) { // 12 parameters = 8 fixed + 5 free
            dataContainer.fittings.fitTotal.emplace_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<FullModelPowerLaw>, minMass, maxMass, 13));
            // Apply positive amplitude constraint to power law background only
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(0, 0., TMath::Infinity()); // only accepts non-negative background fits
            // dataContainer.fittings.fitTotal[iHisto]->SetParLimits(1, -1e10, 0.); // only descending inclination background fits
        } else if (modelToUse == FitModelType::FullPoly2) { // 13 parameters = 8 fixed + 6 free
            dataContainer.fittings.fitTotal.emplace_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<FullModelPoly2>, minMass, maxMass, 14));
            dataContainer.fittings.fitTotal[iHisto]->SetParName(13, "Background C (m^2 term)");  // Free
        }
        // Set initial values and fix parameters (original background amplitude guess = 200000)
        double guessBackgroundA = dataContainer.hDataTotal[iHisto].second->GetBinContent(dataContainer.hDataTotal[iHisto].second->FindBin(m_0_reference)) / 5; // initial guess for background amplitude at the signal peak
        double guessBackgroundB = -10.0; // initial guess for background slope (power law exponent or linear term)
        double guessA1Signal = dataContainer.hDataTotal[iHisto].second->GetBinContent(dataContainer.hDataTotal[iHisto].second->FindBin(m_0_reference)) / 2; // initial guess for A1Signal amplitude from MC fit
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
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.35 * sigma_reference); // sigma_1 = 0.0172
        } else if (iHisto == 5) {   // 6 < pT,D < 7 GeV/c
            sigma_reference = 0.0185;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0185
        } else if (iHisto == 6) {   // 7 < pT,D < 8 GeV/c
            sigma_reference = 0.0199;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.35 * sigma_reference); // sigma_1 = 0.0199
        } else if (iHisto == 7) {   // 8 < pT,D < 12 GeV/c
            sigma_reference = 0.0223;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.35 * sigma_reference); // sigma_1 = 0.0223
        } else if (iHisto == 8) {   // 12 < pT,D < 16 GeV/c
            sigma_reference = 0.0240;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.3 * sigma_reference, 1.5 * sigma_reference); // sigma_1 = 0.0240
        } else if (iHisto == 9) {   // 16 < pT,D < 24 GeV/c
            sigma_reference = 0.0240;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.5 * sigma_reference, 1.25 * sigma_reference); // sigma_1 = 0.0240
        } else if (iHisto == 10) {  // 24 < pT,D < 36 GeV/c
            sigma_reference = 0.0373;
            dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.2 * sigma_reference, 1.35 * sigma_reference); // sigma_1 = 0.0373
        }
        // sigma_reference = 0.012;
        // dataContainer.fittings.fitTotal[iHisto]->SetParLimits(5, 0.35 * sigma_reference, 2.0 * sigma_reference);

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
        dataContainer.hDataTotal[iHisto].second->Fit(dataContainer.fittings.fitTotal[iHisto], "RQN");

        // --- Background only fits ---
        if (modelToUse == FitModelType::FullPowerLaw || modelToUse == FitModelType::StandardSideBand) { // backgroundFunctionPowerLaw
            // Getting total parameter values (f(x) = a * x^b)
            double a_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(0); // Get the value of parameter 'a'
            double b_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(1); // Get the value of parameter 'b'
            dataContainer.fittings.fitBackgroundOnly.push_back(new TF1(Form("backgroundOnlyFit_%zu", iHisto), fitWrapper<PurePowerLawModel>, minMass, maxMass, 13));
            dataContainer.fittings.fitBackgroundOnly[iHisto]->FixParameter(0, a_par);
            dataContainer.fittings.fitBackgroundOnly[iHisto]->FixParameter(1, b_par);
        } else if (modelToUse == FitModelType::FullPoly2) { // backgroundFunctionPoly2
            // Getting total parameter values (f(x) = a*x² + b*x + c)
            double a_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(0); // Get the value of parameter 'a'
            double b_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(1); // Get the value of parameter 'b'
            double c_par = dataContainer.fittings.fitTotal[iHisto]->GetParameter(13); // Get the value of parameter 'c'
            dataContainer.fittings.fitBackgroundOnly.push_back(new TF1(Form("backgroundOnlyFit_%zu", iHisto), fitWrapper<PurePoly2Model>, minMass, maxMass, 14));
            dataContainer.fittings.fitBackgroundOnly[iHisto]->FixParameter(0, a_par);
            dataContainer.fittings.fitBackgroundOnly[iHisto]->FixParameter(1, b_par);
            dataContainer.fittings.fitBackgroundOnly[iHisto]->FixParameter(13, c_par);
        }
        dataContainer.fittings.fitBackgroundOnly[iHisto]->SetLineStyle(kDashed);
        dataContainer.fittings.fitBackgroundOnly[iHisto]->SetLineColor(kRed+2);

        // --- Signal only fits ---
        // Extract signal-related parameters from the total fit
        A1Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(2); // A1 signal
        A1toA2MCSignalRatio = dataContainer.fittings.fitTotal[iHisto]->GetParameter(3); // A2 signal
        m0Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(4); // Signal mean m0
        double sigmaSignal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(5); // Signal width sigma
        Sigma1toSigma2MCSignalRatio = dataContainer.fittings.fitTotal[iHisto]->GetParameter(6); // Sigma1/Sigma2 signal
        
        // Perfoming fit with acquired parameters from total fit
        dataContainer.fittings.fitSignalOnly.push_back(new TF1(Form("signalOnlyFit_%zu", iHisto), fitWrapper<PureSignalModel>, minMass, maxMass, 13));
        dataContainer.fittings.fitSignalOnly[iHisto]->FixParameter(2, A1Signal);
        dataContainer.fittings.fitSignalOnly[iHisto]->FixParameter(3, A1toA2MCSignalRatio);
        dataContainer.fittings.fitSignalOnly[iHisto]->FixParameter(4, m0Signal);
        dataContainer.fittings.fitSignalOnly[iHisto]->FixParameter(5, sigmaSignal);
        dataContainer.fittings.fitSignalOnly[iHisto]->FixParameter(6, Sigma1toSigma2MCSignalRatio);
        dataContainer.fittings.fitSignalOnly[iHisto]->SetLineStyle(kDashed);
        dataContainer.fittings.fitSignalOnly[iHisto]->SetLineColor(kBlue+2);

        // --- Reflections only fits ---
        // Extract reflection-related parameters from the total fit
        double A1SignalToA1ReflectionMCRatios = dataContainer.fittings.fitTotal[iHisto]->GetParameter(7); // A1 signal / A1 reflection
        A1toA2MCReflectionsRatio = dataContainer.fittings.fitTotal[iHisto]->GetParameter(8); // A1 reflection / A2 reflection
        double m0_1Reflection = dataContainer.fittings.fitTotal[iHisto]->GetParameter(9); // Reflection mean m0_1
        double m0_2Reflection = dataContainer.fittings.fitTotal[iHisto]->GetParameter(10); // Reflection mean m0_2
        double sigma1Reflection = dataContainer.fittings.fitTotal[iHisto]->GetParameter(11); // Reflection sigma_1
        double sigma2Reflection = dataContainer.fittings.fitTotal[iHisto]->GetParameter(12); // Reflection sigma_2
        
        // Perfoming fit with acquired parameters from total fit
        dataContainer.fittings.fitReflectionsOnly.push_back(new TF1(Form("reflectionOnlyFit_%zu", iHisto), fitWrapper<PureReflectionsModel>, minMass, maxMass, 13));
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(7, A1SignalToA1ReflectionMCRatios);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(2, A1Signal);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(8, A1toA2MCReflectionsRatio);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(9, m0_1Reflection);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(10, m0_2Reflection);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(11, sigma1Reflection);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->FixParameter(12, sigma2Reflection);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->SetLineStyle(kDashed);
        dataContainer.fittings.fitReflectionsOnly[iHisto]->SetLineColor(kGreen+3);
    }
    std::cout << "Data fits performed." << std::endl;
}

void calculateSignificance(const BinningStruct& binning, BDTScanContainer& dataContainer, BDTScanResult& bdtResults, const double& signalSigmas, const int& startingBackSigma, const int& backgroundSigmas, const int& iCoarseScore) {
    
    // Mark failed fits
    for (size_t iHisto = 0; iHisto < dataContainer.hDataTotal.size(); iHisto++) {

        // Check if signal gaussian acomodates the literature D0 rest mass
        double m_0_reference = 1.86484; // D0 mass in GeV/c^2
        if (didFitFailed(dataContainer.hDataTotal[iHisto].second, dataContainer.fittings.fitTotal[iHisto], m_0_reference, signalSigmas, startingBackSigma)) {            
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
    
    // Loop over pT,HF bins
    for (size_t iHisto = 0; iHisto < dataContainer.hDataTotal.size(); iHisto++) {
        // Calculate significance
        double sig = 0., S = 0., B = 0.;

        // Only compute for working fits that have both total and background fits
        bool hasSignalFit = (iHisto < dataContainer.fittings.fitSignalOnly.size()) && dataContainer.fittings.fitSignalOnly[iHisto];
        bool hasBkgFit = (iHisto < dataContainer.fittings.fitBackgroundOnly.size()) && dataContainer.fittings.fitBackgroundOnly[iHisto];
        bool isWorking = (iHisto < dataContainer.fittings.workingFits.size()) && dataContainer.fittings.workingFits[iHisto];

        if (isWorking && hasSignalFit && hasBkgFit) {
            TF1* fSignal = dataContainer.fittings.fitSignalOnly[iHisto];
            TF1* fBkg = dataContainer.fittings.fitBackgroundOnly[iHisto];
            double m0 = fSignal->GetParameter(4);
            double sigma = fSignal->GetParameter(5);
            std::cout << "m0 = " << m0 << "\t sigma = " << sigma << std::endl;

            // Calculate areas
            double totalInSR = fSignal->Integral(m0 - signalSigmas*sigma, m0 + signalSigmas*sigma);
            double bkgInSR = fBkg->Integral(m0 - signalSigmas*sigma, m0 + signalSigmas*sigma);
            std::cout << "totalInSR = " << totalInSR << "\t bkgInSR = " << bkgInSR << std::endl;
            S = totalInSR;
            B = bkgInSR;
            if (S > 0. && (S + B) > 0.) {
                sig = S / std::sqrt(S + B);
            }
            std::cout << "  [Significance] pT,D bin " << iHisto << ": S=" << S << "  B=" << B << "  sig=" << sig << "\n";

            // Create shaded side-band histograms for visualization
            double binW  = dataContainer.hDataTotal[iHisto].second->GetBinWidth(1);
            double yMax  = dataContainer.hDataTotal[iHisto].second->GetMaximum();

            // Signal region shading (blue)
            double srLow  = m0 - signalSigmas * sigma;
            double srHigh = m0 + signalSigmas * sigma;
            int srLowBin  = safeLowBin(dataContainer.hDataTotal[iHisto].second, srLow);
            int srHighBin = safeHighBin(dataContainer.hDataTotal[iHisto].second, srHigh);
            TH1D* hSigShade = (TH1D*)dataContainer.hDataTotal[iHisto].second->Clone(Form("hSigShade_%zu_%d", iHisto, iCoarseScore));
            hSigShade->Reset();
            for (int iBin = srLowBin; iBin <= srHighBin; ++iBin) {
                hSigShade->SetBinContent(iBin, dataContainer.hDataTotal[iHisto].second->GetBinContent(iBin));
            }
            hSigShade->SetFillColorAlpha(kBlue, 0.2);
            hSigShade->SetLineColor(kBlue);
            dataContainer.signalHistDrawing.emplace_back(hSigShade);

            // Sideband region shading (red)
            // Left sideband: [m0 - sbEnd*sigma, m0 - startingBackSigma*sigma]
            // Right sideband: [m0 + startingBackSigma*sigma, m0 + sbEnd*sigma]
            // Clamp to histogram range
            double sbLow1  = std::max(m0 - (startingBackSigma+backgroundSigmas)*sigma, dataContainer.hDataTotal[iHisto].second->GetXaxis()->GetXmin());
            double sbHigh1 = m0 - startingBackSigma * sigma;
            double sbLow2  = m0 + startingBackSigma * sigma;
            double sbHigh2 = std::min(m0 + (startingBackSigma+backgroundSigmas)*sigma, dataContainer.hDataTotal[iHisto].second->GetXaxis()->GetXmax());

            TH1D* hSBShade = (TH1D*)dataContainer.hDataTotal[iHisto].second->Clone(Form("hSBShade_%zu_%d", iHisto, iCoarseScore));
            hSBShade->Reset();
            int sb1Low  = safeLowBin(dataContainer.hDataTotal[iHisto].second, sbLow1);
            int sb1High = safeHighBin(dataContainer.hDataTotal[iHisto].second, sbHigh1);
            int sb2Low  = safeLowBin(dataContainer.hDataTotal[iHisto].second, sbLow2);
            int sb2High = safeHighBin(dataContainer.hDataTotal[iHisto].second, sbHigh2);
            for (int iBin = sb1Low; iBin <= sb1High; ++iBin) {
                hSBShade->SetBinContent(iBin, dataContainer.hDataTotal[iHisto].second->GetBinContent(iBin));
            }
            for (int iBin = sb2Low; iBin <= sb2High; ++iBin) {
                hSBShade->SetBinContent(iBin, dataContainer.hDataTotal[iHisto].second->GetBinContent(iBin));
            }
            hSBShade->SetFillColorAlpha(kRed, 0.2);
            hSBShade->SetLineColor(kRed);
            dataContainer.sidebandHistDrawing.emplace_back(hSBShade);
        } else {
            std::cout << "  [Significance] pT,D bin " << iHisto << ": skipped (fit failed or missing)\n";
            dataContainer.signalHistDrawing.emplace_back(nullptr);   // ← always one entry per bin
            dataContainer.sidebandHistDrawing.emplace_back(nullptr); // ← always one entry per bin
        }

        dataContainer.fittings.significance.push_back(sig);
        dataContainer.fittings.signalYield.push_back(S);
        dataContainer.fittings.backgroundYield.push_back(B);

        
    }

    bdtResults.significance.emplace_back(dataContainer.fittings.significance);
    std::cout << "Calculated significance." << std::endl;
}

void calculateRun2PromptEfficiency(TFile* fInputMC, const BinningStruct& binning, BDTScanResult& bdtResults, 
                                   const std::vector<std::pair<double, double>>& customBdtCuts, const int& iCoarseScore) {
    std::cout << "Calculating efficiency..." << std::endl;
    // Create pT spectra histograms
    bdtResults.hMcpPt.emplace_back(new TH1D(Form("mcp_pt_%d",iCoarseScore), "Denominator (particle level);#it{p}_{T,D^{0}}^{truth};dN", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data()));
    bdtResults.hMcpPt.back()->Sumw2();
    bdtResults.hMcpPt.back()->SetMarkerColor(kBlue+2);
    bdtResults.hMcpPt.back()->SetLineColor(kBlue+2);
    bdtResults.hMcpPt.back()->SetMarkerStyle(kFullCircle);
    bdtResults.hMcdPt.emplace_back(new TH1D(Form("mcd_pt_%d",iCoarseScore), "Numerator (detector level);#it{p}_{T,D^{0}}^{reco};dN", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data()));
    bdtResults.hMcdPt.back()->Sumw2();
    bdtResults.hMcdPt.back()->SetMarkerColor(kGreen+2);
    bdtResults.hMcdPt.back()->SetLineColor(kGreen+2);
    bdtResults.hMcdPt.back()->SetMarkerStyle(kFullCircle);

    // Fill histograms
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = MCPetaCut;
    const double MCPyCut = 0.8; // on detector level D0
    const double MCDyCut = MCPyCut;
    const double MCPDeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double MCDDeltaRcut = MCPDeltaRcut;
    double hfptMin = binning.ptHFBinEdges_detector[0]; //ptHFBinEdges[0] - should start from 0 or from the lowest pT,HF value?
    double hfptMax = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1];
    const double hfPtMincut = hfptMin; // on particle level HF
    const double hfPtMaxcut = hfptMax; // on particle level HF
    const double jetptMin = binning.ptjetBinEdges_detector[0];
    const double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

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

        // Generator level selection cuts
        bool genJetPtRange = (MCPjetPt >= jetptMin); //  && (MCPjetPt < jetptMax), remove upper bound for particle level
        // bool genHfPtRange = ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)); // remove entirely for particle level
        bool genDeltaRRange = ((MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < MCPDeltaRcut));
        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && genJetPtRange && genDeltaRRange; // && genHfPtRange

        // Reconstruction level selection cuts
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, customBdtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb);
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = ((MCDhfPt >= hfPtMincut) && (MCDhfPt < hfPtMaxcut));
        bool recoDeltaRRange = ((MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < MCDDeltaRcut));
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && recoJetPtRange && recoHfPtRange && recoDeltaRRange; // separated passBDTcut for checking
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isRealD0 = isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs);

        if (genLevelRange) {
            // Matched entries must be of real D0s, not reflections or combinatorial background
            if (MCDhfmatch) {
                if (isRealD0) {
                    if (MCPhfprompt) {
                        bdtResults.hMcpPt.back()->Fill(MCPhfPt); // [iCoarseScore]
                    } else{
                        // dataContainer.hMcpPt[2]->Fill(MCPhfPt);// fill non-prompt efficiency histogram
                    }
                    // Fill detector level entry
                    if (recoLevelRange && passBDTcut) {
                        // fill prompt efficiency histogram
                        if (MCDhfprompt) {
                            bdtResults.hMcdPt.back()->Fill(MCDhfPt); // [iCoarseScore]
                        } else{
                            // fill non-prompt efficiency histogram
                            //dataContainer.hMcdPt[2]->Fill(MCDhfPt);
                        }
                    } else if (recoLevelRange && !passBDTcut) {
                        //dataContainer.hBDTBackgroundScore->Fill(MCDhfMlScore0);
                    }
                }
            } else {// fill particle level even if not matched
                if (MCPhfprompt) {
                    bdtResults.hMcpPt.back()->Fill(MCPhfPt); // [iCoarseScore]
                } else{
                    // dataContainer.hMcpPt[2]->Fill(MCPhfPt);// fill non-prompt efficiency histogram
                }
            }   
        }
    }
    std::cout << "Dividing hMcdPt / hMcpPt..." << std::endl;
    // Calculate efficiency
    bdtResults.efficiencyHist.emplace_back((TH1D*) bdtResults.hMcdPt.back()->Clone(Form("hRun2PromptEfficiency_score%d", iCoarseScore)));
    bdtResults.efficiencyHist.back()->SetTitle("Run 2 style prompt efficiency");
    bdtResults.efficiencyHist.back()->Divide(bdtResults.hMcpPt.back());
    std::cout << "Efficiency distribution calculated." << std::endl;
}

void plotIteration(const BDTScanContainer& dataContainer, const BDTScanResult& bdtResults, const BinningStruct& binning,
                   const std::vector<std::pair<double, double>>& customBdtCuts, const int& scanIteration, const TString& granularity = "") {
    std::cout << "Plotting iteration..." << std::endl;
    const double jetptMin = binning.ptjetBinEdges_detector[0];
    const double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

    // Defining canvases before plotting
    int nHistos = dataContainer.hMCSignals.size();
    // Start with a square layout (or close to it)
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));

    // --- MC templates ---
    TCanvas* cMCPureReflections = new TCanvas("cMCPureReflections","MC pure reflections template", 1800,1000);
    TCanvas* cMCPureSignals = new TCanvas("cMCPureSignals","MC pure signals template", 1800,1000);
    TCanvas* cMCSignalAndReflections = new TCanvas("cMCSignalAndReflections","MC signal and reflections template", 1800,1000);
    TCanvas* cDataTotal = new TCanvas("cDataTotal","Data fits", 1800,1000);
    cMCPureSignals->Divide(nCols,nRows);
    cMCPureReflections->Divide(nCols,nRows);
    cMCSignalAndReflections->Divide(nCols,nRows);
    cDataTotal->Divide(nCols,nRows);

    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.04); // default = 0.05
    for (size_t iHisto = 0; iHisto < nHistos; iHisto++) {

        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box

        // --- MC templates ---
        cMCPureReflections->cd(iHisto+1);//h1->GetYaxis()->SetRangeUser(0.0, 50.0); // Sets Y-axis range
        dataContainer.hMCReflections[iHisto].second->GetYaxis()->SetRangeUser(0.0, dataContainer.hMCReflections[iHisto].second->GetMaximum() * 1.2);
        dataContainer.hMCReflections[iHisto].second->Draw();
        dataContainer.fittings.fMCPureReflectionsFit[iHisto]->SetLineColor(kGreen+2);
        dataContainer.fittings.fMCPureReflectionsFit[iHisto]->SetLineStyle(kSolid);
        dataContainer.fittings.fMCPureReflectionsFit[iHisto]->Draw("same");
        dataContainer.fittings.fIndividualMCPureReflectionsFit[iHisto].first->SetLineColor(kGreen+2);
        dataContainer.fittings.fIndividualMCPureReflectionsFit[iHisto].first->SetLineStyle(kDashed);
        dataContainer.fittings.fIndividualMCPureReflectionsFit[iHisto].first->Draw("same");
        dataContainer.fittings.fIndividualMCPureReflectionsFit[iHisto].second->SetLineColor(kGreen+2);
        dataContainer.fittings.fIndividualMCPureReflectionsFit[iHisto].second->SetLineStyle(kDashed);
        dataContainer.fittings.fIndividualMCPureReflectionsFit[iHisto].second->Draw("same");

        cMCPureSignals->cd(iHisto+1);
        dataContainer.hMCSignals[iHisto].second->GetYaxis()->SetRangeUser(0.0, dataContainer.hMCSignals[iHisto].second->GetMaximum() * 1.2);
        dataContainer.hMCSignals[iHisto].second->Draw();
        dataContainer.fittings.fMCPureSignalsFit[iHisto]->SetLineColor(kBlue+2);
        dataContainer.fittings.fMCPureSignalsFit[iHisto]->SetLineStyle(kSolid);
        dataContainer.fittings.fMCPureSignalsFit[iHisto]->Draw("same");
        dataContainer.fittings.fIndividualMCPureSignalsFit[iHisto].first->SetLineColor(kBlue+2);
        dataContainer.fittings.fIndividualMCPureSignalsFit[iHisto].first->SetLineStyle(kDashed);
        dataContainer.fittings.fIndividualMCPureSignalsFit[iHisto].first->Draw("same");
        dataContainer.fittings.fIndividualMCPureSignalsFit[iHisto].second->SetLineColor(kBlue);
        dataContainer.fittings.fIndividualMCPureSignalsFit[iHisto].second->SetLineStyle(kDashed);
        dataContainer.fittings.fIndividualMCPureSignalsFit[iHisto].second->Draw("same");

        cMCSignalAndReflections->cd(iHisto+1);
        dataContainer.hMCSignalAndReflections[iHisto].second->GetYaxis()->SetRangeUser(0.0, dataContainer.hMCSignalAndReflections[iHisto].second->GetMaximum() * 1.2);
        dataContainer.hMCSignalAndReflections[iHisto].second->Draw();
        dataContainer.fittings.fMCSignalAndReflectionsFit[iHisto]->SetLineColor(kBlack);
        dataContainer.fittings.fMCSignalAndReflectionsFit[iHisto]->SetLineStyle(kSolid);
        dataContainer.fittings.fMCPureReflectionsFit[iHisto]->Draw("same");
        dataContainer.fittings.fMCPureSignalsFit[iHisto]->Draw("same");
        dataContainer.fittings.fMCSignalAndReflectionsFit[iHisto]->Draw("same");
        TLegend* legMCTemplate = new TLegend(0.6,0.57,0.85,0.77);
        legMCTemplate->AddEntry(dataContainer.fittings.fMCPureSignalsFit[iHisto],"Pure signal", "l");
        legMCTemplate->AddEntry(dataContainer.fittings.fMCPureReflectionsFit[iHisto],"Reflections", "l");
        legMCTemplate->AddEntry(dataContainer.fittings.fMCSignalAndReflectionsFit[iHisto],"Total fit", "l");
        legMCTemplate->Draw();

        // Data
        cDataTotal->cd(iHisto+1);
        dataContainer.hDataTotal[iHisto].second->GetYaxis()->SetRangeUser(0.0, dataContainer.hDataTotal[iHisto].second->GetMaximum() * 1.2);
        dataContainer.hDataTotal[iHisto].second->Draw();
        // Only draw shaded regions if they exist for this histogram
        if (iHisto < dataContainer.signalHistDrawing.size() && dataContainer.signalHistDrawing[iHisto]) {
            dataContainer.signalHistDrawing[iHisto]->Draw("same");
        }
        if (iHisto < dataContainer.sidebandHistDrawing.size() && dataContainer.sidebandHistDrawing[iHisto]) {
            dataContainer.sidebandHistDrawing[iHisto]->Draw("same");
        }
        // Only draw fits if this histogram has a working fit
        bool hasWorkingFit = (iHisto < dataContainer.fittings.workingFits.size())
                        && dataContainer.fittings.workingFits[iHisto];
        if (hasWorkingFit) {
            dataContainer.fittings.fitTotal[iHisto]->SetLineColor(kBlack);
            dataContainer.fittings.fitTotal[iHisto]->Draw("same");
            dataContainer.fittings.fitBackgroundOnly[iHisto]->Draw("same");
            dataContainer.fittings.fitSignalOnly[iHisto]->Draw("same");
            dataContainer.fittings.fitReflectionsOnly[iHisto]->Draw("same");
            // ... annotations ...
            double sigma1Signal = dataContainer.fittings.fitTotal[iHisto]->GetParameter(5); // Get the value of parameter 'sigma' from primary signal gaussian
            double sigma2Signal = sigma1Signal / dataContainer.fittings.fitTotal[iHisto]->GetParameter(6); // Get the value of parameter 'sigma' from secondary signal gaussian
            double chi2 = dataContainer.fittings.fitTotal[iHisto]->GetChisquare();
            double degOfFreedom = dataContainer.fittings.fitTotal[iHisto]->GetNDF();
            // latex->DrawLatex(statBoxPos-0.87, 0.23, "Data: LHC23_pass4_Thin_2P3PDstar_D0CJ_4_D0_1"); // JE_HF_LHC23_pass4_Thin_2P3PDstar_D0CJ_4_D0_1
            // latex->DrawLatex(statBoxPos-0.87, 0.28, "MC: LHC24h1c_All_D0"); // HF_LHC24h1c_All_D0
            latex->DrawLatex(statBoxPos-0.27, 0.34, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax)); // Display jet pT cut applied
            latex->SetTextColor(kBlack);
            latex->DrawLatex(0.65, 0.85, Form("#Chi^{2}_{red} = %.3f", degOfFreedom > 0 ? chi2/degOfFreedom : -1.));
            latex->DrawLatex(0.65, 0.8, Form("#sigma_{primary} = %.4f", sigma1Signal));
            latex->DrawLatex(0.65, 0.75, Form("#sigma_{secondary} = %.4f", sigma2Signal));
            double cut = (iHisto < (int)customBdtCuts.size()) ? customBdtCuts[iHisto].second : 0.;
            double sig = (iHisto < (int)dataContainer.fittings.significance.size()) ? dataContainer.fittings.significance[iHisto] : 0.;
            double S = (iHisto < (int)dataContainer.fittings.signalYield.size()) ? dataContainer.fittings.signalYield[iHisto] : 0.;
            double B = (iHisto < (int)dataContainer.fittings.backgroundYield.size()) ? dataContainer.fittings.backgroundYield[iHisto] : 0.;
            latex->SetTextColor(kBlue);
            latex->DrawLatex(0.65, 0.69, Form("BDT cut = %.4f", cut));
            latex->DrawLatex(0.65, 0.64, Form("S = %.0f, B = %.0f", S, B));
            latex->DrawLatex(0.65, 0.59, Form("Sig = %.2f", sig));
        } else {
            // Mark failed fit clearly
            TLatex* texFail = new TLatex();
            texFail->SetNDC();
            texFail->SetTextColor(kRed);
            texFail->SetTextSize(0.05);
            texFail->DrawLatex(0.3, 0.5, "Fit failed");
        }
    }

    // Run 2 style efficiency
    TCanvas* cRun2PromptEfficiency = new TCanvas("cRun2PromptEfficiency","Run 2 style prompt D0s efficiency with custom BDT cuts",1800,1000);
    cRun2PromptEfficiency->Divide(2,2);
    cRun2PromptEfficiency->cd(1);
    bdtResults.hMcpPt[scanIteration-1]->Draw();
    cRun2PromptEfficiency->cd(2);
    bdtResults.hMcdPt[scanIteration-1]->Draw();
    cRun2PromptEfficiency->cd(3);
    bdtResults.efficiencyHist[scanIteration-1]->Draw();
    double statBoxPos = gPad->GetUxmax(); // Height of the stat box
    latex->SetTextColor(kBlack);
    latex->DrawLatex(statBoxPos-0.87, 0.34, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax)); // Display jet pT cut applied

    //
    // Storing images
    //
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString sIteration = Form("%s_%d", granularity.Data(), scanIteration);

    TString imagePath = "../../Images/1-SignalTreatment/BDTOptimization/" + binning.dataPeriod + "/";
    TString jetPtRange = Form("%.0f_to_%.0fGeV", jetptMin, jetptMax);
    cMCPureReflections->Print(Form(imagePath + sIteration + "_%s.pdf(",jetPtRange.Data()));
    cMCPureSignals->Print(Form(imagePath + sIteration + "_%s.pdf",jetPtRange.Data()));
    cMCSignalAndReflections->Print(Form(imagePath + sIteration + "_%s.pdf",jetPtRange.Data()));
    cDataTotal->Print(Form(imagePath + sIteration + "_%s.pdf",jetPtRange.Data()));
    cRun2PromptEfficiency->Print(Form(imagePath + sIteration + "_%s.pdf)",jetPtRange.Data()));

    std::cout << "Iteration plots stored." << std::endl;

}

void plotSignificance(const BDTScanContainer& dataContainer, const BDTScanResult& bdtResults,  const BinningStruct& binning) {
    std::cout << "Plotting significance..." << std::endl;
    const double jetptMin = binning.ptjetBinEdges_detector[0];
    const double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

    const int nPtHFBins = dataContainer.hDataTotal.size();
    const int nScans   = bdtResults.bdtThresholds.size();

    // One pad per pT,D bin
    int nCols = std::min(4, nPtHFBins);
    int nRows = std::ceil((double)nPtHFBins / nCols);
    TCanvas* cOptimalBdtCuts = new TCanvas("cBDTScan", "BDT Significance Scan", 400*nCols, 350*nRows);
    cOptimalBdtCuts->Divide(nCols, nRows);
    // Plot efficiency as a function of the BDT score for each pT,HF bin
    TCanvas* cEfficiencyVsBDT = new TCanvas("cEfficiencyVsBDT","Run 2 style prompt efficiency as a function of the BDT score",1800,1000);
    cEfficiencyVsBDT->Divide(nCols, nRows);

    for (int iPtD = 0; iPtD < nPtHFBins; ++iPtD) {
        cOptimalBdtCuts->cd(iPtD + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

        // Build TGraph: x = BDT threshold, y = significance
        TGraph* graph = new TGraph(nScans);
        graph->SetName(Form("gBDTScan_ptD%d", iPtD));
        graph->SetTitle(Form("%.0f < p_{T,D^{0}} < %.0f GeV/c;" "BDT background score threshold;" "S/#sqrt{S+B}", binning.ptHFBinEdges_detector[iPtD], binning.ptHFBinEdges_detector[iPtD+1]));

        for (int iScan = 0; iScan < nScans; ++iScan) {
            double threshold = bdtResults.bdtThresholds[iScan][iPtD+1];// not to use 0-1 GeV/c interval
            double sig = (iPtD < (int)bdtResults.significance[iScan].size()) ? bdtResults.significance[iScan][iPtD] : 0.;
            graph->SetPoint(iScan, threshold, sig);
        }

        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.2);
        graph->SetMarkerColor(kBlue+1);
        graph->SetLineColor(kBlue+1);
        graph->SetLineWidth(1);
        graph->Draw("APL");
        gPad->SetLogx(1);

        // Mark the optimal cut with a vertical line
        double optCut = (iPtD < (int)bdtResults.optimalCuts.size()) ? bdtResults.optimalCuts[iPtD] : 0.;
        double maxSig = (iPtD < (int)bdtResults.maxSignificance.size()) ? bdtResults.maxSignificance[iPtD] : 0.;
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
        int nFine = bdtResults.bdtThresholdsFine.size();
        if (nFine > 0) {
            TGraph* gFine = new TGraph(nFine);
            for (int j = 0; j < nFine; ++j) {
                gFine->SetPoint(j,
                    bdtResults.bdtThresholdsFine[j][iPtD],
                    bdtResults.significanceFine[j][iPtD]);
            }
            gFine->SetMarkerStyle(21);
            gFine->SetMarkerColor(kRed+1);
            gFine->SetLineColor(kRed+1);
            gFine->Draw("P same"); // points only, no connecting line
            gPad->SetLogx(1);
        }

        // Build TGraph: x = BDT threshold, y = prompt efficiency
        TGraph* gEffVsBdt = new TGraph(nScans);
        gEffVsBdt->SetName(Form("gEffBDTScan_ptD%d", iPtD));
        gEffVsBdt->SetTitle(Form("%.0f < p_{T,D^{0}} < %.0f GeV/c;" "Efficiency as a function of BDT background score threshold;" "Efficiency#times Acceptance", binning.ptHFBinEdges_detector[iPtD], binning.ptHFBinEdges_detector[iPtD+1]));

        for (int iScan = 0; iScan < nScans; ++iScan) {
            double threshold = bdtResults.bdtThresholds[iScan][iPtD+1];
            int HFBin = bdtResults.efficiencyHist[iScan]->FindBin(binning.ptHFBinEdges_detector[iPtD]);
            double efficiency = bdtResults.efficiencyHist[iScan]->GetBinContent(HFBin);
            gEffVsBdt->SetPoint(iScan, threshold, efficiency);
        }
        cEfficiencyVsBDT->cd(iPtD+1);
        gEffVsBdt->Draw();
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
    TString imagePath = "../../Images/1-SignalTreatment/BDTOptimization/" + binning.dataPeriod + "/";
    TString jetPtRange = Form("%.0f_to_%.0fGeV", jetptMin, jetptMax);
    cOptimalBdtCuts->Print(Form(imagePath + "Significance_%s.pdf(",jetPtRange.Data()));
    cEfficiencyVsBDT->Print(Form(imagePath + "Significance_%s.pdf)",jetPtRange.Data()));

    std::cout << "Significance plot stored." << std::endl;
}

void ScanIterator(const double jetptMin, const double jetptMax, BinningStruct& binning, const FitModelType& modelToUse) {

    // Opening files
    TFile* fInputMC = new TFile("../../" + binning.inputMC.second + "/AO2D_mergedDFs.root","read");
    if (!fInputMC || fInputMC->IsZombie()) {
        std::cerr << "Error: Unable to open the input ROOT file." << std::endl;
    }
    TFile* fInputData = new TFile("../../" + binning.inputDATA.second + "/AO2D_mergedDFs.root","read");
    if (!fInputData || fInputData->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT data file." << std::endl;
    }
    
    // --- Configuration ---
    const int nPtHFBins = binning.ptHFBinEdges_detector.size() - 1;
    const int nCoarseScanSteps = 20;  // coarse scan first: 10 steps
    const int nFineScan = 5;  // fine scan: 20 steps around the best region
    double signalSigmas = 2; // default delta = 2
    int startingBackSigma = 4; // default position = 4
    int backgroundSigmas = 4; // default delta = 4

    // Create data container
    BDTScanContainer dataContainer;
    BDTScanResult bdtResults;
    // Loop over coarse scan
    for (size_t iCoarseScore = 1; iCoarseScore <= nCoarseScanSteps; iCoarseScore++) {
        std::cout << "============================================= " << iCoarseScore << "th coarse scan" << " =============================================" << std::endl;

        // Clear all containers before this iteration
        dataContainer.ClearData();

        // Create custom iteration bdt cuts objects
        std::vector<std::pair<double, double>> customBdtCuts = createCustomBdtCuts(bdtResults, binning, iCoarseScore, nCoarseScanSteps);

        // --- Obtain MC templates
        // Create
        dataContainer.hMCReflections = createHistograms(binning, Form("MCreflections_coarseScan%zu", iCoarseScore));
        dataContainer.hMCSignals = createHistograms(binning, Form("MCsignal_coarseScan%zu", iCoarseScore));
        dataContainer.hMCSignalAndReflections = createHistograms(binning, Form("MCsignal_and_reflections_coarseScan%zu", iCoarseScore));

        // Fill histograms
        fillMCHistograms(fInputMC, binning, dataContainer, customBdtCuts);

        // Perform fits
        performMCFits(binning, dataContainer);

        // --- Side-band subtraction
        // Create
        dataContainer.hDataTotal = createHistograms(binning, Form("DataTotal_coarseScan%zu", iCoarseScore));

        // Fill histograms
        fillDataHistograms(fInputData, binning, dataContainer, customBdtCuts);

        // Perform fits
        performDataFits(binning, dataContainer, modelToUse);

        // Calculate significance and store
        calculateSignificance(binning, dataContainer, bdtResults, signalSigmas, startingBackSigma, backgroundSigmas, iCoarseScore);

        // Calculate run 2 style prompt efficiency
        calculateRun2PromptEfficiency(fInputMC, binning, bdtResults, customBdtCuts, iCoarseScore);

        // Plot all histograms with respective fits
        plotIteration(dataContainer, bdtResults, binning, customBdtCuts, iCoarseScore, "CoarseScan");
    }
    std::cout << "Finding best scores" << std::endl;
    // Find the best scores
    bdtResults.optimalCuts = std::vector<double>(nPtHFBins, 0.);
    bdtResults.maxSignificance = std::vector<double>(nPtHFBins, 0.);
    std::vector<int> bestCoarseStep(nPtHFBins, 0);

    int nCompletedScans = (int)bdtResults.significance.size();
    std::cout << "[BDT] Completed " << nCompletedScans << " coarse scans.\n";
    std::cout << "[BDT] bdtResults.significance[0].size() = " 
            << (!bdtResults.significance.empty() ? bdtResults.significance[0].size() : 0) << "\n";

    for (int iPtD = 0; iPtD < nPtHFBins; ++iPtD) {
        // Start with scan 0 as best
        bestCoarseStep[iPtD] = 0;
        if (!bdtResults.significance.empty() && iPtD < (int)bdtResults.significance[0].size()) {
            bdtResults.maxSignificance[iPtD] = bdtResults.significance[0][iPtD];
            bdtResults.optimalCuts[iPtD]     = bdtResults.bdtThresholds[0][iPtD+1];
        }

        // Compare against all other completed scans
        for (int iScan = 1; iScan < nCompletedScans; ++iScan) {
            if (iPtD >= (int)bdtResults.significance[iScan].size()) {
                std::cout << "Warning: significance[" << iScan << "] has only "
                        << bdtResults.significance[iScan].size() << " entries, skipping iPtD=" << iPtD << "\n";
                continue;
            }
            double sig = bdtResults.significance[iScan][iPtD];
            if (sig > bdtResults.maxSignificance[iPtD]) {
                bdtResults.maxSignificance[iPtD] = sig;
                bestCoarseStep[iPtD]             = iScan;
                bdtResults.optimalCuts[iPtD]     = bdtResults.bdtThresholds[iScan][iPtD+1];
            }
        }
        std::cout << "[Coarse best] pT,D bin " << iPtD
                << " [" << binning.ptHFBinEdges_detector[iPtD]
                << "-"  << binning.ptHFBinEdges_detector[iPtD+1] << "]"
                << "  best step=" << bestCoarseStep[iPtD]
                << "  cut="       << bdtResults.optimalCuts[iPtD]
                << "  sig="       << bdtResults.maxSignificance[iPtD] << "\n";
    }
    
    // Loop over fine scan (over the range [best-1poit;best+1point])
    for (int iFine = 0; iFine <= nFineScan; ++iFine) {
        std::cout << "=== Fine scan step " << iFine+1 << "/" << nFineScan+1 << " ===\n";

        dataContainer.ClearData();
        
        // Build per-bin fine cuts
        std::vector<std::pair<double,double>> customBdtCuts = binning.bdtPtCuts;
        std::vector<double> thresholdsThisStep(nPtHFBins);

        for (int iPtD = 0; iPtD < nPtHFBins; ++iPtD) {
            
            // Fine range: one coarse step below and above the best
            double fineMin = (bestCoarseStep[iPtD] > 0) ? bdtResults.bdtThresholds[bestCoarseStep[iPtD]-1][iPtD+1] : 0.5 * bdtResults.bdtThresholds[0][iPtD+1]; // avoid 0
            //double fineMax = (bestCoarseStep[iPtD] < nCoarseScanSteps-1) ? bdtResults.bdtThresholds[bestCoarseStep[iPtD]+1][iPtD] : binning.bdtPtCuts[iPtD+1].second; // use next bin's max as upper bound
            double fineMax = (bestCoarseStep[iPtD] < nCompletedScans - 1) ? bdtResults.bdtThresholds[bestCoarseStep[iPtD]+1][iPtD+1] : bdtResults.bdtThresholds[nCompletedScans-1][iPtD+1];
            
            // Exclude endpoints (don't repeat coarse scan points)
            double cut = fineMin + (iFine+1) * (fineMax - fineMin) / (nFineScan+2);
            customBdtCuts[iPtD+1].second = cut; // +1 because bdtPtCuts[0] = {0, 0} pre-range
            
            thresholdsThisStep[iPtD] = cut;
        }

        bdtResults.bdtThresholdsFine.push_back(thresholdsThisStep);
        std::cout << "MC templates" << std::endl;
        // MC templates
        dataContainer.hMCReflections = createHistograms(binning, Form("MCrefl_f%d",    iFine));
        dataContainer.hMCSignals = createHistograms(binning, Form("MCsig_f%d",     iFine));
        dataContainer.hMCSignalAndReflections = createHistograms(binning, Form("MCsigref_f%d", iFine));
        fillMCHistograms(fInputMC, binning, dataContainer, customBdtCuts);
        performMCFits(binning, dataContainer);
        std::cout << "Data histograms and fits" << std::endl;
        // Data
        dataContainer.hDataTotal = createHistograms(binning, Form("Data_f%d", iFine));
        fillDataHistograms(fInputData, binning, dataContainer, customBdtCuts);
        performDataFits(binning, dataContainer, modelToUse);
        calculateSignificance(binning, dataContainer, bdtResults, signalSigmas, startingBackSigma, backgroundSigmas, iFine+1);
        
        // Calculate run 2 style prompt efficiency
        calculateRun2PromptEfficiency(fInputMC, binning, bdtResults, customBdtCuts, nCoarseScanSteps+1+iFine);

        // Store fine significance for plotting
        std::vector<double> sigThisStep(nPtHFBins, 0.);
        for (int iPtD = 0; iPtD < nPtHFBins; ++iPtD) {
            
            sigThisStep[iPtD] = (iPtD < (int)dataContainer.fittings.significance.size()) ? dataContainer.fittings.significance[iPtD] : 0.;
            
        }
        bdtResults.significanceFine.push_back(sigThisStep);

        // Update optimal cut if fine step gives better significance
        for (int iPtD = 0; iPtD < nPtHFBins; ++iPtD) {
            
            double sig = sigThisStep[iPtD];
            double cut = thresholdsThisStep[iPtD];
            std::cout << "  pT,D bin " << iPtD << "  cut=" << cut << "  sig=" << sig << "\n";
            if (sig > bdtResults.maxSignificance[iPtD]) {
                
                bdtResults.maxSignificance[iPtD] = sig;
                
                bdtResults.optimalCuts[iPtD] = cut;
                
            }
        }

        plotIteration(dataContainer, bdtResults, binning, customBdtCuts, nCoarseScanSteps+iFine, "FineScan");
    }

    // Plot significance as a function of score threshold
    plotSignificance(dataContainer, bdtResults, binning);

    // Store the optimal values to the bdt cuts main object
    for (size_t iPtHF = 1; iPtHF < binning.ptHFBinEdges_detector.size(); iPtHF++) {
        binning.bdtPtCuts[iPtHF].second = bdtResults.optimalCuts[iPtHF-1];
        std::cout << binning.bdtPtCuts[iPtHF].first << "-" << binning.bdtPtCuts[iPtHF+1].first << "\t" << binning.bdtPtCuts[iPtHF].second << std::endl;
    }
    
}

void BDTOptimization() {
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    TString dataPeriod = "2023";
    TString sBinningFileName = "binningInfo_" + dataPeriod + ".root";
    if (!gSystem->AccessPathName(sBinningFileName)) {
        std::cout << "Optimal cuts already calculated. Skipping step..." << std::endl;
        return;
    }
    TFile* fBinning = new TFile(sBinningFileName, "recreate");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the binning info ROOT file." << std::endl;
    }

    // --- Binning objects
    BinningStruct binning;

    // Emma Yeats reported analysis bins:
    binning.useEmmaYeatsBins = false;
    if (binning.useEmmaYeatsBins) {
        // pT,jet cuts
        binning.ptjetBinEdges_detector = {5., 7., 10., 20., 50.};
        // DeltaR bins
        binning.deltaRBinEdges_detector = {0., 0.01, 0.03, 0.05, 0.12, 0.2};
    } else {
        // pT,jet cuts
        binning.ptjetBinEdges_detector = {5., 7., 10., 16., 36., 50.}; // default = {5., 7., 10., 16., 36., 50.}
        // DeltaR bins
        binning.deltaRBinEdges_detector = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5};
    }
    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.05 (from first of pair)
    binning.bdtPtCuts = {
        {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.6}, {8, 0.8}, {10, 0.8}, {12, 1}, {16, 1}, {50, -1}
    }; // on dataset JE_HF_LHC23_pass4_Thin_2P3PDstar_D0CJ_4_D0_1, std::vector<std::pair<double, double>>
    efficiencyBinEdges(binning); // get ptHFEfficiencyBinEdges from the BDT intervals
    binning.ptjetBinEdges_particle = binning.ptjetBinEdges_detector;
    binning.deltaRBinEdges_particle = binning.deltaRBinEdges_detector;
    // the HF binning depends on the ones used in the BDT models, copy binning.ptHFEfficiencyBinEdges_particle minus last edge
    binning.ptHFBinEdges_particle = std::vector<double>(binning.ptHFEfficiencyBinEdges_particle.begin(), binning.ptHFEfficiencyBinEdges_particle.end() - 1);
    binning.ptHFBinEdges_particle.push_back(36.); // add last edge to cover the entire range, as done on Nima's AN, defalut = 36 GeV/c
    //binning.ptHFBinEdges_particle = {1., 2., 3., 4., 5., 6., 7., 8., 12., 16., 24., 36.}; // override with custom binning for HF pT, same as used for BDT models
    binning.ptHFBinEdges_detector = binning.ptHFBinEdges_particle;
    // Define mass range here — this is the ONLY place you need to change it
    binning.massBinDensity = 50.0 / (2.06 - 1.72); // = 147.06 bins/GeV/c²
    // Default: all bins use 2.06 GeV/c²
    binning.minMass.assign(binning.ptHFBinEdges_detector.size() - 1, 1.72);
    binning.maxMass.assign(binning.ptHFBinEdges_detector.size() - 1, 2.06);
    binning.maxMass = {2.02, 2.045, 2.06, 2.08, 2.1, 2.14, 2.20, 2.28, 2.47, 2.47, 2.47};

    // Choose file period
    binning.dataPeriod = dataPeriod;
    if (binning.dataPeriod == "2023") {
        // DATA
        binning.inputDATA.first = "JE_HF_LHC23_pass4_Thin_2P3PDstar_D0CJ_4_D0_1";
        binning.inputDATA.second = "Data/Experimental/Train_643652";
        // Anchored MC
        binning.inputMC.first = "HF_LHC24h1c_All_D0";
        binning.inputMC.second = "Data/MonteCarlo/Train_671273";
    } else if (binning.dataPeriod == "2022") {
        // DATA
        binning.inputDATA.first = "JE_HF_LHC22o_pass7_minBias_2P3PDstar_D0CJ_4_D0_1";
        binning.inputDATA.second = "Data/Experimental/Train_659513";
        // Anchored MC
        binning.inputMC.first = "HF_LHC24g5_All_D0";
        binning.inputMC.second = "Data/MonteCarlo/Train_669231";
    }

    // Choose total fit model
    FitModelType modelToUse = FitModelType::FullPowerLaw;
    // Compute the entire range of pT,jet
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    ScanIterator(jetptMin, jetptMax, binning, modelToUse);

    // Store binning information in a separate file for later use
    storeBinningInFile(fBinning, binning);
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

int main() {
    BDTOptimization();
    return 0;
}



