/*
 * .h file containing common utilities and functions used throughout the analysis
 * 
 * 
 * 
 * Author: Christian Reckziegel
**/


using namespace std;

enum DecayChannelMain : int8_t {
  // D0
  D0ToPiK = 1,      // π+ K−
  D0ToPiKPi0 = 2,   // π+ K− π0
  D0ToPiPi = 3,     // π+ π−
  D0ToPiPiPi0 = 4,  // π+ π− π0
  D0ToKK = 5,       // K+ K−
  // J/ψ
  JpsiToEE = 6,     // e+ e−
  JpsiToMuMu = 7,   // μ+ μ−
  //
  NChannelsMain = JpsiToMuMu // last channel
};
// Selected as...:D0 = +1, D0bar = -1, neither = 0
enum class D0Species : int {
    D0BAR = -1,
    NEITHER = 0,
    D0 = +1
};
// Conversion function
D0Species intToD0Species(int value) {
    switch(value) {
        case -1: return D0Species::D0BAR;
        case 0:  return D0Species::NEITHER;
        case 1:  return D0Species::D0;
        default: return D0Species::NEITHER; // Handle unexpected values
    }
}
bool isReflection(const int& MCDhfMatchedFrom, const int& MCDhfSelectedAs) {
    bool isReflectedD0 = (MCDhfMatchedFrom == static_cast<int>(DecayChannelMain::D0ToPiK)) && (MCDhfSelectedAs == static_cast<int>(D0Species::D0BAR));
    bool isReflectedD0bar = (MCDhfMatchedFrom == -static_cast<int>(DecayChannelMain::D0ToPiK)) && (MCDhfSelectedAs == static_cast<int>(D0Species::D0));
    bool isReflection = isReflectedD0 || isReflectedD0bar;
    
    return isReflection;
}
bool isTrueSignal(const int& MCDhfMatchedFrom, const int& MCDhfSelectedAs) {
    bool isD0 = (MCDhfMatchedFrom == static_cast<int>(DecayChannelMain::D0ToPiK)) && (MCDhfSelectedAs == static_cast<int>(D0Species::D0));
    bool isD0bar = (MCDhfMatchedFrom == -static_cast<int>(DecayChannelMain::D0ToPiK)) && (MCDhfSelectedAs == static_cast<int>(D0Species::D0BAR));
    bool isTrueSignal = isD0 || isD0bar;

    return isTrueSignal;
}
bool isD0ToKPiPi0(const int& MCDhfMatchedFrom, const int& MCDhfSelectedAs) {
    bool cameFromD0ToPiKPi0 = (MCDhfMatchedFrom == static_cast<int>(DecayChannelMain::D0ToPiKPi0));

    return cameFromD0ToPiKPi0;
}
// Calculate delta phi such that 0 < delta phi < 2*pi
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
        if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) { // if pT \in [ptmin; ptmax]
            return bdtPtCuts[i].second;
        }
    }
    return 1.0; // Default: accept all if out of range
}
/*___________________________________________________Binning information (start)______________________________________________________________________________________*/
// Binning object creation, storage and retrieval
struct BinningStruct {
    // Detector level binning
    std::vector<double> ptjetBinEdges_detector;
    std::vector<double> deltaRBinEdges_detector;
    std::vector<double> ptHFBinEdges_detector;
    
    // Particle level binning
    std::vector<double> ptjetBinEdges_particle;
    std::vector<double> deltaRBinEdges_particle;
    std::vector<double> ptHFBinEdges_particle;

    // BDT cut values for each pT,HF bin
    std::vector<std::pair<double, double>> bdtPtCuts;
    std::vector<double> ptHfBfDTBinEdges;
    // Efficiency binning (comes from the BDT cuts binning)
    std::vector<double> ptHFEfficiencyBinEdges_particle;
    std::vector<double> ptHFEfficiencyBinEdges_detector;

    // Flag to indicate whether Emma Yeats binning is used
    bool useEmmaYeatsBins = false;

    // Mass histogram configuration
    double massBinDensity = 50.0 / (2.06 - 1.72);       // bins per GeV/c² → keeps bin width constant
    std::vector<double> minMass;                        // global lower mass edge (same for all pT,D bins)
    std::vector<double> maxMass;                        // per pT,D upper mass edge

    TString dataPeriod;                                 // year
    std::pair<TString, TString> inputMC;                // first = name, second = path
    std::pair<TString, TString> inputDATA;              // first = name, second = path
};
void efficiencyBinEdges(BinningStruct& binning) {

    // Reserve space, "0-1" and "16-50" intervals are excluded
    binning.ptHFEfficiencyBinEdges_particle.reserve(binning.bdtPtCuts.size()-2);
    binning.ptHFEfficiencyBinEdges_detector.reserve(binning.bdtPtCuts.size()-2);

    // Loop through intervals  of pT,HF
    for (size_t iInterval = 1; iInterval < binning.bdtPtCuts.size()-1; iInterval++) {
        binning.ptHFEfficiencyBinEdges_particle.emplace_back(binning.bdtPtCuts[iInterval].first);
    }

    // Add padding bin for folding/unfolding operations
    if (!binning.useEmmaYeatsBins) {
        binning.ptHFEfficiencyBinEdges_particle.emplace_back(36.);
        binning.ptHFEfficiencyBinEdges_particle.emplace_back(50.);
    } else {
        binning.ptHFEfficiencyBinEdges_particle.emplace_back(20.);
        binning.ptHFEfficiencyBinEdges_particle.emplace_back(50.);
    }
    
    // Detector binning is the same as particle binning
    binning.ptHFEfficiencyBinEdges_detector = binning.ptHFEfficiencyBinEdges_particle;
}
std::vector<double> LoadBinning(TFile* fInput, const char* pathInFile) {
    auto* vec = (TVectorD*)fInput->Get(pathInFile);
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", pathInFile));
    }
    return std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
}
void storeBinningInFile(TFile* fOutput, const BinningStruct& binning) {
    
    // Create a directory for axes
    fOutput->mkdir("axes");
    fOutput->cd("axes");
    // Create TVectorD with same content
    // (i) detector level
    TVectorD vecDeltaR_detector(binning.deltaRBinEdges_detector.size());
    for (size_t i = 0; i < binning.deltaRBinEdges_detector.size(); ++i) {
        vecDeltaR_detector[i] = binning.deltaRBinEdges_detector[i];
    }
    vecDeltaR_detector.Write("deltaRBinEdges_detector");
    TVectorD vecPtJet_detector(binning.ptjetBinEdges_detector.size());
    for (size_t i = 0; i < binning.ptjetBinEdges_detector.size(); ++i) {
        vecPtJet_detector[i] = binning.ptjetBinEdges_detector[i];
    }
    vecPtJet_detector.Write("ptjetBinEdges_detector");
    TVectorD vecPtHF_detector(binning.ptHFBinEdges_detector.size());
    for (size_t i = 0; i < binning.ptHFBinEdges_detector.size(); ++i) {
        vecPtHF_detector[i] = binning.ptHFBinEdges_detector[i];
    }
    vecPtHF_detector.Write("ptHFBinEdges_detector");
    // (ii) particle level
    // Create TVectorD with same content
    TVectorD vecDeltaR_particle(binning.deltaRBinEdges_particle.size());
    for (size_t i = 0; i < binning.deltaRBinEdges_particle.size(); ++i) {
        vecDeltaR_particle[i] = binning.deltaRBinEdges_particle[i];
    }
    vecDeltaR_particle.Write("deltaRBinEdges_particle");
    TVectorD vecPtJet_particle(binning.ptjetBinEdges_particle.size());
    for (size_t i = 0; i < binning.ptjetBinEdges_particle.size(); ++i) {
        vecPtJet_particle[i] = binning.ptjetBinEdges_particle[i];
    }
    vecPtJet_particle.Write("ptjetBinEdges_particle");
    TVectorD vecPtHF_particle(binning.ptHFBinEdges_particle.size());
    for (size_t i = 0; i < binning.ptHFBinEdges_particle.size(); ++i) {
        vecPtHF_particle[i] = binning.ptHFBinEdges_particle[i];
    }
    vecPtHF_particle.Write("ptHFBinEdges_particle");
    for (size_t i = 0; i < binning.ptHFBinEdges_particle.size(); ++i) {
        vecPtHF_particle[i] = binning.ptHFBinEdges_particle[i];
    }
    vecPtHF_particle.Write("ptHFBinEdges_particle");
    // Return to root directory
    fOutput->cd();

    // Create a directory for BDT cuts
    fOutput->mkdir("bdt");
    fOutput->cd("bdt");
    TVectorD ptEdges(binning.bdtPtCuts.size());
    TVectorD bdtCuts(binning.bdtPtCuts.size());

    for (size_t i = 0; i < binning.bdtPtCuts.size(); ++i) {
        ptEdges[i] = binning.bdtPtCuts[i].first;
        bdtCuts[i] = binning.bdtPtCuts[i].second;
    }
    ptEdges.Write("ptHfBDTBinEdges");
    bdtCuts.Write("bdtCutValues");

    // Return to root directory
    fOutput->cd();
    TVectorD vecPtHFEff(binning.ptHFEfficiencyBinEdges_particle.size());
    for (size_t i = 0; i < binning.ptHFEfficiencyBinEdges_particle.size(); ++i) {
        vecPtHFEff[i] = binning.ptHFEfficiencyBinEdges_particle[i];
    }
    vecPtHFEff.Write("ptHFEfficiencyBinEdges");

    // Return to root directory
    fOutput->cd();
    TParameter<bool> useEmmaYeatsBins("useEmmaYeatsBins", binning.useEmmaYeatsBins);
    useEmmaYeatsBins.Write();

    // Return to root directory
    fOutput->cd();

    // Store mass range configuration
    fOutput->mkdir("massBins");
    fOutput->cd("massBins");
    

    TVectorD vMassBinDensity(1);
    vMassBinDensity[0] = binning.massBinDensity;
    vMassBinDensity.Write("massBinDensity");

    TVectorD vMinMass(binning.minMass.size());
    for (size_t i = 0; i < binning.minMass.size(); ++i) {
        vMinMass[i] = binning.minMass[i];
    }
    vMinMass.Write("minMass");

    TVectorD vMaxMass(binning.maxMass.size());
    for (size_t i = 0; i < binning.maxMass.size(); ++i) {
        vMaxMass[i] = binning.maxMass[i];
    }
    vMaxMass.Write("maxMass");

    // Return to root directory
    fOutput->cd();

    // Store data period used and local storage address
    fOutput->mkdir("DataPeriod");
    fOutput->cd("DataPeriod");
    TObjString(binning.dataPeriod).Write("dataPeriod");
    fOutput->mkdir("DataPeriod/DATA");
    fOutput->cd("DataPeriod/DATA");
    TObjString(binning.inputDATA.first).Write("name");
    TObjString(binning.inputDATA.second).Write("path");
    fOutput->cd("DataPeriod");
    fOutput->mkdir("DataPeriod/MC");
    fOutput->cd("DataPeriod/MC");
    TObjString(binning.inputMC.first).Write("name");
    TObjString(binning.inputMC.second).Write("path");

    // Return to root directory
    fOutput->cd();
}
BinningStruct retrieveBinningFromFile(TFile* fInput) {
    
    // Create struct to hold all binning information
    BinningStruct binning;

    //
    // Load detector level binning
    //
    TVectorD* vec = nullptr;
    vec = (TVectorD*)fInput->Get("axes/ptjetBinEdges_detector");
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", "axes/ptjetBinEdges_detector"));
    }
    binning.ptjetBinEdges_detector = std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
    vec = (TVectorD*)fInput->Get("axes/deltaRBinEdges_detector");
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", "axes/deltaRBinEdges_detector"));
    }
    binning.deltaRBinEdges_detector = std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
    vec = (TVectorD*)fInput->Get("axes/ptHFBinEdges_detector");
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", "axes/ptHFBinEdges_detector"));
    }
    binning.ptHFBinEdges_detector = std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
    
    //
    // Load particle level binning
    //
    vec = (TVectorD*)fInput->Get("axes/ptjetBinEdges_particle");
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", "axes/ptjetBinEdges_particle"));
    }
    binning.ptjetBinEdges_particle = std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
    vec = (TVectorD*)fInput->Get("axes/deltaRBinEdges_particle");
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", "axes/deltaRBinEdges_particle"));
    }
    binning.deltaRBinEdges_particle = std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
    vec = (TVectorD*)fInput->Get("axes/ptHFBinEdges_particle");
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", "axes/ptHFBinEdges_particle"));
    }
    binning.ptHFBinEdges_particle = std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());

    //
    // --- Load BDT cut values
    //
    fInput->cd(); // return to root directory before accessing bdt directory
    TVectorD* vecPt = (TVectorD*)fInput->Get("bdt/ptHfBDTBinEdges");
    TVectorD* vecCut = (TVectorD*)fInput->Get("bdt/bdtCutValues");
    if (!vecPt || !vecCut) {
        throw std::runtime_error("Could not find BDT cut vectors in file");
    }
    if (vecPt->GetNoElements() != vecCut->GetNoElements()) {
        throw std::runtime_error("BDT pT bins and cut vectors have different sizes");
    }
    
    binning.bdtPtCuts.clear();
    binning.ptHfBfDTBinEdges.clear();
    for (int i = 0; i < vecPt->GetNoElements(); ++i) {
        binning.bdtPtCuts.emplace_back((*vecPt)[i], (*vecCut)[i]);
        binning.ptHfBfDTBinEdges.emplace_back((*vecPt)[i]);
    }

    //
    // --- Load efficiency bin edges, which are the same as the BDT cut pT,HF bin edges, except for the last one which is the upper edge of the last BDT interval
    //
    TVectorD* vecPtHFEff = (TVectorD*)fInput->Get("ptHFEfficiencyBinEdges");
    for (int i = 0; i < vecPtHFEff->GetNoElements(); ++i) {
        binning.ptHFEfficiencyBinEdges_particle.emplace_back((*vecPtHFEff)[i]);
    }
    binning.ptHFEfficiencyBinEdges_detector = binning.ptHFEfficiencyBinEdges_particle;

    // Retrieve boolean flag
    TParameter<bool>* fileUseEmmaYeatsBins = (TParameter<bool>*)fInput->Get("useEmmaYeatsBins");
    if (fileUseEmmaYeatsBins) {
        binning.useEmmaYeatsBins = fileUseEmmaYeatsBins->GetVal();
    } else {
        std::cerr << "Warning: 'useEmmaYeatsBins' boolean flag not found in file, setting to false" << std::endl;
        binning.useEmmaYeatsBins = false;
    }

    //
    // --- Retrieve mass range configuration
    //
    TVectorD* vMassBinDensity = (TVectorD*)fInput->Get("massBins/massBinDensity");
    if (vMassBinDensity) binning.massBinDensity = (*vMassBinDensity)[0];

    TVectorD* vMinMass = (TVectorD*)fInput->Get("massBins/minMass");
    if (vMinMass) {
        binning.minMass.resize(vMinMass->GetNrows());
        for (int i = 0; i < vMinMass->GetNrows(); ++i) {
            binning.minMass[i] = (*vMinMass)[i];
        }
    }

    TVectorD* vMaxMass = (TVectorD*)fInput->Get("massBins/maxMass");
    if (vMaxMass) {
        binning.maxMass.resize(vMaxMass->GetNrows());
        for (int i = 0; i < vMaxMass->GetNrows(); ++i) {
            binning.maxMass[i] = (*vMaxMass)[i];
        }
    }

    //
    // --- Load data period and path
    //
    TObjString* objDataPeriod = (TObjString*)fInput->Get("DataPeriod/dataPeriod");
    binning.dataPeriod = objDataPeriod->GetString();

    TObjString *objName = (TObjString*)fInput->Get("DataPeriod/DATA/name");
    TObjString *objPath = (TObjString*)fInput->Get("DataPeriod/DATA/path");
    binning.inputDATA.first = objName->GetString();
    binning.inputDATA.second = objPath->GetString();

    TObjString *objMCName = (TObjString*)fInput->Get("DataPeriod/MC/name");
    TObjString *objMCPath = (TObjString*)fInput->Get("DataPeriod/MC/path");
    binning.inputMC.first = objMCName->GetString();
    binning.inputMC.second = objMCPath->GetString();

    return binning;
}
/*___________________________________________________Binning information (end)______________________________________________________________________________________________*/

// Check for discontinuities in the histogram
// ToDo: include options for
// - number of bins defining the discontinuity (e.g., 2 empty bins surrounded by 2 filled bins)
// - region of the histogram to check (e.g., only a number of sigmas around the mean, or the entire histogram)
bool checkHistoDiscontinuity(const TH1D* histogram) {
    for (int iBin = 1; iBin < histogram->GetNbinsX() - 1; ++iBin) {
        double countsCurrentBin = histogram->GetBinContent(iBin + 1);
        double countsNextBin = histogram->GetBinContent(iBin + 2);
        double countsNextNextBin = histogram->GetBinContent(iBin + 3);
        if (countsNextBin == 0 && countsCurrentBin > 0 && countsNextNextBin > 0) {
            return true; // Discontinuity found
        }
    }
    return false;
}
/*___________________________________________________Background subtraction (start)______________________________________________________________________________________*/
// Define the enum class with descriptive names
enum class FitModelType {
    StandardSideBand,               // Standard side-band subtraction method
    FullPowerLaw,                   // Signal + Reflections + power law Background,
    FullPoly2,                      // Signal + Reflections + 2nd order polynominal Background,
    FullSingleGaussian,             // Signal as single Gaussian + Reflections + power law background,
    SignalReflectionsOnly,          // Signal + Reflections only
    SignalOnly,                     // Only Signal
    ReflectionsOnly                 // Only Reflections
};
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
// Module struct to perform and store histograms fits
struct FitContainer {
    std::vector<TF1*> fitBackgroundOnly;        // background only fits
    std::vector<TF1*> fitSignalOnly;            // Pure signal only fits
    std::vector<TF1*> fitReflectionsOnly;       // Pure reflections only fits (obtained from constrains from MC templates)
    std::vector<TF1*> fitSignalsAndReflections; // Signal + Reflections only fits
    std::vector<TF1*> fitTotal;                 // Total fit: signal + reflections + background

    // list of working fits (PDG mass inside of signal range)
    std::vector<bool> workingFits;

    // Individual components (primary and secondary Gaussians)
    std::vector<std::pair<TF1*,TF1*>> individualSignals;
    std::vector<std::pair<TF1*,TF1*>> individualReflections;

    // Optimal BDT evaluation
    std::vector<double> significance;    // S/sqrt(S+B) per pT,D bin
    std::vector<double> signalYield;     // S per pT,D bin
    std::vector<double> backgroundYield; // B per pT,D bin
};
template <typename HistT>
int safeLowBin(HistT* histogram, double x) {
    int bin = histogram->GetXaxis()->FindFixBin(x);
    return std::max(bin, 1);
}

template <typename HistT>
int safeHighBin(HistT* histogram, double x) {
    int bin = histogram->GetXaxis()->FindFixBin(x);
    return std::min(bin, histogram->GetXaxis()->GetNbins());
}
// Calculate sidebands' region limits
std::pair<std::array<double, 2>, std::array<double, 2>> calculateSidebandRegions(const size_t iHisto, const FitContainer& fittings, const TH1D* hInvMass, const int startingBackSigma, const int backgroundSigmas) {
    // Define output arrays
    std::array<double, 2> leftRange = {0., 0.}; // remains zero if no sigmas fit inside the left range
    std::array<double, 2> rightRange = {0., 0.};
    
    // Getting fit parameters
    double m_0 = fittings.fitTotal[iHisto]->GetParameter(4); // Get the value of parameter signal m0 of primary gassian
    double signalSigma1 = fittings.fitTotal[iHisto]->GetParameter(5); // Get the value of parameter 'sigma1'
    double signalSigma2 = fittings.fitTotal[iHisto]->GetParameter(5) / fittings.fitTotal[iHisto]->GetParameter(6); // signalSigma2 = sigma1 / sigmaRatio12
    double sigma = signalSigma1;

    // Check which sides should be used for the total side-band distribution (if at least 1 sigma fit inside left range)
    bool useLeftSide = (m_0 - (startingBackSigma+1)*sigma) > hInvMass->GetBinLowEdge(1);
    bool useRightSide = (m_0 + (startingBackSigma+1)*sigma) < hInvMass->GetBinLowEdge(hInvMass->GetNbinsX()+1);
    std::cout << "Check if right side can be used: \t start point = " << m_0 + (startingBackSigma+1)*sigma << " ,\t Right histogram edge = " << hInvMass->GetBinLowEdge(hInvMass->GetNbinsX()+1) << std::endl;
    // if (!useLeftSide && useRightSide) {//useLeftSide && !useRightSide - previous version: m_0 - (startingBackSigma+1)*sigma <  hInvMass->GetBinLowEdge(1)
    //     std::cout << "Using only right sideband." << std::endl;

    //     // Count for how many sigmas is there room inside the left side range
    //     double leftSidebandRange = (m_0 - startingBackSigma * sigma) - hInvMass->GetBinLowEdge(1);
    //     int leftSigmas = static_cast<int>(leftSidebandRange / sigma);// get the integral number
    //     std::cout << leftSigmas << " sigmas fit inside the left side range for the background distribution estimation." << std::endl;

    //     // Count for how many sigmas is there room inside the right side range
    //     double rightSidebandRange = hInvMass->GetBinLowEdge(hInvMass->GetNbinsX()+1) - (m_0 + startingBackSigma*sigma);
    //     int rightSigmas = static_cast<int>(rightSidebandRange / sigma);// get the integral number
    //     std::cout << rightSigmas << " sigmas fit inside the right side range for the background distribution estimation." << std::endl;

    //     // Ensure only a maximum of 4 sigmas (backgroundSigmas) are used
    //     if (leftSigmas > backgroundSigmas) {
    //         leftSigmas = backgroundSigmas;
    //     }
    //     if (rightSigmas > backgroundSigmas) {
    //         rightSigmas = backgroundSigmas;
    //     }

    //     std::cout << "Left extreme = " << m_0 - (startingBackSigma + leftSigmas) * sigma << endl;
    //     std::cout << "Right extreme = " << m_0 + (startingBackSigma + rightSigmas) * sigma << endl;
        

    //     // Calculate right range limits
    //     rightRange[0] = m_0 + startingBackSigma * sigma;
    //     rightRange[1] = m_0 + (startingBackSigma + rightSigmas) * sigma;

    // } else { // else if (useLeftSide && useRightSide)
    //     std::cout << "Using both left and right sidebands.\n";

    //     // Count for how many sigmas is there room inside the left side range
    //     double leftSidebandRange = (m_0 - startingBackSigma * sigma) - hInvMass->GetBinLowEdge(1);
    //     int leftSigmas = static_cast<int>(leftSidebandRange / sigma);// get the integral number
    //     std::cout << leftSigmas << " sigmas fit inside the left side range for the background distribution estimation." << std::endl;

    //     // Count for how many sigmas is there room inside the right side range
    //     double rightSidebandRange = hInvMass->GetBinLowEdge(hInvMass->GetNbinsX()+1) - (m_0 + startingBackSigma*sigma);
    //     int rightSigmas = static_cast<int>(rightSidebandRange / sigma);// get the integral number
    //     std::cout << rightSigmas << " sigmas fit inside the right side range for the background distribution estimation." << std::endl;

    //     // Ensure only a maximum of 4 sigmas (backgroundSigmas) are used
    //     if (leftSigmas > backgroundSigmas) {
    //         leftSigmas = backgroundSigmas;
    //     }
    //     if (rightSigmas > backgroundSigmas) {
    //         rightSigmas = backgroundSigmas;
    //     }

    //     // Calculate left range limits
    //     leftRange[0] = m_0 - (startingBackSigma + leftSigmas) * sigma;
    //     leftRange[1] = m_0 - startingBackSigma * sigma;

    //     // Calculate right range limits
    //     rightRange[0] = m_0 + startingBackSigma * sigma;
    //     rightRange[1] = m_0 + (startingBackSigma + rightSigmas) * sigma;
        
    // }

    bool useIntegerSidebands = false;
    if (useIntegerSidebands) {
        // Count for how many sigmas is there room inside the left side range
        double leftSidebandRange = (m_0 - startingBackSigma * sigma) - hInvMass->GetBinLowEdge(1);
        int leftSigmas = static_cast<int>(leftSidebandRange / sigma);// get the integral number
        std::cout << leftSigmas << " sigmas fit inside the left side range for the background distribution estimation." << std::endl;

        // Count for how many sigmas is there room inside the right side range
        double rightSidebandRange = hInvMass->GetBinLowEdge(hInvMass->GetNbinsX()+1) - (m_0 + startingBackSigma*sigma);
        int rightSigmas = static_cast<int>(rightSidebandRange / sigma);// get the integral number
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
    } else {
        // Count for how many sigmas is there room inside the left side range
        double leftSidebandRange = (m_0 - startingBackSigma * sigma) - hInvMass->GetBinLowEdge(1);
        leftSidebandRange = std::max(leftSidebandRange, 0.0); // Ensure non-negative range
        double leftSigmas = leftSidebandRange / sigma;// get the fractional number
        std::cout << leftSigmas << " sigmas fit inside the left side range for the background distribution estimation." << std::endl;

        // Count for how many sigmas is there room inside the right side range
        double rightSidebandRange = hInvMass->GetBinLowEdge(hInvMass->GetNbinsX()+1) - (m_0 + startingBackSigma*sigma);
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
    }
    
    std::cout << "Left sideband: [" << leftRange[0] << ", " << leftRange[1] << "],\t " << abs(leftRange[1] - leftRange[0]) / sigma << " sigmas used." << std::endl;
    std::cout << "Right sideband: [" << rightRange[0] << ", " << rightRange[1] << "],\t " << abs(rightRange[0] - rightRange[1]) / sigma << " sigmas used." << std::endl;
    std::cout << "Sideband regions calculated.\n";

    return std::make_pair(leftRange, rightRange);
}
// Calculate scalling correction factors
std::array<double, 2> calculateScalingFactor(const size_t iHisto, const FitContainer& fittings,const std::array<double, 2> leftRange, const std::array<double, 2> rightRange, const double signalSigmas) {
    // Define output array
    std::array<double, 2> scallingFactors = {0., 0.};

    // Get total fit parameters
    double m_0 = fittings.fitTotal[iHisto]->GetParameter(4); // Get the value of parameter 'm_0'
    double signalSigma1 = fittings.fitTotal[iHisto]->GetParameter(5); // Get the value of parameter 'sigma1'
    double signalSigma2 = fittings.fitTotal[iHisto]->GetParameter(5) / fittings.fitTotal[iHisto]->GetParameter(6); // signalSigma2 = sigma1 / sigmaRatio12
    double sigma = signalSigma1;
    
    // Define area variables
    double B, Bs, R, Rs, S, Ys, Sb;

    if (leftRange[0] == 0. && leftRange[1] == 0.) {
        // Only right sideband used

        // Calculate background regions' areas
        B = fittings.fitBackgroundOnly[iHisto]->Integral(rightRange[0], rightRange[1]);
        Bs = fittings.fitBackgroundOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of background in signal region

        // Calculate reflection regions' areas
        R = fittings.fitReflectionsOnly[iHisto]->Integral(rightRange[0], rightRange[1]);
        Rs = fittings.fitReflectionsOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of reflections in signal region

        // Calculate signal regions' areas
        S = fittings.fitSignalOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of signal in signal region

        // Calculate total fit area in signal region
        Ys = fittings.fitTotal[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

        // Calculate signal contribution in the right side-band region
        Sb = fittings.fitSignalOnly[iHisto]->Integral(rightRange[0], rightRange[1]);
    } else {
        // Both sidebands used

        // Calculate background regions' areas
        B = fittings.fitBackgroundOnly[iHisto]->Integral(leftRange[0], leftRange[1]) + fittings.fitBackgroundOnly[iHisto]->Integral(rightRange[0], rightRange[1]);
        std::cout << "B_left = " << fittings.fitBackgroundOnly[iHisto]->Integral(leftRange[0], leftRange[1]) << "\t B_right = " << fittings.fitBackgroundOnly[iHisto]->Integral(rightRange[0], rightRange[1]) << std::endl;
        Bs = fittings.fitBackgroundOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of background in signal region

        // Calculate reflection regions' areas
        R = fittings.fitReflectionsOnly[iHisto]->Integral(leftRange[0], leftRange[1]) + fittings.fitReflectionsOnly[iHisto]->Integral(rightRange[0], rightRange[1]);
        std::cout << "R_left = " << fittings.fitReflectionsOnly[iHisto]->Integral(leftRange[0], leftRange[1]) << "\t R_right = " << fittings.fitReflectionsOnly[iHisto]->Integral(rightRange[0], rightRange[1]) << std::endl;
        Rs = fittings.fitReflectionsOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of reflections in signal region

        // Calculate signal regions' areas
        S = fittings.fitSignalOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of signal in signal region

        // Calculate total fit area in signal region
        Ys = fittings.fitTotal[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

        // Calculate signal contribution in both side-band regions
        Sb = fittings.fitSignalOnly[iHisto]->Integral(rightRange[0], rightRange[1]); // right
        Sb = Sb + fittings.fitSignalOnly[iHisto]->Integral(leftRange[0], leftRange[1]); // right
    }

    // Calculate scaling factor
    double alpha = S / (S + Rs - (Bs/B) * (R + Sb)); // ---> Addition of Sb is new!
    if (std::isnan(alpha)) {
        std::cout << "Error: Invalid scaling factor alpha = " << alpha << " for histogram index " << iHisto << ".\n";
        std::cout << "S = " << S << ", Rs = " << Rs << ", Bs = " << Bs << ", B = " << B << ", R = " << R << ", (S + Rs - (Bs/B) * R) = " << (S + Rs - (Bs/B) * R) << ".\n";
        alpha = 0.;
    }

    // Reflections correction
    scallingFactors[0] = alpha;
    // Background correction
    scallingFactors[1] = Bs / B;
    // std::cout << "alpha (signal + reflections correction) = " << scallingFactors[0] << std::endl;
    // std::cout << "Bs / B (background fit area ratios) = " << Bs << " / " << B << " = " << scallingFactors[1] << std::endl;
    std::cout << "fitYs / B (background fit area ratios, calculated from fits) = " << Ys << " / " << B << " = " << Ys / B << std::endl;
    std::cout << "S = " << S << ", Rs = " << Rs << ", Bs = " << Bs << ", B = " << B << ", R = " << R << std::endl;

    std::cout << "Scaling factors calculated.\n";

    return scallingFactors;
}
// Check if fit failed
// Criteria for fit failure could be:
// - PDG mass value outside of signal range (e.g., m_0 < mPDG - n*sigma GeV/c^2 or m_0 > mPDG + n*sigma GeV/c^2)
// - Fit did not converge (e.g., chi2/ndf > 5)
// - ToDo: if mean value is outside of the histogram range (e.g., m_0 < histogram min or m_0 > histogram max)
bool didFitFailed(TH1D* hInvMass, TF1* fitTotal, const double mPDG, const double signalSigmas, const int& startingBackSigma) {
    // The fit can be considered failed for many reasons
    
    double m_0 = fitTotal->GetParameter(4); // Get the value of parameter 'm_0'
    std::cout << "m_0 = " << m_0 << std::endl;
    double signalSigma1 = fitTotal->GetParameter(5); // Get the value of parameter 'sigma1'
    double signalSigma2 = fitTotal->GetParameter(5) / fitTotal->GetParameter(6); // signalSigma2 = sigma1 / sigmaRatio12
    double sigma = signalSigma1;
    
    // 1 - Check if m_0 is within the signal range defined by mPDG +/- signalSigmas * sigma
    if ((m_0 < mPDG - signalSigmas * sigma) || (m_0 > mPDG + signalSigmas * sigma)) {
        std::cout << "Fit failed (reason 1): m_0 = " << m_0 << " is outside of the signal range [" << mPDG - signalSigmas * sigma << ", " << mPDG + signalSigmas * sigma << "].\n";
        hInvMass->SetTitle(TString(hInvMass->GetTitle()) + " (Fit failed)");
        return true;
    }
    
    // 2 - Check if fit converged by looking at chi2/ndf    double chi2 = fitTotal->GetChisquare();
    int ndf = fitTotal->GetNDF();
    if (ndf > 0) {
        double chi2PerNdf = fitTotal->GetChisquare() / ndf;
        if (chi2PerNdf > 320.) {
            std::cout << "Fit failed (reason 2): chi2/ndf = " << chi2PerNdf << " is greater than 320.\n";
            hInvMass->SetTitle(TString(hInvMass->GetTitle()) + " (Fit failed)");
            return true;
        }
    } else {
        std::cout << "Fit failed (reason 2): ndf = " << ndf << " is not greater than 0.\n";
        hInvMass->SetTitle(TString(hInvMass->GetTitle()) + " (Fit failed)");
        return true;
    }
    
    // 3 - Check if mean value is within the histogram range
    double histoMin = hInvMass->GetBinLowEdge(1);
    double histoMax = hInvMass->GetBinLowEdge(hInvMass->GetNbinsX() + 1);
    if (m_0 < histoMin || m_0 > histoMax) {
        std::cout << "Fit failed (reason 3): m_0 = " << m_0 << " is outside of the histogram range [" << histoMin << ", " << histoMax << "].\n";
        hInvMass->SetTitle(TString(hInvMass->GetTitle()) + " (Fit failed)");
        return true;
    }
    
    // 4 - Check if background fit and total fit are too similar on the signal region (e.g., if the ratio of their integrals in the signal region is close to 1, signal quantity is negligible compared to the background)
    // double backgroundIntegral = fitBackgroundOnly->Integral(m_0 - signalSigmas * sigma, m_0 + signalSigmas * sigma);
    // double yieldIntegral = fitTotal->Integral(m_0 - signalSigmas * sigma, m_0 + signalSigmas * sigma);
    // double signal = yieldIntegral - backgroundIntegral;
    // if (yieldIntegral > 0) {
    //     double signalFraction = signal / yieldIntegral;
    //     if (signalFraction < 0.05) {  // <5% signal
    //         std::cout << "Fit failed (reason 4): signal fraction too small (" << signalFraction << ")\n";
    //         hInvMass->SetTitle(TString(hInvMass->GetTitle()) + " (Fit failed)");
    //         return true;
    //     }
    // }
    
    // 5 - Check if at least 3 bins with non-zero content exist in any sideband
    bool leftOutsideOfInvRange = (m_0 - startingBackSigma * signalSigma1) < hInvMass->GetBinLowEdge(1);
    bool rightOutsideOfInvRange = (m_0 + startingBackSigma * signalSigma1) > hInvMass->GetBinLowEdge(hInvMass->GetNbinsX() + 1);

    if (leftOutsideOfInvRange && rightOutsideOfInvRange) {
        std::cout << "Fit failed (reason 5): both sidebands are outside histogram range.\n";
        hInvMass->SetTitle(TString(hInvMass->GetTitle()) + " (Fit failed)");
        return true;
    }

    // Helper lambda: count bins with non-zero content between two x values
    auto countFilledBins = [&](double xLow, double xHigh) -> int {
        int lowBin  = safeLowBin(hInvMass, xLow);
        int highBin = safeHighBin(hInvMass, xHigh);
        int count = 0;
        for (int iBin = lowBin; iBin <= highBin; ++iBin) {
            if (hInvMass->GetBinContent(iBin) > 0) {
                ++count;
            }
        }
        return count;
    };

    // Check sidebands - at least ONE must have >= 3 filled bins
    int filledLeft  = 0;
    int filledRight = 0;

    if (!leftOutsideOfInvRange) {
        double leftSBlow  = hInvMass->GetBinLowEdge(1);
        double leftSBhigh = m_0 - startingBackSigma * signalSigma1;
        filledLeft = countFilledBins(leftSBlow, leftSBhigh);
        std::cout << "Left sideband filled bins: " << filledLeft << "\n";
    }

    if (!rightOutsideOfInvRange) {
        double rightSBlow  = m_0 + startingBackSigma * signalSigma1;
        double rightSBhigh = hInvMass->GetBinLowEdge(hInvMass->GetNbinsX() + 1);
        filledRight = countFilledBins(rightSBlow, rightSBhigh);
        std::cout << "Right sideband filled bins: " << filledRight << "\n";
    }

    // Fail only if NEITHER sideband has enough bins
    bool leftOK  = (!leftOutsideOfInvRange)  && (filledLeft  >= 3);
    bool rightOK = (!rightOutsideOfInvRange) && (filledRight >= 3);

    if (!leftOK && !rightOK) {
        std::cout << "Fit failed (reason 5): neither sideband has >= 3 filled bins "
                << "(left=" << filledLeft << ", right=" << filledRight << ").\n";
        hInvMass->SetTitle(TString(hInvMass->GetTitle()) + " (Fit failed)");
        return true;
    }
    
    return false;
}
// Check if histogram should be erased before storing (in case the fit failed)
bool eraseHistogram(std::vector<bool>& workingFits, const size_t& histoIndex) {
    
    // --- Find largest contiguous block of "true"
    size_t bestStart = 0, bestEnd = 0;
    size_t maxLength = 0;

    for (size_t i = 0; i < workingFits.size(); ++i) {
        if (workingFits[i]) {
            size_t j = i;

            while (j < workingFits.size() && workingFits[j]) {
                ++j;
            }

            size_t length = j - i;
            // Prioritize longer blocks, and in case of ties, the one starting at a higher index
            if (length > maxLength || (length == maxLength && i > bestStart)) {
                maxLength = length;
                bestStart = i;
                bestEnd = j - 1;
            }

            i = j; // skip checked region
        }
    }

    // --- If no valid block exists → erase everything
    if (maxLength == 0) {
        return true;
    }

    // --- Keep only histograms inside the best block

    bool shouldErase = !((histoIndex >= bestStart) && (histoIndex <= bestEnd));
    // Make sure to mark the workingFits vector
    if (shouldErase) {
        workingFits[histoIndex] = false;
    }
    
    return shouldErase;
}
// Emma's work applied specific cuts to certain pT,jet bins
bool passEmmaCut(double jetPt, double hfPt) {
    
    // Define bins: {jetPt_min, jetPt_max, hfPt_min, hfPt_max}
    const std::vector<std::tuple<double,double,double,double>> cuts = {
        {5., 7., 2., 7.},
        {7., 10., 3., 10.},
        {10., 20., 5., 20.},
        {20., 50., 12., 20.},
        // {5., 50., 1., 36.} // allow full range for the inclusive distribution (5 < jetPt < 50)
    };

    for (const auto& [jetMin, jetMax, hfMin, hfMax] : cuts) {
        if (jetPt >= jetMin && jetPt < jetMax) {
            return (hfPt >= hfMin) && (hfPt < hfMax);
        }
    }

    return false; // outside defined ranges → reject
}
bool passChrisCut(double jetPt, double hfPt) {
    
    // Define bins: {jetPt_min, jetPt_max, hfPt_min, hfPt_max}
    const std::vector<std::tuple<double,double,double,double>> cuts = {
        {5., 7., 1., 7.},
        {7., 10., 1., 10.},
        {10., 16., 1., 16.},
        {16., 36., 2., 36.},
        {36., 50., 12., 36.},
        {5., 50., 1., 36.} // allow full range for the inclusive distribution (5 < jetPt < 50)
    };

    for (const auto& [jetMin, jetMax, hfMin, hfMax] : cuts) {
        if (jetPt >= jetMin && jetPt < jetMax) {
            return (hfPt >= hfMin) && (hfPt < hfMax);
        }
    }

    return false; // outside defined ranges → reject
}
void cleanNaNs(TH3D* histogram) {
    int xBins = histogram->GetXaxis()->GetNbins();
    int yBins = histogram->GetYaxis()->GetNbins();
    int zBins = histogram->GetZaxis()->GetNbins();

    for (int x = 1; x <= xBins; x++) {
        for (int y = 1; y <= yBins; y++) {
            for (int z = 1; z <= zBins; z++) {

                double content = histogram->GetBinContent(x, y, z);
                double error   = histogram->GetBinError(x, y, z);

                // Check content
                if (std::isnan(content)) {
                    std::cerr << "Cleaning NaN content at (" << x << "," << y << "," << z << ")\n";

                    histogram->SetBinContent(x, y, z, 0.0);
                    histogram->SetBinError(x, y, z, 0.0);
                }
            }
        }
    }
}
/*___________________________________________________Background subtraction (end)______________________________________________________________________________________________*/


/*___________________________________________________Efficiency (start)______________________________________________________________________________________*/
// calculate number of background subtracted histograms inside file
int HistogramCounter(TFile* file) {
    TList* keys = file->GetListOfKeys();
    std::set<TString> uniqueHistograms;
    int numHistograms = 0;

    for (int i = 0; i < keys->GetSize(); ++i) {
        TKey* key = (TKey*)keys->At(i);
        TString name = key->GetName();

        if (name.BeginsWith("h_back_subtracted_")) {
            // Only consider the latest cycle (ignore duplicates)
            if (uniqueHistograms.find(name) == uniqueHistograms.end()) {
                TObject* obj = file->Get(name);
                if (obj && obj->IsA()->InheritsFrom(TH1::Class())) {
                    uniqueHistograms.insert(name);
                    numHistograms++;
                }
            }
        }
    }

    return numHistograms;
}
/*___________________________________________________Efficiency (end)______________________________________________________________________________________*/

/*___________________________________________________Matching (start)______________________________________________________________________________________*/
// Multiply hTruth and response matrix, The result will the the folded histogram
TH2D* manualFoldingOld(RooUnfoldResponse response, TH2D* hTruth, TH2D* hMeasured) {
    
    // Create empty histogram for the folded data
    TH2D* hFolded = (TH2D*)hMeasured->Clone("hFolded");
    hFolded->Reset();
    hFolded->Sumw2();

    // Get the number of bins in measured and truth histograms
    int nBinsXMeasured = hFolded->GetNbinsX();
    int nBinsYMeasured = hFolded->GetNbinsY();
    int nBinsXTruth = hTruth->GetNbinsX();
    int nBinsYTruth = hTruth->GetNbinsY();
    
    // Debug: print a few values from the response matrix
    bool debug = false;
    if (debug) {
        std::cout << "Measured histogram: " << nBinsXMeasured << " x " << nBinsYMeasured << std::endl;
        std::cout << "Response matrix dimensions: " 
                << response.GetNbinsMeasured() << " x " << response.GetNbinsTruth() << std::endl;
    }
    

    // loop through detector level bins
    for (int iMeasured = 0; iMeasured < nBinsXMeasured; iMeasured++) {
        for (int jMeasured = 0; jMeasured < nBinsYMeasured; jMeasured++) {
            double foldedValue = 0;
            double foldedError2 = 0;

            // obtaining flattened 1D index through row-major ordering
            int index_x_measured = iMeasured + nBinsXMeasured*jMeasured;

            // calculating element iMeasured,jMeasured of folded 2D matrix
            for (int iTruth = 0; iTruth < nBinsXTruth; iTruth++) {
                for (int jTruth = 0; jTruth < nBinsYTruth; jTruth++) {
                    // obtaining flattened 1D index through row-major ordering
                    int index_x_truth = iTruth + nBinsXTruth*jTruth;
                    // calculating matrix element product
                    double truthValue = hTruth->GetBinContent(iTruth + 1,jTruth + 1);
                    double responseValue = response(index_x_measured, index_x_truth);
                    foldedValue += truthValue*responseValue;
                    foldedError2 += std::pow(hTruth->GetBinError(iTruth + 1,jTruth + 1),2) * std::pow(response(index_x_measured, index_x_truth),2);
                }
            }
            hFolded->SetBinContent(iMeasured + 1, jMeasured + 1, foldedValue);
            hFolded->SetBinError(iMeasured + 1, jMeasured + 1, std::sqrt(foldedError2));
        }
    }
    
    std::cout << "Manual folding performed." << std::endl;

    return hFolded;
}
template<typename T>
T manualFolding(RooUnfoldResponse response, T hTruth, T hMeasured) {

    // Check if T is TH1 or TH2 or TH3
    if constexpr (std::is_same_v<T, TH1D*>) {
        TH1D* hFolded = (TH1D*)hMeasured->Clone("hFolded");
        hFolded->Reset();
        hFolded->Sumw2();

        int nBinsXMeasured = hFolded->GetNbinsX();
        int nBinsXTruth = hTruth->GetNbinsX();

        for (int iMeasured = 0; iMeasured < nBinsXMeasured; iMeasured++) {
            double foldedValue = 0;
            double foldedError2 = 0;

            for (int iTruth = 0; iTruth < nBinsXTruth; iTruth++) {
                double truthValue = hTruth->GetBinContent(iTruth + 1);
                double responseValue = response(iMeasured, iTruth);
                
                foldedValue += truthValue * responseValue;
                foldedError2 += std::pow(hTruth->GetBinError(iTruth + 1), 2) * std::pow(responseValue, 2);
            }
            hFolded->SetBinContent(iMeasured + 1, foldedValue);
            hFolded->SetBinError(iMeasured + 1, std::sqrt(foldedError2));
        }
        std::cout << "Manual 1D folding performed." << std::endl;
        return hFolded;
    } else if constexpr (std::is_same_v<T, TH2D*>) {
        // Create empty histogram for the folded data
        TH2D* hFolded = (TH2D*)hMeasured->Clone("hFolded");
        hFolded->Reset();
        hFolded->Sumw2();

        // Get the number of bins in measured and truth histograms
        int nBinsXMeasured = hFolded->GetNbinsX();
        int nBinsYMeasured = hFolded->GetNbinsY();
        int nBinsXTruth = hTruth->GetNbinsX();
        int nBinsYTruth = hTruth->GetNbinsY();
        
        // Debug: print a few values from the response matrix
        bool debug = false;
        if (debug) {
            std::cout << "Measured histogram: " << nBinsXMeasured << " x " << nBinsYMeasured << std::endl;
            std::cout << "Response matrix dimensions: " << response.GetNbinsMeasured() << " x " << response.GetNbinsTruth() << std::endl;
        }

        // loop through detector level bins
        for (int iMeasured = 0; iMeasured < nBinsXMeasured; iMeasured++) {
            for (int jMeasured = 0; jMeasured < nBinsYMeasured; jMeasured++) {
                double foldedValue = 0;
                double foldedError2 = 0;

                // obtaining flattened 1D index through row-major ordering
                int index_x_measured = iMeasured + nBinsXMeasured*jMeasured;

                // calculating element iMeasured,jMeasured of folded 2D matrix
                for (int iTruth = 0; iTruth < nBinsXTruth; iTruth++) {
                    for (int jTruth = 0; jTruth < nBinsYTruth; jTruth++) {
                        // obtaining flattened 1D index through row-major ordering
                        int index_x_truth = iTruth + nBinsXTruth*jTruth;
                        // calculating matrix element product
                        double truthValue = hTruth->GetBinContent(iTruth + 1,jTruth + 1);
                        double responseValue = response(index_x_measured, index_x_truth);
                        foldedValue += truthValue*responseValue;
                        foldedError2 += std::pow(hTruth->GetBinError(iTruth + 1,jTruth + 1),2) * std::pow(response(index_x_measured, index_x_truth),2);
                    }
                }
                hFolded->SetBinContent(iMeasured + 1, jMeasured + 1, foldedValue);
                hFolded->SetBinError(iMeasured + 1, jMeasured + 1, std::sqrt(foldedError2));
            }
        }
        
        std::cout << "Manual 2D folding performed." << std::endl;

        return hFolded;
    }
    else if constexpr (std::is_same_v<T, TH3D*>) {
        // Create empty histogram for the folded data
        TH3D* hFolded = (TH3D*)hMeasured->Clone("hFolded");
        hFolded->Reset();
        hFolded->Sumw2();

        // Get the number of bins in measured and truth histograms
        int nBinsXMeasured = hFolded->GetNbinsX();
        int nBinsYMeasured = hFolded->GetNbinsY();
        int nBinsZMeasured = hFolded->GetNbinsZ();
        int nBinsXTruth = hTruth->GetNbinsX();
        int nBinsYTruth = hTruth->GetNbinsY();
        int nBinsZTruth = hTruth->GetNbinsZ();

        // Debug: print a few values from the response matrix
        bool debug = false;
        if (debug) {
            std::cout << "Measured histogram: " << nBinsXMeasured << " x " << nBinsYMeasured << " x " << nBinsZMeasured << std::endl;
            std::cout << "Response matrix dimensions: " << response.GetNbinsMeasured() << " x " << response.GetNbinsTruth() << std::endl;
        }

        // loop through detector level bins
        for (int iMeasured = 0; iMeasured < nBinsXMeasured; iMeasured++) {
            for (int jMeasured = 0; jMeasured < nBinsYMeasured; jMeasured++) {
                for (int kMeasured = 0; kMeasured < nBinsZMeasured; kMeasured++) {
                    double foldedValue = 0;
                    double foldedError2 = 0;

                    // obtaining flattened 1D index through row-major ordering
                    int index_x_measured = iMeasured + (nBinsXMeasured * jMeasured) + (nBinsXMeasured * nBinsYMeasured * kMeasured);

                    // calculating element iMeasured,jMeasured of folded 2D matrix
                    for (int iTruth = 0; iTruth < nBinsXTruth; iTruth++) {
                        for (int jTruth = 0; jTruth < nBinsYTruth; jTruth++) {
                            for (int kTruth = 0; kTruth < nBinsZTruth; kTruth++) {
                                // obtaining flattened 1D index through row-major ordering
                                int index_x_truth = iTruth + (nBinsXTruth * jTruth) + (nBinsXTruth * nBinsYTruth * kTruth);
                                // calculating matrix element product
                                double truthValue = hTruth->GetBinContent(iTruth + 1,jTruth + 1, kTruth + 1);
                                double responseValue = response(index_x_measured, index_x_truth);
                                foldedValue += truthValue*responseValue;
                                foldedError2 += std::pow(hTruth->GetBinError(iTruth + 1,jTruth + 1, kTruth + 1),2) * std::pow(response(index_x_measured, index_x_truth),2);
                            }
                        }
                    }
                    hFolded->SetBinContent(iMeasured + 1, jMeasured + 1, kMeasured + 1, foldedValue);
                    hFolded->SetBinError(iMeasured + 1, jMeasured + 1, kMeasured + 1, std::sqrt(foldedError2));
                }
            }
        }
        
        std::cout << "Manual 3D folding performed." << std::endl;

        return hFolded;
    } else {
        static_assert(sizeof(T) == 0, "Unsupported histogram type! Only TH1D*, TH2D*, and TH3D* are accepted.");
    }
}
/*___________________________________________________Matching (end)______________________________________________________________________________________*/

/*___________________________________________________Feed-down (end)______________________________________________________________________________________*/
// Get POWHEG and data luminosities
// std::vector<double> getLuminosities(TFile* fPowheg, TFile* fData) {
//     // Accessing total cross section value stored in first bin (in mb)
//     TH1D* xSection_powheg = dynamic_cast<TH1D*>(fPowheg->Get("fHistXsection"));
//     double crossSecPowheg = xSection_powheg->GetBinContent(1);

//     // Accessing number of events
//     TTree* tree_D0 = dynamic_cast<TTree*>(fPowheg->Get("tree_D0"));
//     double numOfEventsPowheg = tree_D0->GetEntries();

//     // integrated POWHEG luminosity 
//     double luminosity_powheg = numOfEventsPowheg/crossSecPowheg;
//     double lumiMC = luminosity_powheg; // Store in dataContainer for later use
//     std::cout << "POWHEG luminosity = " << luminosity_powheg << " mb^-1" << std::endl;

//     // Calculating measured luminosity
//     // folder: jet-luminosity-calculator
//     // Histogram name: counter
//     // bin labels: "BC+TVX", "Coll+TVX", "Coll+TVX+VtxZ+Sel8"
//     TH1D* hDataLumi = (TH1D*)fData->Get("jet-luminosity-calculator/counter");
//     int bin = hDataLumi->GetXaxis()->FindBin("BC+TVX");
//     double bcTVX = hDataLumi->GetBinContent(bin); // BC+TVX
//     std::cout << "On bin number " << bin <<", BC+TVX = " << bcTVX << std::endl;
//     bin = hDataLumi->GetXaxis()->FindBin("Coll+TVX+VtxZ+Sel8");
//     double selection = hDataLumi->GetBinContent(bin); // Coll+TVX+VtxZ+Sel8
//     std::cout << "On bin number " << bin <<", Coll+TVX+VtxZ+Sel8 = " << selection << std::endl;
//     bin = hDataLumi->GetXaxis()->FindBin("Coll+TVX");
//     double collTVX = hDataLumi->GetBinContent(bin); // Coll+TVX
//     std::cout << "On bin number " << bin <<", Coll+TVX = " << collTVX << std::endl;
//     double triggered = bcTVX * selection / collTVX; // number of TVX triggered BC that correspond to your selections and your train efficiencies
//     double runLuminosity = 1.0/0.0594e6; // luminosity value for the runs (// in mb⁻¹?)
//     double dataLuminosity = triggered * runLuminosity;
//     std::cout << "triggered = " << triggered << "; runLuminosity = " << runLuminosity << std::endl;
//     std::cout << "BC+TVX = " << bcTVX << ", Coll+TVX = " << collTVX << ", Selection = " << selection << ", runLuminosity = " << runLuminosity << std::endl;
//     double lumiData = dataLuminosity; // Store in dataContainer for later use
//     std::cout << "Measured luminosity = " << dataLuminosity << " mb^-1" << std::endl;

//     std::vector<double> luminosities = {luminosity_powheg, dataLuminosity};
//     return luminosities;
// }
/*___________________________________________________Feed-down (end)______________________________________________________________________________________*/
/*___________________________________________________Debugging (start)______________________________________________________________________________________*/
void checkCanvas(TCanvas* canvas) {
    if (!canvas) {
        std::cout << "canvas is NULL!" << std::endl;
        return;
    }

    TList* prims = canvas->GetListOfPrimitives();
    if (!prims || prims->GetSize() == 0) {
        std::cout << "canvas has NO primitives!" << std::endl;
        return;
    }

    std::cout << "Canvas: " << canvas->GetName() << std::endl;

    int ipad = 0;
    for (TObject* obj : *prims) {

        if (obj->InheritsFrom(TPad::Class())) {
            TPad* pad = (TPad*)obj;
            ipad++;

            TList* padPrims = pad->GetListOfPrimitives();

            if (!padPrims || padPrims->GetSize() == 0) {
                std::cout << "  Pad " << ipad << " (" << pad->GetName() << ") is EMPTY\n";
            } else {
                std::cout << "  Pad " << ipad << " (" << pad->GetName() << ") has "
                          << padPrims->GetSize() << " objects:\n";

                for (TObject* pobj : *padPrims) {
                    std::cout << "    - " << pobj->GetName()
                              << " (" << pobj->ClassName() << ")\n";
                }
            }
        }
    }
}
bool usingEmmasBinning(TH1D* hIsEmmaYeatsBinning) {
    if (hIsEmmaYeatsBinning->GetBinContent(1) > 0.) {
        return false; // using standard binning
    } else if (hIsEmmaYeatsBinning->GetBinContent(2) > 0.) {
        return true; // using Emma Yeats binning
    }
    return false; // default case
}
/*___________________________________________________Debugging (end)______________________________________________________________________________________*/
