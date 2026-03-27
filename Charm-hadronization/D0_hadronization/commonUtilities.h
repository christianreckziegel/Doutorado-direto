/*
 * .h file containing common utilities and functions used in the closure tests
 * 
 * 
 * 
 * Author: Christian Reckziegel
**/


using namespace std;

// (Matched from) / (selected as)...:Lc = +1, Lcbar = -1, neither = 0
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
        if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) {
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
    // BDT cut values for each pT,HF bin
    std::vector<std::pair<double, double>> bdtPtCuts;
    std::vector<std::pair<double, double>> bdtPtCutsMC;

    // Particle level binning
    std::vector<double> ptjetBinEdges_particle;
    std::vector<double> deltaRBinEdges_particle;
    std::vector<double> ptHFBinEdges_particle;
};
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
    ptEdges.Write("ptBinEdges");
    bdtCuts.Write("bdtCutValues");

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
    TVectorD* vecPt = (TVectorD*)fInput->Get("bdt/ptBinEdges");
    TVectorD* vecCut = (TVectorD*)fInput->Get("bdt/bdtCutValues");
    if (!vecPt || !vecCut) {
        throw std::runtime_error("Could not find BDT cut vectors in file");
    }
    if (vecPt->GetNoElements() != vecCut->GetNoElements()) {
        throw std::runtime_error("BDT pT bins and cut vectors have different sizes");
    }
    binning.bdtPtCuts.clear();
    for (int i = 0; i < vecPt->GetNoElements(); ++i) {
        binning.bdtPtCuts.emplace_back((*vecPt)[i], (*vecCut)[i]);
    }

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
    std::vector<TF1*> fitTotal;                 // Total fit: signal + reflections + background

    // list of working fits (PDG mass inside of signal range)
    std::vector<bool> workingFits;

    // Individual components (primary and secondary Gaussians)
    std::vector<std::pair<TF1*,TF1*>> individualSignals;
    std::vector<std::pair<TF1*,TF1*>> individualReflections;
};
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
    
    std::cout << "Left sideband: [" << leftRange[0] << ", " << leftRange[1] << "]\n";
    std::cout << "Right sideband: [" << rightRange[0] << ", " << rightRange[1] << "]\n";
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
    double B, Bs, R, Rs, S, Ys;

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
    }

    // Calculate scaling factor
    double alpha = S / (S + Rs - (Bs/B) * R);
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
bool didFitFailed() {
    //
    return false;
}
// Check if histogram should be erased before storing (in case the fit failed)
bool eraseHistogram(const std::vector<bool>& workingFits, size_t& histoIndex) {
    bool doEraseHistogram = false;

    // was it stored as sucessfull fit?
    if (workingFits[histoIndex]) {

        // are upper indices all true?
        for (size_t iHisto = histoIndex; iHisto < workingFits.size(); iHisto++) {
            if (!workingFits[iHisto]) {
                doEraseHistogram = true;
                break;
            }
            
        }
        
    } else {
        doEraseHistogram = true;
    }
    
    return doEraseHistogram;
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
TH2D* manualFolding(RooUnfoldResponse response, TH2D* hTruth, TH2D* hMeasured) {
    
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
/*___________________________________________________Matching (end)______________________________________________________________________________________*/