/*
 * .h file containing common utilities and functions used in the closure tests
 * 
 * 
 * 
 * Author: Christian Reckziegel
**/


using namespace std;

// (Matched from) / (selected as)...:Lc = +1, Lcbar = -1, neither = 0
enum class LcSpecies : int {
    LCBAR = -1,
    NEITHER = 0,
    LC = +1
};
// Conversion function
LcSpecies intToLcSpecies(int value) {
    switch(value) {
        case -1: return LcSpecies::LCBAR;
        case 0:  return LcSpecies::NEITHER;
        case 1:  return LcSpecies::LC;
        default: return LcSpecies::NEITHER; // Handle unexpected values
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

// Binning object creation, storage and retrieval
struct BinningStruct {
    // Detector level binning
    std::vector<double> ptjetBinEdges_detector;
    std::vector<double> deltaRBinEdges_detector;
    std::vector<double> ptHFBinEdges_detector;
    // BDT cut values for each pT,HF bin
    std::vector<std::pair<double, double>> bdtPtCuts;

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
