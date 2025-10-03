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
// Define the enum class with descriptive names
enum class FitModelType {
    Full,           // Signal + Reflections + Background
    SignalReflectionsOnly,     // Signal + Reflections only
    SignalOnly,     // Only Signal
    ReflectionsOnly         // Only Reflections
};
//_________________________________Template for on-the-run chosen model (policy-based strategy)________________________________________________
struct BackgroundPolicy {
    static double eval(double* x, double* par) {
        Double_t m = x[0];
        Double_t a = par[0];
        Double_t b = par[1];

        // Defining the custom function
        Double_t result = a * TMath::Power(m, b);
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
using FullModel = FitModel<SignalPolicy, ReflectionPolicy, BackgroundPolicy>;   // Signal + Reflection + Background
using SigRefModel = FitModel<SignalPolicy, ReflectionPolicy>;                   // Signal + Reflection only
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

std::vector<double> LoadBinning(TFile* fInput, const char* pathInFile) {
    auto* vec = (TVectorD*)fInput->Get(pathInFile);
    if (!vec) {
        throw std::runtime_error(Form("Could not find TVectorD at '%s'", pathInFile));
    }
    return std::vector<double>(vec->GetMatrixArray(), vec->GetMatrixArray() + vec->GetNoElements());
}

//__________________________________________________________________________________________________________________________
// Module to create TH2D histograms including interest variable: UNIFORM bin sizes
std::vector<TH2D*> createHistograms(const std::vector<double>& ptDBinEdges, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax, double jetptMin) {
    std::vector<TH2D*> histograms;
    for (size_t i = 0; i < ptDBinEdges.size() - 1; ++i) {
        histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, ybins, ymin, ymax));
        // So that the title adapts to fractional binning title
        if (std::fmod(ptDBinEdges[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(ptDBinEdges[i+1], 1.0) != 0) {
                histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.1f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, ybins, ymin, ymax));
        
            } else {
                histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.1f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, ybins, ymin, ymax));
        
            }
        } else {
            if (std::fmod(ptDBinEdges[i+1], 1.0) != 0) {
                histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, ybins, ymin, ymax));
        
            } else {
                histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, ybins, ymin, ymax));
        
            }
        }
        histograms[i]->Sumw2();
        
    }
    return histograms;
}
// Overloading function
// Module to create TH2D histograms including interest variable: VARIABLE bin sizes
std::vector<TH2D*> createHistograms(const std::vector<double>& ptDBinEdges, int xbins, double xmin, double xmax, const std::vector<double>& yBinEdges, const double& jetptMin) {
    // Change binning to low statistics range [30,50] GeV/c
    if (jetptMin >= 30.) {
        xbins = 25;
    }
    
    std::vector<TH2D*> histograms;
    for (size_t i = 0; i < ptDBinEdges.size() - 1; ++i) {
        // So that the title adapts to fractional binning title
        if (std::fmod(ptDBinEdges[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(ptDBinEdges[i+1], 1.0) != 0) {
                histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.1f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            } else {
                histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.1f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            }
        } else {
            if (std::fmod(ptDBinEdges[i+1], 1.0) != 0) {
                histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            } else {
                histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            }
        }
        histograms[i]->Sumw2();
    }
    return histograms;
}
//__________________________________________________________________________________________________________________________
// Find the appropriate cut value based on pT
double GetBkgProbabilityCut(double pT, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    for (size_t i = 0; i < bdtPtCuts.size() - 1; ++i) {
        if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) {
            return bdtPtCuts[i].second;
        }
    }
    return 1.0; // Default: accept all if out of range
}
// Module to fill 2D histograms from TTree data
void fillHistograms(TFile* fDist, const std::vector<TH2D*>& histograms, double jetptMin, double jetptMax, std::vector<double>& ptDBinEdges, std::vector<double>& deltaRBinEdges, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = deltaRBinEdges[deltaRBinEdges.size() - 1]; // on particle level delta R

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
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= deltaRBinEdges[0]) && (deltaR < DeltaRcut))) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < ptDBinEdges.size() - 1 && !filled; iEdge++) {
                if ((hfPt >= ptDBinEdges[iEdge]) && (hfPt < ptDBinEdges[iEdge + 1])) {
                    // Get the threshold for this pT range
                    double maxBkgProb = GetBkgProbabilityCut(hfPt, bdtPtCuts);

                    // Fill histogram only if the cut is passed
                    if (hfMlScore0 < maxBkgProb) {
                        histograms[iEdge]->Fill(hfMass, deltaR);
                    }
                    filled = true; // Exit the loop once the correct histogram is found (alternative: break)
                }
                
            }
            
        }

        
        
        
    }
    cout << "Histograms filled.\n";
}

// Module to perform fits to histograms
struct FitContainer {
    std::vector<TF1*> fitBackgroundOnly;        // background only fits
    std::vector<TF1*> fitSignalOnly;            // background only fits
    std::vector<TF1*> fitReflectionsOnly;       // background only fits
    std::vector<TF1*> fitTotal;                 // background only fits
    std::vector<bool> workingFits;              // list of working fits (PDG mass inside of signal range)
};
FitContainer performFit(TFile* fReflectionsMC, const std::vector<TH2D*>& histograms2d, double& minMass, double& maxMass, const FitModelType& modelToUse, const double& jetptMin, const double& jetptMax) {
    
    double m_0_reference = 1.86484; // D0 mass in GeV/c^2
    double sigma_reference = 0.012;

    // Create fits container
    FitContainer fittings;
    
    // Create temporary histogram for collecting data
    TH1D* tempHist;
    // Creating 1D mass projection histograms
    std::vector<TH1D*> histograms;

    for (size_t iHisto = 0; iHisto < histograms2d.size(); ++iHisto) {
        /* code */
        tempHist = histograms2d[iHisto]->ProjectionX(Form("h_mass_proj_%zu", iHisto));
        histograms.push_back(tempHist);

        // ----Obtain MC fit parameters----

        // Get TF1 objects from MC file
        TF1* fSignal = (TF1*)fReflectionsMC->Get(Form("signalFit_%zu", iHisto));
        //fSignal->Print("V");
        TF1* fReflections = (TF1*)fReflectionsMC->Get(Form("reflectionsFit_%zu", iHisto));

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
        //fittings.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu", iHisto + 1), customFitFunction, minMass, maxMass, 13)); // 12 parameters = 8 fixed + 5 free
        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        if (modelToUse == FitModelType::Full) {
            fittings.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<FullModel>, minMass, maxMass, 13));
        } else if (modelToUse == FitModelType::SignalReflectionsOnly) {
            fittings.fitTotal.push_back(new TF1(Form("totalFit_histMass%zu_%0.f_to_%0.fGeV", iHisto + 1, jetptMin, jetptMax), fitWrapper<SigRefModel>, minMass, maxMass, 13));
        }

        // Set initial values and fix parameters
        Double_t params[13] = {200000, -10.0, 1000, 1.5, m_0_reference, 0.02, 1.2, 2.0, 1.3, 1.83, 1.85, 0.02, 0.03};
        fittings.fitTotal[iHisto]->SetParameters(params);

        // Set parameter names
        fittings.fitTotal[iHisto]->SetParName(0, "Background A");
        fittings.fitTotal[iHisto]->SetParName(1, "Background B");
        fittings.fitTotal[iHisto]->SetParName(2, "A1 Signal");
        fittings.fitTotal[iHisto]->SetParName(3, "A1/A2 Signal Ratio");
        fittings.fitTotal[iHisto]->SetParName(4, "Signal Mean m0");
        fittings.fitTotal[iHisto]->SetParName(5, "Signal Width Sigma1");
        fittings.fitTotal[iHisto]->SetParName(6, "Sigma Ratio (Sig1/Sig2)");
        fittings.fitTotal[iHisto]->SetParName(7, "A1Signal/A1Reflection Ratio");
        fittings.fitTotal[iHisto]->SetParName(8, "A1/A2 Reflection Ratio");
        fittings.fitTotal[iHisto]->SetParName(9, "Reflection Mean m0_1");
        fittings.fitTotal[iHisto]->SetParName(10, "Reflection Mean m0_2");
        fittings.fitTotal[iHisto]->SetParName(11, "Reflection Width Sigma_1");
        fittings.fitTotal[iHisto]->SetParName(12, "Reflection Width Sigma_2");

        // Fix the parameters from MC fits
        fittings.fitTotal[iHisto]->FixParameter(3, A1toA2MCSignalRatio);   // A1toA2MCSignalRatio
        fittings.fitTotal[iHisto]->FixParameter(6, Sigma1toSigma2MCSignalRatio);   // Sigma1toSigma2MCSignalRatio
        fittings.fitTotal[iHisto]->FixParameter(7, A1SignalToA1ReflectionMCRatio);   // A1SignalToA1ReflectionMCRatio
        fittings.fitTotal[iHisto]->FixParameter(8, A1toA2MCReflectionsRatio);   // A1toA2MCReflectionsRatio
        fittings.fitTotal[iHisto]->FixParameter(9, m0_1Reflections);  // m0_1
        fittings.fitTotal[iHisto]->FixParameter(10, m0_2Reflections); // m0_2
        fittings.fitTotal[iHisto]->FixParameter(11, sigma1Reflections); // sigma_1
        fittings.fitTotal[iHisto]->FixParameter(12, sigma2Reflections); // sigma_2

        // Apply range limits to the parameters
        fittings.fitTotal[iHisto]->SetParLimits(0, 0., TMath::Infinity()); // only accepts non-negative background fits
        fittings.fitTotal[iHisto]->SetParLimits(4, 0.95 * m_0_reference, 1.05 * m_0_reference);
        fittings.fitTotal[iHisto]->SetParLimits(5, 0.35 * sigma_reference, 3.0 * sigma_reference);

        // Choose line color
        fittings.fitTotal[iHisto]->SetLineColor(kBlack);

        // Perform fit with "Q" (quiet) option: no drawing of the fit function
        histograms[iHisto]->Fit(fittings.fitTotal[iHisto], "Q");

        // Check if signal gaussian acomodates the literature D0 rest mass
        double primaryMean = fittings.fitTotal[iHisto]->GetParameter(4);
        double primarySigma = fittings.fitTotal[iHisto]->GetParameter(5);
        double D0LiteratureMass = 1.86484; // GeV/c²
        if ((D0LiteratureMass < (primaryMean - 2*primarySigma)) || (D0LiteratureMass > (primaryMean + 2*primarySigma))) { // fit failed: PDG mass outside of signal region
            // if literature mass is outside, clear histogram as it should be excluded
            //histograms[iHisto]->Reset();
            histograms[iHisto]->SetTitle(TString(histograms[iHisto]->GetTitle()) + " (m_{D^{0}^{literature}} #notin [#bar{m}_{primary}-2#sigma_{primary};#bar{m}_{primary}+2#sigma_{primary}])");
            //histograms2d[iHisto]->Reset("ICES");
            //histograms2d[iHisto]->Reset();
            histograms2d[iHisto]->SetTitle(TString(histograms[iHisto]->GetTitle()) + " (m_{D^{0}^{literature}} #notin [#bar{m}_{primary}-2#sigma_{primary};#bar{m}_{primary}+2#sigma_{primary}])");

            fittings.workingFits.push_back(false);
        } else { // fit worked: PDG mass inside of signal region
            fittings.workingFits.push_back(true);
        }
        

        // Print the fit results
        //fittings.fitTotal[iHisto]->Print("V");

        //std::cout << "Fit " << iHisto << " performed.\n";
    }
    std::cout << "Number of fits performed: " << fittings.workingFits.size() << std::endl;
    for (size_t iFit = 0; iFit < fittings.workingFits.size(); iFit++) {
        std::cout << "Fit number " << iFit << " works: " << (fittings.workingFits[iFit] ? "true" : "false") << std::endl;
    }
    


    // Perform background only fit to each histogram
    int numberOfFits = fittings.fitTotal.size();
    for (size_t iHisto = 0; iHisto < numberOfFits; iHisto++) {
        // Getting total parameter values
        double a_par = fittings.fitTotal[iHisto]->GetParameter(0); // Get the value of parameter 'a'
        double b_par = fittings.fitTotal[iHisto]->GetParameter(1); // Get the value of parameter 'b'
        fittings.fitBackgroundOnly.push_back(new TF1(Form("backgroundOnlyFit_%zu", iHisto), backgroundFunction, minMass, maxMass, 2));
        fittings.fitBackgroundOnly[iHisto]->SetParameters(a_par,b_par); // fittings.size()-1 = the latest added to the vector
        fittings.fitBackgroundOnly[iHisto]->SetLineStyle(kDashed);
        fittings.fitBackgroundOnly[iHisto]->SetLineColor(kRed+1);

    }
    std::cout << "Background only fits performed.\n";

    // Perform signal only fit to each histogram
    for (size_t iHisto = 0; iHisto < numberOfFits; iHisto++) {
        // Extract signal-related parameters from the total fit
        double a1Signal = fittings.fitTotal[iHisto]->GetParameter(2); // A1 signal
        double A1toA2MCSignalRatio = fittings.fitTotal[iHisto]->GetParameter(3); // A2 signal
        double m0Signal = fittings.fitTotal[iHisto]->GetParameter(4); // Signal mean m0
        double sigmaSignal = fittings.fitTotal[iHisto]->GetParameter(5); // Signal width sigma
        double Sigma1toSigma2MCSignalRatio = fittings.fitTotal[iHisto]->GetParameter(6); // Sigma1/Sigma2 signal
        //std::cout << "A1Signal = " << a1Signal << ", A1toA2MCSignalRatio = " << A1toA2MCSignalRatio << ", m0Signal = " << m0Signal << ", sigmaSignal = " << sigmaSignal << ", Sigma1toSigma2MCSignalRatio = " << Sigma1toSigma2MCSignalRatio << std::endl;

        // Perfoming fit with acquired parameters from total fit
        fittings.fitSignalOnly.push_back(new TF1(Form("signalOnlyFit_%zu", iHisto), signalOnlyFunction, minMass, maxMass, 5));
        fittings.fitSignalOnly[iHisto]->SetParameters(a1Signal, A1toA2MCSignalRatio, m0Signal, sigmaSignal, Sigma1toSigma2MCSignalRatio); // fittings.size()-1 = the latest added to the vector
        fittings.fitSignalOnly[iHisto]->SetLineStyle(kDashed);
        fittings.fitSignalOnly[iHisto]->SetLineColor(kBlue+1);

        //fittings.fitSignalOnly[iHisto]->Print("V");
    }
    std::cout << "Signal only fits performed.\n";

    // Perform reflections only fit to each histogram
    for (size_t iHisto = 0; iHisto < numberOfFits; iHisto++) {
        // Extract reflection-related parameters from the total fit
        double A1SignalToA1ReflectionMCRatios = fittings.fitTotal[iHisto]->GetParameter(7); // A1 signal / A1 reflection
        double A1Signal = fittings.fitTotal[iHisto]->GetParameter(2); // A1 signal
        double A1toA2MCReflectionsRatio = fittings.fitTotal[iHisto]->GetParameter(8); // A1 reflection / A2 reflection
        double m0_1Reflection = fittings.fitTotal[iHisto]->GetParameter(9); // Reflection mean m0_1
        double m0_2Reflection = fittings.fitTotal[iHisto]->GetParameter(10); // Reflection mean m0_2
        double sigma1Reflection = fittings.fitTotal[iHisto]->GetParameter(11); // Reflection sigma_1
        double sigma2Reflection = fittings.fitTotal[iHisto]->GetParameter(12); // Reflection sigma_2
        //std:cout << "A1SignalToA1ReflectionMCRatios = " << A1SignalToA1ReflectionMCRatios << ", A1Signal = " << A1Signal << ", A1toA2MCReflectionsRatio = " << A1toA2MCReflectionsRatio << ", m0_1Reflection = " << m0_1Reflection << ", m0_2Reflection = " << m0_2Reflection << ", sigma1Reflection = " << sigma1Reflection << ", sigma2Reflection = " << sigma2Reflection << std::endl;

        // Perfoming fit with acquired parameters from total fit
        fittings.fitReflectionsOnly.push_back(new TF1(Form("reflectionOnlyFit_%zu", iHisto), reflectionOnlyFunction, minMass, maxMass, 7));
        fittings.fitReflectionsOnly[iHisto]->SetParameters(A1SignalToA1ReflectionMCRatios, A1Signal, A1toA2MCReflectionsRatio, m0_1Reflection, m0_2Reflection, sigma1Reflection, sigma2Reflection); // fittings.size()-1 = the latest added to the vector
        fittings.fitReflectionsOnly[iHisto]->SetLineStyle(kDashed);
        fittings.fitReflectionsOnly[iHisto]->SetLineColor(kGreen+3);
        //fittings.fitReflectionsOnly[iHisto]->Print("V");
        //std::cout << "A1 = " << A1Signal / A1SignalToA1ReflectionMCRatios << ", A2 = " << (A1Signal / A1SignalToA1ReflectionMCRatios) / A1toA2MCReflectionsRatio << std::endl;
        //std::cout << "A1/A2 = " << A1toA2MCReflectionsRatio << std::endl;

    }
   //std::cout << "Reflections only fits performed.\n";
    
    cout << "Fits performed.\n";

    return fittings;
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
    double B, Bs, R, Rs, S;

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
    } else {
        // Both sidebands used

        // Calculate background regions' areas
        B = fittings.fitBackgroundOnly[iHisto]->Integral(leftRange[0], leftRange[1]) + fittings.fitBackgroundOnly[iHisto]->Integral(rightRange[0], rightRange[1]);
        Bs = fittings.fitBackgroundOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of background in signal region

        // Calculate reflection regions' areas
        R = fittings.fitReflectionsOnly[iHisto]->Integral(leftRange[0], leftRange[1]) + fittings.fitReflectionsOnly[iHisto]->Integral(rightRange[0], rightRange[1]);
        Rs = fittings.fitReflectionsOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of reflections in signal region

        // Calculate signal regions' areas
        S = fittings.fitSignalOnly[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of signal in signal region
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
    scallingFactors[1] = Bs/B;

    std::cout << "Scaling factors calculated.\n";

    return scallingFactors;
}

// Check if histogram should be erased before storing
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

struct SubtractionResult {
    std::vector<TH1D*> histograms;      // 1D mass projection histograms vector
    std::vector<TH1D*> sidebandHist;    // 1D sideband background histograms vector
    std::vector<TH1D*> signalHist;      // 1D signal histograms vector
    std::vector<TH1D*> subtractedHist;  // 1D sideband subtracted histograms vector
    TH1D* hSubtracted_allPtSummed;      // final deltaR distribution for all pT,D summed
    TH1D* hSignificance;    // final significance distribution for all pT,D summed
    std::vector<std::pair<double, double>> signal_background_values; // vector of signal and background values for each pT,D bin
};
SubtractionResult SideBand(const std::vector<TH2D*>& histograms2d, const FitContainer& fittings, double signalSigmas, int startingBackSigma, int backgroundSigmas, std::vector<double>& ptDBinEdges){
    
    // Creating histograms for collecting data
    TH1D* tempHist; // temporary histogram for collecting data
    TH1D* h_sideBand;
    TH1D* h_signal;
    TH1D* h_back_subtracted;

    // Creating output struct object with histogram vectors
    SubtractionResult outputStruct;

    // Calculating scaling parameter
    for (size_t iHisto = 0; iHisto < histograms2d.size(); ++iHisto) {
        std::cout << "\nPerforming side-band subtraction for fit number " << iHisto << std::endl;
        // Get total fit parameters
        double m_0 = fittings.fitTotal[iHisto]->GetParameter(4); // Get the value of parameter 'm_0'
        double signalSigma1 = fittings.fitTotal[iHisto]->GetParameter(5); // Get the value of parameter 'sigma1'
        double signalSigma2 = fittings.fitTotal[iHisto]->GetParameter(5) / fittings.fitTotal[iHisto]->GetParameter(6); // signalSigma2 = sigma1 / sigmaRatio12
        double sigma = signalSigma1;
        
        // Obtain 1D invariant mass histograms from the projection
        tempHist = histograms2d[iHisto]->ProjectionX(Form("h_mass_proj_%zu", iHisto));
        outputStruct.histograms.push_back(tempHist);

        // Calculate sideband regions
        std::pair<std::array<double, 2>, std::array<double, 2>> sidebandRanges = calculateSidebandRegions(iHisto, fittings, outputStruct.histograms[iHisto], startingBackSigma, backgroundSigmas);

        // Calculate scaling factor to apply to the sideband subtraction original method
        std::array<double, 2> scallingFactors = calculateScalingFactor(iHisto, fittings, sidebandRanges.first, sidebandRanges.second, signalSigmas);

        // Create signal histogram
        int lowBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 - signalSigmas * sigma);
        int highBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 + signalSigmas * sigma);
        h_signal = histograms2d[iHisto]->ProjectionY(Form("h_signal_proj_%zu",iHisto), lowBin, highBin);
        outputStruct.signalHist.push_back(h_signal);

        // Verify match between fit and histogram areas
        bool proportionsPrintCheck = false;

        // Create side-band histogram
        if (sidebandRanges.first[0] == 0. && sidebandRanges.first[1] == 0.) {
            // Use only the right sideband
            lowBin = histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[0]);
            highBin = histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[1]);
            h_sideBand = histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_%zu",iHisto), lowBin, highBin); // sum the right sideband
            double rightSBHistogram = outputStruct.histograms[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram

            double rightSB = fittings.fitTotal[iHisto]->Integral(sidebandRanges.second[0],sidebandRanges.second[1]);
            double fitSigIntegral = fittings.fitTotal[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

            lowBin = outputStruct.histograms[iHisto]->FindBin(m_0 - signalSigmas * sigma);
            highBin = outputStruct.histograms[iHisto]->FindBin(m_0 + signalSigmas * sigma);
            double signalHistogram = outputStruct.histograms[iHisto]->Integral(lowBin, highBin);

            if (proportionsPrintCheck) {
                cout << "Distribution ratio 1 of sideband over raw signal: " << rightSBHistogram / signalHistogram << " = " << rightSBHistogram << "/" << signalHistogram << endl;
                cout << "Distribution ratio 2 of sideband over raw signal: " << h_sideBand->Integral() / h_signal->Integral() << " = " << h_sideBand->Integral() << "/" << h_signal->Integral() <<  endl;
                cout << "Fit function ratio of sideband over raw signal: " << rightSB / fitSigIntegral << " = " << rightSB << "/" << fitSigIntegral <<  endl;
                cout << "---------------" << endl;
            }
            
        } else {
            // Use both sidebands
            lowBin = histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.first[0]);
            highBin = histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.first[1]);
            h_sideBand = histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_left_%zu",iHisto), lowBin, highBin); // sum the left sideband
            double leftSBHistogram = outputStruct.histograms[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram

            lowBin = histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[0]);
            highBin = histograms2d[iHisto]->GetXaxis()->FindBin(sidebandRanges.second[1]);
            h_sideBand->Add(histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_%zu",iHisto), lowBin, highBin)); // sum the right sideband
            double rightSBHistogram = outputStruct.histograms[iHisto]->Integral(lowBin, highBin); // integral of sideband regions of mass histogram

            double leftSB = fittings.fitTotal[iHisto]->Integral(sidebandRanges.first[0],sidebandRanges.first[1]); // in the opposite order because the histogram is decreasing
            double rightSB = fittings.fitTotal[iHisto]->Integral(sidebandRanges.second[0],sidebandRanges.second[1]);
            double fitSigIntegral = fittings.fitTotal[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

            lowBin = outputStruct.histograms[iHisto]->FindBin(m_0 - signalSigmas * sigma);
            highBin = outputStruct.histograms[iHisto]->FindBin(m_0 + signalSigmas * sigma);
            double signalHistogram = outputStruct.histograms[iHisto]->Integral(lowBin, highBin);
            
            if (proportionsPrintCheck) {
                cout << "Distribution ratio 1 of sideband over raw signal: " << (leftSBHistogram + rightSBHistogram) / signalHistogram << " = " << leftSBHistogram + rightSBHistogram << "/" << signalHistogram << endl;
                cout << "Distribution ratio 2 of sideband over raw signal: " << h_sideBand->Integral() / h_signal->Integral() << " = " << h_sideBand->Integral() << "/" << h_signal->Integral() <<  endl;
                cout << "Fit function ratio of sideband over raw signal: " << (leftSB+rightSB)/ fitSigIntegral << " = " << leftSB+rightSB << "/" << fitSigIntegral <<  endl;
                cout << "---------------" << endl;
            }
            
        }
        
        // Scale the sideband histogram by ratio of it in the signal region Bs/(B1+B2)
        h_sideBand->Scale(scallingFactors[1]);
        outputStruct.sidebandHist.push_back(h_sideBand);

        // Subtract background histogram from signal histogram
        h_back_subtracted = (TH1D*)h_signal->Clone(Form("h_back_subtracted_%zu",iHisto));
        h_back_subtracted->Add(h_sideBand,-1.0);

        // Scale by reflections correction
        h_back_subtracted->Scale(scallingFactors[0]);

        // Account for two sigma only area used for signal region
        double coverage = TMath::Erf(signalSigmas / sqrt(2)); // coverage of a Gaussian in a +/- signalSigmas window
        h_back_subtracted->Scale(1 / coverage);
        // If histogram isn't empty, store it in the container
        bool isHistoEmpty = (h_back_subtracted->GetEntries() != 0) ? true : false;
        if (true) {
            outputStruct.subtractedHist.push_back(h_back_subtracted);

            // Calculating significance
            // Get signal yield by integrating Gaussian in a 2σ window
            double mean_signalOnly = fittings.fitSignalOnly[iHisto]->GetParameter(2);
            double sigma1_signalOnly = fittings.fitSignalOnly[iHisto]->GetParameter(3);
            //fittings.fitSignalOnly[iHisto]->Print("V");
            double S = fittings.fitSignalOnly[iHisto]->Integral(mean_signalOnly - 2*sigma1_signalOnly, mean_signalOnly + 2*sigma1_signalOnly);
            double B = fittings.fitBackgroundOnly[iHisto]->Integral(mean_signalOnly - 2*sigma1_signalOnly, mean_signalOnly + 2*sigma1_signalOnly);
            double significance = S / sqrt(S+B); // Binomial-like Approximation: when S is comparable to B, consider both signal and background fluctuations
            // Store signal and background values for each pT,D bin
            outputStruct.signal_background_values.push_back(std::make_pair(S, B));
            // Storing significance in the output struct object
            if (iHisto == 0) {
                // Set a significance value for each corresponding pT,D bin
                outputStruct.hSignificance = new TH1D("hSignificance", "Estimated significance for each m_{invariant} distribution bin;p_{T, D^{0}} (GeV/#it{c});Significance = #frac{S}{#sqrt{S+B}}", ptDBinEdges.size() - 1, ptDBinEdges.data());
                outputStruct.hSignificance->SetBinContent(iHisto+1, significance);
                outputStruct.hSignificance->SetBinError(iHisto+1, 0);
            } else {
                outputStruct.hSignificance->SetBinContent(iHisto+1, significance);
                outputStruct.hSignificance->SetBinError(iHisto+1, 0);
            }
        }
        
        
    }

    // Final adjustments to subtracted histograms
    for (size_t iHisto = 0; iHisto < outputStruct.subtractedHist.size(); iHisto++) {
        
        // Check if histogram should be erased before storing based on fit performed
        bool eraseHist = eraseHistogram(fittings.workingFits, iHisto);

        if (eraseHist) {
            outputStruct.subtractedHist[iHisto]->Reset("ICES");
        } else {
            //
            // Set negative count bin entries to 0
            for (int iBin = 1; iBin <= outputStruct.subtractedHist[iHisto]->GetNbinsX(); iBin++) {
                if (outputStruct.subtractedHist[iHisto]->GetBinContent(iBin) < 0) {
                    outputStruct.subtractedHist[iHisto]->SetBinContent(iBin,0);
                    outputStruct.subtractedHist[iHisto]->SetBinError(iBin,0);
                }
            }

            
        }
        
        // Obtain summed histogram for all pT,D bins after sideband subtraction
        if (iHisto == 0) {
            outputStruct.hSubtracted_allPtSummed = (TH1D*)outputStruct.subtractedHist[iHisto]->Clone("hSubtracted_allPtSummed");
        } else {
            outputStruct.hSubtracted_allPtSummed->Add(outputStruct.subtractedHist[iHisto]);
        }
        
    }

    

    // Return the output struct object containing filled histogram vectors
    return outputStruct;
    
}

TCanvas* plotMCnet_template_mass(std::vector<TH1D*> originalHistograms, const FitContainer fittings, const std::vector<std::pair<double, double>> signal_background_values, int index, double jetptMin, double jetptMax) {
    
    TCanvas* cMCnetFormat = new TCanvas("cMCnetFormat", "MCnet_format", 800, 600);
    cMCnetFormat->cd();
    double statBoxPos = gPad->GetUxmax(); // Height of the stat box
    gStyle->SetOptStat(0); // Turn off the default stats box
    gStyle->SetTextFont(42); // Helvetica-like, cleaner than default

    TH1D* histogram = (TH1D*)originalHistograms[index]->Clone("hMCnetFormat");
    histogram->SetMarkerStyle(kFullCircle); // kFullCircle
    histogram->SetMarkerColor(kBlack);
    histogram->SetLineColor(kBlack);
    TString title = histogram->GetTitle();
    histogram->SetTitle("");
    histogram->GetYaxis()->SetTitle(Form("Counts per %.2f MeV/#it{c}^{2}", 1000*histogram->GetBinWidth(1)));
    histogram->GetXaxis()->SetTitleSize(0.04);
    histogram->GetYaxis()->SetTitleSize(0.04);
    histogram->GetYaxis()->SetTitleOffset(1.25); // default = 1.5
    histogram->SetMinimum(0);
    histogram->Sumw2(); // Enable error calculation
    histogram->Draw();
    histogram->SetMaximum(histogram->GetMaximum() * 1.6); // More vertical space
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    fittings.fitTotal[index]->SetLineStyle(1); // Solid
    fittings.fitTotal[index]->Draw("same");
    fittings.fitSignalOnly[index]->SetLineStyle(2); // Dashed
    fittings.fitSignalOnly[index]->Draw("same");
    fittings.fitBackgroundOnly[index]->SetLineStyle(3); // Dotted
    fittings.fitBackgroundOnly[index]->Draw("same");
    fittings.fitReflectionsOnly[index]->SetLineStyle(4); // Dash-dot
    fittings.fitReflectionsOnly[index]->Draw("same");
    double A1Signal = fittings.fitTotal[index]->GetParameter(2); // Get the value of parameter 'A' from primary signal gaussian
    double m0_1Signal = fittings.fitTotal[index]->GetParameter(4); // Get the value of parameter 'm0' from primary signal gaussian
    double sigma1Signal = fittings.fitTotal[index]->GetParameter(5); // Get the value of parameter 'sigma' from primary signal gaussian
    double A2Signal = A1Signal / fittings.fitTotal[index]->GetParameter(3); // Get the value of parameter 'A' from secondary signal gaussian
    double sigma2Signal = sigma1Signal / fittings.fitTotal[index]->GetParameter(6); // Get the value of parameter 'sigma' from secondary signal gaussian
    double chi2 = fittings.fitTotal[index]->GetChisquare();
    double degOfFreedom = fittings.fitTotal[index]->GetNDF();
    //std::cout << "signal_background_values vector size: " << signal_background_values.size() << std::endl;
    //double S = signal_background_values[index].first; // Signal value
    //double B = signal_background_values[index].second; // Background value
    //double significance = S / sqrt(S + B); // Significance value
    //std::cout << "S = " << S << ", B = " << B << ", significance = " << significance << std::endl;

    float topY = 0.84; // Top Y position for the text
    float deltaY = 0.06; // Vertical spacing between lines
    float middleX = 0.52; // Middle X position for the text
    float rightX = 0.87; // 0.895 or 0.90, depending on your layout
    float textSize = 0.04; // Text size for the labels

    // Allign text to the left
    TLatex* latexLeft = new TLatex();
    latexLeft->SetTextSizePixels(24);
    latexLeft->SetNDC(); // Enables normalized coordinates (0 to 1)
    latexLeft->SetTextSize(textSize); // default = 0.05
    latexLeft->DrawLatex(0.13, topY, "ALICE Performance");
    latexLeft->DrawLatex(0.13, topY-deltaY, "pp, #sqrt{#it{s}} = 13.6 TeV");

    // Allign text to the right
    TLatex* latexRight = new TLatex();
    latexRight->SetTextSizePixels(24);
    latexRight->SetNDC(); // Enables normalized coordinates (0 to 1)
    latexRight->SetTextSize(textSize); // default = 0.05
    latexRight->SetTextAlign(31); // Right-align text to that X position
    latexRight->DrawLatex(rightX, topY, "D^{0}#rightarrow K^{#minus}+#pi^{+} and ch. conj.");
    latexRight->DrawLatex(rightX, topY - deltaY, "in ch-particle jets, anti-#it{k}_{T}, #it{R} = 0.4");
    latexRight->DrawLatex(rightX, topY - 2*deltaY, title+", |#it{y}_{D^{0}}| #leq 0.8");
    //latexRight->DrawLatex(rightX, topY - 3*deltaY, Form("%.0f < #it{p}_{T, ch. jet} (GeV/#it{c}) < %.0f, |#it{#eta}_{jet}| #leq 0.5", jetptMin, jetptMax));// middle GeV/#it{c} - Vit
    latexRight->DrawLatex(rightX, topY - 3*deltaY, Form("%.0f < #it{p}_{T, ch. jet} < %.0f GeV/#it{c}, |#it{#eta}_{jet}| #leq 0.5", jetptMin, jetptMax)); // right GeV/#it{c} - Raymond
    //latexRight->DrawLatex(rightX, topY - 4*deltaY, Form("#mu_{1} = (%.2f #pm %.2f) MeV/#it{c}^{2}", 1000*m0_1Signal, 1000*fittings.fitTotal[index]->GetParError(4)));
    //latexRight->DrawLatex(rightX, topY - 5*deltaY, Form("#sigma_{1} = (%.2f #pm %.2f) MeV/#it{c}^{2}", 1000*sigma1Signal, 1000*fittings.fitTotal[index]->GetParError(5)));

    TLegend* legend = new TLegend(rightX - 0.3, (topY - 3.5*deltaY), rightX, (topY - 3.5*deltaY) - 0.21);
    legend->SetTextSize(0.04);
    //legend->SetTextAlign(12); // Left align text in the legend=12, Right align text in the legend=32
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    //legend->AddEntry(histogram, "Data", "lpe");
    legend->AddEntry(fittings.fitTotal[index], "Total fit funct.", "l");
    legend->AddEntry(fittings.fitSignalOnly[index], "Signal", "l");
    legend->AddEntry(fittings.fitBackgroundOnly[index], "Combinatorial bkg.", "l");
    legend->AddEntry(fittings.fitReflectionsOnly[index], "Reflection contrib.", "l");
    legend->Draw();

    
    cMCnetFormat->Update();

    return cMCnetFormat;
}

void PlotHistograms(const SubtractionResult& outputStruct, const std::vector<TH2D*>& histograms2d, const FitContainer& fittings, double jetptMin, double jetptMax) {
    std::cout << "Plotting histograms..." << std::endl;
    // creating 1D mass projection histograms
    TH1D* tempHist;
    std::vector<TH1D*> histograms;

    // obtaining 1D invariant mass histograms from the projection
    for (size_t iHist = 0; iHist < histograms2d.size(); iHist++) {
        //
        tempHist = histograms2d[iHist]->ProjectionX(Form("h_mass_proj_%zu", iHist));
        histograms.push_back(tempHist);
    }

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.04); // default = 0.05

    // Create a canvas for plotting
    int nHistos = histograms.size();
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
    for(size_t iHisto = 0; iHisto < histograms.size(); ++iHisto) {
        //
        c1d_fit->cd(iHisto+1);
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        histograms[iHisto]->SetMarkerStyle(kDot); //kFullDotMedium
        histograms[iHisto]->SetMarkerColor(kBlack);
        histograms[iHisto]->SetLineColor(kBlack);
        histograms[iHisto]->GetYaxis()->SetTitle("counts");
        histograms[iHisto]->SetMinimum(0);
        histograms[iHisto]->Draw();
        fittings.fitTotal[iHisto]->Draw("same");
        fittings.fitSignalOnly[iHisto]->Draw("same");
        fittings.fitBackgroundOnly[iHisto]->Draw("same");
        fittings.fitReflectionsOnly[iHisto]->Draw("same");

        double A1Signal = fittings.fitTotal[iHisto]->GetParameter(2); // Get the value of parameter 'A' from primary signal gaussian
        double m0_1Signal = fittings.fitTotal[iHisto]->GetParameter(4); // Get the value of parameter 'm0' from primary signal gaussian
        double sigma1Signal = fittings.fitTotal[iHisto]->GetParameter(5); // Get the value of parameter 'sigma' from primary signal gaussian
        double A2Signal = A1Signal / fittings.fitTotal[iHisto]->GetParameter(3); // Get the value of parameter 'A' from secondary signal gaussian
        double sigma2Signal = sigma1Signal / fittings.fitTotal[iHisto]->GetParameter(6); // Get the value of parameter 'sigma' from secondary signal gaussian
        double chi2 = fittings.fitTotal[iHisto]->GetChisquare();
        double degOfFreedom = fittings.fitTotal[iHisto]->GetNDF();
        // original position at statBoxPos-0.35, 0.70 with 0.03 of size
        //latex->DrawLatex(statBoxPos-0.3, 0.70, Form("m_{0} = %.3f #pm %.3f GeV/#it{c}^{2}", m_0,sigma)); // Display parameter 'm_0' value plus width of primary gaussian
        latex->DrawLatex(statBoxPos-0.52, 0.80, Form("A_{1}^{signal} = %.2f, #bar{m_{1}} = %.2f, #sigma_{1} = %.2f GeV/#it{c}^{2}", A1Signal, m0_1Signal,sigma1Signal));
        latex->DrawLatex(statBoxPos-0.52, 0.75, Form("A_{2}^{signal} = %.2f, #bar{m_{2}} = %.2f, #sigma_{2} = %.2f GeV/#it{c}^{2}", A2Signal, m0_1Signal,sigma2Signal));
        latex->DrawLatex(statBoxPos-0.3, 0.70, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax)); // Display jet pT cut applied
        latex->DrawLatex(statBoxPos-0.25, 0.63, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));
        latex->DrawLatex(statBoxPos-0.32, 0.58, "MC: JE_HF_LHC24g5_All_D0");
        latex->DrawLatex(statBoxPos-0.4, 0.53, "Data: LHC23_pass4_Thin_small_2P3PDstar");

        // Drawing 2D histograms
        c_2d->cd(iHisto+1);
        histograms2d[iHisto]->Draw("colz");

    }

    /*TCanvas* cMCnet1 = new TCanvas("cMCnet1", "MCnet1", 800, 600);
    cMCnet1->SetCanvasSize(1800,1000);
    cMCnet1->cd();
    double statBoxPos = gPad->GetUxmax(); // Height of the stat box
    gStyle->SetOptStat(0); // Turn off the default stats box
    histograms[7]->SetMarkerStyle(kDot); //kFullDotMedium
    histograms[7]->SetMarkerColor(kBlack);
    histograms[7]->SetLineColor(kBlack);
    histograms[7]->GetYaxis()->SetTitle("counts");
    histograms[7]->SetMinimum(0);
    histograms[7]->Draw();
    fittings.fitTotal[7]->Draw("same");
    fittings.fitSignalOnly[7]->Draw("same");
    fittings.fitBackgroundOnly[7]->Draw("same");
    fittings.fitReflectionsOnly[7]->Draw("same");
    double A1Signal = fittings.fitTotal[7]->GetParameter(2); // Get the value of parameter 'A' from primary signal gaussian
    double m0_1Signal = fittings.fitTotal[7]->GetParameter(4); // Get the value of parameter 'm0' from primary signal gaussian
    double sigma1Signal = fittings.fitTotal[7]->GetParameter(5); // Get the value of parameter 'sigma' from primary signal gaussian
    double A2Signal = A1Signal / fittings.fitTotal[7]->GetParameter(3); // Get the value of parameter 'A' from secondary signal gaussian
    double sigma2Signal = sigma1Signal / fittings.fitTotal[7]->GetParameter(6); // Get the value of parameter 'sigma' from secondary signal gaussian
    double chi2 = fittings.fitTotal[7]->GetChisquare();
    double degOfFreedom = fittings.fitTotal[7]->GetNDF();
    latex->DrawLatex(statBoxPos-0.52, 0.80, Form("A_{1}^{signal} = %.2f, #bar{m_{1}} = %.2f, #sigma_{1} = %.2f GeV/#it{c}^{2}", A1Signal, m0_1Signal,sigma1Signal));
    latex->DrawLatex(statBoxPos-0.52, 0.75, Form("A_{2}^{signal} = %.2f, #bar{m_{2}} = %.2f, #sigma_{2} = %.2f GeV/#it{c}^{2}", A2Signal, m0_1Signal,sigma2Signal));
    latex->DrawLatex(statBoxPos-0.3, 0.70, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax)); // Display jet pT cut applied
    latex->DrawLatex(statBoxPos-0.25, 0.63, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));
    latex->DrawLatex(statBoxPos-0.32, 0.58, "MC: JE_HF_LHC24g5_All_D0");
    latex->DrawLatex(statBoxPos-0.4, 0.53, "Data: LHC23_pass4_Thin_small_2P3PDstar");*/

    /*TCanvas* cMCnet2 = new TCanvas("cMCnet2", "MCnet2", 800, 600);
    cMCnet2->SetCanvasSize(1800,1000);
    cMCnet2->cd();
    statBoxPos = gPad->GetUxmax(); // Height of the stat box
    gStyle->SetOptStat(0); // Turn off the default stats box
    histograms[10]->SetMarkerStyle(kDot); //kFullDotMedium
    histograms[10]->SetMarkerColor(kBlack);
    histograms[10]->SetLineColor(kBlack);
    histograms[10]->GetYaxis()->SetTitle("counts");
    histograms[10]->SetMinimum(0);
    histograms[10]->Draw();
    fittings.fitTotal[10]->Draw("same");
    fittings.fitSignalOnly[10]->Draw("same");
    fittings.fitBackgroundOnly[10]->Draw("same");
    fittings.fitReflectionsOnly[10]->Draw("same");
    A1Signal = fittings.fitTotal[10]->GetParameter(2); // Get the value of parameter 'A' from primary signal gaussian
    m0_1Signal = fittings.fitTotal[10]->GetParameter(4); // Get the value of parameter 'm0' from primary signal gaussian
    sigma1Signal = fittings.fitTotal[10]->GetParameter(5); // Get the value of parameter 'sigma' from primary signal gaussian
    A2Signal = A1Signal / fittings.fitTotal[10]->GetParameter(3); // Get the value of parameter 'A' from secondary signal gaussian
    sigma2Signal = sigma1Signal / fittings.fitTotal[10]->GetParameter(6); // Get the value of parameter 'sigma' from secondary signal gaussian
    chi2 = fittings.fitTotal[10]->GetChisquare();
    degOfFreedom = fittings.fitTotal[10]->GetNDF();
    latex->DrawLatex(statBoxPos-0.52, 0.80, Form("A_{1}^{signal} = %.2f, #bar{m_{1}} = %.2f, #sigma_{1} = %.2f GeV/#it{c}^{2}", A1Signal, m0_1Signal,sigma1Signal));
    latex->DrawLatex(statBoxPos-0.52, 0.75, Form("A_{2}^{signal} = %.2f, #bar{m_{2}} = %.2f, #sigma_{2} = %.2f GeV/#it{c}^{2}", A2Signal, m0_1Signal,sigma2Signal));
    latex->DrawLatex(statBoxPos-0.3, 0.70, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax)); // Display jet pT cut applied
    latex->DrawLatex(statBoxPos-0.25, 0.63, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));
    latex->DrawLatex(statBoxPos-0.32, 0.58, "MC: JE_HF_LHC24g5_All_D0"); // previous dataset: LHC24d3a_All
    latex->DrawLatex(statBoxPos-0.4, 0.53, "Data: LHC23_pass4_Thin_small_2P3PDstar");*/

    //TCanvas* cMCnet_template_mass = plotMCnet_template_mass(histograms, fittings, outputStruct.signal_background_values, 3, jetptMin, jetptMax);
    //cMCnet_template_mass->Draw();

    // Plotting output observable
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
    legend->AddEntry(outputStruct.sidebandHist[0],"Sideband", "lpe");
    legend->AddEntry(outputStruct.signalHist[0],"Signal", "lpe");
    legend->AddEntry(outputStruct.subtractedHist[0],"Signal (minus background)", "lpe");
    std::cout << "Starting subtracted histograms plotting..." << std::endl;
    for (size_t iHisto = 0; iHisto < outputStruct.subtractedHist.size(); iHisto++) {
        std::cout << "Plotting side-band histogram number " << iHisto << std::endl;
        cSideBand->cd(iHisto+1);
        //outputStruct.sidebandHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        outputStruct.sidebandHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        outputStruct.sidebandHist[iHisto]->GetYaxis()->SetTitle("yields");
        outputStruct.sidebandHist[iHisto]->SetMarkerStyle(kFullTriangleUp);
        outputStruct.sidebandHist[iHisto]->SetMarkerColor(kAzure);
        outputStruct.sidebandHist[iHisto]->SetLineColor(kAzure);
        double statBoxPos = gPad->GetUxmax();
        outputStruct.sidebandHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));

        cSignal->cd(iHisto+1);
        //outputStruct.signalHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        outputStruct.signalHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        outputStruct.signalHist[iHisto]->GetYaxis()->SetTitle("yields");
        outputStruct.signalHist[iHisto]->SetMarkerStyle(kFullSquare);
        outputStruct.signalHist[iHisto]->SetMarkerColor(kViolet);
        outputStruct.signalHist[iHisto]->SetLineColor(kViolet);
        outputStruct.signalHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));
        std::cout << "Plotting subtracted histogram number " << iHisto << std::endl;
        cSubtracted->cd(iHisto+1);
        //outputStruct.subtractedHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        outputStruct.subtractedHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        outputStruct.subtractedHist[iHisto]->GetYaxis()->SetTitle("yields");
        outputStruct.subtractedHist[iHisto]->SetMarkerStyle(kFullCircle);
        outputStruct.subtractedHist[iHisto]->SetMarkerColor(kPink);
        outputStruct.subtractedHist[iHisto]->SetLineColor(kPink);
        outputStruct.subtractedHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));

        cSigPlusBack->cd(iHisto+1);
        //outputStruct.signalHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        outputStruct.signalHist[iHisto]->Draw();
        outputStruct.sidebandHist[iHisto]->Draw("same");
        outputStruct.subtractedHist[iHisto]->Draw("same");
        legend->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.5, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));
    }
    
    TCanvas* cTesting = new TCanvas("cTesting","cTesting");
    cTesting->cd();
    //outputStruct.subtractedHist[8]->Draw();
    std::cout << "Subtracted hist number 8 entries = " << outputStruct.subtractedHist[8]->GetEntries() << std::endl;

    TCanvas* cAllPt = new TCanvas("cAllPt","All pT,D final deltaR distribution");
    cAllPt->SetCanvasSize(1800,1000);
    cAllPt->Divide(1,2);
    cAllPt->cd(1);
    //outputStruct.hSubtracted_allPtSummed->GetXaxis()->SetRangeUser(0.0, 0.5);
    outputStruct.hSubtracted_allPtSummed->SetTitle(";#DeltaR;yields");
    outputStruct.hSubtracted_allPtSummed->SetMarkerStyle(kFullCircle);
    outputStruct.hSubtracted_allPtSummed->SetMarkerColor(kRed+1);
    outputStruct.hSubtracted_allPtSummed->SetLineColor(kRed+1);
    outputStruct.hSubtracted_allPtSummed->Draw();
    double statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.5, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax));
    cAllPt->cd(2);
    outputStruct.hSignificance->Draw();
    latex->DrawLatex(statBoxPos-0.35, 0.5, "Obtained through signal and background fits");
    latex->DrawLatex(statBoxPos-0.35, 0.6, "Calculation region: |m_{inv} - m_{D^{0}}| < 2#sigma");

    //
    // Storing images
    //
    TString imagePath = "../../Images/1-SignalTreatment/SideBand/";
    c1d_fit->Update();
    c1d_fit->SaveAs(imagePath + Form("SB_BackSub_invariant_mass_withReflections_%.0f_to_%.0fGeV_withReflections.png",jetptMin,jetptMax));
    //cMCnet_template_mass->Update();
    //cMCnet_template_mass->SaveAs(imagePath + "SB_BackSub_invariant_mass_withReflections_MCnet_template.eps");
    c_2d->Update();
    c_2d->SaveAs(imagePath + Form("SB_BackSub_2d_deltaR_vs_invmass_withReflections_%.0f_to_%.0fGeV_withReflections.png",jetptMin,jetptMax));
    cSigPlusBack->Update();
    cSigPlusBack->SaveAs(imagePath + Form("SB_BackSub_yield_pT_bins_withReflections_%.0f_to_%.0fGeV_withReflections.png",jetptMin,jetptMax));
    cAllPt->Update();
    cAllPt->SaveAs(imagePath + Form("SB_BackSub_yield_pT_summed_withReflections_%.0f_to_%.0fGeV_withReflections.png",jetptMin,jetptMax));

    //
    // Storing in a single pdf file
    //
    c1d_fit->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf(",jetptMin,jetptMax));
    c_2d->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf",jetptMin,jetptMax));
    //cMCnet_template_mass->Print(Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf",jetptMin,jetptMax));
    cSigPlusBack->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf",jetptMin,jetptMax));
    cAllPt->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf)",jetptMin,jetptMax));

    std::cout << "Histograms and fits plotted and stored to png and pdf files.\n";
}

void SaveData(SubtractionResult& outputStruct, double jetptMin, double jetptMax,
              const std::vector<double>& deltaRBinEdges, const std::vector<double>& ptjetBinEdges, const std::vector<double>& ptDBinEdges) {
    // Open output file
    TFile* fOutput = new TFile(Form("backSub_%d_to_%d_jetpt_with_reflections.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    if (!fOutput || fOutput->IsZombie()) {
        std::cerr << "Error: Unable to open the output ROOT file." << std::endl;
    }

    // Loop over background subtracted histograms
    for (size_t iHisto = 0; iHisto < outputStruct.subtractedHist.size(); iHisto++) {
        // store each histogram in file
        outputStruct.subtractedHist[iHisto]->Write();
    }
    // Significance as a function of pT,D mass histogram
    outputStruct.hSignificance->Write();
    // Final summed pT,D histogram
    outputStruct.hSubtracted_allPtSummed->Write();
    
    // Also store the axes used for the histograms
    // Create a directory for axes
    fOutput->mkdir("axes");
    fOutput->cd("axes");
    // Create TVectorD with same content
    TVectorD vecDeltaR(deltaRBinEdges.size());
    for (size_t i = 0; i < deltaRBinEdges.size(); ++i) {
        vecDeltaR[i] = deltaRBinEdges[i];
    }
    vecDeltaR.Write("deltaRBinEdges_detector");
    TVectorD vecPtJet(ptjetBinEdges.size());
    for (size_t i = 0; i < ptjetBinEdges.size(); ++i) {
        vecPtJet[i] = ptjetBinEdges[i];
    }
    vecPtJet.Write("ptjetBinEdges_detector");
    TVectorD vecPtD(ptDBinEdges.size());
    for (size_t i = 0; i < ptDBinEdges.size(); ++i) {
        vecPtD[i] = ptDBinEdges[i];
    }
    vecPtD.Write("ptDBinEdges_detector");
    // Return to root directory (optional)
    fOutput->cd();
    
    std::cout << "Data saved to file: " << fOutput->GetName() << std::endl;

    // Close output file
    fOutput->Close();
    delete fOutput;
}

void JetPtIterator(const double jetptMin, const double jetptMax, const std::vector<double>& ptjetBinEdges) {

    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.03 (from first of pair)
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {1, 0.03}, {2, 0.03}, {3, 0.05}, {4, 0.05}, {5, 0.08}, {6, 0.15}, {8, 0.22}, {12, 0.35}, {16, 0.47}, {24, 0.47}
    };
    // Dataset: JE_HF_LHC24g5_All_D0
    //std::vector<std::pair<double, double>> bdtPtCuts = {
    //    {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.4}, {8, 0.6}, {10, 0.8}, {12, 0.8}, {16, 1.0}, {50, 1.0}
    //};

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;

    // mass histogram
    int massBins = 50; // default=100 
    double minMass = 1.72; // use from 1.72, used to use 1.67
    double maxMass = 2.06; // old = 2.1, new = 2.05 + 0.01

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
    // Load ΔR bin edges
    std::vector<double> deltaRBinEdges = LoadBinning(fReflectionsMC, "axes/deltaRBinEdges");
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
    // Load pTD bin edges
    std::vector<double> ptDBinEdges = LoadBinning(fReflectionsMC, "axes/ptDBinEdges");

    // Create multiple histograms
    std::vector<TH2D*> histograms2d = createHistograms(ptDBinEdges,                                 // the pT,D edges will determine the number of mass histograms
                                                     massBins, minMass, maxMass,                  // mass histograms binning
                                                     deltaRBinEdges, jetptMin);                             // deltaR histograms with asymmetrical bin widths
                                                     //deltaRbins, minDeltaR, maxDeltaR);         // deltaR histograms
    
    // bin sizes
    cout << "Mass bin width = " << (maxMass-minMass)/massBins << endl;

    // Fill histograms
    fillHistograms(fDist, histograms2d, jetptMin, jetptMax, ptDBinEdges, deltaRBinEdges, bdtPtCuts);

    // Perform fits
    FitModelType modelToUse = FitModelType::Full; // options: Full, SignalReflectionsOnly, SignalOnly, ReflectionsOnly
    FitContainer fittings = performFit(fReflectionsMC, histograms2d, minMass, maxMass, modelToUse, jetptMin, jetptMax);

    // signal/side-band region parameters
    double signalSigmas = 2; // 2
    int startingBackSigma = 4; // 4
    int backgroundSigmas = 4; // 4
    cout << "Signal: |m - m_0| < " << signalSigmas << "sigmas" << endl;
    cout << "Left side-band: " << -(startingBackSigma+backgroundSigmas) << "sigmas << |m - m_0| << " << -startingBackSigma << "sigmas\n";
    cout << "Right side-band: " << startingBackSigma << "sigmas << |m - m_0| << " << (startingBackSigma+backgroundSigmas) << "sigmas\n";

    // Subtract side-band from signal
    SubtractionResult outputStruct = SideBand(histograms2d, fittings, signalSigmas, startingBackSigma, backgroundSigmas, ptDBinEdges);
    
    // Plot histograms
    PlotHistograms(outputStruct, histograms2d, fittings, jetptMin, jetptMax);

    // Storing final histograms to output file
    SaveData(outputStruct,jetptMin,jetptMax, deltaRBinEdges, ptjetBinEdges, ptDBinEdges);

}

void create3DBackgroundSubtracted(const std::vector<double>& ptjetBinEdges, const std::vector<double>& deltaRBinEdges, const std::vector<double>& ptDBinEdges){

    // jet pT cuts
    double jetptMin = ptjetBinEdges[0]; // GeV
    double jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 1]; // GeV
    // deltaR histogram
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];

    TH3D* h3D = new TH3D("h3DBackgroundSubtracted", "Background subtracted; p_{T,jet} (GeV/c); #DeltaR; p_{T,D^{0}} (GeV/c)",
                     ptjetBinEdges.size()-1, ptjetBinEdges.data(),
                     deltaRBinEdges.size()-1, deltaRBinEdges.data(),
                     ptDBinEdges.size()-1, ptDBinEdges.data());
    h3D->Sumw2();

    for (size_t iJetBin = 0; iJetBin < ptjetBinEdges.size() - 1; iJetBin++) {
        TFile* fJetRange = new TFile(Form("backSub_%0.f_to_%0.f_jetpt_with_reflections.root",ptjetBinEdges[iJetBin],ptjetBinEdges[iJetBin+1]),"read");
        if (!fJetRange || fJetRange->IsZombie()) {
            std::cerr << "Error opening file " << Form("backSub_%0.f_to_%0.f_jetpt_with_reflections.root",ptjetBinEdges[iJetBin],ptjetBinEdges[iJetBin+1]) << std::endl;
            continue;
        }

        for (int iHist = 0; iHist < 100; ++iHist) {
            TString histName = Form("h_back_subtracted_%d", iHist);
            TH1D* hDeltaR = (TH1D*)fJetRange->Get(histName);
            if (!hDeltaR) break; // No more histograms
            //std::cout << "Histogram x axis range: " << hDeltaR->GetXaxis()->GetXmin() << " to " << hDeltaR->GetXaxis()->GetXmax() << std::endl;

            // Get bin centers
            double ptJetCenter = 0.5 * (ptjetBinEdges[iJetBin] + ptjetBinEdges[iJetBin+1]);
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

            double ptDlow = -1, ptDhigh = -1;
            if (sscanf(title.Data(), "%lf < #it{p}_{T, D^{0}} < %lf GeV/#it{c}", &ptDlow, &ptDhigh) != 2) {
                std::cerr << "Could not parse pT,D range from histogram title: " << title << std::endl;
                std::cout << "Trying to parse title: '" << title << "'" << std::endl;
                std::cout << "hDeltaR histogram -> " << hDeltaR->GetName() << std::endl;
                std::cout << "In file " << fJetRange->GetName() << std::endl;
                continue;
            }
            double ptDCenter = 0.5 * (ptDlow + ptDhigh);


            for (int iBin = 1; iBin <= hDeltaR->GetNbinsX(); ++iBin) {
                double deltaRcenter = hDeltaR->GetBinCenter(iBin);
                double content = hDeltaR->GetBinContent(iBin);
                double error = hDeltaR->GetBinError(iBin);

                // Fill the TH3D using centers
                h3D->Fill(ptJetCenter, deltaRcenter, ptDCenter, content);
                // Optionally, if you want to preserve error propagation:
                int binX = h3D->GetXaxis()->FindBin(ptJetCenter);
                int binY = h3D->GetYaxis()->FindBin(deltaRcenter);
                int binZ = h3D->GetZaxis()->FindBin(ptDCenter);
                h3D->SetBinError(binX, binY, binZ, error); // only if needed
            }
        }

    }
    
    h3D->Draw("colz");

    TFile* fOut = new TFile(Form("full_merged_ranges_back_sub.root"), "RECREATE");
    h3D->Write();

    // Also store the axes used for the histograms
    // Create a directory for axes
    fOut->mkdir("axes");
    fOut->cd("axes");
    // Create TVectorD with same content
    TVectorD vecDeltaR(deltaRBinEdges.size());
    for (size_t i = 0; i < deltaRBinEdges.size(); ++i) {
        vecDeltaR[i] = deltaRBinEdges[i];
    }
    vecDeltaR.Write("deltaRBinEdges_detector");
    TVectorD vecPtJet(ptjetBinEdges.size());
    for (size_t i = 0; i < ptjetBinEdges.size(); ++i) {
        vecPtJet[i] = ptjetBinEdges[i];
    }
    vecPtJet.Write("ptjetBinEdges_detector");
    TVectorD vecPtD(ptDBinEdges.size());
    for (size_t i = 0; i < ptDBinEdges.size(); ++i) {
        vecPtD[i] = ptDBinEdges[i];
    }
    vecPtD.Write("ptDBinEdges_detector");

    // Return to root directory (optional)
    fOut->cd();

    //fOut->Close();
    std::cout << "3D histogram created and saved to " << fOut->GetName() << " with 3D histogram." << std::endl;

}

void BackgroundSubtraction_with_reflections() {
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // Opening
    double jetptMin = 5.; // GeV
    double jetptMax = 7.; // GeV
    TFile* fReflectionsMC = new TFile(Form("../Reflections/reflections_%.0f_to_%.0fGeV.root",jetptMin,jetptMax),"read");
    if (!fReflectionsMC || fReflectionsMC->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT reflections file." << std::endl;
    }
    // Load pTjet bin edges
    std::vector<double> ptjetBinEdges = LoadBinning(fReflectionsMC, "axes/ptjetBinEdges");
    //std::vector<double> ptjetBinEdges = {30., 50.};
    // Load ΔR bin edges
    std::vector<double> deltaRBinEdges = LoadBinning(fReflectionsMC, "axes/deltaRBinEdges");
    // Load pTD bin edges
    std::vector<double> ptDBinEdges = LoadBinning(fReflectionsMC, "axes/ptDBinEdges");
    fReflectionsMC->Close();
    
    

    for (size_t iJetPt = 0; iJetPt < ptjetBinEdges.size() - 1; iJetPt++) {
        // Apply side-band method to each pT,jet bin
        std::cout << "Processing pT,jet bin: " << ptjetBinEdges[iJetPt] << " to " << ptjetBinEdges[iJetPt+1] << " GeV/c" << std::endl;
        jetptMin = ptjetBinEdges[iJetPt];
        jetptMax = ptjetBinEdges[iJetPt+1];

        JetPtIterator(jetptMin, jetptMax, ptjetBinEdges);
    }
    // Compute the entire range too
    jetptMin = ptjetBinEdges[0];
    jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 1];
    JetPtIterator(jetptMin, jetptMax, ptjetBinEdges);

    // Create 3D final histogram with pT,jet vs DeltaR vs pT,D
    create3DBackgroundSubtracted(ptjetBinEdges,deltaRBinEdges,ptDBinEdges);

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl;
}

int main(){
    BackgroundSubtraction_with_reflections();
    return 0;
}