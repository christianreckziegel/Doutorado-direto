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

//__________________________________________________________________________________________________________________________
// Fit functions

// Custom background fit function
Double_t backgroundFunction(Double_t* x, Double_t* par) {
    Double_t m = x[0];
    Double_t a = par[0];
    Double_t b = par[1];

    // Defining the custom function
    Double_t result = a * TMath::Power(m, b);
    return result;
}
// Custom signal fit function
Double_t signalFunction(Double_t* x, Double_t* par) {
    Double_t m = x[0];
    Double_t C = par[0];
    Double_t m0 = par[1];
    Double_t sigma = par[2];

    // Defining the custom function
    Double_t result = C * TMath::Exp(-TMath::Power((m - m0) / (2 * sigma), 2));
    return result;
}
// Custom sum of background and signal fit function
Double_t customFitFunction(Double_t* x, Double_t* par) {
    return backgroundFunction(x,par) + signalFunction(x,&par[2]);
}

//__________________________________________________________________________________________________________________________
// Module to create TH2D histograms including interest variable: UNIFORM bin sizes
std::vector<TH2D*> createHistograms(const std::vector<double>& ptDBinEdges, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax) {
    std::vector<TH2D*> histograms;
    for (size_t i = 0; i < ptDBinEdges.size() - 1; ++i) {
        histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, ybins, ymin, ymax));
        histograms[i]->Sumw2();
        
    }
    return histograms;
}
// Overloading function
// Module to create TH2D histograms including interest variable: VARIABLE bin sizes
std::vector<TH2D*> createHistograms(const std::vector<double>& ptDBinEdges, int xbins, double xmin, double xmax, const std::vector<double>& yBinEdges) {
    std::vector<TH2D*> histograms;
    for (size_t i = 0; i < ptDBinEdges.size() - 1; ++i) {
        histograms.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < p_{T,D} < %.0f GeV/c;m(K#pi) GeV/c^{2};#DeltaR",ptDBinEdges[i],ptDBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        histograms[i]->Sumw2();
        
    }
    return histograms;
}
//__________________________________________________________________________________________________________________________

// Module to fill 2D histograms from TTree data
void fillHistograms(TFile* fDist, const std::vector<TH2D*>& histograms, double jetptMin, double jetptMax, std::vector<double>& ptDBinEdges) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = 0.4; // on particle level delta R

    // Accessing TTree
    TTree* tree = (TTree*)fDist->Get("DF_2261906078621696/O2jetdisttable");

    // Assuming histograms and tree data correspond in some way
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing the TTree
    float axisDistance, jetPt, jetEta, jetPhi;
    float hfPt, hfEta, hfPhi, hfMass, hfY;

    tree->SetBranchAddress("fJetHfDist",&axisDistance);
    tree->SetBranchAddress("fJetPt",&jetPt);
    tree->SetBranchAddress("fJetEta",&jetEta);
    tree->SetBranchAddress("fJetPhi",&jetPhi);
    tree->SetBranchAddress("fHfPt",&hfPt);
    tree->SetBranchAddress("fHfEta",&hfEta);
    tree->SetBranchAddress("fHfPhi",&hfPhi);
    tree->SetBranchAddress("fHfMass",&hfMass);
    tree->SetBranchAddress("fHfY",&hfY);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(DeltaPhi(jetPhi,hfPhi),2));
        
        // Fill each histogram with their respective pT intervals
        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= 0.) && (deltaR < DeltaRcut))) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < ptDBinEdges.size() - 1 && !filled; iEdge++) {
                if ((hfPt >= ptDBinEdges[iEdge]) && (hfPt < ptDBinEdges[iEdge + 1])) {
                    histograms[iEdge]->Fill(hfMass, deltaR);
                    filled = true; // Exit the loop once the correct histogram is found
                }
                
            }
            
        }

        
        
        
    }
    cout << "Histograms filled.\n";
}

struct InitialParam {
    std::vector<double> paramA; // initial parameter values for parameter a
    std::vector<double> paramB; // initial parameter values for parameter b
    std::vector<double> paramC; // initial parameter values for parameter c
    std::vector<double> paramM0; // initial parameter values for parameter m_0
    std::vector<double> paramSigma; // initial parameter values for parameter sigma
};

// Find best combination of initial parameter values for each histogram fit
std::vector<double> bestFit(TH1D* histogram, double minMass, double maxMass, int numParameters){
    // Create array of best parameter initial values
    std::vector<double> optimalParameters(numParameters);

    // Create TF1 pointer object to be re-used
    TF1* fTestFit = new TF1("fTestFit",customFitFunction, minMass, maxMass, numParameters);

    // Initialize with infinity, for comparison
    double bestChiSquare = std::numeric_limits<double>::infinity();

    // calculating initial parameters
    double b_parameter;
    double a_parameter;
    double sigma_parameter = 0.012;
    double m_0_parameter = 1.86484; // D0 mass in GeV/c^2
    double C_parameter;

    // variables necessary for calculating parameters
    double K = histogram->GetBinContent(1); // Content of the first bin
    double J = histogram->GetBinContent(histogram->GetNbinsX()); // Content of the last bin
    double m_min = histogram->GetBinLowEdge(1); // Lower limit of the first bin
    double m_max = histogram->GetBinLowEdge(histogram->GetNbinsX() + 1); // Upper limit of the last bin
    double I_tot = histogram->Integral();
    double I_3sigma = histogram->Integral(histogram->FindBin(m_0_parameter - 3 * sigma_parameter), histogram->FindBin(m_0_parameter + 3 * sigma_parameter));
    
    // Calculating parameter initial values
    b_parameter = (K > 0 && J > 0) ? TMath::Log(K / J) / (m_max - m_min) : -1;
    a_parameter = I_tot / (TMath::Power(m_max, b_parameter + 1) - TMath::Power(m_min, b_parameter + 1)) / (b_parameter + 1); // a(b)
    C_parameter = (I_3sigma - a_parameter * (TMath::Power(m_0_parameter + 3 * sigma_parameter, b_parameter + 1) - TMath::Power(m_0_parameter - 3 * sigma_parameter, b_parameter + 1)) / (b_parameter + 1)) / (TMath::Sqrt(2 * TMath::Pi()) * sigma_parameter); // C(a,b)

    double rangeFactor = 0.3; // 0.1=10% around the initially calculated values, default = 0.3
    int stepsNumber = 10; // default = 10
    // loop over b parameter range
    for (double iB = b_parameter*(1-rangeFactor); iB < b_parameter*(1+rangeFactor); iB+=(b_parameter*2*rangeFactor/stepsNumber)) {
        //cout << "iB = " << iB << endl;
        // loop over a parameter range
        for (double iA = a_parameter*(1-rangeFactor); iA < a_parameter*(1+rangeFactor); iA+=(a_parameter*2*rangeFactor/stepsNumber)) {
            //cout << "iA = " << iA << endl;
            // loop over C parameter range
            for (double iC = C_parameter*(1-rangeFactor); iC < C_parameter*(1+rangeFactor); iC+=(C_parameter*2*rangeFactor/stepsNumber)) {
                //cout << "iC = " << iC << endl;
                for (double iS = sigma_parameter*(1-rangeFactor); iS < sigma_parameter*(1+rangeFactor*10); iS+=(sigma_parameter*2*rangeFactor/(stepsNumber))) {
                    // Perform the fit
                    fTestFit = new TF1("fTestFit",customFitFunction, minMass, maxMass, numParameters);
                    double m_0_lower_limit = 0.95 * m_0_parameter;
                    double m_0_upper_limit = 1.95 * m_0_parameter;
                    //fTestFit->SetParLimits(2, 0., DBL_MAX); // Set lower boundary of parameter iC to 0, and higher to maximum representable value for a double -> no negative gaussians
                    //fTestFit->SetParLimits(2, 0., std::numeric_limits<double>::infinity());
                    fTestFit->SetParLimits(3, m_0_lower_limit, m_0_upper_limit); // constraints on m_0_parameter: [0.95*m0, 1.05*m0]
                    fTestFit->SetParLimits(4, 0, 1); // constraints on sigma: [0, 1]
                    fTestFit->SetParameters(iA, iB, iC, m_0_parameter, iS);
                    fTestFit->SetParNames("a","b", "C", "m_0", "sigma");
                    fTestFit->SetLineColor(kBlue);
                    histogram->Fit(fTestFit, "Q");// "Q" option performs quiet fit without drawing the fit function

                    // Verify fit quality with current parameter values
                    if ((fTestFit->GetChisquare() < bestChiSquare) && (fTestFit->GetParameter(2) > 0)) { // (fTestFit->GetChisquare() < bestChiSquare) && (fTestFit->GetParameter(2) > 0)
                        bestChiSquare = fTestFit->GetChisquare();
                        optimalParameters[0] = iA;
                        optimalParameters[1] = iB;
                        optimalParameters[2] = iC;
                        optimalParameters[3] = m_0_parameter;
                        optimalParameters[4] = iS;
                    }

                    delete fTestFit;
                }
                
                


            }
        }
    }
    
    cout << "Optimal initial fit parameters values found = [" << optimalParameters[0] << "," << optimalParameters[1] << "," << optimalParameters[2] << "," << optimalParameters[3] << "," << optimalParameters[4] << "]\n";


    return optimalParameters;

}


/**
 * @brief Perform fits.
 *
 * This function takes the 2D histograms, perfom the fit of the signal+background model
 * and store the n TF1 signal+background model fit objects and n TF1 background model 
 * fit objects .
 * 
 * @param[in] histograms Histograms vector
 * @return Vector of fit objects TF1.
 */
// Module to perform fits to histograms
std::vector<TF1*> performFit(const std::vector<TH2D*>& histograms2d, InitialParam parametersVectors, double minMass, double maxMass) {
    
    TH1D* tempHist;
    // creating 1D mass projection histograms
    std::vector<TH1D*> histograms;

    // obtaining 1D invariant mass histograms from the projection
    for (size_t iHist = 0; iHist < histograms2d.size(); iHist++) {
        //
        tempHist = histograms2d[iHist]->ProjectionX(Form("h_mass_proj_%zu", iHist));
        histograms.push_back(tempHist);
    }
    
    
    //
    std::vector<TF1*> fittings;
    
    // calculating initial parameters
    double b_parameter;
    double a_parameter;
    double sigma_parameter = 0.012;
    double m_0_parameter = 1.86484; // D0 mass in GeV/c^2
    double C_parameter;
    
    // Perform total fit to each histogram
    for (size_t iHisto = 0; iHisto < histograms.size(); ++iHisto) {
        

        // variables necessary for calculating parameters
        double K = histograms[iHisto]->GetBinContent(1); // Content of the first bin
        double J = histograms[iHisto]->GetBinContent(histograms[iHisto]->GetNbinsX()); // Content of the last bin
        double m_min = histograms[iHisto]->GetBinLowEdge(1); // Lower limit of the first bin
        double m_max = histograms[iHisto]->GetBinLowEdge(histograms[iHisto]->GetNbinsX() + 1); // Upper limit of the last bin
        double I_tot = histograms[iHisto]->Integral();
        double I_3sigma = histograms[iHisto]->Integral(histograms[iHisto]->FindBin(m_0_parameter - 3 * sigma_parameter), histograms[iHisto]->FindBin(m_0_parameter + 3 * sigma_parameter));
        
        // Calculating parameter initial values
        b_parameter = (K > 0 && J > 0) ? TMath::Log(K / J) / (m_max - m_min) : -1;
        a_parameter = I_tot / (TMath::Power(m_max, b_parameter + 1) - TMath::Power(m_min, b_parameter + 1)) / (b_parameter + 1); // a(b)
        C_parameter = (I_3sigma - a_parameter * (TMath::Power(m_0_parameter + 3 * sigma_parameter, b_parameter + 1) - TMath::Power(m_0_parameter - 3 * sigma_parameter, b_parameter + 1)) / (b_parameter + 1)) / (TMath::Sqrt(2 * TMath::Pi()) * sigma_parameter); // C(a,b)

        bool printParam = false;
        if (printParam) {
            cout << "Initial parameter values:\n"
                 << "b_parameter = " << b_parameter << endl
                 << "a_parameter = " << a_parameter << endl
                 << "C_parameter = " << C_parameter << endl << endl;
        }
        
        

        // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
        fittings.push_back(new TF1(Form("totalFit_histMass%zu", iHisto + 1), customFitFunction, minMass, maxMass, 5)); // 5 parameters

        bool usePreSetParam = false; // use given pre-set initial parameter value
        if (usePreSetParam) {
            a_parameter = parametersVectors.paramA[iHisto];
            b_parameter = parametersVectors.paramB[iHisto];
            C_parameter = parametersVectors.paramC[iHisto];
            m_0_parameter = parametersVectors.paramM0[iHisto];
            sigma_parameter = parametersVectors.paramSigma[iHisto];
            cout << "Pre set parameters used = [" << parametersVectors.paramA[iHisto] << "," << parametersVectors.paramB[iHisto] << "," << parametersVectors.paramC[iHisto] << "," << parametersVectors.paramM0[iHisto] << "," << parametersVectors.paramSigma[iHisto] << "]\n";

        }

        // Calculating optimal parameters
        std::vector<double> optimalParameters = bestFit(histograms[iHisto], minMass, maxMass, 5);
        double m_0_lower_limit = 0.95 * m_0_parameter;
        double m_0_upper_limit = 1.95 * m_0_parameter;
        //fittings[iHisto]->SetParLimits(2, 0., DBL_MAX); // Set lower boundary of parameter iC to 0, and higher to maximum representable value for a double -> no negative gaussians
        fittings[iHisto]->SetParLimits(3, m_0_lower_limit, m_0_upper_limit); // constraints on m_0_parameter: [0.95*m0, 1.05*m0]
        fittings[iHisto]->SetParLimits(4, 0, 1); // constraints on sigma: [0, 1]
        //fittings[iHisto]->SetParameters(a_parameter, b_parameter, C_parameter, m_0_parameter, sigma_parameter);
        fittings[iHisto]->SetParameters(optimalParameters[0], optimalParameters[1], optimalParameters[2], optimalParameters[3], optimalParameters[4]);
        fittings[iHisto]->SetParNames("a","b", "C", "m_0", "sigma");
        fittings[iHisto]->SetLineColor(kBlue);
        histograms[iHisto]->Fit(fittings[iHisto], "Q");// "Q" option performs quiet fit without drawing the fit function

    }

    // Perform background only fit to each histogram
    int numberOfFits = fittings.size();
    for (size_t iHisto = 0; iHisto < numberOfFits; iHisto++) {
        // Getting total parameter values
        double a_par = fittings[iHisto]->GetParameter(0); // Get the value of parameter 'a'
        double b_par = fittings[iHisto]->GetParameter(1); // Get the value of parameter 'b'
        fittings.push_back(new TF1(Form("backFit_histMass%zu", iHisto), backgroundFunction, minMass, maxMass, 2));
        fittings[fittings.size()-1]->SetParameters(a_par,b_par); // fittings.size()-1 = the latest added to the vector
        fittings[fittings.size()-1]->SetLineStyle(kDashed);
        fittings[fittings.size()-1]->SetLineColor(kGreen);

    }
    
    cout << "Fit performed.\n";
    return fittings;
    // first n fittings are the total fit functions (from 0 to n-1)
    // second n fittings are the background fit functions (n to 2n-1)
}

struct SubtractionResult {
    std::vector<TH1D*> histograms;      // 1D mass projection histograms vector
    std::vector<TH1D*> sidebandHist;    // 1D sideband background histograms vector
    std::vector<TH1D*> signalHist;      // 1D signal histograms vector
    std::vector<TH1D*> subtractedHist;  // 1D sideband subtracted histograms vector
    TH1D* hSubtracted_allPtSummed;      // final deltaR distribution for all pT,D summed
};

SubtractionResult SideBand(const std::vector<TH2D*>& histograms2d, const std::vector<TF1*>& fittings, int signalSigmas, int startingBackSigma, int backgroundSigmas){
    
    // Creating histograms for collecting data
    TH1D* tempHist; // temporary histogram for collecting data
    TH1D* h_sideBand;
    TH1D* h_signal;
    TH1D* h_back_subtracted;

    // Creating output struct object with histogram vectors
    SubtractionResult vectorOutputs;

    // obtaining 1D invariant mass histograms from the projection
    for (size_t iHisto = 0; iHisto < histograms2d.size(); iHisto++) {
        //
        tempHist = histograms2d[iHisto]->ProjectionX(Form("h_mass_proj_%zu", iHisto));
        vectorOutputs.histograms.push_back(tempHist);
    }

    // Calculating scaling parameter
    for (size_t iHisto = 0; iHisto < vectorOutputs.histograms.size(); ++iHisto) {
        
        // Corresponding background fit index
        int backgroundHist = vectorOutputs.histograms.size() + iHisto;

        // Getting fit parameters
        double m_0 = fittings[iHisto]->GetParameter(3); // Get the value of parameter 'm_0'
        double sigma = fittings[iHisto]->GetParameter(4); // Get the value of parameter 'sigma'

        // Check which sides should be used for the total side-band distribution (if at least 1 sigma fit inside left range)
        if (m_0 - (startingBackSigma+1)*sigma <  vectorOutputs.histograms[iHisto]->GetBinLowEdge(1)) {
            std::cout << "Using only right sideband." << std::endl;

            // Count for how many sigmas is there room inside the left side range
            double leftSidebandRange = (m_0 - startingBackSigma * sigma) - vectorOutputs.histograms[iHisto]->GetBinLowEdge(1);
            int leftSigmas = static_cast<int>(leftSidebandRange / sigma);// get the integral number
            std::cout << leftSigmas << " sigmas fit inside the left side range for the background distribution estimation." << std::endl;
            // Count for how many sigmas is there room inside the right side range
            double rightSidebandRange = vectorOutputs.histograms[iHisto]->GetBinLowEdge(vectorOutputs.histograms[iHisto]->GetNbinsX()+1) - (m_0 + startingBackSigma*sigma);
            int rightSigmas = static_cast<int>(rightSidebandRange / sigma);// get the integral number
            std::cout << rightSigmas << " sigmas fit inside the right side range for the background distribution estimation." << std::endl;

            // Ensure only a maximum of 4 sigmas (backgroundSigmas) are used
            if (leftSigmas > backgroundSigmas) {
                leftSigmas = backgroundSigmas;
            }
            if (rightSigmas > backgroundSigmas) {
                rightSigmas = backgroundSigmas;
            }

            cout << "Left extreme = " << m_0 - (startingBackSigma + leftSigmas) * sigma << endl;
            cout << "Right extreme = " << m_0 + (startingBackSigma + rightSigmas) * sigma << endl;
            
            // Finding side-band region: 4*sigma < side-band < 8*sigma
            double rightSB = fittings[backgroundHist]->Integral(m_0 + startingBackSigma * sigma,m_0 + (startingBackSigma + rightSigmas) * sigma);
            // Finding side-band region: signal < |2*sigma|
            double backSigIntegral = fittings[backgroundHist]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of background in signal region
            double fitSigIntegral = fittings[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

            // Calculating scaling parameter
            double alpha = backSigIntegral / rightSB;
            std::cout << "alpha = " << alpha << std::endl;

            // Create side-band histogram: use only the right sideband
            int lowBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 + startingBackSigma * sigma);
            int highBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 + (startingBackSigma + rightSigmas) * sigma);
            h_sideBand = histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_%zu",iHisto), lowBin, highBin); // sum the right sideband 

            // Create signal histogram
            lowBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 - signalSigmas * sigma);
            highBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 + signalSigmas * sigma);
            h_signal = histograms2d[iHisto]->ProjectionY(Form("h_signal_proj_%zu",iHisto), lowBin, highBin);
            vectorOutputs.signalHist.push_back(h_signal); 

            cout << "Distribution ratio of sideband over raw signal: " << h_sideBand->Integral() / h_signal->Integral() << endl;
            cout << "Fit function ratio of sideband over raw signal: " << (rightSB)/ fitSigIntegral << endl;

            // Apply scaling to background function
            h_sideBand->Scale(alpha);
            vectorOutputs.sidebandHist.push_back(h_sideBand);

            // Subtract background histogram from signal histogram
            h_back_subtracted = (TH1D*)h_signal->Clone(Form("h_back_subtracted_%zu",iHisto));
            h_back_subtracted->Add(h_sideBand,-1.0);

            // Account for two sigma only area used for signal region
            h_back_subtracted->Scale(1/0.9545);
            vectorOutputs.subtractedHist.push_back(h_back_subtracted);

        } else {
            std::cout << "Using both left and right sidebands.\n";

            // Count for how many sigmas is there room inside the left side range
            double leftSidebandRange = (m_0 - startingBackSigma * sigma) - vectorOutputs.histograms[iHisto]->GetBinLowEdge(1);
            int leftSigmas = static_cast<int>(leftSidebandRange / sigma);// get the integral number
            std::cout << leftSigmas << " sigmas fit inside the left side range for the background distribution estimation." << std::endl;
            // Count for how many sigmas is there room inside the right side range
            double rightSidebandRange = vectorOutputs.histograms[iHisto]->GetBinLowEdge(vectorOutputs.histograms[iHisto]->GetNbinsX()+1) - (m_0 + startingBackSigma*sigma);
            int rightSigmas = static_cast<int>(rightSidebandRange / sigma);// get the integral number
            std::cout << rightSigmas << " sigmas fit inside the right side range for the background distribution estimation." << std::endl;

            // Ensure only a maximum of 4 sigmas (backgroundSigmas) are used
            if (leftSigmas > backgroundSigmas) {
                leftSigmas = backgroundSigmas;
            }
            if (rightSigmas > backgroundSigmas) {
                rightSigmas = backgroundSigmas;
            }

            cout << "Left extreme = " << m_0 - (startingBackSigma + leftSigmas) * sigma << endl;
            cout << "Right extreme = " << m_0 + (startingBackSigma + rightSigmas) * sigma << endl;

            // Finding side-band region: 4*sigma < side-band < 8*sigma
            double leftSB = fittings[backgroundHist]->Integral(m_0 - (startingBackSigma + leftSigmas) * sigma,m_0 - startingBackSigma * sigma);
            double rightSB = fittings[backgroundHist]->Integral(m_0 + startingBackSigma * sigma,m_0 + (startingBackSigma + rightSigmas) * sigma);
            // Finding side-band region: signal < |2*sigma|
            double backSigIntegral = fittings[backgroundHist]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of background in signal region
            double fitSigIntegral = fittings[iHisto]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

            // Calculating scaling parameter
            double alpha = backSigIntegral/(leftSB+rightSB);
            std::cout << "alpha = " << alpha << std::endl;
        
            // Create side-band histogram
            int lowBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 - (startingBackSigma + leftSigmas) * sigma);
            int highBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 - startingBackSigma * sigma);
            h_sideBand = histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_%zu",iHisto),lowBin, highBin); // start with left sideband
            lowBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 + startingBackSigma * sigma);
            highBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 + (startingBackSigma + rightSigmas) * sigma);
            TH1D* tempHist = histograms2d[iHisto]->ProjectionY(Form("h_sideband_proj_temp_%zu",iHisto), lowBin, highBin); // sum the right sideband
            h_sideBand->Add(tempHist);

            // Apply scaling to background function
            //h_sideBand->Scale(alpha);
            //vectorOutputs.sidebandHist.push_back(h_sideBand);

            // Create signal histogram
            lowBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 - signalSigmas * sigma);
            highBin = histograms2d[iHisto]->GetXaxis()->FindBin(m_0 + signalSigmas * sigma);
            h_signal = histograms2d[iHisto]->ProjectionY(Form("h_signal_proj_%zu",iHisto), lowBin, highBin);
            vectorOutputs.signalHist.push_back(h_signal); 

            cout << "Distribution ratio of sideband over raw signal: " << h_sideBand->Integral() / h_signal->Integral() << endl;
            cout << "Fit function ratio of sideband over raw signal: " << (leftSB+rightSB)/ fitSigIntegral << endl;

            // Apply scaling to background function
            h_sideBand->Scale(alpha);
            vectorOutputs.sidebandHist.push_back(h_sideBand);

            // Subtract background histogram from signal histogram
            h_back_subtracted = (TH1D*)h_signal->Clone(Form("h_back_subtracted_%zu",iHisto));
            h_back_subtracted->Add(h_sideBand,-1.0);

            // Account for two sigma only area used for signal region
            h_back_subtracted->Scale(1/0.9545);
            vectorOutputs.subtractedHist.push_back(h_back_subtracted);
            
        }

        // Add to final all pT,D summed histogram
        if (iHisto == 0) {
            vectorOutputs.hSubtracted_allPtSummed = (TH1D*)h_back_subtracted->Clone("hSubtracted_allPtSummed_SB");
        } else {
            vectorOutputs.hSubtracted_allPtSummed->Add(h_back_subtracted);
        }
        
        
    }

    // Return the output struct object containing filled histogram vectors
    return vectorOutputs;
    
}

void PlotHistograms(const std::vector<TH2D*>& histograms2d, const std::vector<TF1*>& fittings, SubtractionResult outputStruct, double jetptMin, double jetptMax) {
    
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
    latex->SetTextSize(0.05);

    // Create a canvas for plotting
    TCanvas* c1d_fit = new TCanvas("c1d_fit", "1D histograms with Fit", 800, 600);
    c1d_fit->SetCanvasSize(1800,1000);
    c1d_fit->Divide(3,static_cast<int>(histograms.size() / 3)); // columns, lines
    TCanvas* c_2d = new TCanvas("c_2d", "2D histograms", 800, 600);
    c_2d->SetCanvasSize(1800,1000);
    c_2d->Divide(3,static_cast<int>(histograms.size() / 3)); // columns, lines

    // Loop through all histograms and fitting functions
    for(size_t iHisto = 0; iHisto < histograms.size(); ++iHisto) {
        //
        c1d_fit->cd(iHisto+1);
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        histograms[iHisto]->SetMarkerStyle(kFullDotMedium);
        histograms[iHisto]->SetMarkerColor(kBlack);
        histograms[iHisto]->SetLineColor(kGray);
        histograms[iHisto]->GetYaxis()->SetTitle("counts");
        histograms[iHisto]->Draw();
        fittings[iHisto]->Draw("same");
        int backgroundHist = histograms.size() + iHisto;
        fittings[backgroundHist]->Draw("same");
        double m_0 = fittings[iHisto]->GetParameter(3); // Get the value of parameter 'm_0'
        double sigma = fittings[iHisto]->GetParameter(4); // Get the value of parameter 'sigma'
        double chi2 = fittings[iHisto]->GetChisquare();
        double degOfFreedom = fittings[iHisto]->GetNDF();
        // original position at statBoxPos-0.35, 0.70 with 0.03 of size
        latex->DrawLatex(statBoxPos-0.3, 0.70, Form("m_{0} = %.3f #pm %.3f GeV/c^{2}", m_0,sigma)); // Display parameter 'm_0' value
        latex->DrawLatex(statBoxPos-0.3, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax)); // Display jet pT cut applied
        latex->DrawLatex(statBoxPos-0.3, 0.58, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));

        // Drawing 2D histograms
        c_2d->cd(iHisto+1);
        histograms2d[iHisto]->GetYaxis()->SetTitle("#DeltaR");
        histograms2d[iHisto]->Draw("colz");

    }

    // Plot fitting functions only
    TCanvas* cFitsOnly = new TCanvas("cFitsOnly","Fit functions");
    cFitsOnly->SetCanvasSize(1800,1000);
    cFitsOnly->Divide(3,static_cast<int>(histograms.size() / 3)); // columns, lines
    for (size_t iHisto = 0; iHisto < histograms.size(); iHisto++) {
        cFitsOnly->cd(iHisto+1);
        int backgroundHist = histograms.size() + iHisto;
        fittings[iHisto]->Draw();
    }
    

    // Plotting output observable
    TCanvas* cSideBand = new TCanvas("cSideBand", "delta R for side-band", 800, 600);
    cSideBand->SetCanvasSize(1800,1000);
    cSideBand->Divide(3,static_cast<int>(outputStruct.histograms.size() / 3)); // columns, lines
    TCanvas* cSignal = new TCanvas("cSignal", "delta R for signal", 800, 600);
    cSignal->SetCanvasSize(1800,1000);
    cSignal->Divide(3,static_cast<int>(outputStruct.histograms.size() / 3)); // columns, lines
    TCanvas* cSubtracted = new TCanvas("cSubtracted", "delta R for side-band subtracted signal", 800, 600);
    cSubtracted->SetCanvasSize(1800,1000);
    cSubtracted->Divide(3,static_cast<int>(outputStruct.histograms.size() / 3)); // columns, lines
    TCanvas* cSigPlusBack = new TCanvas("cSigPlusBack", "delta R for side-band and signal in the same plot", 800, 600);
    cSigPlusBack->SetCanvasSize(1800,1000);
    cSigPlusBack->Divide(3,static_cast<int>(outputStruct.histograms.size() / 3)); // columns, lines

    TLegend* legend = new TLegend(0.6,0.57,0.9,0.77);
    legend->AddEntry(outputStruct.sidebandHist[0],"Sideband", "lpe");
    legend->AddEntry(outputStruct.signalHist[0],"Signal", "lpe");
    legend->AddEntry(outputStruct.subtractedHist[0],"Signal (minus background)", "lpe");

    for (size_t iHisto = 0; iHisto < outputStruct.histograms.size(); iHisto++) {
        cSideBand->cd(iHisto+1);
        //outputStruct.sidebandHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        outputStruct.sidebandHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        outputStruct.sidebandHist[iHisto]->GetYaxis()->SetTitle("yields");
        outputStruct.sidebandHist[iHisto]->SetMarkerStyle(kFullTriangleUp);
        outputStruct.sidebandHist[iHisto]->SetMarkerColor(kAzure);
        outputStruct.sidebandHist[iHisto]->SetLineColor(kAzure);
        double statBoxPos = gPad->GetUxmax();
        outputStruct.sidebandHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));

        cSignal->cd(iHisto+1);
        //outputStruct.signalHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        outputStruct.signalHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        outputStruct.signalHist[iHisto]->GetYaxis()->SetTitle("yields");
        outputStruct.signalHist[iHisto]->SetMarkerStyle(kFullSquare);
        outputStruct.signalHist[iHisto]->SetMarkerColor(kViolet);
        outputStruct.signalHist[iHisto]->SetLineColor(kViolet);
        outputStruct.signalHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));

        cSubtracted->cd(iHisto+1);
        //outputStruct.subtractedHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        outputStruct.subtractedHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        outputStruct.subtractedHist[iHisto]->GetYaxis()->SetTitle("yields");
        outputStruct.subtractedHist[iHisto]->SetMarkerStyle(kFullCircle);
        outputStruct.subtractedHist[iHisto]->SetMarkerColor(kPink);
        outputStruct.subtractedHist[iHisto]->SetLineColor(kPink);
        outputStruct.subtractedHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));

        cSigPlusBack->cd(iHisto+1);
        //outputStruct.signalHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        outputStruct.signalHist[iHisto]->Draw();
        outputStruct.sidebandHist[iHisto]->Draw("same");
        outputStruct.subtractedHist[iHisto]->Draw("same");
        legend->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.5, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    }
    
    TCanvas* cAllPt = new TCanvas("cAllPt","All pT,D final deltaR distribution");
    cAllPt->SetCanvasSize(1800,1000);
    //outputStruct.hSubtracted_allPtSummed->GetXaxis()->SetRangeUser(0.0, 0.5);
    outputStruct.hSubtracted_allPtSummed->SetTitle(";#DeltaR;yields");
    outputStruct.hSubtracted_allPtSummed->SetMarkerStyle(kFullCircle);
    outputStruct.hSubtracted_allPtSummed->SetMarkerColor(kRed);
    outputStruct.hSubtracted_allPtSummed->SetLineColor(kRed);
    outputStruct.hSubtracted_allPtSummed->Draw();
    double statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.5, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));

    cout << "Plotting...\n";

    //
    // Storing images
    //
    TString imagePath = "../../Images/1-SignalTreatment/SideBand/";
    c1d_fit->Update();
    c1d_fit->SaveAs(imagePath + "SB_BackSub_invariant_mass.png");
    c_2d->Update();
    c_2d->SaveAs(imagePath + "SB_BackSub_2d_deltaR_vs_invmass.png");
    cSigPlusBack->Update();
    cSigPlusBack->SaveAs(imagePath + "SB_BackSub_yield_pT_bins.png");
    cAllPt->Update();
    cAllPt->SaveAs(imagePath + "SB_BackSub_yield_pT_summed.png");

    //
    // Storing in a single pdf file
    //
    c1d_fit->Print(Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    c_2d->Print(Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cSigPlusBack->Print(Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cAllPt->Print(Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
}

void SaveData(SubtractionResult outputStruct, double jetptMin, double jetptMax){
    // Open output file
    TFile* outFile = new TFile(Form("backSub_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    // Loop over background subtracted histograms
    for (size_t iHisto = 0; iHisto < outputStruct.signalHist.size(); iHisto++) {
        // store each histogram in file
        outputStruct.subtractedHist[iHisto]->Write();
    }
    outFile->Close();
    delete outFile;
    
    
}

void BackgroundSubtraction(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;

    // jet pT cuts
    std::vector<double> ptjetBinEdges = {5., 7., 15., 30.};
    double jetptMin = ptjetBinEdges[0]; // GeV
    double jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 1]; // GeV

    // pT,D bins
    std::vector<double> ptDBinEdges = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};

    // deltaR histogram
    int deltaRbins = 10000; // deltaRbins = numberOfPoints, default=10 bins for [0. 0.4]
    //std::vector<double> deltaRBinEdges = {0.,0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.4}; // TODO: investigate structure before 0.005
    std::vector<double> deltaRBinEdges = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4}; // chosen by Nima
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];

    // mass histogram
    int massBins = 100; // default=100 
    double minMass = 1.72; // use from 1.72, used to use 1.67
    double maxMass = 2.1;
    
    // Initial parameter values
    InitialParam parametersVectors;
    parametersVectors.paramA = {6047.36, 3334.83, 3819.97, 1988.98, 1.1*1049.22, 1.1*2048.5, 944.77, 1460.87, 1448.5,};
    parametersVectors.paramB = {3.10233, 3.56564, 2.971, 4.11542, 4.19286, 3.56594, 3.36828, 3.21522, 5.28297e-316};
    parametersVectors.paramC = {1.41624e+06, 2.23897e+06, 1.56683e+06, 1.52035e+06, 916326, 1.3886e+06, 626330, 858739, 673854};
    parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
    parametersVectors.paramSigma = { sigmaInitial, 2*sigmaInitial, 1.4*sigmaInitial, 1.25*sigmaInitial, 1.5*sigmaInitial, 1.5*sigmaInitial, 2*sigmaInitial, 3*sigmaInitial, 2*sigmaInitial };
    

    // Opening data file
    TFile* fDist = new TFile("../../ExperimentalData/Hyperloop_output/AO2D.root","read");
    if (!fDist || fDist->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT file." << std::endl;
    }

    // Create multiple histograms
    std::vector<TH2D*> histograms = createHistograms(ptDBinEdges,                                 // the pT,D edges will determine the number of mass histograms
                                                     massBins, minMass, maxMass,                  // mass histograms binning
                                                     deltaRBinEdges);                             // deltaR histograms with asymmetrical bin widths
                                                     //deltaRbins, minDeltaR, maxDeltaR);         // deltaR histograms
    
    // bin sizes
    cout << "Mass bin width = " << (maxMass-minMass)/massBins << endl;


    // Fill histograms
    fillHistograms(fDist, histograms, jetptMin, jetptMax, ptDBinEdges);

    // Perform fits
    std::vector<TF1*> fittings = performFit(histograms, parametersVectors, minMass, maxMass);

    // signal/side-band region parameters
    int signalSigmas = 2; // 2
    int startingBackSigma = 4; // 4
    int backgroundSigmas = 4; // 4
    cout << "Signal: |m - m_0| < " << signalSigmas << "sigmas" << endl;
    cout << "Left side-band: " << -(startingBackSigma+backgroundSigmas) << "sigmas << |m - m_0| << " << -startingBackSigma << "sigmas\n";
    cout << "Right side-band: " << startingBackSigma << "sigmas << |m - m_0| << " << (startingBackSigma+backgroundSigmas) << "sigmas\n";

    // Subtract side-band from signal
    SubtractionResult finalDeltaR = SideBand(histograms, fittings, signalSigmas, startingBackSigma, backgroundSigmas);
    
    // Plot histograms
    PlotHistograms(histograms, fittings, finalDeltaR, jetptMin, jetptMax);

    // Storing final histograms to output file
    SaveData(finalDeltaR,jetptMin,jetptMax);

    // Cleanup
    for (auto hist : histograms) {
        delete hist;
    }
    histograms.clear();

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    BackgroundSubtraction();
    return 0;
}



//
// Fitting parameter values
//
// set alpha: 5-15 GeV/c jet pT initial fit parameters (deltaRbins = 100, massBins = 100)
/*parametersVectors.paramA = {5354.55, 4825.02, 3236.23, 1600.03, 926.955, 1684.15, 737.055, 572.857, 1.};
parametersVectors.paramB = {2.8335, 3.11021, 3.40724, 3.58429, 3.72211, 3.44699, 3.40815+0.5, 3.47049-4.5, 1.};
parametersVectors.paramC = {1.98347e+06, 2.20959e+06, 1.37818e+06, 920285, 1.18674e+06, 502167, 394194, 1.};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { 2*sigmaInitial, sigmaInitial, 2.5*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial };*/
// set beta: 15-30 GeV/c jet pT initial fit parameters (deltaRbins = 100, massBins = 100)
/*parametersVectors.paramA = {151.495, 160.988, 107.586, 95.7014, 82.6835, 221.485, 207.169, 891.023, 1448.5};
parametersVectors.paramB = {2.65389+0.5, 2.883, 2.96528, 3.19736, 3.23468, 3.3857, 3.29321, 3.21522};
parametersVectors.paramC = {21040.2, 30964.5, 29364.8, 29796.1, 35084.3, 98945.7, 124168, 464545, 673854};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { sigmaInitial, sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial };*/
// set gamma: 5-30 GeV/c jet pT initial fit parameters (deltaRbins = 100, massBins = 100)
/*parametersVectors.paramA = {5497.6, 4963.91, 3334.83, 1681.7, 1006.26, 1904.17, 944.77, 1460.87, 1448.5};
parametersVectors.paramB = {3.10233, 3.39654, 3.56564, 3.69124, 3.42524, 3.40274, 3.36828, 3.21522, };
parametersVectors.paramC = {1.41624e+06, 2.01454e+06, 2.23897e+06, 1.40804e+06, 955357, 1.28567e+06, 626330, 858739, 673854};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { 3*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 1.5*sigmaInitial, 1.5*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial };*/
// set delta: 5-30 GeV/c jet pT initial fit parameters (deltaRbins = 10, massBins = 50)
/*parametersVectors.paramA = {6047.36, 3334.83, 3819.97, 1988.98, 1.1*1049.22, 1.1*2048.5, 944.77, 1460.87, 1448.5,};
parametersVectors.paramB = {3.10233, 3.56564, 2.971, 4.11542, 4.19286, 3.56594, 3.36828, 3.21522, 5.28297e-316};
parametersVectors.paramC = {1.41624e+06, 2.23897e+06, 1.56683e+06, 1.52035e+06, 916326, 1.3886e+06, 626330, 858739, 673854};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { sigmaInitial, 2*sigmaInitial, 1.4*sigmaInitial, 1.25*sigmaInitial, 1.5*sigmaInitial, 1.5*sigmaInitial, 2*sigmaInitial, 3*sigmaInitial, 2*sigmaInitial };
1st fit worked with = [6047.36,3.10233,1.41624e+06,1.86484,0.018]
2nd fit worked with = [3334.83,3.56564,2.23897e+06,1.86484,0.024]
3rd fit worked with = [3819.97,2.971,1.56683e+06,1.86484,0.0084]
4th fit worked with = [1988.98,4.11542,1.52035e+06,1.86484,0.0096]
5th fit worked with = [1049.22,4.19286,916326,1.86484,0.0096]
6th fit worked with = [2048.5,3.56594,1.3886e+06,1.86484,0.0096]
7th fit worked with = [944.77,3.36828,626330,1.86484,0.024]
8th fit worked with = [1460.87,3.21522,858739,1.86484,0.036]
9th fit worked with = [1448.5,5.28297e-316,673854,1.86484,0.024]
*/