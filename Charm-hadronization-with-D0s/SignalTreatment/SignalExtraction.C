/**
 * @file SignalExtraction.C
 * @author Christian Reckziegel
 * @brief Macro for executing the signal extraction procedure of D0 candidates
**/


#include <bits/stdc++.h> // for measuring execution time

using namespace std;

// calculate delta phi such that 0 < delta phi < pi
double DeltaPhi(double phi1, double phi2) {
    // Compute the absolute difference between phi1 and phi2
    double dphi = std::abs(phi1 - phi2); 
    if (dphi > M_PI) {
        // subtract 2pi if the difference if bigger than pi
        dphi = dphi - 2*M_PI;
    }

    return dphi;
}

//
//
// Fit functions
//
//
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



// Module to create TH2D histograms including interest variable
std::vector<TH2D*> createHistograms(const std::vector<const char*>& names, const std::vector<const char*>& titles,
                                    int xbins, double xmin, double xmax, int ybins, double ymin, double ymax) {
    std::vector<TH2D*> histograms;
    for (size_t i = 0; i < names.size(); ++i) {
        histograms.push_back(new TH2D(names[i], titles[i], xbins, xmin, xmax, ybins, ymin, ymax));
        histograms[i]->Sumw2();
        
    }
    return histograms;
}

// Module to fill 2D histograms from TTree data
void fillHistograms(TTree* tree, const std::vector<TH2D*>& histograms, double jetptMin, double jetptMax) {
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
        //double deltaR = sqrt(pow(jetEta-hfEta,2) + pow(jetPhi-hfPhi,2));

        // Fill each histogram with their respective pT intervals
        if ((hfPt > 3 && hfPt < 4) && (histograms.size() > 0) && (jetPt > jetptMin && jetPt < jetptMax)) {
            histograms[0]->Fill(hfMass, deltaR);
        } else if ((hfPt > 4 && hfPt < 5) && (histograms.size() > 1) && (jetPt > jetptMin && jetPt < jetptMax)) {
            histograms[1]->Fill(hfMass, deltaR);
        } else if ((hfPt > 5 && hfPt < 6) && (histograms.size() > 2) && (jetPt > jetptMin && jetPt < jetptMax)) {
            histograms[2]->Fill(hfMass, deltaR);
        } else if ((hfPt > 6 && hfPt < 7) && (histograms.size() > 3) && (jetPt > jetptMin && jetPt < jetptMax)) {
            histograms[3]->Fill(hfMass, deltaR);
        } else if ((hfPt > 7 && hfPt < 8) && (histograms.size() > 4) && (jetPt > jetptMin && jetPt < jetptMax)) {
            histograms[4]->Fill(hfMass, deltaR);
        } else if ((hfPt > 8 && hfPt < 10) && (histograms.size() > 5) && (jetPt > jetptMin && jetPt < jetptMax)) {
            histograms[5]->Fill(hfMass, deltaR);
        } else if ((hfPt > 10 && hfPt < 12) && (histograms.size() > 6) && (jetPt > jetptMin && jetPt < jetptMax)) {
            histograms[6]->Fill(hfMass, deltaR);
        } else if ((hfPt > 12 && hfPt < 15) && (histograms.size() > 7) && (jetPt > jetptMin && jetPt < jetptMax)) {
            histograms[7]->Fill(hfMass, deltaR);
        } else if ((hfPt > 15 && hfPt < 30) && (histograms.size() > 8) && (jetPt > jetptMin && jetPt < jetptMax)) {
            histograms[8]->Fill(hfMass, deltaR);
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

    double rangeFactor = 0.1; // 0.1=10% around the initially calculated values
    int stepsNumber = 10;
    // loop over b parameter range
    for (double iB = b_parameter*(1-rangeFactor); iB < b_parameter*(1+rangeFactor); iB+=(b_parameter*2*rangeFactor/stepsNumber)) {
        //cout << "iB = " << iB << endl;
        // loop over a parameter range
        for (double iA = a_parameter*(1-rangeFactor); iA < a_parameter*(1+rangeFactor); iA+=(a_parameter*2*rangeFactor/stepsNumber)) {
            //cout << "iA = " << iA << endl;
            // loop over C parameter range
            for (double iC = C_parameter*(1-rangeFactor); iC < C_parameter*(1+rangeFactor); iC+=(C_parameter*2*rangeFactor/stepsNumber)) {
                //cout << "iC = " << iC << endl;
                for (double iS = sigma_parameter*(1-rangeFactor); iS < sigma_parameter*(1+rangeFactor); iS+=(sigma_parameter*2*rangeFactor/stepsNumber)) {
                //for (double iS = sigma_parameter*(1-0.5); iS < sigma_parameter*(1+3); iS+=(sigma_parameter*3.5/stepsNumber)) {
                    // Perform the fit
                    fTestFit = new TF1("fTestFit",customFitFunction, minMass, maxMass, numParameters);
                    double m_0_lower_limit = 0.95 * m_0_parameter;
                    double m_0_upper_limit = 1.95 * m_0_parameter;
                    //fTestFit->SetParLimits(2, 0., DBL_MAX); // Set lower boundary of parameter iC to 0, and higher to maximum representable value for a double -> no negative gaussians
                    //fTestFit->SetParLimits(2, 0., std::numeric_limits<double>::infinity());
                    fTestFit->SetParLimits(3, m_0_lower_limit, m_0_upper_limit); // constraints on m_0_parameter: [0.95*m0, 1.05*m0]
                    fTestFit->SetParLimits(4, 0, 1); // constraints on sigma: [0, 1]
                    fTestFit->SetParameters(iA, iB, iC, m_0_parameter, iS);
                    //fTestFit->SetParameters(iA, iB, C_parameter, m_0_parameter, iS);
                    fTestFit->SetParNames("a","b", "C", "m_0", "sigma");
                    fTestFit->SetLineColor(kBlue);
                    histogram->Fit(fTestFit, "Q");// "Q" option performs quiet fit without drawing the fit function

                    // Verify fit quality with current parameter values
                    if ((fTestFit->GetChisquare() < bestChiSquare) // select best fit
                         && (fTestFit->GetParameter(2) > 0) // avoiding negative gaussians
                         && (fTestFit->GetParameter(4) > sigma_parameter/2) // avoiding point delta distributions
                         && (abs(fTestFit->GetParameter(3)-m_0_parameter) < 4*sigma_parameter)) { // avoid mean values too far from literature value for HF invariant mass
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
 * @param[in] names Histogram names vector
 * @param[in] histograms Histograms vector
 * @return Vector of fit objects TF1.
 */
// Module to perform fits to histograms
std::vector<TF1*> performFit(const std::vector<const char*>& names, const std::vector<TH2D*>& histograms2d, InitialParam parametersVectors, 
                             double minMass, double maxMass, double minDeltaR, double maxDeltaR, int numberOfPoints, std::vector<TFitResultPtr>& fitresultptr) {
    // creating temporary histogram for storing to the vector
    TH1D* tempHist;

    // creating 1D mass projection histograms
    std::vector<TH1D*> histograms;

    double deltaRinterval = (maxDeltaR-minDeltaR)/numberOfPoints;

    // obtaining 1D invariant mass histograms from the projection
    for (size_t iHist = 0; iHist < histograms2d.size(); iHist++) {
        // create numberOfPoints histograms for each selection cut (jet pT and HF pT)
        for (size_t iPoint = 0; iPoint < numberOfPoints; iPoint++) {
            double minBin = histograms2d[iHist]->GetYaxis()->FindBin(deltaRinterval*iPoint);
            double maxBin = histograms2d[iHist]->GetYaxis()->FindBin(deltaRinterval*(iPoint+1));
            tempHist = histograms2d[iHist]->ProjectionX(Form("h_mass_proj_%zu_%zu", iHist, iPoint), minBin, maxBin);
            histograms.push_back(tempHist);
        }
        // histograms are filled in vector as
        // histograms = {histMass_0_0, histMass_0_1, histMass_0_2, ..., histMass_0_(numberOfPoints-1), histMass_1_0, histMass_1_1, histMass_1_2, ...}
    }
    
    
    // Creating vector for storing total fits (and/or empty fits) and signal fits
    std::vector<TF1*> fittings;
    
    // calculating initial parameters
    double b_parameter;
    double a_parameter;
    double sigma_parameter = 0.012;
    double m_0_parameter = 1.86484; // D0 mass in GeV/c^2
    double C_parameter;
    
    // Perform total fit to each histogram
    for (size_t iHisto = 0; iHisto < (histograms.size()/numberOfPoints); ++iHisto) {
        //cout << "iHisto = " << iHisto << endl;
        for (size_t iPoint = 0; iPoint < numberOfPoints; iPoint++) {
           // cout << "iPoint = " << iPoint << endl;
            // variables necessary for calculating parameters
            double K = histograms[iHisto*numberOfPoints+iPoint]->GetBinContent(1); // Content of the first bin
            double J = histograms[iHisto*numberOfPoints+iPoint]->GetBinContent(histograms[iHisto*numberOfPoints+iPoint]->GetNbinsX()); // Content of the last bin
            double m_min = histograms[iHisto*numberOfPoints+iPoint]->GetBinLowEdge(1); // Lower limit of the first bin
            double m_max = histograms[iHisto*numberOfPoints+iPoint]->GetBinLowEdge(histograms[iHisto*numberOfPoints+iPoint]->GetNbinsX() + 1); // Upper limit of the last bin
            double I_tot = histograms[iHisto*numberOfPoints+iPoint]->Integral();
            double I_3sigma = histograms[iHisto*numberOfPoints+iPoint]->Integral(histograms[iHisto*numberOfPoints+iPoint]->FindBin(m_0_parameter - 3 * sigma_parameter), histograms[iHisto*numberOfPoints+iPoint]->FindBin(m_0_parameter + 3 * sigma_parameter));
            
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
            
            // Store a fit object regardless of the histogram being empty or not
            if (histograms[iHisto*numberOfPoints+iPoint]->GetEntries() == 0) {
                // Store empty TF1 object for associated empty histogram
                fittings.push_back(new TF1(Form("totalFit_%s_%zu",names[iHisto],iPoint),"")); // 0 parameters
                //cout << histograms[iHisto*numberOfPoints+iPoint]->GetTitle() << endl;
            } else{
                // Create the fit function, enforce constraints on some parameters, set initial parameter values and performing fit
                fittings.push_back(new TF1(Form("totalFit_%s_%zu",names[iHisto],iPoint), customFitFunction, minMass, maxMass, 5)); // 5 parameters
            }
            
            // fittings are filled in vector as
            // fittings = {totalFit_0_0, totalFit_0_1, ..., totalFit_0_(numberOfPoints-1), totalFit_1_0, totalFit_1_1, totalFit_1_2, ...}

            bool usePreSetParam = false; // use given pre-set initial parameter value
            if (usePreSetParam) {
                a_parameter = parametersVectors.paramA[iHisto];
                b_parameter = parametersVectors.paramB[iHisto];
                C_parameter = parametersVectors.paramC[iHisto];
                m_0_parameter = parametersVectors.paramM0[iHisto];
                sigma_parameter = parametersVectors.paramSigma[iHisto];
            }
            

            double m_0_lower_limit = 0.95 * m_0_parameter;
            double m_0_upper_limit = 1.95 * m_0_parameter;
            
            // Perform fit only for not empty TF1 objects
            if (fittings[iHisto*numberOfPoints+iPoint]->GetNpar() > 0 && histograms[iHisto*numberOfPoints+iPoint]->GetEntries() != 0) { // histograms[iHisto*numberOfPoints+iPoint]->GetEntries() != 0, fittings[iHisto*numberOfPoints+iPoint]->GetNpar() > 0
                // Calculating optimal parameters
                std::vector<double> optimalParameters = bestFit(histograms[iHisto*numberOfPoints+iPoint], minMass, maxMass, 5);
                
                //fittings[iHisto*numberOfPoints+iPoint]->SetParLimits(2, 0., DBL_MAX); // Set lower boundary of parameter iC to 0, and higher to maximum representable value for a double -> no negative gaussians
                fittings[iHisto*numberOfPoints+iPoint]->SetParLimits(3, m_0_lower_limit, m_0_upper_limit); // constraints on m_0_parameter: [0.95*m0, 1.05*m0]
                fittings[iHisto*numberOfPoints+iPoint]->SetParLimits(4, 0, 1); // constraints on sigma: [0, 1]
                //fittings[iHisto*numberOfPoints+iPoint]->SetParameters(a_parameter, b_parameter, C_parameter, m_0_parameter, sigma_parameter);
                fittings[iHisto*numberOfPoints+iPoint]->SetParameters(optimalParameters[0], optimalParameters[1], optimalParameters[2], optimalParameters[3], optimalParameters[4]);
                fittings[iHisto*numberOfPoints+iPoint]->SetParNames("a","b", "C", "m_0", "sigma");
                fittings[iHisto*numberOfPoints+iPoint]->SetLineColor(kBlue);
                histograms[iHisto*numberOfPoints+iPoint]->Fit(fittings[iHisto*numberOfPoints+iPoint], "Q");// "Q" option performs quiet fit without drawing the fit function
            } else{
                cout << "No fit performed (empty histogram)\n";
            }
            
        }
        
        

    }


    // Create and store signal only fit to each histogram
    for (size_t iHisto = 0; iHisto < (histograms.size()/numberOfPoints); iHisto++) {
        for (size_t iPoint = 0; iPoint < numberOfPoints; iPoint++) {
            // For non empty fits
            if (fittings[iHisto*numberOfPoints+iPoint]->GetNpar() > 0 && histograms[iHisto*numberOfPoints+iPoint]->GetEntries() != 0) { // histograms[iHisto*numberOfPoints+iPoint]->GetEntries() != 0, fittings[iHisto*numberOfPoints+iPoint]->GetNpar() > 0
                // Getting total parameter values
                double c_par = fittings[iHisto*numberOfPoints+iPoint]->GetParameter(2); // Get the value of parameter 'c'
                double m_par = fittings[iHisto*numberOfPoints+iPoint]->GetParameter(3); // Get the value of parameter 'm_0'
                double sigma_par = fittings[iHisto*numberOfPoints+iPoint]->GetParameter(4); // Get the value of parameter 'sigma'
                fittings.push_back(new TF1(Form("sigFit_%s_%zu", names[iHisto],iPoint), signalFunction, minMass, maxMass, 3));
                fittings[fittings.size()-1]->SetParameters(c_par,m_par,sigma_par); // fittings.size()-1 = index from the latest fitting added
                fittings[fittings.size()-1]->SetLineStyle(kDashed);
                fittings[fittings.size()-1]->SetLineColor(kGreen);
            } else{
                fittings.push_back(new TF1(Form("sigFit_%s_%zu", names[iHisto],iPoint),""));
            }

            
            // signal fits filled in vector the same way as total fits
        }
        
        

    }
    
    cout << "Fits performed.\n";
    return fittings;
    // first n fittings are the total fit functions (from 0 to n-1)
    // second n fittings are the signal fit functions (n to 2n-1)
}

struct ExtractionResult {
    std::vector<TH1D*> histograms; // 1D mass projection histograms vector
    std::vector<TH1D*> signalHist; // 1D signal histograms vector
};

ExtractionResult PointsExtraction(const std::vector<TH2D*>& histograms2d, const std::vector<TF1*>& fittings, int signalSigmas, int startingBackSigma, int backgroundSigmas,
                                  double minDeltaR, double maxDeltaR, int numberOfPoints, double minMass, double maxMass){
    
    // Creating temporary histogram for collecting data
    TH1D* tempHist;

    // Creating output struct object with histogram vectors
    ExtractionResult vectorOutputs;

    double deltaRinterval = (maxDeltaR-minDeltaR)/numberOfPoints;

    // obtaining 1D invariant mass histograms from the projection
    for (size_t iHisto = 0; iHisto <  histograms2d.size(); iHisto++) {
        cout << "iHisto = " << iHisto << endl;
        // create numberOfPoints histograms for each selection cut (jet pT and HF pT)
        for (size_t iPoint = 0; iPoint < numberOfPoints; iPoint++) {
            double minBin = histograms2d[iHisto]->GetYaxis()->FindBin(deltaRinterval*iPoint);
            double maxBin = histograms2d[iHisto]->GetYaxis()->FindBin(deltaRinterval*(iPoint+1));
            tempHist = histograms2d[iHisto]->ProjectionX(Form("h_mass_proj_%zu_%zu", iHisto, iPoint), minBin,maxBin);
            cout << "Number of entries on point histogram " << iPoint << " = " << tempHist->GetEntries() << endl;
            vectorOutputs.histograms.push_back(tempHist);
        }
        // histograms are filled in vector as
        // histograms = {histMass_0_0, histMass_0_1, histMass_0_2, ..., histMass_0_(numberOfPoints-1), histMass_1_0, histMass_1_1, histMass_1_2, ...}
    }

    TF1* tempFit;
    // Calculating delta R histograms for signal area only
    for (size_t iHisto = 0; iHisto < histograms2d.size(); iHisto++) { // loop through selection cuts: HF pT intervals
        // define asymmetrical bin widths manually chosen
        Double_t binEdges[] = {0., 0.05, 0.075, 0.1, 0.125, 0.15, 0.2};
        // Creating histogram delta R
        tempHist = new TH1D(Form("h_sig_extracted_%zu",iHisto),";#DeltaR;dN/d(#DeltaR)",numberOfPoints,minDeltaR, maxDeltaR);

        // Construct each histogram bin-by-bin
        for (size_t iPoint = 0; iPoint < numberOfPoints; iPoint++) { // loop through delta R intervals (for each HF pT interval)
            
            //int signalFitIndex = vectorOutputs.histograms.size() + iHisto + iPoint;
            int signalFitIndex = vectorOutputs.histograms.size() + iHisto*numberOfPoints + iPoint; 
            
            // Calculating signal area
            double binWidth;
            double m_0; // Get the value of parameter 'm_0'
            double sigma; // Get the value of parameter 'sigma'
            double fitSigIntegral;
            double fitSigIntegralError;

            if (fittings[signalFitIndex]->GetNpar() > 0) { // if fit is not empty
                // Calculating signal area
                binWidth = vectorOutputs.histograms[iHisto*numberOfPoints + iPoint]->GetXaxis()->GetBinWidth(1);
                m_0 = fittings[signalFitIndex]->GetParameter(1); // Get the value of parameter 'm_0'
                sigma = fittings[signalFitIndex]->GetParameter(2); // Get the value of parameter 'sigma'
                fitSigIntegral = fittings[signalFitIndex]->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma) / binWidth;
                //fitSigIntegralError = fittings[signalFitIndex]->IntegralError(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma,fittings[signalFitIndex]->GetParams(), fittings[signalFitIndex]->GetCovarianceMatrix()->GetMatrixArray() );
                cout << "Fitting number " << signalFitIndex << " has integral of " << fitSigIntegral << endl;
                //cout << "Chi2 = " << fittings[iHisto*numberOfPoints + iPoint]->GetChisquare() << ", NDF = " << fittings[iHisto*numberOfPoints + iPoint]->GetNDF() << endl;
            } else{ // if fit is empty
                fitSigIntegral = 0.;
                fitSigIntegralError = 0.;
                cout << "Fitting number " << signalFitIndex << " did not happen.\n\n";

            }

            bool test = false;
            if (fittings[iHisto*numberOfPoints+iPoint]->GetNpar() > 0 && test) {
                double c_par = fittings[iHisto*numberOfPoints+iPoint]->GetParameter(2); // Get the value of parameter 'c'
                double m_par = fittings[iHisto*numberOfPoints+iPoint]->GetParameter(3); // Get the value of parameter 'm_0'
                double sigma_par = fittings[iHisto*numberOfPoints+iPoint]->GetParameter(4); // Get the value of parameter 'sigma'
                tempFit = new TF1(Form("newsigFit_%s_%zu", histograms2d[iHisto]->GetName(),iPoint), signalFunction, minMass, maxMass, 3);
                tempFit->SetParameters(c_par,m_par,sigma_par); // fittings.size()-1 = index from the latest fitting added

                binWidth = vectorOutputs.histograms[iHisto*numberOfPoints + iPoint]->GetXaxis()->GetBinWidth(1);
                m_0 = tempFit->GetParameter(1); // Get the value of parameter 'm_0'
                sigma = tempFit->GetParameter(2); // Get the value of parameter 'sigma'
                fitSigIntegral = tempFit->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma) / binWidth;
                cout << "Fitting number " << iHisto*numberOfPoints+iPoint << " has integral of " << fitSigIntegral << endl;
                delete tempFit;
            }
            
            

            // Creating delta R bin entry
            tempHist->SetBinContent(iPoint+1,fitSigIntegral);
            //tempHist->SetBinError(iPoint,fitSigIntegralError);
            

        }
        
        // Storing delta R distribution
        vectorOutputs.signalHist.push_back(tempHist);

    }

    
    

    cout << "Points extracted.\n";
    // Return the output struct object containing filled histogram vectors
    return vectorOutputs;
    
}

void PlotHistograms(const std::vector<TH2D*>& histograms2d, const std::vector<TF1*>& fittings, ExtractionResult outputStruct, double jetptMin, double jetptMax,
                    int numberOfPoints) {
    
    // creating 1D mass projection histograms
    TH1D* tempHist;
    std::vector<TH1D*> histograms;

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    // obtaining 1D invariant mass histograms from the projection
    for (size_t iHist = 0; iHist < histograms2d.size(); iHist++) {
        //
        tempHist = histograms2d[iHist]->ProjectionX(Form("h_mass_proj_%zu", iHist));
        histograms.push_back(tempHist);
    }

    // Plotting point mass histograms
    std::vector<TCanvas*> cPointMass;
    //TCanvas* cPointMass = new TCanvas("cPointMass", "Point mass histograms", 800, 600);
    //cPointMass->Divide(2+1,numberOfPoints / 2); // columns, lines

    for (size_t iHisto = 0; iHisto < histograms2d.size(); iHisto++) {
        //
        cPointMass.push_back(new TCanvas(Form("histoMass_%zu",iHisto),Form("Points histograms for interval number %zu",iHisto)));
        if(numberOfPoints%2 == 1){
            cPointMass[iHisto]->Divide(2+1,numberOfPoints / 2); // columns, lines
        } else{
            cPointMass[iHisto]->Divide(2,numberOfPoints / 2); // columns, lines
        }
        for (size_t iPoint = 0; iPoint < numberOfPoints; iPoint++) {
            cPointMass[iHisto]->cd(iPoint+1);
            double statBoxPos = gPad->GetUxmax();
            outputStruct.histograms[iHisto*numberOfPoints + iPoint]->SetMarkerStyle(kFullDotMedium);
            outputStruct.histograms[iHisto*numberOfPoints + iPoint]->SetMarkerColor(kBlack);
            outputStruct.histograms[iHisto*numberOfPoints + iPoint]->SetLineColor(kGray);
            outputStruct.histograms[iHisto*numberOfPoints + iPoint]->GetYaxis()->SetTitle("counts");
            outputStruct.histograms[iHisto*numberOfPoints + iPoint]->Draw();
            fittings[iHisto*numberOfPoints + iPoint]->Draw("same");
            //cout << "Number of entries on point histogram " << iPoint << " = " << outputStruct.histograms[iPoint]->GetEntries() << endl;
            latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
        }
    }
    
    
    
    //
    // Plotting output observable
    //
    TCanvas* cSignal = new TCanvas("cSignal", "delta R for extracted signal", 800, 600);
    cSignal->Divide(3,static_cast<int>(outputStruct.signalHist.size() / 3)); // columns, lines

    for (size_t iHisto = 0; iHisto < outputStruct.signalHist.size(); iHisto++) {
        cSignal->cd(iHisto+1);
        outputStruct.signalHist[iHisto]->GetXaxis()->SetRangeUser(0.0, 0.5);
        outputStruct.signalHist[iHisto]->GetXaxis()->SetTitle("#DeltaR");
        outputStruct.signalHist[iHisto]->GetYaxis()->SetTitle("yields");
        outputStruct.signalHist[iHisto]->SetMarkerStyle(kOpenCircle);
        outputStruct.signalHist[iHisto]->SetMarkerColor(kRed);
        outputStruct.signalHist[iHisto]->SetLineColor(kRed);
        outputStruct.signalHist[iHisto]->SetTitle(histograms2d[iHisto]->GetTitle());
        double statBoxPos = gPad->GetUxmax();
        outputStruct.signalHist[iHisto]->Sumw2();
        outputStruct.signalHist[iHisto]->Draw();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    }
    
    

    cout << "Plotting...\n";

    // Save the canvas as an image or display it
    for (size_t iHisto = 0; iHisto < histograms2d.size(); iHisto++) {
        //cPointMass[iHisto]->SetCanvasSize(1800,1800);
        if (iHisto == 0) {
            cPointMass[iHisto]->Print(Form("sig_extraction_deltaR_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
        } else{
            cPointMass[iHisto]->Print(Form("sig_extraction_deltaR_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
        }
        
        
    }
    //cSignal->SetCanvasSize(1800,1800);
    cSignal->Print(Form("sig_extraction_deltaR_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));
}

void SaveData(ExtractionResult outputStruct, double jetptMin, double jetptMax){
    // Open output file
    TFile* outFile = new TFile(Form("sigExt_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    // Loop over signal extracted histograms
    for (size_t iHisto = 0; iHisto < outputStruct.signalHist.size(); iHisto++) {
        // store each histogram in file
        outputStruct.signalHist[iHisto]->Write();
    }
    outFile->Close();
    delete outFile;
    
    
}

void SignalExtraction(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;
    // jet pT cuts
    double jetptMin = 5; // GeV
    double jetptMax = 30; // GeV
    // deltaR histogram
    int deltaRbins = 50; // 20 bins
    double minDeltaR = 0.;
    double maxDeltaR = 0.4;
    int numberOfPoints = 10;
    // mass histogram
    int massBins = 50; 
    double minMass = 1.67;
    double maxMass = 2.1;
    // Initial parameter values
    InitialParam parametersVectors;
    parametersVectors.paramA = {5497.6, 4963.91, 3334.83, 1681.7, 1006.26, 1904.17, 944.77, 1460.87, 1448.5};
    parametersVectors.paramB = {3.10233, 3.39654, 3.56564, 3.69124, 3.42524, 3.40274, 3.36828, 3.21522, };
    parametersVectors.paramC = {1.41624e+06, 2.01454e+06, 2.23897e+06, 1.40804e+06, 955357, 1.28567e+06, 626330, 858739, 673854};
    parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
    parametersVectors.paramSigma = { 3*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 1.5*sigmaInitial, 1.5*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial };
    

    // opening file
    TFile* fDist = new TFile("../ExperimentalData/Hyperloop_output/AO2D.root","read");
    if (!fDist || fDist->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT file." << std::endl;
    }

    // Create multiple histograms
    std::vector<const char*> names = {"histMass1", "histMass2", "histMass3", 
                                      "histMass4", "histMass5", "histMass6",
                                      "histMass7", "histMass8", "histMass9"}; // Names of histograms
    std::vector<const char*> titles = {"3 < p_{T,D} < 4 GeV/c;m(K#pi) GeV/c^{2};#DeltaR", "4 < p_{T,D} < 5 GeV/c;m(K#pi) GeV/c^{2};#DeltaR", "5 < p_{T,D} < 6 GeV/c;m(K#pi) GeV/c^{2};#DeltaR",
                                       "6 < p_{T,D} < 7 GeV/c;m(K#pi) GeV/c^{2};#DeltaR", "7 < p_{T,D} < 8 GeV/c;m(K#pi) GeV/c^{2};#DeltaR", "8 < p_{T,D} < 10 GeV/c;m(K#pi) GeV/c^{2};#DeltaR",
                                       "10 < p_{T,D} < 12 GeV/c;m(K#pi) GeV/c^{2};#DeltaR", "12 < p_{T,D} < 15 GeV/c;m(K#pi) GeV/c^{2};#DeltaR", "15 < p_{T,D} < 30 GeV/c;m(K#pi) GeV/c^{2};#DeltaR"}; // Titles of histograms
    std::vector<TH2D*> histograms = createHistograms(names, titles, 
                                                       massBins, minMass, maxMass, // mass histograms
                                                       deltaRbins, minDeltaR, maxDeltaR); // deltaR histograms
    
    // bin sizes
    cout << "Mass bin width = " << (maxMass-minMass)/massBins << endl;
    cout << "DeltaR bin width = " << (maxDeltaR-minDeltaR)/deltaRbins << endl;

    // Accessing TTree
    TTree* tDist = (TTree*)fDist->Get("DF_2261906078621696/O2jetdisttable");

    // Fill histograms
    fillHistograms(tDist, histograms, jetptMin, jetptMax);

    std::vector<TFitResultPtr> fitResult;
    // Perform fits
    std::vector<TF1*> fittings = performFit(names,histograms, parametersVectors, minMass, maxMass, minDeltaR, maxDeltaR, numberOfPoints, fitResult);

    // signal/side-band region parameters
    int signalSigmas = 2; // 2
    int startingBackSigma = 4; // 4
    int backgroundSigmas = 4; // 4
    cout << "Signal: |m - m_0| < " << signalSigmas << "sigmas" << endl;
    cout << "Left side-band: " << -(startingBackSigma+backgroundSigmas) << "sigmas < |m - m_0| < " << -startingBackSigma << "sigmas\n";
    cout << "Right side-band: " << startingBackSigma << "sigmas << |m - m_0| < " << (startingBackSigma+backgroundSigmas) << "sigmas\n";

    // Subtract side-band from signal
    ExtractionResult finalDeltaR = PointsExtraction(histograms, fittings, signalSigmas, startingBackSigma, backgroundSigmas, minDeltaR, maxDeltaR, numberOfPoints, minMass, maxMass);
    
    // Plot histograms
    PlotHistograms(histograms, fittings, finalDeltaR, jetptMin, jetptMax, numberOfPoints);

    // Storing final histograms to output file
    SaveData(finalDeltaR,jetptMin,jetptMax);

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 
    
}

int main(){
    SignalExtraction();
    return 0;
}



//
// Fitting parameter values
//
// 5-15 GeV/c jet pT initial fit parameters (deltaRbins = 100, massBins = 100)
/*parametersVectors.paramA = {5354.55, 4825.02, 3236.23, 1600.03, 926.955, 1684.15, 737.055, 572.857, 1.};
parametersVectors.paramB = {2.8335, 3.11021, 3.40724, 3.58429, 3.72211, 3.44699, 3.40815+0.5, 3.47049-4.5, 1.};
parametersVectors.paramC = {1.98347e+06, 2.20959e+06, 1.37818e+06, 920285, 1.18674e+06, 502167, 394194, 1.};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { 2*sigmaInitial, sigmaInitial, 2.5*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial };*/
// 15-30 GeV/c jet pT initial fit parameters (deltaRbins = 100, massBins = 100)
/*parametersVectors.paramA = {151.495, 160.988, 107.586, 95.7014, 82.6835, 221.485, 207.169, 891.023, 1448.5};
parametersVectors.paramB = {2.65389+0.5, 2.883, 2.96528, 3.19736, 3.23468, 3.3857, 3.29321, 3.21522};
parametersVectors.paramC = {21040.2, 30964.5, 29364.8, 29796.1, 35084.3, 98945.7, 124168, 464545, 673854};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { sigmaInitial, sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial };*/
// 5-30 GeV/c jet pT initial fit parameters (deltaRbins = 100, massBins = 100)
/*parametersVectors.paramA = {5497.6, 4963.91, 3334.83, 1681.7, 1006.26, 1904.17, 944.77, 1460.87, 1448.5};
parametersVectors.paramB = {3.10233, 3.39654, 3.56564, 3.69124, 3.42524, 3.40274, 3.36828, 3.21522, };
parametersVectors.paramC = {1.41624e+06, 2.01454e+06, 2.23897e+06, 1.40804e+06, 955357, 1.28567e+06, 626330, 858739, 673854};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { 3*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 1.5*sigmaInitial, 1.5*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial };*/