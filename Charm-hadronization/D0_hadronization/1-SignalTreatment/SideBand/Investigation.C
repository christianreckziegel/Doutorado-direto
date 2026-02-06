/*
 *
 *
 * Macro for executing the background subtraction procedure of D0 candidates
 * using the side-band technique
 * 
 * 
 * 
 * 
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


void Investigation(){
    
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = 0.4; // on particle level delta R

    // signal/side-band region parameters
    int signalSigmas = 2; // 2
    int startingBackSigma = 4; // 4
    int backgroundSigmas = 4; // 4

    // pT,jet
    double jetptMin = 5;
    double jetptMax = 30;

    // deltaR histogram
    std::vector<double> deltaRBinEdges = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4};
    // mass histogram
    int massBins = 100; // default = 50
    double minMass = 1.72;// in the past = 1.67, current default = 1.72
    double maxMass = 2.1;

    // creating histograms
    TH2D* h2_HFmass = new TH2D("histMass", "6 < p_{T,D} < 7 GeV/c;m(K#pi) GeV/c^{2};#DeltaR", massBins, minMass, maxMass, deltaRBinEdges.size() - 1, deltaRBinEdges.data()); 

    // opening file
    TFile* fDist = new TFile("../../ExperimentalData/Hyperloop_output/AO2D.root","read");
    if (!fDist || fDist->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT file." << std::endl;
    }

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

        if ((abs(hfEta) < etaCut) && (abs(hfY) < yCut) && ((jetPt >= jetptMin) && (jetPt < jetptMax)) && ((deltaR >= 0.) && (deltaR < DeltaRcut))) {
            if (hfPt > 6. && hfPt < 7.) {
                h2_HFmass->Fill(hfMass, deltaR);
            }
        }
        
    }
    cout << "Histograms filled.\n";

    //
    // Perform fits
    //
    TH1D* hInvMass = h2_HFmass->ProjectionX("h_mass_proj");
    // calculating initial parameters
    double b_parameter;
    double a_parameter;
    double sigma_parameter = 0.012;
    double m_0_parameter = 1.86484; // D0 mass in GeV/c^2
    double C_parameter;
    // variables necessary for calculating parameters
    double K = hInvMass->GetBinContent(1); // Content of the first bin
    double J = hInvMass->GetBinContent(hInvMass->GetNbinsX()); // Content of the last bin
    double m_min = hInvMass->GetBinLowEdge(1); // Lower limit of the first bin
    double m_max = hInvMass->GetBinLowEdge(hInvMass->GetNbinsX() + 1); // Upper limit of the last bin
    double I_tot = hInvMass->Integral();
    double I_3sigma = hInvMass->Integral(hInvMass->FindBin(m_0_parameter - 3 * sigma_parameter), hInvMass->FindBin(m_0_parameter + 3 * sigma_parameter));
    
    // Calculating parameter initial values
    b_parameter = (K > 0 && J > 0) ? TMath::Log(K / J) / (m_max - m_min) : -1;
    a_parameter = I_tot / (TMath::Power(m_max, b_parameter + 1) - TMath::Power(m_min, b_parameter + 1)) / (b_parameter + 1); // a(b)
    C_parameter = (I_3sigma - a_parameter * (TMath::Power(m_0_parameter + 3 * sigma_parameter, b_parameter + 1) - TMath::Power(m_0_parameter - 3 * sigma_parameter, b_parameter + 1)) / (b_parameter + 1)) / (TMath::Sqrt(2 * TMath::Pi()) * sigma_parameter); // C(a,b)

    
    TF1* fitSigBack = new TF1("fitSigBack", customFitFunction, minMass, maxMass, 5);
    // Calculating optimal parameters
    std::vector<double> optimalParameters = bestFit(hInvMass, minMass, maxMass, 5);
    double m_0_lower_limit = 0.95 * m_0_parameter;
    double m_0_upper_limit = 1.95 * m_0_parameter;
    //fitSigBack->SetParLimits(2, 0., DBL_MAX); // Set lower boundary of parameter iC to 0, and higher to maximum representable value for a double -> no negative gaussians
    fitSigBack->SetParLimits(3, m_0_lower_limit, m_0_upper_limit); // constraints on m_0_parameter: [0.95*m0, 1.05*m0]
    fitSigBack->SetParLimits(4, 0, 1); // constraints on sigma: [0, 1]
    //fitSigBack->SetParameters(a_parameter, b_parameter, C_parameter, m_0_parameter, sigma_parameter);
    fitSigBack->SetParameters(optimalParameters[0], optimalParameters[1], optimalParameters[2], optimalParameters[3], optimalParameters[4]);
    fitSigBack->SetParNames("a","b", "C", "m_0", "sigma");
    fitSigBack->SetLineColor(kBlue);
    hInvMass->Fit(fitSigBack, "Q");// "Q" option performs quiet fit without drawing the fit function

    TF1* fitBack = new TF1("fitBack", backgroundFunction, minMass, maxMass, 2);
    double a_par = fitSigBack->GetParameter(0); // Get the value of parameter 'a'
    double b_par = fitSigBack->GetParameter(1); // Get the value of parameter 'b'
    fitBack->SetParameters(a_par,b_par); // fittings.size()-1 = the latest added to the vector
    fitBack->SetLineStyle(kDashed);
    fitBack->SetLineColor(kGreen);

    std::cout << "Fits performed." << std::endl;

    //
    // Subtract the background
    //
    // Getting fit parameters
    double m_0 = fitSigBack->GetParameter(3); // Get the value of parameter 'm_0'
    double sigma = fitSigBack->GetParameter(4); // Get the value of parameter 'sigma'

    TH1D* h_sideBand;
    TH1D* h_signal;
    TH1D* h_back_subtracted;

    // Check which sides should be used for the total side-band distribution (if at least 1 sigma fit inside left range)
    if (m_0 - (startingBackSigma+1)*sigma <  h2_HFmass->ProjectionX("h_mass_")->GetBinLowEdge(1)) {
        std::cout << "Using only right sideband." << std::endl;

        // Count for how many sigmas is there room inside the left side range
        double leftSidebandRange = (m_0 - startingBackSigma * sigma) - h2_HFmass->ProjectionX("h_mass_")->GetBinLowEdge(1);
        int leftSigmas = static_cast<int>(leftSidebandRange / sigma);// get the integral number
        std::cout << leftSigmas << " sigmas fit inside the left side range for the background distribution estimation." << std::endl;
        // Count for how many sigmas is there room inside the right side range
        double rightSidebandRange = h2_HFmass->ProjectionX("h_mass_")->GetBinLowEdge(h2_HFmass->ProjectionX("h_mass_")->GetNbinsX()+1) - (m_0 + startingBackSigma*sigma);
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
        double rightSB = fitBack->Integral(m_0 + startingBackSigma * sigma,m_0 + (startingBackSigma + rightSigmas) * sigma);
        // Finding side-band region: signal < |2*sigma|
        double backSigIntegral = fitBack->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of background in signal region
        double fitSigIntegral = fitSigBack->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

        // Calculating scaling parameter
        double alpha = backSigIntegral/rightSB;
        std::cout << "alpha = " << alpha << std::endl;

        // Create side-band histogram: use only the right sideband
        int lowBin = h2_HFmass->GetXaxis()->FindBin(m_0 + startingBackSigma * sigma);
        int highBin = h2_HFmass->GetXaxis()->FindBin(m_0 + (startingBackSigma + rightSigmas) * sigma);
        h_sideBand = h2_HFmass->ProjectionY("h_sideband_proj_temp", lowBin, highBin); // sum the right sideband

        // Apply scaling to background function
        //h_sideBand->Scale(alpha);

        // Create signal histogram
        lowBin = h2_HFmass->GetXaxis()->FindBin(m_0 - signalSigmas * sigma);
        highBin = h2_HFmass->GetXaxis()->FindBin(m_0 + signalSigmas * sigma);
        h_signal = h2_HFmass->ProjectionY("h_signal_proj", lowBin, highBin);

        cout << "Distribution ratio of sideband over raw signal: " << h_sideBand->Integral() / h_signal->Integral() << endl;
        cout << "Fit function ratio of sideband over raw signal: " << (rightSB)/ fitSigIntegral << endl;

        h_sideBand->Scale(alpha);

        // Subtract background histogram from signal histogram
        h_back_subtracted = (TH1D*)h_signal->Clone("h_back_subtracted");
        h_back_subtracted->Add(h_sideBand,-1.0);

        h_back_subtracted->Scale(1/0.9545);

    } else {
        std::cout << "Using both left and right sidebands.\n";

        // Count for how many sigmas is there room inside the left side range
        double leftSidebandRange = (m_0 - startingBackSigma * sigma) - h2_HFmass->ProjectionX("h_mass_")->GetBinLowEdge(1);
        int leftSigmas = static_cast<int>(leftSidebandRange / sigma);// get the integral number
        std::cout << leftSigmas << " sigmas fit inside the left side range for the background distribution estimation." << std::endl;
        // Count for how many sigmas is there room inside the right side range
        double rightSidebandRange = h2_HFmass->ProjectionX("h_mass_")->GetBinLowEdge(h2_HFmass->ProjectionX("h_mass_")->GetNbinsX()+1) - (m_0 + startingBackSigma*sigma);
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
        double leftSB = fitBack->Integral(m_0 - (startingBackSigma + leftSigmas) * sigma,m_0 - startingBackSigma * sigma);
        double rightSB = fitBack->Integral(m_0 + startingBackSigma * sigma,m_0 + (startingBackSigma + rightSigmas) * sigma);
        // Finding side-band region: signal < |2*sigma|
        double backSigIntegral = fitBack->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma); // yield of background in signal region
        double fitSigIntegral = fitSigBack->Integral(m_0 - signalSigmas * sigma,m_0 + signalSigmas * sigma);

        // Calculating scaling parameter
        double alpha = backSigIntegral/(leftSB+rightSB);
        std::cout << "alpha = " << alpha << std::endl;
    
        // Create side-band histogram
        int lowBin = h2_HFmass->GetXaxis()->FindBin(m_0 - (startingBackSigma + leftSigmas) * sigma);
        int highBin = h2_HFmass->GetXaxis()->FindBin(m_0 - startingBackSigma * sigma);
        h_sideBand = h2_HFmass->ProjectionY("h_sideband_proj",lowBin, highBin); // start with left sideband
        lowBin = h2_HFmass->GetXaxis()->FindBin(m_0 + startingBackSigma * sigma);
        highBin = h2_HFmass->GetXaxis()->FindBin(m_0 + (startingBackSigma + rightSigmas) * sigma);
        TH1D* tempHist = h2_HFmass->ProjectionY("h_sideband_proj_temp", lowBin, highBin); // sum the right sideband
        h_sideBand->Add(tempHist);

        // Apply scaling to background function
        //h_sideBand->Scale(alpha);

        // Create signal histogram
        lowBin = h2_HFmass->GetXaxis()->FindBin(m_0 - signalSigmas * sigma);
        highBin = h2_HFmass->GetXaxis()->FindBin(m_0 + signalSigmas * sigma);
        h_signal = h2_HFmass->ProjectionY("h_signal_proj", lowBin, highBin);

        cout << "Distribution ratio of sideband over raw signal: " << h_sideBand->Integral() / h_signal->Integral() << endl;
        cout << "Fit function ratio of sideband over raw signal: " << (leftSB + rightSB)/ fitSigIntegral << endl;

        h_sideBand->Scale(alpha);

        // Subtract background histogram from signal histogram
        h_back_subtracted = (TH1D*)h_signal->Clone("h_back_subtracted");
        h_back_subtracted->Add(h_sideBand,-1.0);

        h_back_subtracted->Scale(1/0.9545);
        
    }

    std::cout << "Background subtracted." << std::endl;

    //
    // Plotting data
    //
    TCanvas* cInvMass = new TCanvas("cInvMass","Invariant mass plot with fits");
    cInvMass->Divide(2,2);
    cInvMass->cd(1);
    hInvMass->Draw();
    fitSigBack->Draw("same");
    fitBack->Draw("same");
    // Adding text with fit parameters and fit quality
    double chi2 = fitSigBack->GetChisquare();
    int ndf = fitSigBack->GetNDF();
    double chi2_ndf = chi2 / ndf;
    double fitProb = fitSigBack->GetProb();
    m_0 = fitSigBack->GetParameter(3);
    sigma = fitSigBack->GetParameter(4);
    TLatex latex;
    // Adding chi2/ndf and fit probability
    latex.DrawLatexNDC(0.6, 0.60, Form("#chi^{2}/ndf = %.2f", chi2_ndf));
    latex.DrawLatexNDC(0.6, 0.55, Form("Fit Prob = %f", fitProb));
    latex.DrawLatexNDC(0.6, 0.50, Form("m_{0} = %.2f", m_0));
    latex.DrawLatexNDC(0.6, 0.45, Form("#sigma = %.2f", sigma));

    cInvMass->cd(2);
    h_back_subtracted->Draw();
    cInvMass->cd(3);
    h_sideBand->Draw();
    cInvMass->cd(4);
    h_signal->Draw();
    

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

}

int main(){
    Investigation();
    return 0;
}



//
// Fitting parameter values
//
// 5-15 GeV/c jet pT initial fit parameters (deltaRbins = 500,minDeltaR = 0.,maxDeltaR = 1., massBins = 1500, minMass = 1.67, minMass = 2.1)
/*parametersVectors.paramA = {10455.4, 3442.01, 2922.09, 1147.83, 697.953, 1810.32, 789.441, 453.783, 1.};
parametersVectors.paramB = {2.84081, 3.41307, 3.50045, 3.89475, 3.98976+0.1, 3.38129, 3.34591, 3.68545, -1.};
parametersVectors.paramC = {2.73194e+06, 1.92476e+06, 2.1404e+06, 1.3353e+06, 891027-30, 1.14255e+06, 484831, 380790, 1.};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { 3.5*sigmaInitial, 3*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, sigmaInitial, 0.5*sigmaInitial, 3*sigmaInitial, 3*sigmaInitial, sigmaInitial };*/
// 5-15 GeV/c jet pT initial fit parameters (deltaRbins = 200, massBins = 200)
/*parametersVectors.paramA = {};
parametersVectors.paramB = {};
parametersVectors.paramC = {};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial };*/