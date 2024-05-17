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



void Investigation(){
    double jetptMin = 5;
    double jetptMax = 30;
    // deltaR histogram
    int deltaRbins = 50; // 20 bins
    double minDeltaR = 0.;
    double maxDeltaR = 0.5;
    int numberOfPoints = 10;
    // mass histogram
    int massBins = 100; 
    double minMass = 1.67;
    double maxMass = 2.1;
    double deltaRinterval = (maxDeltaR-minDeltaR)/numberOfPoints;
    cout << "deltaRinterval = " << deltaRinterval << endl;

    // creating histograms
    TH1D* hHfPt = new TH1D("hHfPt",";p_{T,D};dN/dp_{T,D}",1000,3.,30.);
    TH1D* hHfMass = new TH1D("hHfMass",";m_{0} (K#pi) GeV/c^{2};dN/dm_{0}",1000,0,10.);
    TH2D* h2_HFmass = new TH2D("h2_HFmass", ";m(K#pi) GeV/c^{2};#DeltaR", massBins, minMass, maxMass, deltaRbins, minDeltaR, maxDeltaR);
    TH1D* hDeltaR = new TH1D("hDeltaR",";#DeltaR;dN/d(#DeltaR)",1000,minDeltaR,maxDeltaR);
    // mass histograms
    std::vector<const char*> names = {"hDeltaRMass1", "hDeltaRMass2", "hDeltaRMass3", 
                                      "hDeltaRMass4", "hDeltaRMass5", "hDeltaRMass6",
                                      "hDeltaRMass7", "hDeltaRMass8", "hDeltaRMass9"}; // Names of histograms
    std::vector<const char*> titles = {"3 < p_{T,D} < 4 GeV/c;#DeltaR;counts", "4 < p_{T,D} < 5 GeV/c;#DeltaR;counts", "5 < p_{T,D} < 6 GeV/c;#DeltaR;counts",
                                       "6 < p_{T,D} < 7 GeV/c;#DeltaR;counts", "7 < p_{T,D} < 8 GeV/c;#DeltaR;counts", "8 < p_{T,D} < 10 GeV/c;#DeltaR;counts",
                                       "10 < p_{T,D} < 12 GeV/c;#DeltaR;counts", "12 < p_{T,D} < 15 GeV/c;#DeltaR;counts", "15 < p_{T,D} < 30 GeV/c;#DeltaR;counts"}; // Titles of histograms
    std::vector<TH1D*> hDeltaRMass;
    for (size_t i = 0; i < names.size(); ++i) {
        hDeltaRMass.push_back(new TH1D(names[i], titles[i],1000,minDeltaR,maxDeltaR));
        //hDeltaRMass[i]->Sumw2();
        
    }

    // opening file
    TFile* fDist = new TFile("AO2D.root","read");
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

        // Fill histograms for 15 < p_{T,jet} < 30 GeV/c
        if (jetPt > jetptMin && jetPt < jetptMax) {
            hHfPt->Fill(hfPt);
            if (hfPt > 3. && hfPt < 4.) {
                hHfMass->Fill(hfMass);
                h2_HFmass->Fill(hfMass, deltaR);
                hDeltaR->Fill(deltaR);
            }
            
        }

        // Fill each histogram with their respective pT intervals
        if ((hfPt > 3 && hfPt < 4) && (hDeltaRMass.size() > 0) && (jetPt > jetptMin && jetPt < jetptMax)) {
            hDeltaRMass[0]->Fill(deltaR);
        } else if ((hfPt > 4 && hfPt < 5) && (hDeltaRMass.size() > 1) && (jetPt > jetptMin && jetPt < jetptMax)) {
            hDeltaRMass[1]->Fill(deltaR);
        } else if ((hfPt > 5 && hfPt < 6) && (hDeltaRMass.size() > 2) && (jetPt > jetptMin && jetPt < jetptMax)) {
            hDeltaRMass[2]->Fill(deltaR);
        } else if ((hfPt > 6 && hfPt < 7) && (hDeltaRMass.size() > 3) && (jetPt > jetptMin && jetPt < jetptMax)) {
            hDeltaRMass[3]->Fill(deltaR);
        } else if ((hfPt > 7 && hfPt < 8) && (hDeltaRMass.size() > 4) && (jetPt > jetptMin && jetPt < jetptMax)) {
            hDeltaRMass[4]->Fill(deltaR);
        } else if ((hfPt > 8 && hfPt < 10) && (hDeltaRMass.size() > 5) && (jetPt > jetptMin && jetPt < jetptMax)) {
            hDeltaRMass[5]->Fill(deltaR);
        } else if ((hfPt > 10 && hfPt < 12) && (hDeltaRMass.size() > 6) && (jetPt > jetptMin && jetPt < jetptMax)) {
            hDeltaRMass[6]->Fill(deltaR);
        } else if ((hfPt > 12 && hfPt < 15) && (hDeltaRMass.size() > 7) && (jetPt > jetptMin && jetPt < jetptMax)) {
            hDeltaRMass[7]->Fill(deltaR);
        } else if ((hfPt > 15 && hfPt < 30) && (hDeltaRMass.size() > 8) && (jetPt > jetptMin && jetPt < jetptMax)) {
            hDeltaRMass[8]->Fill(deltaR);
        }
        
    }
    cout << "Histograms filled.\n";

    // creating temporary histogram for storing to the vector
    TH1D* tempHist;
    // Creating point histograms
    std::vector<TH1D*> histograms;
    // obtaining 1D invariant mass histograms from the projection
    
    // create numberOfPoints histograms for each selection cut (jet pT and HF pT)
    for (size_t iPoint = 0; iPoint < numberOfPoints; iPoint++) {
        double minBin = h2_HFmass->GetYaxis()->FindBin(deltaRinterval*iPoint);
        double maxBin = h2_HFmass->GetYaxis()->FindBin(deltaRinterval*(iPoint+1));
        cout << "minBin = " << minBin << endl;
        cout << "maxBin = " << maxBin << endl;
        cout << endl;
        tempHist = h2_HFmass->ProjectionX(Form("h_mass_proj_%zu", iPoint), minBin,maxBin);
        histograms.push_back(tempHist);
    }
        // histograms are filled in vector as
        // histograms = {histMass_0_0, histMass_0_1, histMass_0_2, ..., histMass_0_(numberOfPoints-1), histMass_1_0, histMass_1_1, histMass_1_2, ...}

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    // plotting histograms
    TCanvas* cHf = new TCanvas("cHfPt","D0 plots");
    cHf->Divide(2,2);
    cHf->cd(1);
    hHfPt->Draw();
    double statBoxPos = gPad->GetUxmax(); // Height of the stat box
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax)); // Display jet pT cut applied
    cHf->cd(2);
    hHfMass->Draw();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.60, Form("3 < p_{T,D} < 4 GeV/c"));
    cHf->cd(3);
    h2_HFmass->Draw("colz");
    cHf->cd(4);
    hDeltaR->Draw();

    TCanvas* cPointHists = new TCanvas("cPointHists","Point histograms for one set of selection cut");
    cPointHists->Divide(2,static_cast<int>(histograms.size() / 2));
    for (size_t iHisto = 0; iHisto < histograms.size(); iHisto++) {
        //
        cPointHists->cd(iHisto+1);
        histograms[iHisto]->SetTitle(Form("Point %zu",iHisto+1));
        histograms[iHisto]->Draw();
        double statBoxPos = gPad->GetUxmax();
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
        latex->DrawLatex(statBoxPos-0.35, 0.60, Form("3 < p_{T,D} < 4 GeV/c"));
    }
    
    TCanvas* cDeltaRMasses = new TCanvas("cDeltaRMasses","Delta R histograms for different HF mass intervals");
    cDeltaRMasses->Divide(3,static_cast<int>(hDeltaRMass.size() / 3)); // columns, lines

    // Loop through all histograms and fitting functions
    for(size_t iHisto = 0; iHisto < hDeltaRMass.size(); iHisto++) {
        //
        cDeltaRMasses->cd(iHisto+1);
        //hDeltaRMass[iHisto]->SetMarkerStyle(kFullDotMedium);
        //hDeltaRMass[iHisto]->SetMarkerColor(kBlack);
        //hDeltaRMass[iHisto]->SetLineColor(kGray);
        //hDeltaRMass[iHisto]->GetYaxis()->SetTitle("counts");
        hDeltaRMass[iHisto]->Draw();
    }

    // Cleanup (if done, the histograms disappear from the canvas)
    //for (auto hist : histograms) {
    //    delete hist;
    //}
    //fDist->Close();
    //delete fDist;
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