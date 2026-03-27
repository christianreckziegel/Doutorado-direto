
/**
 * @brief Manual folding of TH2D data using a RooUnfoldResponse object.
 *
 * This helper function perfom the manual matrix multiplication of the 2D true distribution by the response matrix
 * It uses row-major bin index flattening, since RooUnfoldResponse is a 2D structure.
 * Bins should be accessed through: binx + nx*(biny + ny*binz).
 *
 * @param response The response matrix.
 * @param hTruth The 2D truth level distribution to be folded.
 *
 * @return The folded TH2D distribution.
 *
 * @note uses 0-based indexing for calculation but 1-based indexing for ROOT histograms.
 *
 * @see smearGeneratorData() [Makes use of Nima's folding method.]
 */
TH2D* nimaFolding(RooUnfoldResponse response, TH2D* hTruth, TH2D* hMeasured) {
    std::cout << "Nima folding:\n";
    // Create empty histogram for the folded data
    TH2D* hFolded = (TH2D*)hMeasured->Clone("hFolded");
    hFolded->Reset();

    // Get the number of bins in measured and truth histograms
    int nBinsXMeasured = hFolded->GetNbinsX();
    int nBinsYMeasured = hFolded->GetNbinsY();
    std::cout << "Measured histogram: " << nBinsXMeasured << " x " << nBinsYMeasured << std::endl;
    int nBinsXTruth = hTruth->GetNbinsX();
    int nBinsYTruth = hTruth->GetNbinsY();
    std::cout << "Measured histogram: " << nBinsXTruth << " x " << nBinsYTruth << std::endl;
    // Debug: print a few values from the response matrix
    bool debug = true;
    if (debug) {
        std::cout << "Response matrix dimensions: " 
                << response.GetNbinsMeasured() << " x " << response.GetNbinsTruth() << std::endl;
    }
    

    // loop through detector level bins
    for (int iMeasured = 0; iMeasured < nBinsXMeasured; iMeasured++) {
        for (int jMeasured = 0; jMeasured < nBinsYMeasured; jMeasured++) {
            double foldedValue = 0;
            double foldedError2 = 0;

            // obtaining flattened 1D index through row-major ordering
            int index_x_measured = iMeasured + nBinsXMeasured*jMeasured; // Nima version
            //int index_x_measured = iMeasured*nBinsYMeasured + jMeasured; // Christian version

            // calculating element iMeasured,jMeasured of folded 2D matrix
            for (int iTruth = 0; iTruth < nBinsXTruth; iTruth++) {
                for (int jTruth = 0; jTruth < nBinsYTruth; jTruth++) {
                    // obtaining flattened 1D index through row-major ordering
                    int index_x_truth = iTruth + nBinsXTruth*jTruth; // Nima version
                    //int index_x_truth = iTruth*nBinsYTruth + jTruth; // Christian version
                    
                    //std::cout << "Accessing true value\n";
                    // calculating matrix element product
                    double truthValue = hTruth->GetBinContent(iTruth + 1,jTruth + 1);
                    //std::cout << "Accessing response value\n";
                    double responseValue = response(index_x_measured, index_x_truth);
                    
                    //double responseValue = response(index_x_truth, index_x_measured); // in case it is transposed
                    foldedValue += truthValue*responseValue;
                    //std::cout << "Flag 2: calculated folded error to the power of 2\n";
                    //foldedValue = hTruth->GetBinContent(iTruth + 1,jTruth + 1) * response(index_x_measured, index_x_truth);
                    foldedError2 = std::pow(hTruth->GetBinError(iTruth + 1,jTruth + 1),2) * std::pow(response(index_x_measured, index_x_truth),2);
                    std::cout << "Response indexes (" << index_x_measured << "," << index_x_truth << ") has value " << response(index_x_measured, index_x_truth) << "\n";
                    std::cout << "Truth value is " << truthValue << std::endl;
                    std::cout << "and folded value of " << foldedValue << std::endl;
                }
            }
            hFolded->SetBinContent(iMeasured + 1, jMeasured + 1, foldedValue);
            std::cout << "bin: (" << iMeasured + 1 << "," << jMeasured +1 << ") = " << foldedValue << std::endl << std::endl;
            //hFolded->SetBinContent(nBinsXMeasured - (iMeasured + 1), nBinsYMeasured - (jMeasured + 1), foldedValue); // axis inversion
            hFolded->SetBinError(iMeasured + 1, jMeasured + 1, std::sqrt(foldedError2));
        }
    }
    
    return hFolded;
}

// calculate delta phi such that 0 < delta phi < 2*pi
double DeltaPhi(double phi1, double phi2) {
    // Compute the absolute difference between phi1 and phi2
    double dphi = std::abs(phi1 - phi2); 
    if (dphi > M_PI) {
        // subtract 2pi if the difference if bigger than pi
        dphi = dphi - 2*M_PI;
    }

    return dphi;
}

void Testing() {
    // Defining cuts
    double jetptMin = 5; // GeV, jet pT cuts
    double jetptMax = 30; // GeV
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on jet
    const double yCut = 0.8; // on D0

    // Defining histogram bin edges
    // x -> delta R
    // y -> pT,D
    // z -> pT,jet
    std::vector<double> zBinEdges_particle = {5., 7., 15., 30.};
    std::vector<double> zBinEdges_detector = {5., 7., 15., 30.};
    std::vector<double> xBinEdges_particle = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5}; // TODO: investigate structure before 0.005: 0.,0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.4
    std::vector<double> xBinEdges_detector = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5};
    std::vector<double> yBinEdges_particle = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};
    std::vector<double> yBinEdges_detector = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};
    // Testing edges
    /*std::vector<double> zBinEdges_particle = {5., 15., 30.};
    std::vector<double> zBinEdges_detector = {5., 15., 30.};
    std::vector<double> xBinEdges_particle = {0., 0.3, 0.5}; // TODO: investigate structure before 0.005: 0.,0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.4
    std::vector<double> xBinEdges_detector = {0., 0.3, 0.5};
    std::vector<double> yBinEdges_particle = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};
    std::vector<double> yBinEdges_detector = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};*/

    // Opening files
    TFile* fPowheg = new TFile("../SimulatedData/POWHEG/trees_powheg_fd_central.root","read");
    TFile* fSimulatedO2 = new TFile("../SimulatedData/Hyperloop_output/McChargedMatched/AO2D_merged_All.root","read");
    TFile* fEfficiency = new TFile(Form("../2-Efficiency/backSubEfficiency_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");

    // Creating histograms
    int xNumBinEdges = xBinEdges_particle.size();
    int yNumBinEdges = yBinEdges_particle.size();
    int zNumBinEdges = zBinEdges_particle.size();
    TH3D* hPowheg = new TH3D("h_deltaR_vs_pt", "POWHEG + PYTHIA;#DeltaR;p_{T,D}^{gen};p_{T,jet}^{ch}", xNumBinEdges-1, xBinEdges_particle.data(), 
                                                                                                            yNumBinEdges-1, yBinEdges_particle.data(), 
                                                                                                            zNumBinEdges-1, zBinEdges_particle.data());
    TH2D* hTruth = new TH2D("hTruth2D", "Truth inside response ranges;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges_particle.data(), zNumBinEdges-1, zBinEdges_particle.data());
    xNumBinEdges = xBinEdges_detector.size();
    yNumBinEdges = yBinEdges_detector.size();
    zNumBinEdges = zBinEdges_detector.size();
    TH2D* hMeasured = new TH2D("hMeasured2D", "Measured inside response ranges;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges_detector.data(), zNumBinEdges-1, zBinEdges_detector.data());
    TH2D* hFolded = new TH2D("hFolded2D", "Folded;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges_detector.data(), zNumBinEdges-1, zBinEdges_detector.data());

    TH2D* hJetPtMatched = new TH2D("hJetPtMatched",";p_{T,jet}^{det};p_{T,jet}^{part}",2000,0.,200,2000,0.,200);
    TH2D* hDeltaRMatched = new TH2D("hDeltaRMatched",";#DeltaR^{det};#DeltaR^{part}",100,0.,0.4,100,0.,0.4);

    // Luminosity (for now arbitrary)
    double luminosity = 1000000;

    //
    // MC generator level tree and histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)fPowheg->Get("tree_D0");

    // Check for correct access
    if (!tree) {
        cout << "Error opening POWHEG tree.\n";
    }
    
    // defining variables for accessing data on TTree
    double pt_cand, eta_cand, phi_cand, y_cand;
    double pt_jet, eta_jet, phi_jet, delta_r_jet;

    tree->SetBranchAddress("pt_cand",&pt_cand);
    tree->SetBranchAddress("eta_cand",&eta_cand);
    tree->SetBranchAddress("phi_cand",&phi_cand);
    tree->SetBranchAddress("y_cand",&y_cand);
    tree->SetBranchAddress("pt_jet",&pt_jet);
    tree->SetBranchAddress("eta_jet",&eta_jet);
    tree->SetBranchAddress("phi_jet",&phi_jet);
    tree->SetBranchAddress("delta_r_jet",&delta_r_jet);

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        //double deltaR = sqrt(pow(eta_jet-eta_cand,2) + pow(DeltaPhi(phi_jet,phi_cand),2));

        // Fill 2D histogram considering jet pT and detector acceptance
        if ((abs(eta_cand) < etaCut) && (abs(y_cand) < yCut) && (pt_jet > jetptMin) && (pt_jet < jetptMax)) {
            
            hPowheg->Fill(delta_r_jet, pt_cand, pt_jet);
            
        }
        
    }
    cout << "Generator level (POWHEG+PYTHIA) histograms filled.\n";


    //
    // O2 matching task measured and truth histograms
    //
    // Accessing TTree
    tree = (TTree*)fSimulatedO2->Get("O2matchtable");
    // Check for correct access
    if (!tree) {
        cout << "Error opening O2 matching tree.\n";
    }

    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;

    // particle level branches
    tree->SetBranchAddress("fMCJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMCJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMCJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMCJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMCHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMCHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMCHfPhi",&MCPhfPhi);
    MCPhfMass = 1.86483; // D0 rest mass in GeV/c^2
    tree->SetBranchAddress("fMCHfY",&MCPhfY);
    tree->SetBranchAddress("fMCHfPrompt",&MCPhfprompt);
    // detector level branches
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fHfPt",&MCDhfPt);
    tree->SetBranchAddress("fHfEta",&MCDhfEta);
    tree->SetBranchAddress("fHfPhi",&MCDhfPhi);
    tree->SetBranchAddress("fHfMass",&MCDhfMass);
    tree->SetBranchAddress("fHfY",&MCDhfY);
    tree->SetBranchAddress("fHfPrompt",&MCDhfprompt);

    // Create response matrix
    RooUnfoldResponse response(hMeasured, hTruth);
    

    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        // Fill histograms considering jet pT and detector acceptance for NON-PROMPT particles
        if ((abs(MCPjetEta) < etaCut) && (abs(MCPhfY) < yCut) && (MCPjetPt > jetptMin) && (MCPjetPt < jetptMax) && !MCPhfprompt &&
            (abs(MCDjetEta) < etaCut) && (abs(MCDhfY) < yCut) && (MCDjetPt > jetptMin) && (MCDjetPt < jetptMax) && !MCDhfprompt) {
            // Filling measured 2D histogram
            hMeasured->Fill(MCDDeltaR, MCDjetPt);

            // Filling truth 2D histogram
            hTruth->Fill(MCPDeltaR, MCPjetPt);
            // Fill 4D RooUnfoldResponse object
            response.Fill(MCDDeltaR, MCDjetPt, MCPDeltaR, MCPjetPt);
        
        }
        
        // Fill matching histograms, sanity check
        hJetPtMatched->Fill(MCDjetPt,MCPjetPt);
        hDeltaRMatched->Fill(MCDDeltaR,MCPDeltaR);
    }

    cout << "Response matched histograms filled.\n";

    //
    // Smear data
    //
    TH2D* hAllptDPow = static_cast<TH2D*>(hPowheg->Project3D("zx"));
    hAllptDPow->Scale(1/luminosity);

    TH1D* hEffPrompt;
    TH1D* hEffNonPrompt;
    hEffPrompt = (TH1D*)fEfficiency->Get("efficiency_prompt");
    hEffNonPrompt = (TH1D*)fEfficiency->Get("efficiency_nonprompt");
    if (!hEffPrompt || !hEffNonPrompt) {
        std::cerr << "Error: could not retrieve efficiency histograms.\n";
    }
    
    // manually scaling over each bin of pT,D dimension of 2D histogram (y-axis)
    int nBinsX = hPowheg->GetNbinsX();
    int nBinsY = hPowheg->GetNbinsY();
    int nBinsZ = hPowheg->GetNbinsZ();
    // the order of nested loop dictates which is the scaled axis
    for (int xBin = 1; xBin <= nBinsX; xBin++) {
        for (int zBin = 1; zBin <= nBinsZ; zBin++) {
            // Inner loop over y-axis bins: this is scaled axis (pT,D bins)
            for (int yBin = 1; yBin <= nBinsY; yBin++) {
                double binContent = hPowheg->GetBinContent(xBin, yBin, zBin);
                double binError = hPowheg->GetBinError(xBin, yBin, zBin);

                // Obtain efficiencies ratio for specific bin
                double effPrompt = hEffPrompt->GetBinContent(yBin);
                double effNonPrompt = hEffNonPrompt->GetBinContent(yBin);

                if (effPrompt == 0) {
                    std::cerr << "Error: null prompt efficiency value for bin (" << xBin << "," << yBin << "," << zBin << ").\n";
                }
                

                // Scale the y-axis content and error by non-prompt efficiency/prompt efficiency
                hPowheg->SetBinContent(xBin, yBin, zBin, binContent * effNonPrompt / effPrompt);
                hPowheg->SetBinError(xBin, yBin, zBin, binError * effNonPrompt / effPrompt);
            }
        }
    }// end of efficiency scaling

    // Folding
    hFolded = nimaFolding(response, hAllptDPow, hMeasured);
    hFolded->SetTitle("Folded with Nima's 1D index function;#frac{1}{L_{int}}#DeltaR^{b #rightarrow D^{0}}_{reco};p_{T,jet}^{ch}");

    TCanvas* cResponse = new TCanvas("cResponse","Response plots");
    cResponse->Divide(2,3);
    cResponse->cd(1);
    const TH2* hResponse2D = response.Hresponse();
    TH2D* hResponse2DClone = static_cast<TH2D*>(hResponse2D->Clone("hResponse2DClone"));
    std::cout << "hResponse2DClone has " << hResponse2DClone->GetEntries() << " bins." << std::endl;
    hResponse2DClone->SetTitle("2D response matrix from method 1;2D Reconstructed;2D Truth");
    hResponse2DClone->Draw("colz");
    cResponse->cd(3);
    hAllptDPow->Draw("colz");
    
    cResponse->cd(4);
    hFolded->Draw("colz");
    cResponse->cd(5);
    hTruth->Draw("colz");
    cResponse->cd(6);
    hMeasured->Draw("colz");

    gStyle->SetPalette(kRainbow);
    TCanvas* cMatching = new TCanvas("cMatching","Sanity matching plots");
    cMatching->Divide(2,2);
    cMatching->cd(1);
    hJetPtMatched->Draw("colz");
    cMatching->cd(2);
    hDeltaRMatched->Draw("colz");

}

int main() {
    Testing();
    return 0;
}
