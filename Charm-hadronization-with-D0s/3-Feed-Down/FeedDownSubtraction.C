/*
 *
 *
 * Macro for performing B feed-down subtraction 
 * and applying correction to BackgroundSubtraction.C
 * and SignalExtraction.C resulting distributions.
 * 
 * 
 * 
 * 
**/


using namespace std;

// calculate number of objects inside file
int HistogramCounter(TFile* file) {
    // Get the list of keys (i.e., objects) in the file
    TList* keys = file->GetListOfKeys();

    int numHistograms = 0;

    // Loop over all keys in the file
    for (int i = 0; i < keys->GetSize(); ++i) {
        TKey* key = (TKey*)keys->At(i);
        TObject* obj = file->Get(key->GetName());

        // Check if the object is a histogram
        if (obj->IsA()->InheritsFrom(TH1::Class())) {
            numHistograms++;
        }
    }

    return numHistograms;
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

struct FeedDownData {
    std::pair<TH2D*, TH1D*> hMeasured;                      // 2D and 1D representation of response measured data
    std::pair<TH2D*, TH1D*> hTruth;                         // 2D and 1D representation of response truth data
    std::pair<TH2D*, TH1D*> hFolded;                        // 2D and 1D representation of folded data from method 1
    RooUnfoldResponse response;                             // response matrix for non-prompt D0s only, for overall pT,D; using RooUnfoldResponse object, method 1
    THnSparseD* hResponse;                                  // response matrix for non-prompt D0s only, for overall pT,D; using manual folding operation (matrix multiplication), method 2
    std::vector<TH3D*> hPowheg;                             // deltaR vs pT,D vs pT,jet; [0] = generator level, [1] = treated but not yet detector level
    std::vector<TH2D*> hAllptDPowheg;                       // deltaR vs pT,jet; [0] = generator level, [1] = method 1 folded, [2] = method 2 folded
    std::vector<TH1D*> hEfficiencies;                       // inclusive = 0, prompt only = 1, non-prompt only = 2
    TH1D* hBackSubCorrected;                                // prompt efficiency corrected Delta R distributions
    TH1D* hSBFeedDownSubtracted;                            // non-prompt subtracted Delta R distribution, all pT,D
};

// Obtain the bin edges of a histogram (useful for asymmetrical bin sizes)
std::vector<double> getBinEdges(const TAxis* axis) {
    
    if (!axis) {
        std::cerr << "Error: null axis pointer.\n";
    }
    
    int nBins = axis->GetNbins();

    // vector for storing bin edges
    std::vector<double> binEdges(nBins + 1);

    for (int iBin = 0; iBin <= nBins; iBin++) {
        binEdges[iBin] = axis->GetBinLowEdge(iBin + 1);
    }

    return binEdges;
}

// Module to create TH2D histograms including interest variable
FeedDownData createHistograms(const std::vector<double>& xBinEdges, 
                              const std::vector<double>& yBinEdges, 
                              const double& jetptMin, const double& jetptMax) {
    // Create struct to store data
    FeedDownData dataContainer;

    int xNumBinEdges = xBinEdges.size();
    int yNumBinEdges = yBinEdges.size();
    // Obtain fixed width zBinEdges
    int zNumBinEdges = 26;
    std::vector<double> zBinEdges(zNumBinEdges); // Define bin edges for the z-axis
    double binWidth = (jetptMax - jetptMin) / (zNumBinEdges-1);
    for (int iBin = 0; iBin <= zNumBinEdges; ++iBin) {
        zBinEdges[iBin] = jetptMin + iBin * binWidth;
    }

    // Create 2D histogram
    dataContainer.hPowheg.emplace_back(new TH3D("h_deltaR_vs_pt", "POWHEG + PYTHIA;#DeltaR;p_{T,D}^{gen};p_{T,jet}^{ch}", xNumBinEdges-1, xBinEdges.data(), 
                                                                                                            yNumBinEdges-1, yBinEdges.data(), 
                                                                                                            zNumBinEdges-1, zBinEdges.data()));
    dataContainer.hPowheg[0]->SetMarkerColor(kGreen);
    dataContainer.hPowheg[0]->SetLineColor(kGreen);
    dataContainer.hPowheg[0]->SetMarkerStyle(kCircle);
    dataContainer.hPowheg[0]->Sumw2();
    dataContainer.hPowheg[0]->SetStats(0);

    cout << "POWHEG histograms created.\n";

    //
    // Matching histograms for folding process
    //
    // 2D original information
    dataContainer.hMeasured.first = new TH2D("hMeasured2D", "Measured;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges.data(), zNumBinEdges-1, zBinEdges.data());
    dataContainer.hTruth.first = new TH2D("hTruth2D", "Truth;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges.data(), zNumBinEdges-1, zBinEdges.data());
    dataContainer.hFolded.first = new TH2D("hFolded2D", "Folded;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges.data(), zNumBinEdges-1, zBinEdges.data());

    // 1D flattened (binLow and High need correction)
    int totalBins = (xBinEdges.size() - 1) * (zBinEdges.size() - 1);
    dataContainer.hMeasured.second = new TH1D("hMeasured1D", "Measured flattened;#DeltaR,p_{T,jet}", totalBins, 0, totalBins);
    dataContainer.hTruth.second = new TH1D("hTruth1D", "Truth flattened;#DeltaR,p_{T,jet}", totalBins, 0, totalBins);
    dataContainer.hFolded.second = new TH1D("hFolded1D", "Folded flattened;#DeltaR,p_{T,jet}", totalBins, 0, totalBins);

    // Set bin labels for 1D histograms
    for (int i = 0; i < xNumBinEdges - 1; ++i) {
        for (int j = 0; j < zNumBinEdges - 1; ++j) {
            int flatBin = i * (zNumBinEdges - 1) + j;
            std::string label = Form("#DeltaR[%.2f,%.2f],p_{T,jet}[%.1f,%.1f]", 
                                     xBinEdges[i], xBinEdges[i+1], 
                                     zBinEdges[j], zBinEdges[j+1]);
            dataContainer.hMeasured.second->GetXaxis()->SetBinLabel(flatBin + 1, label.c_str());
            dataContainer.hTruth.second->GetXaxis()->SetBinLabel(flatBin + 1, label.c_str());
            dataContainer.hFolded.second->GetXaxis()->SetBinLabel(flatBin + 1, label.c_str());
        }
    }
    cout << "Matching histograms created.\n";

    //
    // Response matrix as a 4D histogram object (method 2)
    //
    std::vector<double> xBinsVec = getBinEdges(dataContainer.hMeasured.first->GetXaxis());
    std::vector<double> yBinsVec = getBinEdges(dataContainer.hMeasured.first->GetYaxis());
    int nBins[4] = {
        static_cast<int>(xBinsVec.size() - 1),
        static_cast<int>(yBinsVec.size() - 1),
        static_cast<int>(xBinsVec.size() - 1),
        static_cast<int>(yBinsVec.size() - 1)
    };

    // Set the bin edges for each dimension
    int numDim = 4; // number of dimensions of histogram
    dataContainer.hResponse = new THnSparseD("hResponse", "4D response matrix", numDim, nBins, nullptr, nullptr);
    dataContainer.hResponse->GetAxis(0)->Set(nBins[0], xBinsVec.data());
    dataContainer.hResponse->GetAxis(1)->Set(nBins[1], yBinsVec.data());
    dataContainer.hResponse->GetAxis(2)->Set(nBins[2], xBinsVec.data());
    dataContainer.hResponse->GetAxis(3)->Set(nBins[3], yBinsVec.data());

    // Set axis titles
    dataContainer.hResponse->GetAxis(0)->SetTitle("#DeltaR_{measured}");
    dataContainer.hResponse->GetAxis(1)->SetTitle("p_{T,jet}^{measured}");
    dataContainer.hResponse->GetAxis(2)->SetTitle("#DeltaR_{truth}");
    dataContainer.hResponse->GetAxis(3)->SetTitle("p_{T,jet}^{truth}");


    return dataContainer;
}

// Module to fill histograms from 
// POWHEG+PYTHIA TFile data (all data in file is non-prompt and on particle level)
// and
// O2 matching task simulation
void fillHistograms(TFile* fPowheg, TFile* fSimulatedO2, FeedDownData& dataContainer, double jetptMin, double jetptMax) {
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on jet
    const double yCut = 0.8; // on D0
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
            
            dataContainer.hPowheg[0]->Fill(delta_r_jet, pt_cand, pt_jet);
            
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

    nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        // Fill histograms considering jet pT and detector acceptance
        if ((abs(MCPhfEta) < etaCut) && (abs(MCPhfY) < yCut) && (MCPjetPt > jetptMin) && (MCPjetPt < jetptMax)) {
            // Filling measured 2D histogram
            dataContainer.hMeasured.first->Fill(MCDDeltaR, MCDjetPt);

            // Filling truth 2D histogram
            dataContainer.hTruth.first->Fill(MCPDeltaR, MCPjetPt);

            // Fill 4D response matrix
            dataContainer.hResponse->Fill(MCDDeltaR, MCDjetPt, MCPDeltaR, MCPjetPt);
        
        }
        
    }

    cout << "Response matched histograms filled.\n";

}

/**
 * @brief Flattens 2D histogram into a 1D histogram.
 *
 * This helper function converts 2D bin indice, content and error into a single 1D index, content and error respectivelly.
 * It uses row-major flattening, which means the 2D histogram is flattened row by row (line by line).
 * This is useful when working with response matrices that can only handle 1D entries.
 *
 * @param hist2D The 2D histogram object to be flattened.
 * @param hist1D The resulting flattened 1D histogram object.
 *
 * @return The bin index of the new flattened 1D distribution.
 *
 * @note uses 0-based indexing for calculation but 1-based indexing for ROOT histograms.
 *
 * @see unflattenHistogram() [Corresponding function to convert back to 2D] and buildResponseMatrix() [Uses flattened histograms]
 */
void flattenHistogram(TH2D* hist2D, TH1D* hist1D) {
    
    std::cout << "Flattening 2D histogram.\n";
    
    // Get number of bins
    int nBinsX = hist2D->GetNbinsX();
    int nBinsY = hist2D->GetNbinsY();
    int totalBins = nBinsX * nBinsY;

    std::cout << "Original 2D histogram had " << nBinsX << " bins in x axis and " << nBinsY << " bins in y axis.\n";
    std::cout << "New 1D histogram has " << totalBins << " bins in one single axis." << std::endl << std::endl;

    // Loop through 2D histogram bins
    for (int binX = 1; binX <= nBinsX; binX++) {
        for (int binY = 1; binY <= nBinsY; binY++) {
            
            // maps a 2D coordinate to a 1D index in a row-major order
            int flatBin = (binX - 1) * nBinsY + (binY - 1);

            // Get content of bin from 2D histogram
            double content = hist2D->GetBinContent(binX, binY);
            double error = hist2D->GetBinError(binX, binY);

            // Set content of bin from 1D histogram
            hist1D->SetBinContent(flatBin + 1,content);
            hist1D->SetBinError(flatBin + 1, error);

            // Already done in histogram creation: set bin label of 1D flattened histogram based on 2D histogram
            /*double xLow = hist2D->GetXaxis()->GetBinLowEdge(binX);
            double xUp = hist2D->GetXaxis()->GetBinUpEdge(binX);
            double yLow = hist2D->GetYaxis()->GetBinLowEdge(binY);
            double yUp = hist2D->GetYaxis()->GetBinUpEdge(binY);
            std::string label = Form("x[%.2f,%.2f],y[%.2f,%.2f]",xLow,xUp,yLow,yUp);
            hist1D->GetXaxis()->SetBinLabel(flatBin + 1, label.c_str());*/
        }
        
    }
    
}

void unflattenHistogram(TH1D* hist1D, TH2D* hist2D) {
    std::cout << "Unflattening histogram.\n";
    
    // Get number of bins
    int nBinsX = hist2D->GetNbinsX();
    int nBinsY = hist2D->GetNbinsY();

    if (hist1D->GetNbinsX() != nBinsX*nBinsY) {
        std::cerr << "Error: 1D histogram bin count does not match 2D histogram total bins.\n";
        std::cout << "Original 1D histogram had " << hist1D->GetNbinsX() << " bins in one single axis.\n";
        std::cout << "New 2D histogram has " << nBinsX << " bins in x axis and " << nBinsY << " bins in y axis." << std::endl << std::endl;
    }
    
    

    for (int binX = 1; binX <= nBinsX; binX++) {
        for (int binY = 1; binY <= nBinsY; binY++) {
            
            // maps a 2D coordinate to a 1D index in a row-major order
            int flatBin = (binX - 1) * nBinsY + (binY - 1);

            // Get content of bin from 1D histogram
            double content = hist1D->GetBinContent(flatBin + 1);
            double error = hist1D->GetBinError(flatBin + 1);

            // Set content of bin from 2D histogram
            hist2D->SetBinContent(binX, binY, content);
            hist2D->SetBinError(binX, binY, error);

            
        }
        
    }
    
}

// Module to build 2D response matrix out of flattened 1D input data
void buildResponseMatrix(FeedDownData& dataContainer) {
    
    // Flattening 2D histograms (RooUnfold only supports 1D distributions)
    flattenHistogram(dataContainer.hMeasured.first, dataContainer.hMeasured.second); // input = 2D, output = 1D
    flattenHistogram(dataContainer.hTruth.first, dataContainer.hTruth.second);

    std::cout << "Measured and truth histograms flattened.\n";

    // method 1: create 2D response matrix of non-prompt flattened Delta R and pT,jet, for overall pT,D
    //(DeltaR_detector, pTjet_detector, DeltaR_particle, pTjet_particle) -> flattened to (detector, particle)
    dataContainer.response = RooUnfoldResponse(dataContainer.hMeasured.second, dataContainer.hTruth.second);

    std::cout << "Response matrix created.\n";

    // Fill response matrix correlating measured and truth flattened distributions
    int totalBins = dataContainer.hTruth.second->GetNbinsX();
    if (dataContainer.hTruth.second->GetNbinsX() != dataContainer.hMeasured.second->GetNbinsX()) {
        std::cout << "WARNING: measured and truth flattened histograms have different number of bins.\n";
    }
    
    // Get direct access to bin contents for efficiency
    const double* measuredArray = dataContainer.hMeasured.second->GetArray();
    const double* truthArray = dataContainer.hTruth.second->GetArray();

    // Fill correlated values for method 1
    for (int iBin = 1; iBin <= totalBins; iBin++) {
        // Get flattened values
        double measured = measuredArray[iBin];
        double truth = truthArray[iBin]; // dataContainer.hTruth.second->GetBinContent(iBin);

        dataContainer.response.Fill(measured, truth);
        //std::cout << "Filling response matrix: measured = " << measured << ", truth = " << truth << std::endl;
    }    
    
    std::cout << "Response matrix filled.\n";
}

// Perform manual matrix multiplication of the 4D response matrix by the 2D truth level distribution used by smearGeneratorData()
TH2D* manualFolding(THnSparseD* hResponse, TH2D* hTruth) {

    // Create empty histogram for the folded data
    TH2D* hFolded = (TH2D*)hTruth->Clone("hFolded");
    hFolded->Reset();

    // Get the number of bins in measured and truth histograms
    int nBinsXMeasured = hFolded->GetNbinsX();
    int nBinsYMeasured = hFolded->GetNbinsY();
    int nBinsXTruth = hTruth->GetNbinsX();
    int nBinsYTruth = hTruth->GetNbinsY();

    // Check bin consistency
    if (nBinsXTruth != hResponse->GetAxis(2)->GetNbins() || nBinsYTruth != hResponse->GetAxis(3)->GetNbins()) {
        std::cerr << "Error: truth histogram dimensions do not match response matrix dimensions.\n";
        std::cout << "Truth X histogram has " << nBinsXTruth << " while response matrix 3rd axis has " << hResponse->GetAxis(2)->GetNbins() << " bins.\n";
        std::cout << "Truth Y histogram has " << nBinsYTruth << " while response matrix 3rd axis has " << hResponse->GetAxis(3)->GetNbins() << " bins.\n";
    }
    
    // Prepare arrays for faster bin access
    double truthContent[nBinsXTruth + 2][nBinsYTruth + 2];
    for (int xBin = 0; xBin <= nBinsXTruth + 1; xBin++) {
        for (int yBin = 0; yBin <= nBinsYTruth + 1; yBin++) {
            truthContent[xBin][yBin] = hTruth->GetBinContent(xBin,yBin);
        }
        
    }
    

    // Perform the folding operation manually
    for (int iMeasured = 1; iMeasured <= nBinsXMeasured; iMeasured++) {
        for (int jMeasured = 1; jMeasured <= nBinsYMeasured; jMeasured++) {
            // variable for calculating truth folded value
            double foldedValue = 0.;
            double foldedError2 = 0.;

            // calculating element iMeasured,jMeasured of folded 2D matrix
            for (int iTruth = 1; iTruth <= nBinsXTruth; iTruth++) {
                for (int jTruth = 1; jTruth <= nBinsYTruth; jTruth++) {
                    //double truthValue = hTruth->GetBinContent(iTruth, jTruth);
                    double truthValue = truthContent[iTruth][jTruth];
                    int bins[4] = {iMeasured, jMeasured, iTruth, jTruth};
                    double responseValue = hResponse->GetBinContent(bins);
                    foldedValue += truthValue * responseValue;

                    // Error propagation (assuming Poisson statistics)
                    foldedError2 += truthValue* responseValue * responseValue;
                }
                
            }

            // store it in the folded matrix as the bin content
            hFolded->SetBinContent(iMeasured, jMeasured, foldedValue);
            hFolded->SetBinError(iMeasured, jMeasured, sqrt(foldedError2));
        }
        
    }

    return hFolded;
}

// Module for folding particle level data from POWHEG simulation
void smearGeneratorData(FeedDownData& dataContainer, double& luminosity, TFile* fEfficiency) {
    //
    // 0th step: clone 3D histogram in order to save the smeared at the end (while keeping the original)
    //
    TH3D* hTreatedPowheg = static_cast<TH3D*>(dataContainer.hPowheg[0]->Clone("h_subtraction_Bs"));
    hTreatedPowheg->SetTitle("POWHEG + PYTHIA efficiency scaled;#frac{1}{L_{int}}#DeltaR^{b #rightarrow D^{0}};#frac{#epsilon_{non-prompt}}{prompt}#frac{1}{L_{int}}p_{T,D}^{b #rightarrow D^{0}};p_{T,jet}^{ch}");

    //
    // 1st step: obtain the Delta R vs pT,jet projection for all pT,D bins before folding
    //
    // using a smart pointer provided by the C++ Standard Library
    int nBinsY = hTreatedPowheg->GetNbinsY();
    TH2D* hAllptDPow = static_cast<TH2D*>(hTreatedPowheg->Project3D("zx"));
    hAllptDPow->SetName("h_true_nonprompt");
    //int binYMin = hTreatedPowheg->GetYaxis()->FindBin(1);
    //int binYMax = hTreatedPowheg->GetYaxis()->FindBin(nBinsY);
    //hTreatedPowheg->GetYaxis()->SetRange(binYMin, binYMax);
    //hAllptDPow->SetTitle(";#frac{1}{L_{int}}#DeltaR^{b #rightarrow D^{0}}_{truth};p_{T,jet}^{ch}");
    hAllptDPow->SetTitle("Before folding;#DeltaR^{b #rightarrow D^{0}}_{truth};p_{T,jet}^{ch}");
    // transfer ownership of the histogram from the unique_ptr to dataContainer object
    dataContainer.hAllptDPowheg.emplace_back(hAllptDPow); // [0] = not folded

    //
    // 2nd step: scale by integrated luminosity
    // (skip this step for now)
    //hTreatedPowheg->Scale(1/luminosity);

    //
    // 3rd step: scale by efficiency ratio
    //
    TH1D* hEffPrompt;
    TH1D* hEffNonPrompt;
    hEffPrompt = (TH1D*)fEfficiency->Get("efficiency_prompt");
    hEffNonPrompt = (TH1D*)fEfficiency->Get("efficiency_nonprompt");
    if (!hEffPrompt || !hEffNonPrompt) {
        std::cerr << "Error: could not retrieve efficiency histograms.\n";
    }
    
    // manually scaling over each bin of pT,D dimension of 2D histogram (y-axis)
    int nBinsX = hTreatedPowheg->GetNbinsX();
    //int nBinsY = hTreatedPowheg->GetNbinsY();
    int nBinsZ = hTreatedPowheg->GetNbinsZ();
    // the order of nested loop dictates which is the scaled axis
    for (int xBin = 1; xBin <= nBinsX; xBin++) {
        for (int zBin = 1; zBin <= nBinsZ; zBin++) {
            // Inner loop over y-axis bins: this is scaled axis (pT,D bins)
            for (int yBin = 1; yBin <= nBinsY; yBin++) {
                double binContent = hTreatedPowheg->GetBinContent(xBin, yBin, zBin);
                double binError = hTreatedPowheg->GetBinError(xBin, yBin, zBin);

                // Obtain efficiencies ratio for specific bin
                double effPrompt = hEffPrompt->GetBinContent(yBin);
                double effNonPrompt = hEffNonPrompt->GetBinContent(yBin);

                if (effPrompt == 0) {
                    std::cerr << "Error: null prompt efficiency value for bin (" << xBin << "," << yBin << "," << zBin << ").\n";
                }
                

                // Scale the y-axis content and error by non-prompt efficiency/prompt efficiency
                hTreatedPowheg->SetBinContent(xBin, yBin, zBin, binContent * effNonPrompt / effPrompt);
                hTreatedPowheg->SetBinError(xBin, yBin, zBin, binError * effNonPrompt / effPrompt);
            }
        }
        
        
        
    }
    // Save modified 3D histogram
    dataContainer.hPowheg.emplace_back(hTreatedPowheg);
    
    
    

    //
    // 4th step: fold Delta R vs pT,jet distribution using detector response matrix of non-prompt D0 jets
    //
    
    // a) method 1: AppluTruth()
    TH1D* hAllptDPowheg_flat = dynamic_cast<TH1D*>(dataContainer.hTruth.second->Clone("hAllptDPowheg_flat"));
    hAllptDPowheg_flat->Reset();
    flattenHistogram(dataContainer.hAllptDPowheg[0], hAllptDPowheg_flat); // input 2D, output 1D
    dataContainer.hFolded.second = (TH1D*)dataContainer.response.ApplyToTruth(hAllptDPowheg_flat);
    unflattenHistogram(dataContainer.hFolded.second, dataContainer.hFolded.first); // input 1D, output 2D
    dataContainer.hFolded.first->SetName("hFoldedNonPrompt_1");
    dataContainer.hFolded.first->SetTitle("Folded with ApplyToTruth() method;#frac{1}{L_{int}}#DeltaR^{b #rightarrow D^{0}}_{reco};p_{T,jet}^{ch}");
    dataContainer.hAllptDPowheg.emplace_back(dataContainer.hFolded.first); // [1] = folded method 1
    
    // b) method 2: manual matrix multiplication
    // since TH2D objects can't be used for 4D object, each bin needs to be treated directly
    TH2D* hFoldedNonPrompt_2 = manualFolding(dataContainer.hResponse, dataContainer.hAllptDPowheg[0]);
    hFoldedNonPrompt_2->SetName("hFoldedNonPrompt_2");
    hFoldedNonPrompt_2->SetTitle("Folded with manual matrix multiplication method;#frac{1}{L_{int}}#DeltaR^{b #rightarrow D^{0}}_{reco};p_{T,jet}^{ch}");
    dataContainer.hAllptDPowheg.emplace_back(hFoldedNonPrompt_2); // [2] = folded method 2

    
    std::cout << "Generator data smeared.\n";
}

// Module to subtract non-prompt D0 jets from prompt efficiency corrected distribution
void feedDown(FeedDownData& dataContainer, const double jetptMin, const double jetptMax, TFile* fEfficiency) {
    
    // Temporary histogram for data treatment and storage
    TH1D* h_back_FDsubtracted = (TH1D*)fEfficiency->Get("hBackSubCorrected_allpt");

    if (!h_back_FDsubtracted) {
        std::cerr << "Error: could not retrieve hBackSubCorrected_allpt histogram" << std::endl;
    }
    

    // Obtain non-prompt distribution within required pT,jet interval by projecting it along x axis
    int minBin = dataContainer.hAllptDPowheg[2]->GetYaxis()->FindBin(jetptMin);
    int maxBin = dataContainer.hAllptDPowheg[2]->GetYaxis()->FindBin(jetptMax);
    TH1D* h_feed_down = dataContainer.hAllptDPowheg[2]->ProjectionX("hBackSub_FDsubtracted", minBin, maxBin);
    h_back_FDsubtracted->Add(h_feed_down, -1.0);
    dataContainer.hSBFeedDownSubtracted = h_back_FDsubtracted;
    dataContainer.hSBFeedDownSubtracted->SetLineColor(46); // 46 = pastel red
    dataContainer.hSBFeedDownSubtracted->SetMarkerColor(46);
    dataContainer.hSBFeedDownSubtracted->SetMarkerStyle(kOpenSquare);
    dataContainer.hSBFeedDownSubtracted->SetStats(0);
    
    std::cout << "Feed-down subtracted from experimental data.\n";
}


void plotHistograms(const FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax) {
    cout << "Plotting histograms...\n";

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    //
    // 3D histograms
    //
    TCanvas* cPowheg = new TCanvas("cPowheg","POWHEG data");
    cPowheg->SetCanvasSize(1800,1000);
    cPowheg->Divide(2,2);
    cPowheg->cd(1);
    dataContainer.hPowheg[0]->Draw("colz");
    cPowheg->cd(2);
    dataContainer.hPowheg[1]->Draw("colz");

    //
    // Response matrix representation in 2D histogram
    //
    TCanvas* cResponse = new TCanvas("cResponse","Response matrices for all pT,D bins");
    cResponse->SetCanvasSize(1800,1000);
    cResponse->cd();
    const TH2* hResponse2D = dataContainer.response.Hresponse();
    TH2D* hResponse2DClone = static_cast<TH2D*>(hResponse2D->Clone("hResponse2DClone"));
    std::cout << "hResponse2DClone has " << hResponse2DClone->GetEntries() << " bins." << std::endl;
    hResponse2DClone->SetTitle("2D response matrix from method 1;#DeltaR_{reco}^{b #rightarrow D^{0}};#DeltaR_{truth}^{b #rightarrow D^{0}}");
    hResponse2DClone->Draw("colz");
    
    //
    // Detector and particle level matrices used to build the manual response matrix
    //
    TCanvas* cMatching = new TCanvas("cMatching","Matching detector and particle 2D histograms");
    cMatching->SetCanvasSize(1800,1000);
    cMatching->Divide(2,2);
    cMatching->cd(1);
    dataContainer.hMeasured.first->Draw("colz");
    cMatching->cd(2);
    dataContainer.hMeasured.second->Draw();
    cMatching->cd(3);
    dataContainer.hTruth.first->Draw("colz");
    cMatching->cd(4);
    dataContainer.hTruth.second->Draw();

    //
    // 2D True prompt and non-prompt Delta R distribution for all pT,D: not folded, folded method 1, folded method 2
    //
    //TLegend* legendNonPrompt = new TLegend(0.65,0.49,0.8,0.62);
    //legendNonPrompt->AddEntry(dataContainer.hAllptDPowheg[0],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    //legendNonPrompt->AddEntry(dataContainer.hAllptDPowheg[1],"Reconstructed - method 1", "lpe");
    //legendNonPrompt->AddEntry(dataContainer.hAllptDPowheg[2],"Reconstructed - method 2", "lpe");
    TCanvas* cNonPrompt = new TCanvas("cNonPrompt","Non-prompt Delta R plots");
    cNonPrompt->Divide(2,3);
    cNonPrompt->SetCanvasSize(1800,1000);
    cNonPrompt->cd(1);
    dataContainer.hAllptDPowheg[0]->SetStats(0);
    dataContainer.hAllptDPowheg[0]->Draw("colz");
    cNonPrompt->cd(3);
    dataContainer.hAllptDPowheg[1]->SetStats(0);
    dataContainer.hAllptDPowheg[1]->Draw("colz");
    cNonPrompt->cd(4);
    dataContainer.hAllptDPowheg[2]->SetStats(0);
    dataContainer.hAllptDPowheg[2]->Draw("colz");
    cNonPrompt->cd(5);
    dataContainer.hFolded.first->Draw();
    cNonPrompt->cd(6);
    dataContainer.hFolded.second->Draw();
    //legendNonPrompt->Draw();
    double statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.75, "Folded non-prompt D0 jets");
    

    //
    // Feed-down subtracted Delta R distributions
    //
    TCanvas* cSBFeedDown = new TCanvas("cSBFeedDown","Non-prompt subtracted Delta R plots");
    cSBFeedDown->SetCanvasSize(1800,1000);
    cSBFeedDown->cd();
    dataContainer.hSBFeedDownSubtracted->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.75, "Feed-down subtracted");


    cPowheg->Print(Form("pT_feeddown_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cMatching->Print(Form("pT_feeddown_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cNonPrompt->Print(Form("pT_feeddown_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cResponse->Print(Form("pT_feeddown_%.0f_to_%.0fGeV.pdf",jetptMin,jetptMax));
    cSBFeedDown->Print(Form("pT_feeddown_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

}

void saveData(const FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax){
    // Open output file
    TFile* outFile = new TFile(Form("backSubFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");

    // store each histogram in file
    dataContainer.hSBFeedDownSubtracted->Write();
    
    outFile->Close();
    delete outFile;
    
    cout << "Data stored in file" << Form("backSubFeedDown_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)) << endl;
}

void FeedDownSubtraction(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Luminosity (for now arbitrary)
    double luminosity = 10000;

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;
    // jet pT cuts
    double jetptMin = 5; // GeV
    double jetptMax = 30; // GeV
    std::vector<double> ptjetBinEdges = {5., 10., 15., 20., 25., 30.};
    // deltaR histogram
    int deltaRbins = 10000; // deltaRbins = numberOfPoints, default=10 bins for [0. 0.4]
    double minDeltaR = 0.;
    double maxDeltaR = 0.4;
    std::vector<double> deltaRBinEdges = {0.,0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.4}; // TODO: investigate structure before 0.005
    // mass histogram
    int massBins = 100; 
    double minMass = 1.67;
    double maxMass = 2.1;
    // pT,D histograms
    int ptBins = 100;
    double minPt = 0.;
    double maxPt = 30.;
    std::vector<double> ptDBinEdges = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};

    // opening files
    TFile* fPowheg = new TFile("../SimulatedData/POWHEG/trees_powheg_fd_central.root","read");
    TFile* fSimulatedO2 = new TFile("../SimulatedData/Hyperloop_output/McChargedMatched/AO2D_merged_All.root","read");
    TFile* fEfficiency = new TFile(Form("../2-Efficiency/backSubEfficiency_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fBackSub = new TFile(Form("../1-SignalTreatment/backSub_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fSigExt = new TFile(Form("../1-SignalTreatment/sigExt_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fSimulatedO2 || fSimulatedO2->IsZombie()) {
        std::cerr << "Error: Unable to open simulated data ROOT file." << std::endl;
    }
    if (!fBackSub || fBackSub->IsZombie()) {
        std::cerr << "Error: Unable to open background subtracted data ROOT file." << std::endl;
    }
    if (!fSigExt || fSigExt->IsZombie()) {
        std::cerr << "Error: Unable to open signal extracted ROOT file." << std::endl;
    }
    
    std::vector<const char*> names = {"histPt1", "histPt2", "histPt3", 
                                      "histPt4", "histPt5", "histPt6",
                                      "histPt7", "histPt8", "histPt9"};                                                     // Names of histograms
    std::vector<const char*> titles = {"3 < p_{T,D} < 4 GeV/c", "4 < p_{T,D} < 5 GeV/c", "5 < p_{T,D} < 6 GeV/c",
                                       "6 < p_{T,D} < 7 GeV/c", "7 < p_{T,D} < 8 GeV/c", "8 < p_{T,D} < 10 GeV/c",
                                       "10 < p_{T,D} < 12 GeV/c", "12 < p_{T,D} < 15 GeV/c", "15 < p_{T,D} < 30 GeV/c"};    // Titles of histograms
    FeedDownData dataContainer = createHistograms(deltaRBinEdges, ptDBinEdges, jetptMin, jetptMax);

    // Fill histograms with POWHEG simulation data
    fillHistograms(fPowheg, fSimulatedO2, dataContainer, jetptMin, jetptMax);

    // Create response matrices for all pT,D bins considered
    buildResponseMatrix(dataContainer);

    // Fold data using two methods
    smearGeneratorData(dataContainer, luminosity, fEfficiency);

    // Subtract non-prompt distribution from prompt efficiency corrected ones
    feedDown(dataContainer, jetptMin, jetptMax, fEfficiency);

    // Plot the efficiency histogram and further corrected histograms
    plotHistograms(dataContainer, jetptMin, jetptMax);

    // Save corrected distributions to file
    saveData(dataContainer, jetptMin, jetptMax);


    

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    FeedDownSubtraction();
    return 0;
}
