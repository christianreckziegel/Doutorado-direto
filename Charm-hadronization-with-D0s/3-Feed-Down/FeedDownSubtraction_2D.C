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
    // non-prompt D0s
    TH2D* hMeasured;                                        // 2D representation of response measured data
    TH2D* hMeasuredTotalRange;                              // 2D total measured range data
    TH2D* hDivMeasuredRange;                                // 2D division of (total range) / (inside range) for measured data
    TH2D* hTruth;                                           // 2D representation of response truth data
    TH2D* hTruthTotalRange;                                 // 2D total truth range data
    TH2D* hDivTruthRange;                                   // 2D division of (inside range) / (total range) for truth data

    TH2D* nimaFolded;                                       // 2D folded data with Nima's folding function
    RooUnfoldResponse response;                             // response matrix for non-prompt D0s only, for overall pT,D; using RooUnfoldResponse object, method 1
    std::vector<TH3D*> hPowheg;                             // deltaR vs pT,D vs pT,jet; [0] = generator level, [1] = treated but not yet detector level
    std::vector<TH2D*> hAllptDPowheg;                       // deltaR vs pT,jet; [0] = generator level, [1] = method 1 folded, [2] = method 2 folded
    std::vector<TH1D*> hEfficiencies;                       // inclusive = 0, prompt only = 1, non-prompt only = 2
    TH2D* hSBFeedDownSubtracted;                            // non-prompt subtracted Delta R distribution, all pT,D

    // Testing histograms
    TH2D* MCPoutRespInput;                                  // outside response matrix before folding on particle level data
    TH2D* JetPtOutRespInput;                                // outside response matrix before folding detector vs. particle level jet pT
    TH2D* DeltaROutRespInput;                               // outside response matrix before folding detector vs. particle level delta R

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
FeedDownData createHistograms(const std::vector<double>& xBinEdges_particle, const std::vector<double>& yBinEdges_particle, const std::vector<double>& zBinEdges_particle,
                              const std::vector<double>& xBinEdges_detector, const std::vector<double>& yBinEdges_detector, const std::vector<double>& zBinEdges_detector) {
                              //const double& jetptMin, const double& jetptMax) {
    // Create struct to store data
    FeedDownData dataContainer;

    int xNumBinEdges = xBinEdges_particle.size();
    int yNumBinEdges = yBinEdges_particle.size();
    int zNumBinEdges = zBinEdges_particle.size();

    // Create 2D histogram
    dataContainer.hPowheg.emplace_back(new TH3D("h_deltaR_vs_pt", "POWHEG + PYTHIA;#DeltaR;p_{T,D}^{gen};p_{T,jet}^{ch}", xNumBinEdges-1, xBinEdges_particle.data(), 
                                                                                                            yNumBinEdges-1, yBinEdges_particle.data(), 
                                                                                                            zNumBinEdges-1, zBinEdges_particle.data()));
    dataContainer.hPowheg[0]->SetMarkerColor(30);
    dataContainer.hPowheg[0]->SetLineColor(30); // 30 = pastel green
    dataContainer.hPowheg[0]->SetMarkerStyle(kCircle);
    dataContainer.hPowheg[0]->Sumw2();
    dataContainer.hPowheg[0]->SetStats(0);

    cout << "POWHEG histograms created.\n";

    //
    // Matching histograms for folding process
    //
    // 2D truth and particle level data
    dataContainer.hTruth = new TH2D("hTruth2D", "Truth;#DeltaR^{part};p_{T,jet}^{part}", xNumBinEdges-1, xBinEdges_particle.data(), zNumBinEdges-1, zBinEdges_particle.data());
    dataContainer.hTruthTotalRange = new TH2D("hTruth2D_totalRange", "Truth;#DeltaR^{part};p_{T,jet}^{part}", xNumBinEdges-1, xBinEdges_particle.data(), zNumBinEdges-1, zBinEdges_particle.data());
    
    //dataContainer.MCPoutRespInput = new TH2D("MCPoutRespInput", "Truth outside response, before folding;#DeltaR;p_{T,jet}", 100, 0., 5., 350, 0., 350.);
    dataContainer.MCPoutRespInput = new TH2D("MCPoutRespInput", "Truth outside response, before folding;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges_particle.data(), zNumBinEdges-1, zBinEdges_particle.data()); 

    // 2D measured and detector level data
    xNumBinEdges = xBinEdges_detector.size();
    yNumBinEdges = yBinEdges_detector.size();
    zNumBinEdges = zBinEdges_detector.size();
    dataContainer.hMeasured = new TH2D("hMeasured2D", "Measured;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges_detector.data(), zNumBinEdges-1, zBinEdges_detector.data());
    dataContainer.hMeasuredTotalRange = new TH2D("hMeasured2D_totalRange", "Measured;#DeltaR;p_{T,jet}", xNumBinEdges-1, xBinEdges_detector.data(), zNumBinEdges-1, zBinEdges_detector.data());

    // Particle level outside response range, inside particle level range
    dataContainer.JetPtOutRespInput = new TH2D("JetPtOutRespInput", "Truth outside response, before folding;p_{T,jet}^{det};p_{T,jet}^{part}", 100, 0., 100., 100, 0., 100.);
    dataContainer.DeltaROutRespInput = new TH2D("DeltaROutRespInput", "Truth outside response, before folding;#DeltaR^{det};#DeltaR^{part}", 100, -0.1, 3.0, 100, -0.1, 3.0);

    cout << "Matching histograms created.\n";

    return dataContainer;
}

// Module to fill histograms from 
// POWHEG+PYTHIA TFile data (all data in file is non-prompt and on particle level)
// and
// O2 matching task simulation
/**
 * @brief Fill histograms.
 * 
 * This helper function fill two kinds of histograms:
 * - POWHEG+PYTHIA simulated data on particle level with non-prompt D0 
 * - O2 simulated matched data between particle and detector level
 *
 * @param fPowheg The ROOT file with POWHEG+PYTHIA non-prompt D0 simulated data.
 * @param fSimulatedO2 The ROOT file with matched O2 simulated data.
 * @param dataContainer Container with already instanciated histograms.
 * @param jetptMin Min jet pT cut.
 * @param jetptMax Max jet pT cut.
 *
 * @note uses 0-based indexing for calculation but 1-based indexing for ROOT histograms.
 *
 * @see createHistograms() [Instanciate histograms.]
 */
void fillHistograms(TFile* fPowheg, TFile* fSimulatedO2, FeedDownData& dataContainer, double jetptMin, double jetptMax) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = 0.4; // on particle level delta R
    const double MCDDeltaRcut = 0.4; // on detector level delta R
    const double MCPHfPtMincut = 3.; // on particle level
    const double MCDHfPtMincut = 3.; // on detector level
    const double MCPHfPtMaxcut = 30.; // on particle level
    const double MCDHfPtMaxcut = 30.; // on detector level

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
        double deltaR = sqrt(pow(eta_jet-eta_cand,2) + pow(DeltaPhi(phi_jet,phi_cand),2));
        // Fill 2D histogram considering jet pT and detector acceptance
        if ((abs(eta_jet) < MCPetaCut) && (abs(y_cand) < MCPyCut) && ((pt_jet >= jetptMin) && (pt_jet < jetptMax)) && ((deltaR >= 0.) && (deltaR < MCPDeltaRcut)) && ((pt_cand >= MCPHfPtMincut) && (pt_cand < MCPHfPtMaxcut))) {
            
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

        // Fill histograms considering jet pT and detector acceptance for NON-PROMPT particles, inside response range (truth and measured levels)
        if ((abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= 0.) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)) && !MCPhfprompt
            && (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= 0.) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut)) && !MCDhfprompt) {

            // Filling measured 2D histogram
            dataContainer.hMeasured->Fill(MCDDeltaR, MCDjetPt);

            // Filling truth 2D histogram
            dataContainer.hTruth->Fill(MCPDeltaR, MCPjetPt);

        
        } else {
            // all entries outside response ranges, but inside particle level range (non-prompt particles)
            if ((abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= 0.) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)) && !MCPhfprompt
                //&& (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && !MCDhfprompt) {
                ) {
                dataContainer.MCPoutRespInput->Fill(MCPDeltaR, MCPjetPt);
                dataContainer.JetPtOutRespInput->Fill(MCDjetPt, MCPjetPt);
                dataContainer.DeltaROutRespInput->Fill(MCDDeltaR, MCPDeltaR);

                

            }
        }

        // Fill reference input total truth range only histograms for non-prompt particles
        if ((abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= 0.) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)) && !MCPhfprompt) {
            dataContainer.hTruthTotalRange->Fill(MCPDeltaR, MCPjetPt);
        }
        
        // Fill reference output total measured range histograms for non-prompt particles
        if ((abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= 0.) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut)) && !MCDhfprompt) {
            dataContainer.hMeasuredTotalRange->Fill(MCDDeltaR, MCDjetPt);
        }

    }

    cout << "Response matched histograms filled.\n";

}

// Module to build 2D response matrix out of flattened 1D input data
void buildResponseMatrix(FeedDownData& dataContainer, TFile* fSimulatedO2, TFile* fEfficiency, std::vector<double>& ptjetBinEdges_particle, std::vector<double>& ptjetBinEdges_detector) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = 0.4; // on particle level delta R
    const double MCDDeltaRcut = 0.4; // on detector level delta R
    const double MCPHfPtMincut = 3.; // on particle level D0
    const double MCDHfPtMincut = 3.; // on detector level D0
    const double MCPHfPtMaxcut = 30.; // on particle level D0
    const double MCDHfPtMaxcut = 30.; // on detector level D0
    const double jetptMin = ptjetBinEdges_particle[0];
    const double jetptMax = ptjetBinEdges_particle[ptjetBinEdges_particle.size() - 1];

    // method 1: create 2D response matrix of non-prompt flattened Delta R and pT,jet, for overall pT,D
    //(DeltaR_detector, pTjet_detector, DeltaR_particle, pTjet_particle) -> flattened to (detector, particle)
    dataContainer.response = RooUnfoldResponse(dataContainer.hMeasured, dataContainer.hTruth);

    std::cout << "Response matrix created.\n";

    // Access of correspondent prompt efficiency histogram
    TH1D* hEffPrompt;

    //__________________________________________-
    //
    // O2 matching task measured and truth histograms
    //
    // Accessing TTree
    TTree* tree = (TTree*)fSimulatedO2->Get("O2matchtable");
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

    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        // calculating delta R
        double MCPDeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2));
        double MCDDeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));

        // Fill histograms considering jet pT and detector acceptance
        if ((abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && ((MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax)) && ((MCPDeltaR >= 0.) && (MCPDeltaR < MCPDeltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)) && !MCPhfprompt
            && (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax)) && ((MCDDeltaR >= 0.) && (MCDDeltaR < MCDDeltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut)) && !MCDhfprompt) {
            
            // Access the efficiency histogram correspondent to the correct pT,jet range
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < ptjetBinEdges_particle.size() - 1 && !filled; iEdge++) {
                if ((MCPjetPt >= ptjetBinEdges_particle[iEdge]) && (MCPjetPt < ptjetBinEdges_particle[iEdge + 1])) {
                    hEffPrompt = (TH1D*)fEfficiency->Get(Form("JetPtRange_%.0f_%0.f/efficiency_prompt",ptjetBinEdges_particle[iEdge],ptjetBinEdges_particle[iEdge + 1]));
                    if (!hEffPrompt || hEffPrompt->IsZombie()) {
                        std::cerr << Form("Error: Unable to open efficiency histogram JetPtRange_%.0f_%0.f/efficiency_prompt in ROOT data file.",jetptMin,jetptMax) << std::endl;
                    }
                    filled = true; // Exit the loop once the correct histogram is found

                    // Find the bin corresponding to the given pT,D value
                    int bin = hEffPrompt->FindBin(MCDhfPt);
                    // Get the efficiency value from the bin content
                    double efficiency_prompt = hEffPrompt->GetBinContent(bin);

                    // Fill 4D RooUnfoldResponse object
                    //dataContainer.response.Fill(MCDDeltaR, MCDjetPt, MCPDeltaR, MCPjetPt);
                    dataContainer.response.Fill(MCDDeltaR, MCDjetPt, MCPDeltaR, MCPjetPt, 1./efficiency_prompt); // jet pT shape is influenced by D0 pT efficiency
                }
                
            }
            
            
        }
        
    }
    
    std::cout << "Response matrix filled.\n";
}

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
                    foldedError2 = std::pow(hTruth->GetBinError(iTruth + 1,jTruth + 1),2) * std::pow(response(index_x_measured, index_x_truth),2);
                }
            }
            hFolded->SetBinContent(iMeasured + 1, jMeasured + 1, foldedValue);
            hFolded->SetBinError(iMeasured + 1, jMeasured + 1, std::sqrt(foldedError2));
        }
    }
    
    std::cout << "Folded manually with bin index flattening." << std::endl;

    return hFolded;
}

/**
 * @brief Remove entries in region not treated by response matrix from the folding input data.
 *
 * This helper function
 *
 * @param dataContainer Container with total and intersection range data and the POWHEG data to be corrected.
 *
 * @return Corrected POWHEG truth input 2D distribution
 *
 * @note 
 *
 * @see smearGeneratorData() [Correct truth POWHEG data before folding.]
 */
TH2D* removeOutsideData(FeedDownData& dataContainer) {

    // Copy original POWHEG distribution structure
    TH2D* hPowhegCorrected = (TH2D*)dataContainer.hAllptDPowheg[0]->Clone("hPowhegRangeCorrected");

    // Start with histogram of intersection area only (inside response ranges)
    dataContainer.hDivTruthRange = dynamic_cast<TH2D*>(dataContainer.hTruth->Clone("hInsideOverTotalDivision"));

    // Get kinematic efficiency 2D histogram diving intersection over total range
    dataContainer.hDivTruthRange->Divide(dataContainer.hTruthTotalRange);
    dataContainer.hDivTruthRange->SetTitle("Inside response range / Total truth range (Kinematic efficiency?)");

    // Correct POWHEG data multiplying each bin by corresponding kinematic efficiency
    hPowhegCorrected->Multiply(dataContainer.hDivTruthRange);

    return hPowhegCorrected;

}

/**
 * @brief Add entries of the folding output data from region not treated by response matrix .
 *
 * This helper function
 *
 * @param dataContainer Container with total and intersection range data and the folded data to be corrected.
 *
 * @return Corrected Folded detector level output 2D distribution
 *
 * @note 
 *
 * @see smearGeneratorData() [Correct detector level folded data.]
 */
TH2D* addOutsideData(FeedDownData& dataContainer) {

    // Copy original immediately folded distribution structure
    TH2D* hFoldedCorrected = (TH2D*)dataContainer.hAllptDPowheg[0]->Clone("hPowhegRangeCorrected");

    // Start with histogram of total range area
    dataContainer.hDivMeasuredRange = dynamic_cast<TH2D*>(dataContainer.hMeasured->Clone("hInsideOverTotalDivision"));

    // Get kinematic efficiency 2D histogram diving total range over intersection
    dataContainer.hDivMeasuredRange->Divide(dataContainer.hMeasuredTotalRange);
    dataContainer.hDivMeasuredRange->SetTitle("Inside response range / Total measured range (kinematic efficiency?)");

    // Correct POWHEG data multiplying each bin by corresponding kinematic efficiency
    hFoldedCorrected->Divide(dataContainer.hDivMeasuredRange);

    return hFoldedCorrected;

}

// Get POWHEG and data luminosities
void getLuminosities(TFile* fPowheg, double luminosity_powheg, double luminosity) {
    // Accessing total cross section value stored in first bin (in mb)
    TH1D* xSection_powheg = dynamic_cast<TH1D*>(fPowheg->Get("fHistXsection"));
    double crossSecPowheg = xSection_powheg->GetBinContent(1);

    // Accessing number of events
    TTree* tree_D0 = dynamic_cast<TTree*>(fPowheg->Get("tree_D0"));
    double numOfEventsPowheg = tree_D0->GetEntries();

    // integrated POWHEG luminosity 
    luminosity_powheg = numOfEventsPowheg/crossSecPowheg;
    std::cout << "POWHEG luminosity = " << luminosity_powheg << " mb^-1" << std::endl;

    // Calculating measured luminosity
    double numOfEventsMeasured = 3000000;
    double crossSecMeasured = 57.8; // in mb units
    luminosity = numOfEventsMeasured/crossSecMeasured;
    std::cout << "Measured luminosity = " << luminosity << " mb^-1" << std::endl;

}

// Module for folding particle level data from POWHEG simulation
void smearGeneratorData(FeedDownData& dataContainer, double& luminosity_powheg, TFile* fEfficiency, double& luminosity, double& BR) {
    //
    // 0th step: clone 3D histogram in order to save the smeared at the end (while keeping the original)
    //
    TH3D* hTreatedPowheg = static_cast<TH3D*>(dataContainer.hPowheg[0]->Clone("h_subtraction_Bs"));
    hTreatedPowheg->SetTitle("POWHEG + PYTHIA efficiency scaled;#frac{1}{L_{int}}#DeltaR^{b #rightarrow D^{0}};#frac{#epsilon_{non-prompt}}{prompt}#frac{1}{L_{int}}p_{T,D}^{b #rightarrow D^{0}};p_{T,jet}^{ch}");

    //
    // 1st step: obtain the 2D Delta R vs pT,jet projection for all pT,D bins before folding
    //
    // using a smart pointer provided by the C++ Standard Library
    TH2D* hAllptDPow = static_cast<TH2D*>(hTreatedPowheg->Project3D("zx"));
    hAllptDPow->SetName("h_true_nonprompt");
    //hAllptDPow->SetTitle(";#frac{1}{L_{int}}#DeltaR^{b #rightarrow D^{0}}_{truth};p_{T,jet}^{ch}");
    hAllptDPow->SetTitle("Before folding;#DeltaR^{b #rightarrow D^{0}}_{truth};p_{T,jet}^{ch}");
    // transfer ownership of the histogram from the unique_ptr to dataContainer object
    dataContainer.hAllptDPowheg.emplace_back(hAllptDPow); // [0] = not folded

    //
    // 2nd step: scale by 1 over POWHEG integrated luminosity
    // (skip this step for now)
    dataContainer.hAllptDPowheg[0]->Scale(1/luminosity_powheg); //luminosity/luminosity_powheg

    //
    // 3rd step: scale by efficiency ratio
    //
    TH1D* hEffPrompt;
    TH1D* hEffNonPrompt;
    hEffPrompt = (TH1D*)fEfficiency->Get("hptD0_vs_ptJet_prompt");
    hEffNonPrompt = (TH1D*)fEfficiency->Get("hptD0_vs_ptJet_nonprompt");
    if (!hEffPrompt || !hEffNonPrompt) {
        std::cerr << "Error: could not retrieve efficiency histograms.\n";
    }
    
    // manually scaling over each bin of pT,D dimension of 2D histogram (y-axis)
    int nBinsX = hTreatedPowheg->GetNbinsX(); // Delta R
    int nBinsY = hTreatedPowheg->GetNbinsY(); // pT,D
    int nBinsZ = hTreatedPowheg->GetNbinsZ(); // pT,jet
    // the order of nested loop dictates which is the scaled axis
    for (int xBin = 1; xBin <= nBinsX; xBin++) {
        for (int zBin = 1; zBin <= nBinsZ; zBin++) {
            // Inner loop over y-axis bins (pT,D): this is scaled axis (pT,D bins)
            for (int yBin = 1; yBin <= nBinsY; yBin++) {
                double binContent = hTreatedPowheg->GetBinContent(xBin, yBin, zBin);
                double binError = hTreatedPowheg->GetBinError(xBin, yBin, zBin);

                // Obtain efficiencies ratio for specific bin
                double effPrompt = hEffPrompt->GetBinContent(yBin, zBin);
                double effNonPrompt = hEffNonPrompt->GetBinContent(yBin, zBin);

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
    // 4th step: remove outside of response range data in POWHEG
    //
    TH2D* hPowhegOutRange = removeOutsideData(dataContainer);
    hPowhegOutRange->SetTitle("Before folding, with outside range correction (kinematic efficiency?)");
    dataContainer.hAllptDPowheg.emplace_back(hPowhegOutRange); // [1] = data outside response range removed

    //
    // 5th step: manualy fold Delta R vs pT,jet distribution using detector response matrix of non-prompt D0 jets
    //
    dataContainer.nimaFolded = nimaFolding(dataContainer.response, dataContainer.hAllptDPowheg[1], dataContainer.hMeasured);
    dataContainer.nimaFolded->SetTitle("Folded with Nima's 1D index function;#frac{1}{L_{int}}#DeltaR^{b #rightarrow D^{0}}_{reco};p_{T,jet}^{ch}");

    //
    // 6th step: add outside of response range data in folded data
    //
    TH2D* hFoldedOutRange = addOutsideData(dataContainer);
    dataContainer.hAllptDPowheg.emplace_back(hFoldedOutRange); // [2] = entries outside response range added to folded data
    dataContainer.hAllptDPowheg[2]->SetTitle("Folded data with outside range correction");

    //
    // 7th step::scale by measured integrated luminosity and BR of D0 decay channel
    //
    dataContainer.hAllptDPowheg[2]->Scale(BR*luminosity);

    std::cout << "Generator data smeared.\n";
}

// Module to subtract non-prompt D0 jets from prompt efficiency corrected distribution
void feedDown(FeedDownData& dataContainer, const double jetptMin, const double jetptMax, TFile* fEfficiency) {
    
    // Accessing 2D background subtracted data corrected by prompt efficiency
    TH1D* h_total_signal = (TH1D*)fEfficiency->Get("hDeltaR_vs_ptJet");
    if (!h_total_signal) {
        std::cerr << "Error: could not retrieve hDeltaR_vs_ptJet histogram" << std::endl;
    }
    dataContainer.hSBFeedDownSubtracted = (TH2D*)h_total_signal->Clone("h_DeltaR_vs_ptJet_fd_subtracted");

    // Perform subtraction of B contribution
    dataContainer.hSBFeedDownSubtracted->Add(dataContainer.hAllptDPowheg[2],-1);
    dataContainer.hSBFeedDownSubtracted->SetTitle("Feed-down subtracted");
    dataContainer.hSBFeedDownSubtracted->SetLineColor(46); // 46 = pastel red
    dataContainer.hSBFeedDownSubtracted->SetMarkerColor(46);
    dataContainer.hSBFeedDownSubtracted->SetMarkerStyle(kOpenSquare);
    dataContainer.hSBFeedDownSubtracted->SetStats(0);
    
    std::cout << "Feed-down estimation subtracted from experimental data.\n";
}


void plotHistograms(const FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax) {
    cout << "Plotting histograms...\n";

    // Changing color palette to a bigger one
    gStyle->SetPalette(kRainBow);

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
    hResponse2DClone->SetTitle("2D response matrix from 4D RooUnfoldResponse;2D Reconstructed;2D Truth");
    hResponse2DClone->Draw("colz");
    
    //
    // Detector and particle level matrices used to build the manual response matrix
    //
    TCanvas* cMatching = new TCanvas("cMatching","Matching detector and particle 2D histograms");
    cMatching->SetCanvasSize(1800,1000);
    cMatching->Divide(2,2);
    cMatching->cd(1);
    dataContainer.hMeasured->Draw("colz");
    cMatching->cd(2);
    dataContainer.hTruth->Draw("colz");
    cMatching->cd(3);
    TH1D* hMatchProjX1 = dataContainer.hMeasured->ProjectionX();
    hMatchProjX1->SetMarkerStyle(kCircle);
    hMatchProjX1->SetMarkerColor(38); // 38 = pastel blue
    hMatchProjX1->SetLineColor(38);
    hMatchProjX1->SetStats(0);
    hMatchProjX1->Sumw2();
    hMatchProjX1->Draw();
    cMatching->cd(4);
    TH1D* hMatchProjX2 = dataContainer.hTruth->ProjectionX();
    hMatchProjX2->SetMarkerStyle(kCircle);
    hMatchProjX2->SetMarkerColor(38); // 38 = pastel blue
    hMatchProjX2->SetLineColor(38);
    hMatchProjX2->SetStats(0);
    hMatchProjX2->Sumw2();
    hMatchProjX2->Draw();

    //
    // 2D True prompt and non-prompt Delta R distribution for all pT,D: not folded, folded method 1, folded method 2
    //
    //TLegend* legendNonPrompt = new TLegend(0.65,0.49,0.8,0.62);
    //legendNonPrompt->AddEntry(dataContainer.hAllptDPowheg[0],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    //legendNonPrompt->AddEntry(dataContainer.hAllptDPowheg[1],"Reconstructed - method 1", "lpe");
    //legendNonPrompt->AddEntry(dataContainer.hAllptDPowheg[2],"Reconstructed - method 2", "lpe");
    TCanvas* cNonPrompt = new TCanvas("cNonPrompt","Non-prompt Delta R plots");
    cNonPrompt->Divide(2,4);
    cNonPrompt->SetCanvasSize(1800,1000);
    cNonPrompt->cd(1);
    TH1D* hProjectionX0 = dataContainer.hAllptDPowheg[0]->ProjectionX();
    hProjectionX0->SetMarkerStyle(kCircle);
    hProjectionX0->SetMarkerColor(29); // 30 = pastel green
    hProjectionX0->SetLineColor(30);
    hProjectionX0->SetStats(0);
    hProjectionX0->Draw();
    cNonPrompt->cd(2);
    dataContainer.hAllptDPowheg[0]->SetStats(0);
    dataContainer.hAllptDPowheg[0]->Draw("colz");
    cNonPrompt->cd(3);
    TH1D* hProjectionX1 = dataContainer.hAllptDPowheg[1]->ProjectionX();
    hProjectionX1->SetMarkerStyle(kCircle);
    hProjectionX1->SetMarkerColor(30); // 30 = pastel green
    hProjectionX1->SetLineColor(30);
    hProjectionX1->SetStats(0);
    hProjectionX1->Draw();
    cNonPrompt->cd(4);
    dataContainer.hAllptDPowheg[1]->SetStats(0);
    dataContainer.hAllptDPowheg[1]->Draw("colz");
    cNonPrompt->cd(5);
    TH1D* hProjectionX2 = dataContainer.hAllptDPowheg[2]->ProjectionX();
    hProjectionX2->SetMarkerStyle(kCircle);
    hProjectionX2->SetMarkerColor(31); // 30 = pastel green
    hProjectionX2->SetLineColor(31);
    hProjectionX2->SetStats(0);
    hProjectionX2->Draw();
    cNonPrompt->cd(6);
    dataContainer.hAllptDPowheg[2]->SetStats(0);
    dataContainer.hAllptDPowheg[2]->Draw("colz");
    cNonPrompt->cd(7);
    TH1D* hProjectionX3 = dataContainer.nimaFolded->ProjectionX();
    hProjectionX3->SetMarkerStyle(kCircle);
    hProjectionX3->SetMarkerColor(32); // 30 = pastel green
    hProjectionX3->SetLineColor(32);
    hProjectionX3->SetStats(0);
    hProjectionX3->Draw();
    cNonPrompt->cd(8);
    dataContainer.nimaFolded->Draw("colz");
    //legendNonPrompt->Draw();
    double statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.75, "Folded non-prompt D0 jets");

    //
    // Truth data range plots
    //
    TCanvas* cInputRangeCorrection = new TCanvas("cInputRangeCorrection","Outside response range data removal");
    cInputRangeCorrection->Divide(3,2);
    cInputRangeCorrection->cd(1);
    dataContainer.hTruth->SetStats(0);
    dataContainer.hTruth->SetTitle("Inside response range truth data");
    dataContainer.hTruth->Draw("text");
    cInputRangeCorrection->cd(2);
    dataContainer.hTruthTotalRange->SetStats(0);
    dataContainer.hTruthTotalRange->SetTitle("Total truth range data");
    dataContainer.hTruthTotalRange->Draw("text");
    cInputRangeCorrection->cd(4);
    dataContainer.hAllptDPowheg[0]->SetStats(0);
    dataContainer.hAllptDPowheg[0]->Draw("colz");
    cInputRangeCorrection->cd(5);
    dataContainer.hAllptDPowheg[1]->SetStats(0);
    dataContainer.hAllptDPowheg[1]->Draw("colz");
    cInputRangeCorrection->cd(3);
    dataContainer.hDivTruthRange->SetStats(0);
    dataContainer.hDivTruthRange->Draw("text");

    //
    // Detector level data range plots
    //
    TCanvas* cOutputRangeCorrection = new TCanvas("cOutputRangeCorrection","Outside response range data addition");
    cOutputRangeCorrection->Divide(3,2);
    cOutputRangeCorrection->cd(1);
    dataContainer.hMeasured->SetStats(0);
    dataContainer.hMeasured->SetTitle("Inside response range measured data");
    dataContainer.hMeasured->Draw("text");
    cOutputRangeCorrection->cd(2);
    dataContainer.hMeasuredTotalRange->SetStats(0);
    dataContainer.hMeasuredTotalRange->SetTitle("Total measured range data");
    dataContainer.hMeasuredTotalRange->Draw("text");
    cOutputRangeCorrection->cd(4);
    dataContainer.nimaFolded->SetStats(0);
    dataContainer.nimaFolded->Draw("colz");
    cOutputRangeCorrection->cd(5);
    dataContainer.hAllptDPowheg[2]->SetStats(0);
    dataContainer.hAllptDPowheg[2]->Draw("colz");
    cOutputRangeCorrection->cd(3);
    dataContainer.hDivMeasuredRange->SetStats(0);
    dataContainer.hDivMeasuredRange->Draw("text");

    //
    // Detector vs. particle match
    //
    TCanvas* cDetPartMatch = new TCanvas("cDetPartMatch","Detector vs. particle level match");
    cDetPartMatch->Divide(2,2);
    cDetPartMatch->cd(1);
    dataContainer.JetPtOutRespInput->SetStats(0);
    dataContainer.JetPtOutRespInput->Draw("colz");
    cDetPartMatch->cd(2);
    dataContainer.DeltaROutRespInput->SetStats(0);
    dataContainer.DeltaROutRespInput->Draw("colz");

    //
    // For image printing
    //
    TCanvas* cFoldedData = new TCanvas("cFoldedData","Before and after folding data");
    cFoldedData->SetCanvasSize(1800,1000);
    cFoldedData->Divide(2,2);
    cFoldedData->cd(1);
    dataContainer.hAllptDPowheg[0]->SetStats(0);
    dataContainer.hAllptDPowheg[0]->Draw("colz");
    cFoldedData->cd(2);
    dataContainer.hAllptDPowheg[2]->SetStats(0);
    dataContainer.hAllptDPowheg[2]->Draw("colz");
    cFoldedData->cd(3);
    TH1D* hProjectionX4 = dataContainer.hAllptDPowheg[0]->ProjectionX();
    hProjectionX4->SetMarkerStyle(kCircle);
    hProjectionX4->SetMarkerColor(30); // 30 = pastel green
    hProjectionX4->SetLineColor(30);
    hProjectionX4->SetStats(0);
    hProjectionX4->Draw();
    cFoldedData->cd(4);
    hProjectionX2->Draw();

    


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
    
    //
    // Storing images
    //
    TString imagePath = "../Images/3-Feed-Down/";
    cResponse->Update();
    cResponse->SaveAs(imagePath + "FD_response_matrix_2D.png");
    cMatching->Update();
    cMatching->SaveAs(imagePath + "FD_response_match_histograms_2D.png");
    cSBFeedDown->Update();
    cSBFeedDown->SaveAs(imagePath + "FD_subtracted_2D.png");
    cFoldedData->Update();
    cFoldedData->SaveAs(imagePath + "FD_folded_data_2D.png");

    //
    // Storing in a single pdf file
    //
    cPowheg->Print(Form("pT_feeddown_%.0f_to_%.0fGeV_2D.pdf(",jetptMin,jetptMax));
    cMatching->Print(Form("pT_feeddown_%.0f_to_%.0fGeV_2D.pdf",jetptMin,jetptMax));
    cNonPrompt->Print(Form("pT_feeddown_%.0f_to_%.0fGeV_2D.pdf",jetptMin,jetptMax));
    cResponse->Print(Form("pT_feeddown_%.0f_to_%.0fGeV_2D.pdf",jetptMin,jetptMax));
    cSBFeedDown->Print(Form("pT_feeddown_%.0f_to_%.0fGeV_2D.pdf)",jetptMin,jetptMax));

}

void saveData(const FeedDownData& dataContainer, const double& jetptMin, const double& jetptMax){
    // Open output file
    TFile* outFile = new TFile(Form("backSubFeedDown_%d_to_%d_jetpt_2D.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");

    // store each histogram in file
    dataContainer.hSBFeedDownSubtracted->Write();
    
    outFile->Close();
    delete outFile;
    
    cout << "Data stored in file" << Form("backSubFeedDown_%d_to_%d_jetpt_2D.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)) << endl;
}

void FeedDownSubtraction_2D(){

    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Luminosity (for now arbitrary)
    double luminosity = 0;
    double luminosity_powheg = 0;
    double BR = 0.0393; // D0 -> KPi decay channel branching ratio = (3.93 +- 0.04) %

    // jet pT cuts
    std::vector<double> ptjetBinEdges_particle = {5., 7., 15., 30.};
    std::vector<double> ptjetBinEdges_detector = {5., 7., 15., 30.};
    double jetptMin = ptjetBinEdges_particle[0]; // GeV
    double jetptMax = ptjetBinEdges_particle[ptjetBinEdges_particle.size() - 1]; // GeV

    // deltaR histogram
    std::vector<double> deltaRBinEdges_particle = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4}; // TODO: investigate structure before 0.005: 0.,0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.4
    std::vector<double> deltaRBinEdges_detector = {0.,0.05, 0.1, 0.15, 0.2, 0.3, 0.4};
    double minDeltaR = deltaRBinEdges_particle[0];
    double maxDeltaR = deltaRBinEdges_particle[deltaRBinEdges_particle.size() - 1];
    
    // pT,D histograms
    std::vector<double> ptDBinEdges_particle = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};
    std::vector<double> ptDBinEdges_detector = {3., 4., 5., 6., 7., 8., 10., 12., 15., 30.};
    double minPtD = ptDBinEdges_particle[0];
    double maxPtD = ptDBinEdges_particle[ptDBinEdges_particle.size() - 1];

    // Opening files
    TFile* fPowheg = new TFile("../SimulatedData/POWHEG/trees_powheg_fd_central.root","read");
    TFile* fSimulatedO2 = new TFile("../SimulatedData/Hyperloop_output/McChargedMatched/AO2D_merged_All.root","read");
    TFile* fEfficiency = new TFile(Form("../2-Efficiency/backSubEfficiency_%d_to_%d_jetpt_2D.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fBackSub = new TFile(Form("../1-SignalTreatment/SideBand/backSub_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    TFile* fSigExt = new TFile(Form("../1-SignalTreatment/SignalExtraction/sigExt_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fSimulatedO2 || fSimulatedO2->IsZombie()) {
        std::cerr << "Error: Unable to open simulated data ROOT file." << std::endl;
    }
    if (!fBackSub || fBackSub->IsZombie()) {
        std::cerr << "Error: Unable to open background subtracted data ROOT file." << std::endl;
    }
    if (!fSigExt || fSigExt->IsZombie()) {
        std::cerr << "Error: Unable to open signal extracted ROOT file." << std::endl;
    }
    
    //FeedDownData dataContainer = createHistograms(deltaRBinEdges, ptDBinEdges, jetptMin, jetptMax);
    FeedDownData dataContainer = createHistograms(deltaRBinEdges_particle, ptDBinEdges_particle, ptjetBinEdges_particle,
                                                  deltaRBinEdges_detector, ptDBinEdges_detector, ptjetBinEdges_detector);

    // Fill histograms with POWHEG simulation data
    fillHistograms(fPowheg, fSimulatedO2, dataContainer, jetptMin, jetptMax);

    // Create response matrices for all pT,D bins considered
    buildResponseMatrix(dataContainer, fSimulatedO2, fEfficiency, ptjetBinEdges_particle, ptjetBinEdges_detector);

    // Get POWHEG and data luminosities
    getLuminosities(fPowheg, luminosity_powheg, luminosity);

    // Fold data using two methods
    smearGeneratorData(dataContainer, luminosity_powheg, fEfficiency, luminosity, BR);

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
    FeedDownSubtraction_2D();
    return 0;
}
