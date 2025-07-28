/*
 * Macro for performing side-band subtraction procedure to second closure test
 * 
 * 
 * 
 * Author: Christian Reckziegel
**/


using namespace std;

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

// Get the optimal BDT score cut for the corresponding pT,D of the entry
double GetBkgProbabilityCut(double pT, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    for (size_t i = 0; i < bdtPtCuts.size() - 1; ++i) {
        if (pT >= bdtPtCuts[i].first && pT < bdtPtCuts[i + 1].first) {
            return bdtPtCuts[i].second;
        }
    }
    return 1.0; // Default: accept all if out of range
}

// Module to perform fits to histograms
struct FitContainer {
    std::vector<TF1*> fitBackgroundOnly;      // background only fits
    std::vector<TF1*> fitSignalOnly;      // background only fits
    std::vector<TF1*> fitReflectionsOnly;      // background only fits
    std::vector<TF1*> fitTotal;      // background only fits
};
FitContainer calculateFitTemplates() {
    //
};
// Save images to pdf file
void storeImages(std::vector<TH2D*>& hInvariantMass2D, std::vector<TH1D*>& hInvariantMass1D, const double& jetptMin, const double& jetptMax) {
    
    int nHistos = hInvariantMass2D.size();
    // Start with a square layout (or close to it)
    int nCols = static_cast<int>(std::ceil(std::sqrt(nHistos)));
    int nRows = static_cast<int>(std::ceil(nHistos / static_cast<double>(nCols)));
    
    // Creating canvases
    TCanvas* cInvariantMass2D = new TCanvas("cInvariantMass2D","DeltaR as a function of invariant mass plots");
    cInvariantMass2D->SetCanvasSize(1800,1000);
    cInvariantMass2D->Divide(nCols,nRows); // columns, lines
    TCanvas* cInvariantMass1D = new TCanvas("cInvariantMass1D","Invariant mass plots");
    cInvariantMass1D->SetCanvasSize(1800,1000);
    cInvariantMass1D->Divide(nCols,nRows); // columns, lines

    // Loop through all histograms (and fitting functions in the future)
    for(size_t iHisto = 0; iHisto < hInvariantMass2D.size() - 1; ++iHisto) {
        
        // Skip not filled histograms
        if (hInvariantMass2D[iHisto]->GetEntries() == 0) {
            continue;
        }
        
        // Drawing 2D histograms
        cInvariantMass2D->cd(iHisto+1);
        hInvariantMass2D[iHisto]->Draw("colz");

        // Drawing 1D histograms
        cInvariantMass1D->cd(iHisto+1);
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        hInvariantMass1D[iHisto]->SetMarkerStyle(kDot); //kFullDotMedium
        hInvariantMass1D[iHisto]->SetMarkerColor(kBlack);
        hInvariantMass1D[iHisto]->SetLineColor(kBlack);
        hInvariantMass1D[iHisto]->GetYaxis()->SetTitle("counts");
        //hInvariantMass1D->SetMinimum(0);
        hInvariantMass1D[iHisto]->Draw();
    }

    //
    // Storing images in a single pdf file
    //
    TString imagePath = "../Images/5-ClosureTest/Second/";
    cInvariantMass2D->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cInvariantMass1D->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

    delete cInvariantMass2D;
    delete cInvariantMass1D;

}

TH2D* AnalyzeJetPtRange(TFile* fClosureInput, const double& jetptMin, const double& jetptMax, 
                        const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& ptDBinEdges_detector,
                        const std::vector<std::pair<double, double>>& bdtPtCuts) {

    // mass histogram
    int massBins = 50; // default=100 
    double minMass = 1.72; // use from 1.72, used to use 1.67
    double maxMass = 2.06; // old = 2.1, new = 2.05 + 0.01

    // ----- Create and fill invariant mass distributions with detector level data
    std::vector<TH2D*> hInvariantMass2D;
    for (size_t i = 0; i < ptDBinEdges_detector.size() - 1; ++i) {
        // So that the title adapts to fractional binning title
        if (std::fmod(ptDBinEdges_detector[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(ptDBinEdges_detector[i+1], 1.0) != 0) {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), massBins, minMass, maxMass, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
        
            } else {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.1f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), massBins, minMass, maxMass, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
        
            }
        } else {
            if (std::fmod(ptDBinEdges_detector[i+1], 1.0) != 0) {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < #it{p}_{T, D^{0}} < %.1f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), massBins, minMass, maxMass, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
        
            } else {
                hInvariantMass2D.push_back(new TH2D(Form("histMass%zu_jetpt_from_%.0f_to_%.0fGeV",i+1,jetptMin,jetptMax), Form("%.0f < #it{p}_{T, D^{0}} < %.0f GeV/#it{c};#it{M}(K#pi) (GeV/#it{c}^{2});#DeltaR",ptDBinEdges_detector[i],ptDBinEdges_detector[i+1]), massBins, minMass, maxMass, deltaRBinEdges_detector.size() - 1, deltaRBinEdges_detector.data()));
        
            }
        }
        hInvariantMass2D[i]->Sumw2();
    }
    // ----- Fill data
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCDetaCut = 0.9 - jetRadius; // on detector level jet
    const double MCDyCut = 0.8; // on detector level D0
    const double MCDDeltaRcut = deltaRBinEdges_detector[deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    const double MCDjetptMin = jetptMin;
    const double MCDjetptMax = jetptMax;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;
    // variables for reflection contribution selection
    int MCDhfMatchedFrom, MCDhfSelectedAs;
    // Accessing TTree
    TTree* tree = (TTree*)fClosureInput->Get("InputTree");
    // Check for correct access
    if (!tree) {
        cout << "Error opening correction data tree.\n";
    }
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
    tree->SetBranchAddress("fHfMlScore0",&MCDhfMlScore0);
    tree->SetBranchAddress("fHfMlScore1",&MCDhfMlScore1);
    tree->SetBranchAddress("fHfMlScore2",&MCDhfMlScore2);
    tree->SetBranchAddress("fHfMatchedFrom",&MCDhfMatchedFrom);
    tree->SetBranchAddress("fHfSelectedAs",&MCDhfSelectedAs);
    int nEntries = tree->GetEntries();
    for (int entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // calculating delta R
        double MCDdeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));
        
        // Fill each histogram with their respective pT intervals
        if ((abs(MCDhfEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && ((MCDjetPt >= MCDjetptMin) && (MCDjetPt < MCDjetptMax)) && ((MCDdeltaR >= deltaRBinEdges_detector[0]) && (MCDdeltaR < MCDDeltaRcut))) {
            
            bool filled = false;
            // Loop through pT,D bin edges to find the appropriate histogram and fill it
            for (size_t iEdge = 0; iEdge < ptDBinEdges_detector.size() - 1 && !filled; iEdge++) {
                if ((MCDhfPt >= ptDBinEdges_detector[iEdge]) && (MCDhfPt < ptDBinEdges_detector[iEdge + 1])) {
                    // Get the threshold for this pT range
                    double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, bdtPtCuts);

                    // Fill histogram only if the cut is passed
                    if (MCDhfMlScore0 < maxBkgProb) {
                        hInvariantMass2D[iEdge]->Fill(MCDhfMass, MCDdeltaR);
                    }
                    filled = true; // Exit the loop once the correct histogram is found (alternative: break)
                }
                
            }
            
        }

        
        
        
    }
    // ----- Obtain 1D projections: invariant mass distributions
    std::vector<TH1D*> hInvariantMass1D;
    for (size_t iHisto = 0; iHisto < hInvariantMass2D.size() - 1; iHisto++) {

        hInvariantMass1D.push_back(hInvariantMass2D[iHisto]->ProjectionX(Form("h_mass_proj_%zu", iHisto)));
    }
    
    // ----- Obtain template fits from MC particle level data
    FitContainer fitTemplates = calculateFitTemplates();

    // ----- Fit detector level data

    // ----- Perform side-band subtraction to fitted detector level data

    // ----- Create 2D distribution of DeltaR vs. pT,D0
    TH2D* hDeltaR_vs_PtD;

    // Store images
    storeImages(hInvariantMass2D, hInvariantMass1D, jetptMin, jetptMax);

    // Clean vector of pointers
    for (auto hist : hInvariantMass2D) {
        delete hist;
    }
    hInvariantMass2D.clear();
    for (auto hist : hInvariantMass1D) {
        delete hist;
    }
    hInvariantMass1D.clear();

    // Return resulting ditribution
    return hDeltaR_vs_PtD;
}

TH3D* create3DBackgroundSubtracted(const std::vector<TH2D*>& outputHistograms) {
    //

    TH3D* hBackgroundSubtracted;

    return hBackgroundSubtracted;
}

TH3D* SidebandClosure(TFile* fClosureInput, const std::vector<double>& ptjetBinEdges_detector, const std::vector<double>& deltaRBinEdges_detector, const std::vector<double>& ptDBinEdges_detector, const std::vector<std::pair<double, double>>& bdtPtCuts) {

    // One TH2D for each pT,jet range computed
    std::vector<TH2D*> outputHistograms;
    for (size_t iPtJetRange = 0; iPtJetRange < ptjetBinEdges_detector.size() - 1; iPtJetRange++) {
        TH2D* hDeltaR_vs_PtD = AnalyzeJetPtRange(fClosureInput, ptjetBinEdges_detector[iPtJetRange], ptjetBinEdges_detector[iPtJetRange + 1], deltaRBinEdges_detector, ptDBinEdges_detector, bdtPtCuts);
        outputHistograms.push_back(hDeltaR_vs_PtD);
    }

    TH3D* hBackgroundSubtracted = create3DBackgroundSubtracted(outputHistograms);
    return hBackgroundSubtracted;
}
