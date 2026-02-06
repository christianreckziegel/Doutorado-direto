/**
 * Lambda_c hadron analysis
 * @file BackgroundSubtraction.C
 * @brief Background subtraction using sideband method
 * Input: Reflections.root: contain histograms and fit templates of Lambda_c signal and reflections
 * Outputs: BackgroundSubtraction.root
 * 
 * @author: Christian Reckziegel
 * Date: February 2026
 */

 #include "../../commonUtilities.h"
using namespace std;

struct SidebandData {
    std::vector<TH2D*> histograms2d;
    std::vector<TH1D*> histograms1d;
};

// Module to create TH2D histograms including interest variable: VARIABLE bin sizes
SidebandData createHistograms(const std::vector<double>& ptHFBinEdges, int xbins, double xmin, double xmax, const std::vector<double>& yBinEdges, const double& jetptMin) {
    // Change binning to low statistics range [30,50] GeV/c
    if (jetptMin >= 30.) {
        xbins = xbins / 2;
    }
    
    SidebandData dataContainer;
    for (size_t i = 0; i < ptHFBinEdges.size() - 1; ++i) {
        // So that the title adapts to fractional binning title
        if (std::fmod(ptHFBinEdges[i], 1.0) != 0) { // if the first bin edge is not an integer
            if (std::fmod(ptHFBinEdges[i+1], 1.0) != 0) {
                dataContainer.histograms2d.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.1f < #it{p}_{T, #Lambda_{c}} < %.1f GeV/#it{c};#it{M}(#piKPr) (GeV/#it{c}^{2});#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            } else {
                dataContainer.histograms2d.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.1f < #it{p}_{T, #Lambda_{c}} < %.0f GeV/#it{c};#it{M}(#piKPr) (GeV/#it{c}^{2});#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            }
        } else {
            if (std::fmod(ptHFBinEdges[i+1], 1.0) != 0) {
                dataContainer.histograms2d.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < #it{p}_{T, #Lambda_{c}} < %.1f GeV/#it{c};#it{M}(#piKPr) (GeV/#it{c}^{2});#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            } else {
                dataContainer.histograms2d.push_back(new TH2D(Form("histMass%zu",i+1), Form("%.0f < #it{p}_{T, #Lambda_{c}} < %.0f GeV/#it{c};#it{M}(#piKPr) (GeV/#it{c}^{2});#DeltaR",ptHFBinEdges[i],ptHFBinEdges[i+1]), xbins, xmin, xmax, yBinEdges.size() - 1, yBinEdges.data()));
        
            }
        }
        dataContainer.histograms2d[i]->Sumw2();
    }

    return dataContainer;
}
//__________________________________________________________________________________________________________________________

void fillHistograms(TFile* fDist, SidebandData& dataContainer, double jetptMin, double jetptMax, std::vector<double>& ptHFBinEdges, std::vector<double>& deltaRBinEdges, const std::vector<std::pair<double, double>>& bdtPtCuts) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level D0
    const double DeltaRcut = deltaRBinEdges[deltaRBinEdges.size() - 1]; // on particle level delta R

    // Accessing TTree
    TTree* tree = (TTree*)fDist->Get("DF_2261906155510144/O2jetdisttable");

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
            for (size_t iEdge = 0; iEdge < ptHFBinEdges.size() - 1 && !filled; iEdge++) {
                if ((hfPt >= ptHFBinEdges[iEdge]) && (hfPt < ptHFBinEdges[iEdge + 1])) {
                    // Get the threshold for this pT range
                    double maxBkgProb = GetBkgProbabilityCut(hfPt, bdtPtCuts);

                    // Fill histogram only if the cut is passed
                    if (hfMlScore0 < maxBkgProb) {
                        dataContainer.histograms2d[iEdge]->Fill(hfMass, deltaR);
                    }
                    filled = true; // Exit the loop once the correct histogram is found (alternative: break)
                }
                
            }
            
        } // kinematic range choice

    } // end of TTree entries loop
    cout << "Histograms filled.\n";

    // Creating 1D mass projections
    TH1D* tempHist;
    // Obtaining 1D invariant mass histograms from the projection
    for (size_t iHist = 0; iHist < dataContainer.histograms2d.size(); iHist++) {
        //
        tempHist = dataContainer.histograms2d[iHist]->ProjectionX(Form("h_mass_proj_%zu", iHist));
        dataContainer.histograms1d.push_back(tempHist);
    }
}







void PlotHistograms(const SidebandData& dataContainer, double jetptMin, double jetptMax, const std::vector<double>& ptHFBinEdges) {
    std::cout << "Plotting histograms..." << std::endl;


    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.04); // default = 0.05

    // Create a canvas for plotting
    int nHistos = dataContainer.histograms1d.size();
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
    for(size_t iHisto = 0; iHisto < dataContainer.histograms1d.size(); ++iHisto) {
        //
        c1d_fit->cd(iHisto+1);
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        gStyle->SetOptStat(0); // Turn off the default stats box
        dataContainer.histograms1d[iHisto]->SetMarkerStyle(kDot); //kFullDotMedium
        dataContainer.histograms1d[iHisto]->SetMarkerColor(kBlack);
        dataContainer.histograms1d[iHisto]->SetLineColor(kBlack);
        dataContainer.histograms1d[iHisto]->GetYaxis()->SetTitle("counts");
        dataContainer.histograms1d[iHisto]->SetMinimum(0);
        
        dataContainer.histograms1d[iHisto]->Draw();
        // fittings.fitTotal[iHisto]->Draw("same");
        // fittings.fitSignalOnly[iHisto]->Draw("same");
        // fittings.fitBackgroundOnly[iHisto]->Draw("same");
        // fittings.fitReflectionsOnly[iHisto]->Draw("same");

        // double A1Signal = fittings.fitTotal[iHisto]->GetParameter(2); // Get the value of parameter 'A' from primary signal gaussian
        // double m0_1Signal = fittings.fitTotal[iHisto]->GetParameter(4); // Get the value of parameter 'm0' from primary signal gaussian
        // double sigma1Signal = fittings.fitTotal[iHisto]->GetParameter(5); // Get the value of parameter 'sigma' from primary signal gaussian
        // double A2Signal = A1Signal / fittings.fitTotal[iHisto]->GetParameter(3); // Get the value of parameter 'A' from secondary signal gaussian
        // double sigma2Signal = sigma1Signal / fittings.fitTotal[iHisto]->GetParameter(6); // Get the value of parameter 'sigma' from secondary signal gaussian
        // double chi2 = fittings.fitTotal[iHisto]->GetChisquare();
        // double degOfFreedom = fittings.fitTotal[iHisto]->GetNDF();
        
        // latex->DrawLatex(statBoxPos-0.52, 0.80, Form("A_{1}^{signal} = %.2f, #bar{m_{1}} = %.2f, #sigma_{1} = %.2f GeV/#it{c}^{2}", A1Signal, m0_1Signal,sigma1Signal));
        // latex->DrawLatex(statBoxPos-0.52, 0.75, Form("A_{2}^{signal} = %.2f, #bar{m_{2}} = %.2f, #sigma_{2} = %.2f GeV/#it{c}^{2}", A2Signal, m0_1Signal,sigma2Signal));
        latex->DrawLatex(statBoxPos-0.3, 0.70, Form("%.0f < p_{T, ch. jet} < %.0f GeV/#it{c}",jetptMin,jetptMax)); // Display jet pT cut applied
        // latex->DrawLatex(statBoxPos-0.25, 0.63, Form("#Chi^{2}_{red} = %.3f",chi2/degOfFreedom));
        latex->DrawLatex(statBoxPos-0.32, 0.58, "MC: HF_LHC24g5_All");
        latex->DrawLatex(statBoxPos-0.4, 0.53, "Data: HF_LHC22o_pass7_minBias_small_2P3PDstar");

        // Drawing 2D histograms
        c_2d->cd(iHisto+1);
        dataContainer.histograms2d[iHisto]->Draw("colz");

        //
        // Storing images
        //
        TString imagePath = "../../Images/1-SignalTreatment/SideBand/";

        //
        // Storing in a single pdf file
        //
        c1d_fit->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf(",jetptMin,jetptMax));
        c_2d->Print(imagePath + Form("sb_subtraction_deltaR_%.0f_to_%.0fGeV_withReflections.pdf)",jetptMin,jetptMax));

    }
}


void JetPtIterator(const double jetptMin, const double jetptMax, const std::vector<double>& ptjetBinEdges) {
    std::cout << "============================================= " << jetptMin << " GeV/c < pT,jet < " << jetptMax << " GeV/c =============================================" << std::endl;

    // BDT background probability cuts based on pT,D ranges. Example: 1-2 GeV/c -> 0.05 (from first of pair)
    std::vector<std::pair<double, double>> bdtPtCuts = {
        {1, 0.05}, {2, 0.08}, {3, 0.08}, {4, 0.1}, {5, 0.2}, {6, 0.25}, {7, 0.3}, {8, 0.4}, {12, 0.5}, {24, 0.7}, {100, 0.7}
    }; // on dataset HF_LHC22o_pass7_minBias_small_2P3PDstar

    // Lc mass in GeV/c^2
    double m_0_parameter = 2.28646;
    double sigmaInitial = 0.012;

    // mass histogram
    int massBins = 50; // default=50 
    double minMass = 2.1; // default = 2.1
    double maxMass = 2.49; // default = 2.49

    // Opening data file
    TFile* fDist = new TFile("../../Data/Experimental/Train_606159/AO2D.root","read");
    if (!fDist || fDist->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT data file." << std::endl;
    }
    // Opening
    TFile* fReflectionsMC = new TFile(Form("../Reflections/reflections_%.0f_to_%.0fGeV.root",jetptMin,jetptMax),"read");
    if (!fReflectionsMC || fReflectionsMC->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT reflections file." << std::endl;
    }
    // Load ΔR bin edges
    // std::vector<double> deltaRBinEdges = LoadBinning(fReflectionsMC, "axes/deltaRBinEdges");
    std::vector<double> deltaRBinEdges = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5};
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
    // Load pTD bin edges
    // std::vector<double> ptHFBinEdges = LoadBinning(fReflectionsMC, "axes/ptHFBinEdges");
    std::vector<double> ptHFBinEdges = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.};
    // Load binning and BDT cuts from reflections file
    BinningStruct binning = retrieveBinningFromFile(fReflectionsMC);

    // Create multiple histograms
    SidebandData dataContainer = createHistograms(ptHFBinEdges,                                 // the pT,HF edges will determine the number of mass histograms
                                                     massBins, minMass, maxMass,                  // mass histograms binning
                                                     deltaRBinEdges, jetptMin);                             // deltaR histograms with asymmetrical bin widths
    //
    // Fill histograms
    fillHistograms(fDist, dataContainer, jetptMin, jetptMax, ptHFBinEdges, deltaRBinEdges, bdtPtCuts);

    // Perform fits

    // signal/side-band region parameters
    double signalSigmas = 2; // default delta = 2
    int startingBackSigma = 4; // default position = 4
    int backgroundSigmas = 4; // default delta = 4

    // Subtract side-band from signal

    // Plot histograms
    PlotHistograms(dataContainer, jetptMin, jetptMax, ptHFBinEdges);

    // Storing final histograms to output file

}

void create3DBackgroundSubtracted(const std::vector<double>& ptjetBinEdges, const std::vector<double>& deltaRBinEdges, const std::vector<double>& ptHFBinEdges){

    // jet pT cuts
    double jetptMin = ptjetBinEdges[0]; // GeV
    double jetptMax = ptjetBinEdges[ptjetBinEdges.size() - 1]; // GeV
    // deltaR histogram
    double minDeltaR = deltaRBinEdges[0];
    double maxDeltaR = deltaRBinEdges[deltaRBinEdges.size() - 1];
}


void BackgroundSubtraction() {
    
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // Opening
    double jetptMin = 5.; // default = 5 GeV, for Emma's use 5 GeV
    bool useEmmaYeatsBins = false;
    double jetptMax = useEmmaYeatsBins ? 10. : 7.; // default = 7GeV, for Emma's use 10 GeV
    // TFile* fReflectionsMC = new TFile(Form("../Reflections/reflections_%.0f_to_%.0fGeV.root",jetptMin,jetptMax),"read");
    // if (!fReflectionsMC || fReflectionsMC->IsZombie()) {
    //     std::cerr << "Error: Unable to open the first ROOT reflections file." << std::endl;
    // }
    // Load pTjet bin edges
    // std::vector<double> ptjetBinEdges = LoadBinning(fReflectionsMC, "axes/ptjetBinEdges");
    std::vector<double> ptjetBinEdges = {5., 7., 10., 15., 30., 50.};
    // std::vector<double> ptjetBinEdges = {30., 50.};
    // Load ΔR bin edges
    // std::vector<double> deltaRBinEdges = LoadBinning(fReflectionsMC, "axes/deltaRBinEdges");
    std::vector<double> deltaRBinEdges = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5};
    // Load pTD bin edges
    // std::vector<double> ptHFBinEdges = LoadBinning(fReflectionsMC, "axes/ptHFBinEdges");
    std::vector<double> ptHFBinEdges = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 18., 30.};
    // fReflectionsMC->Close();


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
    create3DBackgroundSubtracted(ptjetBinEdges, deltaRBinEdges, ptHFBinEdges);

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl;
}

int main() {
    BackgroundSubtraction();
    return 0;
}