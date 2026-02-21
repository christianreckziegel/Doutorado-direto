/**
 * Lambda_c hadron analysis
 * @file Efficiency_run2_style.C
 * @brief Calculates pT dependent efficiency of HF selections using methodology from run 2
 * Input: series of backSub_%.jetptMin_to_%.jetptMax_jetpt_with_reflections.root files -> contain histograms such that the background+reflections contribution was removed
 * Outputs: one run2_style_efficiency_%.jetptMin_to_%.jetptMax_jetpt.root -> contain pT,HF dependent run 2 style efficiency, to be used by Efficiency_run3_Style.C
 * 
 * @author: Christian Reckziegel
 * Date: February 2026
 */
#include <set>
#include "../commonUtilities.h"

using namespace std;

struct EfficiencyData {
    std::vector<TH1D*> hMcpPt;              // particle level pT,D distribution: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hMcdPt;              // detector level pT,D distribution: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hEfficiencies;       // efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hBackSubCorrected;   // background subtracted (sideband method) DeltaR distributions corrected by efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2
    TH1D* hBackSubCorrected_allpt;          // resulting summed all pT,D bins distribution from sidebad method
    std::vector<TH1D*> hSigExtCorrected;    // background subtracted (signal extraction method) DeltaR distributions corrected by efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2

    // Investigation histograms
    TH1D* hBDTBackgroundScore;              // rejected entries that didn't pass the BDT background score cuts
};

// Module to create histograms including interest variable
EfficiencyData createHistograms(const BinningStruct& binning) {

    // Create struct to store data
    EfficiencyData dataContainer;

    // Create 3 histogram cases: inclusive = 0, prompt only = 1, non-prompt only = 2
    int histCaseNum = 3;
    for (size_t iCase = 0; iCase < histCaseNum; ++iCase) {
        dataContainer.hMcpPt.push_back(new TH1D(Form("mcp_pt_%zu",iCase), ";#it{p}_{T,#Lambda_{c}^{+}}^{truth};dN/d#it{p}_{T,#Lambda_{c}^{+}}^{truth}", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data()));
        dataContainer.hMcpPt[iCase]->SetMarkerColor(kBlack);
        dataContainer.hMcpPt[iCase]->SetLineColor(kBlack);
        dataContainer.hMcpPt[iCase]->SetMarkerStyle(kOpenCircle);
        dataContainer.hMcpPt[iCase]->Sumw2();
        dataContainer.hMcpPt[iCase]->SetStats(0);
        dataContainer.hMcdPt.push_back(new TH1D(Form("mcd_pt_%zu",iCase), ";#it{p}_{T,D}^{reco};dN/d#it{p}_{T,D}^{reco}", binning.ptHFBinEdges_detector.size() - 1, binning.ptHFBinEdges_detector.data()));
        dataContainer.hMcdPt[iCase]->SetMarkerColor(kBlue);
        dataContainer.hMcdPt[iCase]->SetLineColor(kBlue);
        dataContainer.hMcdPt[iCase]->SetMarkerStyle(kFullCircle);
        dataContainer.hMcdPt[iCase]->Sumw2();
        dataContainer.hMcdPt[iCase]->SetStats(0);
    }

    // Creating investigation histogram
    dataContainer.hBDTBackgroundScore = new TH1D("hBDTBackgroundScore", "Entries that didn't pass the cuts;BDT background score;Counts", 100, 0, 1);

    cout << "Histograms created.\n";

    return dataContainer;
}

// Module to fill histograms from TFile data
void fillHistograms(TFile* fSimulated, EfficiencyData& dataContainer, const BinningStruct& binning, double& jetptMin, double& jetptMax) {
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double etaCut = 0.9 - jetRadius; // on particle level jet
    const double yCut = 0.8; // on detector level Lc
    const double deltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on particle level delta R
    double hfptMin = binning.ptHFBinEdges_detector[0]; //ptHFBinEdges[0] - should start from 0 or from the lowest pT,Lc value?
    double hfptMax = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1];
    const double MCPHfPtMincut = hfptMin; // on particle level HF
    const double MCDHfPtMincut = hfptMin; // on detector level HF
    const double MCPHfPtMaxcut = hfptMax; // on particle level HF
    const double MCDHfPtMaxcut = hfptMax; // on detector level HF

    // Accessing detector level data TTree
    TTree* tree = (TTree*)fSimulated->Get("DF_merged/O2matchtable");

    // Assuming histograms and tree data correspond in some way
    if (!tree) {
        cout << "Error opening tree.\n";
    }
    // defining variables for accessing particle level data on TTree
    float MCPaxisDistance, MCPjetPt, MCPjetEta, MCPjetPhi;
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY, MCPjetNConst;
    bool MCPhfprompt;
    // defining variables for accessing detector level data on TTree
    float MCDaxisDistance, MCDjetPt, MCDjetEta, MCDjetPhi;
    float MCDhfPt, MCDhfEta, MCDhfPhi, MCDhfMass, MCDhfY;
    bool MCDhfprompt;
    int MCDhfMatchedFrom, MCDhfSelectedAs, MCDjetNConst;
    // defining ML score variables for accessing the TTree
    float MCDhfMlScore0, MCDhfMlScore1, MCDhfMlScore2;

    // particle level branches
    tree->SetBranchAddress("fMcJetHfDist",&MCPaxisDistance);
    tree->SetBranchAddress("fMcJetPt",&MCPjetPt);
    tree->SetBranchAddress("fMcJetEta",&MCPjetEta);
    tree->SetBranchAddress("fMcJetPhi",&MCPjetPhi);
    tree->SetBranchAddress("fMcJetNConst",&MCPjetNConst); // double
    tree->SetBranchAddress("fMcHfPt",&MCPhfPt);
    tree->SetBranchAddress("fMcHfEta",&MCPhfEta);
    tree->SetBranchAddress("fMcHfPhi",&MCPhfPhi);
    MCPhfMass = 2.28646; // Lc rest mass in GeV/c^2
    tree->SetBranchAddress("fMcHfY",&MCPhfY);
    tree->SetBranchAddress("fMcHfPrompt",&MCPhfprompt);
    // detector level branches
    tree->SetBranchAddress("fJetHfDist",&MCDaxisDistance);
    tree->SetBranchAddress("fJetPt",&MCDjetPt);
    tree->SetBranchAddress("fJetEta",&MCDjetEta);
    tree->SetBranchAddress("fJetPhi",&MCDjetPhi);
    tree->SetBranchAddress("fJetNConst",&MCDjetNConst); // int
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

        // Particle level histograms
        double MCPdeltaR = sqrt(pow(MCPjetEta-MCPhfEta,2) + pow(DeltaPhi(MCPjetPhi,MCPhfPhi),2)); // one can also use MCPaxisDistance
        bool genLevelRange = (abs(MCPhfEta) < etaCut) && (abs(MCPhfY) < yCut) && (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax) && ((MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < deltaRcut)) && ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut));
        // Fill histograms considering jet pT and ALICE's acceptance
        if (genLevelRange) {
            
            // Fill inclusive histogram
            dataContainer.hMcpPt[0]->Fill(MCPhfPt);
            // fill prompt efficiency histogram
            if (MCPhfprompt) {
                dataContainer.hMcpPt[1]->Fill(MCPhfPt);
            } else{
                // fill non-prompt efficiency histogram
                dataContainer.hMcpPt[2]->Fill(MCPhfPt);
            }
        }

        // Detector level histograms
        // only compute matched detector level candidates, but compute all particle level ones
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isReflection = (MCDhfMatchedFrom != MCDhfSelectedAs) ? true : false;
        if (!MCDhfmatch || isReflection) { // !MCDhfmatch || isReflection
            continue;
        }
        
        // calculating delta R
        double MCDdeltaR = sqrt(pow(MCDjetEta-MCDhfEta,2) + pow(DeltaPhi(MCDjetPhi,MCDhfPhi),2));
        bool recoLevelRange = (abs(MCDhfEta) < etaCut) && (abs(MCDhfY) < yCut) && (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax) && ((MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < deltaRcut)) && ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));

        // Fill histograms considering jet pT and detector acceptance
        if (recoLevelRange) {
            
            // Get the threshold for this pT range
            double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);

            // Fill histogram only if the BDT cut is passed
            if (MCDhfMlScore0 < maxBkgProb) {
                // Fill inclusive histogram
                dataContainer.hMcdPt[0]->Fill(MCDhfPt);
                // fill prompt efficiency histogram
                if (MCDhfprompt) {
                    dataContainer.hMcdPt[1]->Fill(MCDhfPt);
                } else{
                    // fill non-prompt efficiency histogram
                    dataContainer.hMcdPt[2]->Fill(MCDhfPt);
                }
            } else {
                dataContainer.hBDTBackgroundScore->Fill(MCDhfMlScore0);
            }
                        
        }
    }
    std::cout << "Both generator and reconstruction level histograms filled.\n";
}

void CalculateEfficiency(EfficiencyData& dataContainer){
    const char* efficiencyNames[] = {"inclusive", "prompt", "nonprompt"};

    // Loop through histogram cases: inclusive = 0, prompt only = 1, non-prompt only = 2
    int histCaseNum = 3;
    for (int iEff = 0; iEff < histCaseNum; iEff++) {
        // Obtain MC pT distributions
        TH1D* mcdHist = dataContainer.hMcdPt[iEff];
        TH1D* mcpHist = dataContainer.hMcpPt[iEff];

        // Otain efficiency by dividing tempHist by parameter -> mcd/mcp
        TH1D* efficiencyHist = static_cast<TH1D*>(mcdHist->Clone(Form("efficiency_%s", efficiencyNames[iEff]))); // static needed to explicitly cast the cloned histogram to TH1D. ensures that the function returns a TH1D*, matching the return type.
        efficiencyHist->Divide(mcpHist);
        efficiencyHist->SetMarkerStyle(kFullCircle);
        efficiencyHist->SetMarkerColor(30+iEff*10); // 30 = not so bright green
        efficiencyHist->SetLineColor(30+iEff*10);
        efficiencyHist->SetTitle(";#it{p}_{T,#Lambda_{c}^{+}}^{truth};Efficiency#times Acceptance");
        efficiencyHist->SetStats(0);
        
        // Store efficiency distribution in struct
        dataContainer.hEfficiencies.push_back(efficiencyHist);

    }

    std::cout << "Efficiency calculated.\n";
}

void PerformCorrections(EfficiencyData& dataContainer, TFile* fBackSub) {
    
}

void PlotHistograms(const EfficiencyData& dataContainer, double jetptMin, double jetptMax) {

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);

    std::vector<TLegend*> legendPt;
    legendPt.push_back(new TLegend(0.65,0.49,0.8,0.62));
    legendPt[0]->AddEntry(dataContainer.hMcpPt[0],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    legendPt[0]->AddEntry(dataContainer.hMcdPt[0],"Reconstructed", "lpe");
    legendPt.push_back(new TLegend(0.65,0.49,0.8,0.62));
    legendPt[1]->AddEntry(dataContainer.hMcpPt[1],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    legendPt[1]->AddEntry(dataContainer.hMcdPt[1],"Reconstructed", "lpe");
    legendPt.push_back(new TLegend(0.65,0.49,0.8,0.62));
    legendPt[2]->AddEntry(dataContainer.hMcpPt[2],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    legendPt[2]->AddEntry(dataContainer.hMcdPt[2],"Reconstructed", "lpe");
    legendPt.push_back(new TLegend(0.65,0.47,0.85,0.62));
    legendPt[3]->AddEntry(dataContainer.hMcpPt[0],"Truth", "lpe"); // inclusive efficiency only is used for efficiency correction
    legendPt[3]->AddEntry(dataContainer.hMcdPt[0],"Inclusive reconstructed", "lpe");
    legendPt[3]->AddEntry(dataContainer.hMcdPt[1],"Prompt #Lambda_{c}^{+} reconstructed", "lpe");
    legendPt[3]->AddEntry(dataContainer.hMcdPt[2],"Non-prompt #Lambda_{c}^{+} reconstructed", "lpe");

    
    //
    // 1 - plot pT,D of generator and reconstruction level
    //
    TCanvas* cEff_1 = new TCanvas("cEff_1","pT,D spectra");
    cEff_1->Divide(2,2);
    cEff_1->SetCanvasSize(1800,1000);
    cEff_1->cd(1);
    dataContainer.hMcpPt[0]->Draw();
    dataContainer.hMcdPt[0]->Draw("same");
    legendPt[0]->Draw();
    double statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Inclusive");
    cEff_1->cd(2);
    dataContainer.hMcpPt[1]->Draw();
    dataContainer.hMcdPt[1]->Draw("same");
    legendPt[1]->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Prompt #Lambda_{c}^{+}'s");
    cEff_1->cd(3);
    dataContainer.hMcpPt[2]->Draw();
    dataContainer.hMcdPt[2]->Draw("same");
    legendPt[2]->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Non-prompt #Lambda_{c}^{+}'s");
    cEff_1->cd(4);
    dataContainer.hMcdPt[0]->SetMarkerColor(30+0*10); // 30 = not so bright green
    dataContainer.hMcdPt[0]->SetLineColor(30+0*10); // 30 = not so bright green
    dataContainer.hMcdPt[1]->SetMarkerColor(30+1*10); // 30 = not so bright green
    dataContainer.hMcdPt[1]->SetLineColor(30+1*10); // 30 = not so bright green
    dataContainer.hMcdPt[2]->SetMarkerColor(30+2*10); // 30 = not so bright green
    dataContainer.hMcdPt[2]->SetLineColor(30+2*10); // 30 = not so bright green
    dataContainer.hMcdPt[0]->Draw();
    dataContainer.hMcdPt[1]->Draw("same");
    dataContainer.hMcdPt[2]->Draw("same");
    legendPt[3]->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Reconstructed: inclusive, prompt,\n non-prompt");

    //
    // 2 - plot three kinds of efficiencies: [0] = inclusive, [1] = prompt D0, [2] = non-prompt D0
    //
    // Loop through efficiency cases
    TCanvas* cEff_2 = new TCanvas("cEff_2","Prompt/non-prompt fficiency plots");
    cEff_2->SetCanvasSize(1800,1000);
    cEff_2->cd();
    TLegend* legendEff = new TLegend(0.65,0.45,0.8,0.58);
    legendEff->AddEntry(dataContainer.hEfficiencies[0],"Inclusive", "lpe");
    legendEff->AddEntry(dataContainer.hEfficiencies[1],"Prompt #Lambda_{c}^{+}", "lpe");
    legendEff->AddEntry(dataContainer.hEfficiencies[2],"Non-prompt #Lambda_{c}^{+}", "lpe");
    for (size_t iEff = 0; iEff < dataContainer.hEfficiencies.size(); iEff++) {
        // if the current is the first plot in the canvas
        if (iEff == 0) {
            dataContainer.hEfficiencies[iEff]->SetMinimum(0);
            dataContainer.hEfficiencies[iEff]->SetMaximum(0.8);
            dataContainer.hEfficiencies[iEff]->Draw();
        } else {
            dataContainer.hEfficiencies[iEff]->Draw("same");
        }
        
        statBoxPos = gPad->GetUxmax();
        latex->DrawLatex(statBoxPos-0.35, 0.61, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    }
    legendEff->Draw();
    cEff_2->Update();

    //
    //
    //
    TCanvas* cInvestigation = new TCanvas("cInvestigation","Investigation");
    cInvestigation->SetCanvasSize(1800,1000);
    cInvestigation->cd();
    dataContainer.hBDTBackgroundScore->Draw();
    cInvestigation->SetCanvasSize(1800,1000);
    //
    // Storing images
    //
    TString imagePath = "../Images/2-Efficiency/Run2Style/";
    cEff_1->Update();
    cEff_1->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cEff_1->SaveAs(imagePath + "Efficiency_dNdpT_run2_style.png");
    cEff_2->Update();
    cEff_2->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cEff_2->SaveAs(imagePath + "Efficiency_acceptance_run2_style.png");

    //
    // Storing in a single pdf file
    //
    cEff_1->Print(Form(imagePath + "pT_efficiency_%.0f_to_%.0fGeV.pdf(",jetptMin,jetptMax));
    cEff_2->Print(Form(imagePath + "pT_efficiency_%.0f_to_%.0fGeV.pdf)",jetptMin,jetptMax));

    std::cout << "Histograms plotted.\n";
}

void SaveData(const EfficiencyData& dataContainer, const BinningStruct& binning, double& jetptMin, double& jetptMax){
    // Open output file
    TFile* fOutput = new TFile(Form("run2_style_efficiency_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    // Loop over efficiency histograms
    for (size_t iHisto = 0; iHisto < dataContainer.hEfficiencies.size(); iHisto++) {
        // store each histogram in file
        dataContainer.hEfficiencies[iHisto]->Write();
    }

    // Also store the axes used for the histograms
    storeBinningInFile(fOutput, binning);

    fOutput->Close();
    delete fOutput;
    
    std::cout << "Data stored.\n";
}



void Efficiency_run2_style(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // Load binning information
    TFile* fBinning = new TFile(Form("../1-SignalTreatment/Sideband/full_merged_ranges_back_sub.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    // Load binning from background subtracted file
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

    // Opening files
    TFile* fSimulated = new TFile("../Data/MonteCarlo/Train_607573/AO2D_mergedDFs.root","read","read");

    EfficiencyData dataContainer = createHistograms(binning); // pT histograms

    // Fill histograms
    fillHistograms(fSimulated, dataContainer, binning, jetptMin, jetptMax);

    // Calculate efficiency distribution
    CalculateEfficiency(dataContainer);

    // Apply efficiency corrections to distributions
    //PerformCorrections(dataContainer, fBackSub);

    // Plot the efficiency histogram and further corrected histograms
    PlotHistograms(dataContainer, jetptMin, jetptMax);

    // Save corrected distributions to file
    SaveData(dataContainer, binning, jetptMin, jetptMax);

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    Efficiency_run2_style();
    return 0;
}
