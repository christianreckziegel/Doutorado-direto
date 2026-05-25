/**
 * D0 meson analysis
 * @file Efficiency_run2_style.C
 * @brief Calculates pT dependent efficiency of HF selections using methodology from run 2
 * Input: series of backSub_%.jetptMin_to_%.jetptMax_jetpt.root files -> contain histograms such that the background+reflections contribution was removed
 * Outputs: one run2_style_efficiency_%.jetptMin_to_%.jetptMax_jetpt.root -> contain pT,HF dependent run 2 style efficiency, to be used by Efficiency_run3_Style.C
 * 
 * @author: Christian Reckziegel
 * Date: February 2026
 */
#include <set>
#include "../commonUtilities.h"

using namespace std;

struct EfficiencyData {
    std::vector<TH1D*> hMcpPt;              // particle level pT,HF distribution: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hMcdPt;              // detector level pT,HF distribution: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hEfficiencies;       // efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2
    std::vector<TH1D*> hBackSubCorrected;   // background subtracted (sideband method) DeltaR distributions corrected by efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2
    TH1D* hBackSubCorrected_allpt;          // resulting summed all pT,HF bins distribution from sidebad method
    std::vector<TH1D*> hSigExtCorrected;    // background subtracted (signal extraction method) DeltaR distributions corrected by efficiency: inclusive = 0, prompt only = 1, non-prompt only = 2

    // 3D efficiency calculation
    TH3D* hMcpPt_vs_ptJet_vs_deltaR;
    TH3D* hMcdPt_vs_ptJet_vs_deltaR;
    TH3D* hSelEff_run2style3d;
    std::vector<TH1D*> hSelEff_run2style_per_jetpt;
    std::vector<TH1D*> hSelEff_run2style_per_deltaR;

    // 2D efficiency calculation
    std::vector<TH2D*> hSelEff_run2style2d;

    // Run 3 particle level prompt selection efficiency
    TH1D* hDenominator_prompt_part_run3style;
    TH1D* hNumerator_prompt_part_run3style;
    TH1D* hEfficiency_prompt_part_run3style;

    // Run 3 particle level non-prompt selection efficiency
    TH1D* hDenominator_nonprompt_part_run3style;
    TH1D* hNumerator_nonprompt_part_run3style;
    TH1D* hEfficiency_nonprompt_part_run3style;

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
        dataContainer.hMcpPt.emplace_back(new TH1D(Form("mcp_pt_%zu",iCase), ";#it{p}_{T,D^{0}}^{truth};dN/d#it{p}_{T,D^{0}}^{truth}", binning.ptHFEfficiencyBinEdges_detector.size() - 1, binning.ptHFEfficiencyBinEdges_detector.data()));
        dataContainer.hMcpPt[iCase]->SetMarkerColor(kBlack);
        dataContainer.hMcpPt[iCase]->SetLineColor(kBlack);
        dataContainer.hMcpPt[iCase]->SetMarkerStyle(kOpenCircle);
        dataContainer.hMcpPt[iCase]->Sumw2();
        dataContainer.hMcpPt[iCase]->SetStats(0);
        dataContainer.hMcdPt.emplace_back(new TH1D(Form("mcd_pt_%zu",iCase), ";#it{p}_{T,D}^{reco};dN/d#it{p}_{T,D}^{reco}", binning.ptHFEfficiencyBinEdges_detector.size() - 1, binning.ptHFEfficiencyBinEdges_detector.data()));
        dataContainer.hMcdPt[iCase]->SetMarkerColor(kBlue);
        dataContainer.hMcdPt[iCase]->SetLineColor(kBlue);
        dataContainer.hMcdPt[iCase]->SetMarkerStyle(kFullCircle);
        dataContainer.hMcdPt[iCase]->Sumw2();
        dataContainer.hMcdPt[iCase]->SetStats(0);
    }

    // Creating investigation histogram
    dataContainer.hBDTBackgroundScore = new TH1D("hBDTBackgroundScore", "Entries that didn't pass the cuts;BDT background score;Counts", 100, 0, 1);

    // Create 2D version histograms
    dataContainer.hSelEff_run2style2d.emplace_back(new TH2D(Form("mcp2d_pt"), "Prompt D^{0}s denominator (particle level);#it{p}_{T,jet}^{truth};#it{p}_{T,D^{0}}^{truth}", binning.ptjetBinEdges_particle.size() - 1, binning.ptjetBinEdges_particle.data(), binning.ptHFEfficiencyBinEdges_particle.size() - 1, binning.ptHFEfficiencyBinEdges_particle.data()));
    dataContainer.hSelEff_run2style2d.back()->SetStats(0);
    dataContainer.hSelEff_run2style2d.emplace_back(new TH2D(Form("mcd2d_pt"), "Prompt D^{0}s numerator (detector level);#it{p}_{T,jet}^{reco};;#it{p}_{T,D^{0}}^{reco}", binning.ptjetBinEdges_detector.size() - 1, binning.ptjetBinEdges_detector.data(), binning.ptHFEfficiencyBinEdges_detector.size() - 1, binning.ptHFEfficiencyBinEdges_detector.data()));
    dataContainer.hSelEff_run2style2d.back()->SetStats(0);
    
    // Run 3 style particle level prompt selection efficiency histograms
    dataContainer.hDenominator_prompt_part_run3style = new TH1D("mcp_pt_prompt_part_run3style", ";#it{p}_{T,D^{0}}^{truth};dN/d#it{p}_{T,D^{0}}^{truth}", binning.ptHFEfficiencyBinEdges_particle.size() - 1, binning.ptHFEfficiencyBinEdges_particle.data());
    dataContainer.hDenominator_prompt_part_run3style->SetMarkerStyle(kOpenCircle);
    dataContainer.hDenominator_prompt_part_run3style->Sumw2();
    dataContainer.hDenominator_prompt_part_run3style->SetStats(0);
    dataContainer.hNumerator_prompt_part_run3style = new TH1D("mcd_pt_prompt_det_run3style", ";#it{p}_{T,D^{0}}^{truth};dN/d#it{p}_{T,D^{0}}^{truth}", binning.ptHFEfficiencyBinEdges_particle.size() - 1, binning.ptHFEfficiencyBinEdges_particle.data());
    dataContainer.hNumerator_prompt_part_run3style->SetMarkerStyle(kOpenCircle);
    dataContainer.hNumerator_prompt_part_run3style->SetMarkerColor(kBlue);
    dataContainer.hNumerator_prompt_part_run3style->SetLineColor(kBlue);
    dataContainer.hNumerator_prompt_part_run3style->Sumw2();
    dataContainer.hNumerator_prompt_part_run3style->SetStats(0);

    // Run 3 style particle level non-prompt selection efficiency histograms
    dataContainer.hDenominator_nonprompt_part_run3style = new TH1D("mcp_pt_nonprompt_part_run3style", ";#it{p}_{T,D^{0}}^{truth};dN/d#it{p}_{T,D^{0}}^{truth}", binning.ptHFEfficiencyBinEdges_particle.size() - 1, binning.ptHFEfficiencyBinEdges_particle.data());
    dataContainer.hDenominator_nonprompt_part_run3style->SetMarkerStyle(kOpenCircle);
    dataContainer.hDenominator_nonprompt_part_run3style->Sumw2();
    dataContainer.hDenominator_nonprompt_part_run3style->SetStats(0);
    dataContainer.hNumerator_nonprompt_part_run3style = new TH1D("mcd_pt_nonprompt_det_run3style", ";#it{p}_{T,D^{0}}^{truth};dN/d#it{p}_{T,D^{0}}^{truth}", binning.ptHFEfficiencyBinEdges_particle.size() - 1, binning.ptHFEfficiencyBinEdges_particle.data());
    dataContainer.hNumerator_nonprompt_part_run3style->SetMarkerStyle(kOpenCircle);
    dataContainer.hNumerator_nonprompt_part_run3style->SetMarkerColor(kRed);
    dataContainer.hNumerator_nonprompt_part_run3style->SetLineColor(kRed);
    dataContainer.hNumerator_nonprompt_part_run3style->Sumw2();
    dataContainer.hNumerator_nonprompt_part_run3style->SetStats(0);

    cout << "Histograms created.\n";

    return dataContainer;
}

// Module to fill histograms from TFile data
void fillHistograms(TFile* fSimulated, EfficiencyData& dataContainer, const BinningStruct& binning, double& jetptMin, double& jetptMax) {
    
    std::cout << "pT,HF range (GeV)\t\t BDT score" << std::endl;
    for (size_t iHFBin = 0; iHFBin < binning.bdtPtCuts.size() -1; iHFBin++) {
        std::cout << "\t" << binning.bdtPtCuts[iHFBin].first << "-" << binning.bdtPtCuts[iHFBin+1].first << "\t\t\t\t" << binning.bdtPtCuts[iHFBin].second << std::endl;
    }
    
    std::cout << "Efficiency bin edges (GeV):" << std::endl;
    for (size_t iHFBin = 0; iHFBin < binning.ptHFEfficiencyBinEdges_detector.size() - 1; iHFBin++) {
        std::cout << "\t" << binning.ptHFEfficiencyBinEdges_detector[iHFBin] << "-" << binning.ptHFEfficiencyBinEdges_detector[iHFBin+1] << std::endl;
    }
    
    // Defining cuts
    const double jetRadius = 0.4;
    const double MCPetaCut = 0.9 - jetRadius; // on particle level jet
    const double MCDetaCut = MCPetaCut; // on detector level jet
    const double MCPyCut = 0.8; // on particle level D0
    const double MCDyCut = 0.8; // on detector level D0
    const double MCPDeltaRcut = binning.deltaRBinEdges_particle[binning.deltaRBinEdges_particle.size() - 1]; // on particle level delta R
    const double MCDDeltaRcut = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1]; // on detector level delta R
    double hfptMin = binning.ptHFEfficiencyBinEdges_detector[0]; //ptHFBinEdges[0] - should start from 0 or from the lowest pT,HF value?
    double hfptMax = binning.ptHFEfficiencyBinEdges_detector[binning.ptHFEfficiencyBinEdges_detector.size() - 1];
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
    float MCPhfPt, MCPhfEta, MCPhfPhi, MCPhfMass, MCPhfY;
    bool MCPhfprompt;
    int MCPjetNConst;
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
    MCPhfMass = 1.86484; // D0 rest mass in GeV/c^2 (hard-coded)
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

        // Generator level selection cuts
        bool genJetPtRange = (MCPjetPt >= jetptMin) && (MCPjetPt < jetptMax);//, remove upper bound for particle level
        bool genHfPtRange = ((MCPhfPt >= MCPHfPtMincut) && (MCPhfPt < MCPHfPtMaxcut)); // remove entirely for particle level
        bool genDeltaRRange = ((MCPaxisDistance >= binning.deltaRBinEdges_particle[0]) && (MCPaxisDistance < MCPDeltaRcut));
        bool genLevelRange = (abs(MCPjetEta) < MCPetaCut) && (abs(MCPhfY) < MCPyCut) && genJetPtRange && genDeltaRRange && genHfPtRange;
        // Reconstruction level selection cuts
        double maxBkgProb = GetBkgProbabilityCut(MCDhfPt, binning.bdtPtCuts);
        bool passBDTcut = (MCDhfMlScore0 < maxBkgProb);
        bool recoJetPtRange = (MCDjetPt >= jetptMin) && (MCDjetPt < jetptMax);
        bool recoHfPtRange = ((MCDhfPt >= MCDHfPtMincut) && (MCDhfPt < MCDHfPtMaxcut));
        bool recoDeltaRRange = ((MCDaxisDistance >= binning.deltaRBinEdges_detector[0]) && (MCDaxisDistance < MCDDeltaRcut));
        bool recoLevelRange = (abs(MCDjetEta) < MCDetaCut) && (abs(MCDhfY) < MCDyCut) && recoJetPtRange && recoHfPtRange && recoDeltaRRange; // separated passBDTcut for checking
        bool MCDhfmatch = (MCDjetNConst != -2) ? true : false;
        bool isRealD0 = isTrueSignal(MCDhfMatchedFrom, MCDhfSelectedAs);
        if (genLevelRange) {
            // Matched entries must be of real D0s, not reflections or combinatorial background
            if (MCDhfmatch) {

                if (isRealD0) {
                    // Fill particle level entry
                    dataContainer.hMcpPt[0]->Fill(MCPhfPt);// Fill inclusive histogram
                    if (MCPhfprompt) {
                        dataContainer.hMcpPt[1]->Fill(MCPhfPt);// fill prompt efficiency histogram
                        dataContainer.hSelEff_run2style2d[0]->Fill(MCPjetPt, MCPhfPt);
                        // Run 3 style prompt particle level efficiency histogram
                        dataContainer.hDenominator_prompt_part_run3style->Fill(MCPhfPt);
                    } else{
                        dataContainer.hMcpPt[2]->Fill(MCPhfPt);// fill non-prompt efficiency histogram
                        dataContainer.hDenominator_nonprompt_part_run3style->Fill(MCPhfPt);
                    }
                    // Fill detector level entry
                    if (recoLevelRange && passBDTcut) {
                        // Fill inclusive histogram
                        dataContainer.hMcdPt[0]->Fill(MCDhfPt);
                        // fill prompt efficiency histogram
                        if (MCDhfprompt) {
                            dataContainer.hMcdPt[1]->Fill(MCDhfPt);
                            dataContainer.hSelEff_run2style2d[1]->Fill(MCDjetPt, MCDhfPt);
                            // Run 3 style prompt particle level efficiency histogram
                            dataContainer.hNumerator_prompt_part_run3style->Fill(MCPhfPt);
                        } else{
                            // fill non-prompt efficiency histogram
                            dataContainer.hMcdPt[2]->Fill(MCDhfPt);
                            dataContainer.hNumerator_nonprompt_part_run3style->Fill(MCPhfPt);
                        }
                    } else if (recoLevelRange && !passBDTcut) {
                        dataContainer.hBDTBackgroundScore->Fill(MCDhfMlScore0);
                    }
                }
            } else {// fill particle level even if not matched
                dataContainer.hMcpPt[0]->Fill(MCPhfPt);// Fill inclusive histogram
                if (MCPhfprompt) {
                    dataContainer.hMcpPt[1]->Fill(MCPhfPt);// fill prompt efficiency histogram
                    dataContainer.hSelEff_run2style2d[0]->Fill(MCPjetPt, MCPhfPt);
                    // Run 3 style prompt particle level efficiency histogram
                    dataContainer.hDenominator_prompt_part_run3style->Fill(MCPhfPt);
                } else{
                    dataContainer.hMcpPt[2]->Fill(MCPhfPt);// fill non-prompt efficiency histogram
                    dataContainer.hDenominator_nonprompt_part_run3style->Fill(MCPhfPt);
                }
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
        TH1D* efficiencyHist = static_cast<TH1D*>(mcdHist->Clone(Form("efficiency_run2style_%s", efficiencyNames[iEff]))); // static needed to explicitly cast the cloned histogram to TH1D. ensures that the function returns a TH1D*, matching the return type.
        efficiencyHist->Divide(mcpHist);
        efficiencyHist->SetMarkerStyle(kFullCircle);
        efficiencyHist->SetMarkerColor(30+iEff*10); // 30 = not so bright green
        efficiencyHist->SetLineColor(30+iEff*10);
        efficiencyHist->SetTitle(";#it{p}_{T,D^{0}}^{truth};Efficiency#times Acceptance");
        efficiencyHist->SetStats(0);
        
        // Store efficiency distribution in struct
        dataContainer.hEfficiencies.push_back(efficiencyHist);

    }

    // 2D efficiency calculation
    TH2D* hEfficiencyHist2D = static_cast<TH2D*>(dataContainer.hSelEff_run2style2d[1]->Clone(Form("efficiency_prompt_2d")));
    hEfficiencyHist2D->Divide(dataContainer.hSelEff_run2style2d[0]);
    hEfficiencyHist2D->SetTitle("2D prompt run 2 style efficiency;#it{p}_{T,jet}^{truth};#it{p}_{T,D^{0}}^{truth};Efficiency#times Acceptance");
    dataContainer.hSelEff_run2style2d.emplace_back(hEfficiencyHist2D);

    // Run 3 style particle level prompt selection efficiency calculation
    dataContainer.hEfficiency_prompt_part_run3style = static_cast<TH1D*>(dataContainer.hNumerator_prompt_part_run3style->Clone("efficiency_prompt_run3style_particleLevel"));
    dataContainer.hEfficiency_prompt_part_run3style->Divide(dataContainer.hDenominator_prompt_part_run3style);

    // Run 3 style non-prompt particle level selection efficiency calculation
    dataContainer.hEfficiency_nonprompt_part_run3style = static_cast<TH1D*>(dataContainer.hNumerator_nonprompt_part_run3style->Clone("efficiency_nonprompt_run3style_particleLevel"));
    dataContainer.hEfficiency_nonprompt_part_run3style->Divide(dataContainer.hDenominator_nonprompt_part_run3style);

    std::cout << "Efficiency calculated.\n";
}

void PlotHistograms(const EfficiencyData& dataContainer, double jetptMin, double jetptMax, const BinningStruct& binning) {
     // Create a TLatex object to display text on the canvas

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
    legendPt[3]->AddEntry(dataContainer.hMcdPt[1],"Prompt D^{0} reconstructed", "lpe");
    legendPt[3]->AddEntry(dataContainer.hMcdPt[2],"Non-prompt D^{0} reconstructed", "lpe");

    
    //
    // 1 - plot pT,HF of generator and reconstruction level
    //
    TCanvas* cEff_1 = new TCanvas("cEff_1","pT,HF spectra");
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
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Prompt D^{0}'s");
    cEff_1->cd(3);
    dataContainer.hMcpPt[2]->Draw();
    dataContainer.hMcdPt[2]->Draw("same");
    legendPt[2]->Draw();
    statBoxPos = gPad->GetUxmax();
    latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < #it{p}_{T,jet} < %.0f GeV/c",jetptMin,jetptMax));
    latex->DrawLatex(statBoxPos-0.35, 0.70, "Non-prompt D^{0}'s");
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
    legendEff->AddEntry(dataContainer.hEfficiencies[1],"Prompt D^{0}", "lpe");
    legendEff->AddEntry(dataContainer.hEfficiencies[2],"Non-prompt D^{0}", "lpe");
    for (size_t iEff = 0; iEff < dataContainer.hEfficiencies.size(); iEff++) {
        // if the current is the first plot in the canvas
        if (iEff == 0) {
            double maxHeight = 1.2 * std::max(std::max(dataContainer.hEfficiencies[iEff]->GetMaximum(0),dataContainer.hEfficiencies[iEff]->GetMaximum(1)),dataContainer.hEfficiencies[iEff]->GetMaximum(2));
            dataContainer.hEfficiencies[iEff]->SetMinimum(0);
            dataContainer.hEfficiencies[iEff]->SetMaximum(maxHeight);
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
    // 3 - BDT investigation
    //
    TCanvas* cInvestigation = new TCanvas("cInvestigation","Investigation");
    cInvestigation->SetCanvasSize(1800,1000);
    cInvestigation->cd();
    dataContainer.hBDTBackgroundScore->Draw();
    cInvestigation->SetCanvasSize(1800,1000);

    TCanvas* cEfficiency2D = new TCanvas("cEfficiency2D","2D prompt D0s efficiency");
    cEfficiency2D->Divide(2,2);
    cEfficiency2D->cd(1);
    dataContainer.hSelEff_run2style2d[0]->Draw("colz");
    cEfficiency2D->cd(2);
    dataContainer.hSelEff_run2style2d[1]->Draw("colz");
    cEfficiency2D->cd(3);
    dataContainer.hSelEff_run2style2d[2]->Draw("colz");
    cEfficiency2D->cd(4);
    dataContainer.hSelEff_run2style2d[2]->Draw("text");

    //
    // 4 - Compare the three types of efficiencies calculated
    //
    TLegend* legendEffTypes = new TLegend(0.65,0.45,0.8,0.58);
    TCanvas* cEfficiencyTypes = new TCanvas("cEfficiencyTypes","Efficiency types comparison",1800,1000);
    cEfficiencyTypes->cd();
    // Run 3 style particle level
    dataContainer.hEfficiency_prompt_part_run3style->Draw();
    legendEffTypes->AddEntry(dataContainer.hEfficiency_prompt_part_run3style,"Run 3 style prompt particle level efficiency", "lpe");
    // Run 2 style
    dataContainer.hEfficiencies[1]->Draw("same");
    legendEffTypes->AddEntry(dataContainer.hEfficiencies[1],"Run 2 style prompt efficiency", "lpe");
    legendEffTypes->Draw();

    //
    // Storing images
    //
    TString sEmmaBins;
    if (binning.useEmmaYeatsBins) {
        sEmmaBins = "EmmaYeatsBins";
    } else {
        sEmmaBins = "";
    }
    TString imagePath = "../Images/2-Efficiency/Run2Style/" + sEmmaBins + "/" + binning.dataPeriod + "/";
    
    cEff_1->Update();
    cEff_1->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cEff_1->SaveAs(imagePath + "Efficiency_dNdpT_run2_style_" + sEmmaBins + ".png");
    cEff_2->Update();
    cEff_2->SetWindowSize(1920, 1080);  // Optional: Match window size for viewing.
    cEff_2->SaveAs(imagePath + "Efficiency_acceptance_run2_style_" + sEmmaBins + ".png");

    //
    // Storing in a single pdf file
    //
    cEff_1->Print(Form(imagePath + "pT_efficiency_%.0f_to_%.0fGeV_run2_style.pdf(",jetptMin,jetptMax));
    cEff_2->Print(Form(imagePath + "pT_efficiency_%.0f_to_%.0fGeV_run2_style.pdf",jetptMin,jetptMax));
    cEfficiency2D->Print(Form(imagePath + "pT_efficiency_%.0f_to_%.0fGeV_run2_style.pdf",jetptMin,jetptMax));
    cEfficiencyTypes->Print(Form(imagePath + "pT_efficiency_%.0f_to_%.0fGeV_run2_style.pdf)",jetptMin,jetptMax));

    std::cout << "Histograms plotted.\n";
}

void SaveData(const EfficiencyData& dataContainer, const BinningStruct& binning, double& jetptMin, double& jetptMax){
    // Open output file
    TFile* fOutput = new TFile(Form("run2_style_efficiency_%d_to_%d_jetpt_" + binning.dataPeriod + ".root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"recreate");
    // Loop over efficiency histograms
    for (size_t iHisto = 0; iHisto < dataContainer.hEfficiencies.size(); iHisto++) {
        // store each histogram in file
        dataContainer.hEfficiencies[iHisto]->Write();
    }

    // Run 3 particle level prompt selection efficiency
    dataContainer.hEfficiency_prompt_part_run3style->Write();
    dataContainer.hEfficiency_nonprompt_part_run3style->Write();
    

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
    TFile* fBinning = new TFile(Form("../1-SignalTreatment/Reflections/binningInfo.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    // Load binning from background subtracted file
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];

    // Opening files
    TFile* fSimulated = new TFile("../" + binning.inputMC.second + "/AO2D_mergedDFs.root","read","read");

    EfficiencyData dataContainer = createHistograms(binning); // pT histograms

    // Fill histograms
    fillHistograms(fSimulated, dataContainer, binning, jetptMin, jetptMax);

    // Calculate efficiency distribution
    CalculateEfficiency(dataContainer);

    // Plot the efficiency histogram and further corrected histograms
    PlotHistograms(dataContainer, jetptMin, jetptMax, binning);

    // Save corrected distributions to file
    SaveData(dataContainer, binning, jetptMin, jetptMax);

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    std::time_t now = std::time(nullptr);
    std::cout << "Finished at: " << std::ctime(&now);
}

int main(){
    Efficiency_run2_style();
    return 0;
}
