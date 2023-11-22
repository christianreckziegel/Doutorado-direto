
void myPlots(){
    // hard coded variables
    TString sNumEvents = "1m";// 1k, 10k, 100k, 1m

    //access file
    TFile* file1 = new TFile("AnalysisResults_"+sNumEvents+"_events_option_b_Merged.root");
    TH1F* hMergedJetsPt = (TH1F*)file1->Get("hMergedJetsPt");
    TH1F* hMergedHFCandidateJetsPt = (TH1F*)file1->Get("hMergedHFCandidateJetsPt");
    
    // Normalizing inclusive jets pT histogram
    Double_t integralInclusive = hMergedJetsPt->Integral();
    if (integralInclusive > 0.0) {
        hMergedJetsPt->Scale(1.0 / integralInclusive);
        hMergedJetsPt->SetTitle("Normalized jets p_{T}");
        hMergedJetsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hMergedJetsPt->GetYaxis()->SetTitle("Relative frequency");
        hMergedJetsPt->Rebin(8); //divide current size by 8
        hMergedJetsPt->SetStats(0);
        hMergedJetsPt->GetXaxis()->SetRangeUser(0, 20);
        hMergedJetsPt->SetLineColor(kBlack);
        hMergedJetsPt->SetMarkerColor(kBlack);
        hMergedJetsPt->SetMarkerStyle(kFullCircle);
        hMergedJetsPt->SetLineWidth(2);
        
    }

    // Normalizing HF jet candidate pT histogram
    Double_t integralHF = hMergedHFCandidateJetsPt->Integral();
    if (integralHF > 0.0) {
        hMergedHFCandidateJetsPt->Scale(1.0 / integralHF);
        hMergedHFCandidateJetsPt->Rebin(8);
        hMergedHFCandidateJetsPt->SetStats(0);
        hMergedHFCandidateJetsPt->GetXaxis()->SetRangeUser(0, 20);
        hMergedHFCandidateJetsPt->SetLineColor(kRed);
        hMergedHFCandidateJetsPt->SetMarkerColor(kRed);
        hMergedHFCandidateJetsPt->SetMarkerStyle(kCircle);
        hMergedHFCandidateJetsPt->SetLineWidth(2);
    }

    // set statistic box off
    //gStyle->SetOptStat(0);

    TCanvas* cJetsPt = new TCanvas("cJetsPt", "Normalized jets' pT");
    cJetsPt->cd();
    // Draw the histograms on the same canvas
    hMergedJetsPt->Draw();
    hMergedHFCandidateJetsPt->Draw("same");

    double meanPtInclusive = hMergedJetsPt->GetMean();
    double meanPtCand = hMergedHFCandidateJetsPt->GetMean();

    // Add a legend
    TLegend* legIncHF = new TLegend(0.6, 0.8, 0.9, 0.9);
    //legend->AddEntry(hMergedJetsPt, "Inclusive jets", "l");
    //legend->AddEntry(hMergedHFCandidateJetsPt, "HF candidate jets", "l");
    legIncHF->AddEntry(hMergedJetsPt, Form("Inclusive jets, #bar p_{T}^{inc} = %.2f GeV/c", meanPtInclusive), "l");
    legIncHF->AddEntry(hMergedHFCandidateJetsPt, Form("HF candidate jets, #bar p_{T}^{HF cand} = %.2f GeV/c", meanPtCand), "l");
    legIncHF->Draw();

    // Add text labels
    TLatex* textIncHF = new TLatex();
    textIncHF->SetNDC();
    textIncHF->SetTextSize(0.03);
    textIncHF->SetTextFont(42); // Font 42 is the default
    textIncHF->DrawLatex(0.15, 0.85, "PYTHIA 8, pp, #sqrt{s} = 13 TeV");
    textIncHF->DrawLatex(0.15, 0.80, "Anti-k_{T}, R = 0.4");

    //
    //
    //
    TCanvas* cPromptJetsPt = new TCanvas("cPromptJetsPt", "Normalized prompt and non-prompt jets pT");
    cPromptJetsPt->cd();
    TH1F* hMergedPromptD0_JetsPt = (TH1F*)file1->Get("hMergedPromptD0_JetsPt");
    TH1F* hMergedNonPromptD0_JetsPt = (TH1F*)file1->Get("hMergedNonPromptD0_JetsPt");
    // Normalizing prompt D0 tagged jets pT histogram
    Double_t integralPrompt = hMergedPromptD0_JetsPt->Integral();
    if (integralPrompt > 0.0) {
        //hMergedPromptD0_JetsPt->Scale(1.0 / integralPrompt);
        hMergedPromptD0_JetsPt->SetTitle("Prompt and non-prompt true D^{0} jets p_{T}");
        hMergedPromptD0_JetsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hMergedPromptD0_JetsPt->GetYaxis()->SetTitle("Counts");
        hMergedPromptD0_JetsPt->Rebin(8); //divide current size by 8
        hMergedPromptD0_JetsPt->SetStats(0);
        hMergedPromptD0_JetsPt->GetXaxis()->SetRangeUser(0, 52);
        hMergedPromptD0_JetsPt->SetLineColor(kBlue);
        hMergedPromptD0_JetsPt->SetMarkerColor(kBlue);
        hMergedPromptD0_JetsPt->SetMarkerStyle(kFullCircle);
        hMergedPromptD0_JetsPt->SetLineWidth(2);
        
    }
    // Normalizing non-prompt D0 tagged jets pT histogram
    Double_t integralNonPrompt = hMergedNonPromptD0_JetsPt->Integral();
    if (integralNonPrompt > 0.0) {
        //hMergedNonPromptD0_JetsPt->Scale(1.0 / integralNonPrompt);
        hMergedNonPromptD0_JetsPt->Rebin(8); //divide current size by 8
        hMergedNonPromptD0_JetsPt->SetStats(0);
        hMergedNonPromptD0_JetsPt->GetXaxis()->SetRangeUser(0, 52);
        hMergedNonPromptD0_JetsPt->SetLineColor(kGreen+2);
        hMergedNonPromptD0_JetsPt->SetMarkerColor(kGreen);
        hMergedNonPromptD0_JetsPt->SetMarkerStyle(kCircle);
        hMergedNonPromptD0_JetsPt->SetLineWidth(2);
        
    }
    hMergedPromptD0_JetsPt->Draw();
    hMergedNonPromptD0_JetsPt->Draw("same");
    // Add a legend
    TLegend* legPrompt = new TLegend(0.6, 0.8, 0.9, 0.9);
    //legend->AddEntry(hMergedPromptD0_JetsPt, "Inclusive jets", "l");
    //legend->AddEntry(hMergedNonPromptD0_JetsPt, "HF candidate jets", "l");
    double meanPtPrompt = hMergedPromptD0_JetsPt->GetMean();
    double meanPtNonPrompt = hMergedNonPromptD0_JetsPt->GetMean();
    legPrompt->AddEntry(hMergedPromptD0_JetsPt, Form("Prompt D^{0} jets, #bar p_{T}^{prompt} = %.2f GeV/c", meanPtPrompt), "l");
    legPrompt->AddEntry(hMergedNonPromptD0_JetsPt, Form("Non-prompt D^{0} jets, #bar p_{T}^{non-prompt} = %.2f GeV/c", meanPtNonPrompt), "l");
    legPrompt->Draw();
    // Add text labels
    TLatex* textPrompt = new TLatex();
    textPrompt->SetNDC();
    textPrompt->SetTextSize(0.03);
    textPrompt->SetTextFont(42); // Font 42 is the default
    double nonPromptFraction = 100*integralNonPrompt/(integralNonPrompt+integralPrompt);
    textPrompt->DrawLatex(0.15, 0.85, Form("Non-prompt fraction = %.2f %%",nonPromptFraction));
    textPrompt->DrawLatex(0.15, 0.80, "Non-prompt = B^{0}, B^{+}, B^{-}, B_{s}^{0}, B_{c}^{+} and B_{c}^{-} meson decays");

}

int main(){
    myPlots();
    return 0;
}