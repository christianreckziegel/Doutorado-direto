
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
        //hMergedJetsPt->SetBins(250, hMergedJetsPt->GetXaxis()->GetXmin(), hMergedJetsPt->GetXaxis()->GetXmax()); // only use SetBins() function before creating histogram
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
        //hMergedHFCandidateJetsPt->SetBins(250, hMergedHFCandidateJetsPt->GetXaxis()->GetXmin(), hMergedHFCandidateJetsPt->GetXaxis()->GetXmax()); // 500 bins
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
    TLegend* legend = new TLegend(0.6, 0.8, 0.9, 0.9);
    //legend->AddEntry(hMergedJetsPt, "Inclusive jets", "l");
    //legend->AddEntry(hMergedHFCandidateJetsPt, "HF candidate jets", "l");
    legend->AddEntry(hMergedJetsPt, Form("Inclusive jets, #bar p_{T}^{inc} = %.2f GeV/c", meanPtInclusive), "l");
    legend->AddEntry(hMergedHFCandidateJetsPt, Form("HF candidate jets, #bar p_{T}^{HF cand} = %.2f GeV/c", meanPtCand), "l");
    legend->Draw();

    // Add text labels
    TLatex* text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextFont(42); // Font 42 is the default
    text->DrawLatex(0.15, 0.85, "PYTHIA 8, pp, #sqrt{s} = 13 TeV");
    text->DrawLatex(0.15, 0.80, "Anti-k_{T}, R = 0.4");

}

int main(){
    myPlots();
    return 0;
}