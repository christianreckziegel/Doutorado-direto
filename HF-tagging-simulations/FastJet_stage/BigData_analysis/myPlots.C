
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
        //hMergedJetsPt->SetBins(500, hMergedJetsPt->GetXaxis()->GetXmin(), hMergedJetsPt->GetXaxis()->GetXmax());
    }

    // Normalizing HF jet candidate pT histogram
    Double_t integralHF = hMergedHFCandidateJetsPt->Integral();
    if (integralHF > 0.0) {
        hMergedHFCandidateJetsPt->Scale(1.0 / integralHF);
        //hMergedHFCandidateJetsPt->SetBins(500, hMergedHFCandidateJetsPt->GetXaxis()->GetXmin(), hMergedHFCandidateJetsPt->GetXaxis()->GetXmax()); // 500 bins
    }

    TCanvas* cJetsPt = new TCanvas("cJetsPt", "Normalized jets' pT");
    cJetsPt->cd();
    // Draw the histograms on the same canvas
    hMergedJetsPt->SetLineColor(kBlue);
    hMergedJetsPt->Draw();
    hMergedHFCandidateJetsPt->SetLineColor(kRed);
    hMergedHFCandidateJetsPt->Draw("same");

    // Add a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.8, 0.75);
    legend->AddEntry(hMergedJetsPt, "Inclusive jets", "l");
    legend->AddEntry(hMergedHFCandidateJetsPt, "HF candidate jet", "l");
    legend->Draw();
}

int main(){
    myPlots();
    return 0;
}