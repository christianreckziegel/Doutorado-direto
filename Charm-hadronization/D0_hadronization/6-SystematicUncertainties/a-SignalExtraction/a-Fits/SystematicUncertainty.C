



void SystematicUncertainty() {
    // List of your files
    std::vector<std::string> files = {
        "unfolding_5_to_50_jetpt_1_middle_fit_range.root",
        "unfolding_5_to_50_jetpt_2_full_fit_range.root",
        "unfolding_5_to_50_jetpt_3_minor_fit_range.root",
        "unfolding_5_to_50_jetpt_4_left_fit_range.root",
        "unfolding_5_to_50_jetpt_5_right_fit_range.root",
    };
    std::vector<std::string> legendTexts = {
        "Middle range: | --- |",
        "Full range: |-----|",
        "Minor range: |  -  |",
        "Left delocated range: |---  |",
        "Right delocated range: |  ---|"
    }; // 1.6 increase to 1.64, full ranges, change a little bit in the borders

    TCanvas *c = new TCanvas("c","Projection comparison",800,600);
    c->cd();

    int color = 1; // start color index
    TLegend *leg = new TLegend(0.65,0.65,0.88,0.88);

    for (size_t i=0; i < files.size(); i++) {
        TFile *f = TFile::Open(files[i].c_str(),"READ");
        if (!f || f->IsZombie()) {
            std::cout << "Could not open file " << files[i] << std::endl;
            continue;
        }

        // Go inside "Unfolded" directory
        f->cd("Unfolded");
        TIter next(gDirectory->GetListOfKeys());
        TKey *key;
        TH2D *h2 = nullptr;

        int counter = 0;
        while ((key = (TKey*)next())) {
            counter++;
            if (counter == 8) { // <-- the 8th histogram
                h2 = (TH2D*)key->ReadObj();
                break;
            }
        }

        if (!h2) {
            std::cout << "8th histogram not found in " << files[i] << std::endl;
            continue;
        }
        //h2->Sumw2(); // Ensure Sumw2 is called for proper error handling

        // Project Y
        TH1D *hy = h2->ProjectionY(Form("py_%zu",i));
        if (!hy || hy->GetEntries() == 0) {
            std::cout << "⚠️ Empty projection for file: " << files[i] << std::endl;
            continue;
        }
        if (!std::isfinite(hy->GetMaximum())) {
            std::cout << "⚠️ Invalid values in histogram from file: " << files[i] << std::endl;
            continue;
        }
        hy->SetTitle("Uncertainty due to fit range choice & model");
        //hy->Sumw2(); // Ensure Sumw2 is called for proper error handling

        // Style
        auto colors = std::vector<int>{
            kRed+1, kBlue+1, kGreen+2, kMagenta+1, kCyan+2,
            kOrange+7, kViolet+1, kPink+9, kTeal-5, kAzure+7,
            kSpring+9, kGray+2, kBlack
        };
        hy->SetLineColor(colors[i % colors.size()]);
        hy->SetLineWidth(2);
        //hy->SetLineStyle(1 + (i % 4)); // optional: varying style

        if (i==0) hy->Draw("");  // first one creates the axes
        else      hy->Draw("same");

        leg->AddEntry(hy, legendTexts[i].c_str(), "l");

        //color++;
        //if (color==5) color=6; // skip yellow
    }

    leg->Draw();

}

int main() {
    SystematicUncertainty();
    return 0;
}