



void SystematicUncertainty() {
    // List of your files
    std::vector<std::string> files = {
        "unfolding_5_to_50_jetpt_1.6sigma_signal.root",
        "unfolding_5_to_50_jetpt_1.7sigma_signal.root",
        "unfolding_5_to_50_jetpt_1.8sigma_signal.root",
        "unfolding_5_to_50_jetpt_1.9sigma_signal.root",
        "unfolding_5_to_50_jetpt_2.1sigma_signal.root",
        "unfolding_5_to_50_jetpt_2.4sigma_signal.root",
        "unfolding_5_to_50_jetpt_2.5sigma_signal.root",
        "unfolding_5_to_50_jetpt_3-5sigma_sideband.root",
        "unfolding_5_to_50_jetpt_3-6sigma_sideband.root",
        "unfolding_5_to_50_jetpt_4-5sigma_sideband.root",
        "unfolding_5_to_50_jetpt_4-6sigma_sideband.root",
        "unfolding_5_to_50_jetpt_4-7sigma_sideband.root",
        "unfolding_5_to_50_jetpt_default.root"
    };
    std::vector<std::string> legendTexts = {
        "signal #in #bar{m} #pm 1.6#sigma",
        "signal #in #bar{m} #pm 1.7#sigma",
        "signal #in #bar{m} #pm 1.8#sigma",
        "signal #in #bar{m} #pm 1.9#sigma",
        "signal #in #bar{m} #pm 2.1#sigma",
        "signal #in #bar{m} #pm 2.4#sigma",
        "signal #in #bar{m} #pm 2.5#sigma",
        "side-band #in #bar{m} #pm [3#sigma;5#sigma]",
        "side-band #in #bar{m} #pm [3#sigma;6#sigma]",
        "side-band #in #bar{m} #pm [4#sigma;5#sigma]",
        "side-band #in #bar{m} #pm [3#sigma;6#sigma]",
        "side-band #in #bar{m} #pm [4#sigma;7#sigma]",
        "Default: #Delta signal = 2#sigma; #Delta side-band = 4#sigma"
    };

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
        hy->SetTitle("Uncertainty due to signal & sideband region definitions");
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