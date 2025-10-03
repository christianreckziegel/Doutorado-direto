



void SystematicUncertainty() {
    // List of your files
    std::vector<std::string> files = {
        "unfolding_5_to_50_jetpt_1.6sigma_signal.root",
        "unfolding_5_to_50_jetpt_1.7sigma_signal.root",
        "unfolding_5_to_50_jetpt_1.8sigma_signal.root",
        "unfolding_5_to_50_jetpt_1.9sigma_signal.root",
        "unfolding_5_to_50_jetpt_2sigma_signal.root",
        "unfolding_5_to_50_jetpt_2.1sigma_signal.root",
        "unfolding_5_to_50_jetpt_2.4sigma_signal.root",
        "unfolding_5_to_50_jetpt_2.5sigma_signal.root"
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
        //hy->Sumw2(); // Ensure Sumw2 is called for proper error handling

        // Style
        hy->SetLineColor(color);
        hy->SetLineWidth(2);

        if (i==0) hy->Draw("");  // first one creates the axes
        else      hy->Draw("same");

        leg->AddEntry(hy, files[i].c_str(), "l");

        color++;
        if (color==5) color=6; // skip yellow
    }

    leg->Draw();

}

int main() {
    SystematicUncertainty();
    return 0;
}