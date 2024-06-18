/*
 *
 *
 * Macro for plotting result histograms from 
 * BackgroundSubtraction.C and SignalExtraction.C macros
 * 
 * 
 * 
 * 
**/


using namespace std;

int HistogramCounter(TFile* file) {
    // Get the list of keys (i.e., objects) in the file
    TList* keys = file->GetListOfKeys();

    int numHistograms = 0;

    // Loop over all keys in the file
    for (int i = 0; i < keys->GetSize(); ++i) {
        TKey* key = (TKey*)keys->At(i);
        TObject* obj = file->Get(key->GetName());

        // Check if the object is a histogram
        if (obj->IsA()->InheritsFrom(TH1::Class())) {
            numHistograms++;
        }
    }

    return numHistograms;
}

void PlotRatios(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    // D0 mass in GeV/c^2
    double m_0_parameter = 1.86484;
    double sigmaInitial = 0.012;
    // jet pT cuts
    double jetptMin = 5; // GeV
    double jetptMax = 30; // GeV
    // deltaR histogram
    int deltaRbins = 10; // deltaRbins = numberOfPoints, default=100 bins for [0. 1.0]
    double minDeltaR = 0.;
    double maxDeltaR = 0.4;
    // mass histogram
    int massBins = 100; 
    double minMass = 1.67;
    double maxMass = 2.1;
    

    // opening files
    TFile* fSigExt = new TFile(Form("sigExt_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fSigExt || fSigExt->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT file." << std::endl;
    }
    TFile* fBackSub = new TFile(Form("backSub_%d_to_%d_jetpt.root",static_cast<int>(jetptMin),static_cast<int>(jetptMax)),"read");
    if (!fBackSub || fBackSub->IsZombie()) {
        std::cerr << "Error: Unable to open the ROOT file." << std::endl;
    }

    // Creating temporary histograms for managing data
    TH1D* hTempSigExt;
    TH1D* hTempSBackSub;
    TH1D* hRatio;

    // Count number of histograms stored in file
    int numHistograms = HistogramCounter(fSigExt);

    // Creating canvases
    // same pad plots canvas
    TCanvas* cOverPlots = new TCanvas("cOverPlots","Plots one over the other");
    cOverPlots->Divide(3,numHistograms/3);
    // ratio canvas, (signal extracted method)/(background subtracted method)
    TCanvas* cRatio = new TCanvas("cRatio","Ratio plots");
    cRatio->Divide(3,numHistograms/3);

    // Create a TLatex object to display text on the canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(); // Set the coordinates to be normalized device coordinates
    latex->SetTextSize(0.03);


    // Looping through histograms
    for (size_t iHisto = 0; iHisto < numHistograms; iHisto++) {
        // Accessing histograms
        hTempSigExt = (TH1D*)fSigExt->Get(Form("h_sig_extracted_%zu",iHisto));
        if (!hTempSigExt) {
            cout << "Unable to find " << iHisto << " signal extracted histogram.\n";
        }
        hTempSBackSub = (TH1D*)fBackSub->Get(Form("h_back_subtracted_%zu",iHisto));
        if (!hTempSBackSub) {
            cout << "Unable to find " << iHisto << " background subtracted histogram.\n";
        }

        // Creating legend entries
        TLegend* legend = new TLegend(0.7,0.47,0.8,0.57);
        legend->AddEntry(hTempSigExt,"SigExt", "p");
        legend->AddEntry(hTempSBackSub,"SideBand", "p");
        
        // Plot one distribution over the other
        cOverPlots->cd(iHisto+1);
        hTempSigExt->Draw();
        hTempSBackSub->Draw("same");
        legend->Draw();
        double statBoxPos = gPad->GetUxmax(); // Height of the stat box
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax)); // Display jet pT cut applied

        // Create and plot ratio histograms
        //hRatio = hTempSigExt->Clone(Form("Ratio_%zu",iHisto));
        hRatio = dynamic_cast<TH1D*>(hTempSigExt->Clone(Form("Ratio_%zu", iHisto)));
        hRatio->Divide(hTempSBackSub);
        hRatio->SetMarkerStyle(kFullDotLarge);
        hRatio->SetMarkerColor(8); // not that bright green
        hRatio->SetLineColor(kBlack);
        hRatio->GetYaxis()->SetTitle("Ratio = Sig/SB");
        cRatio->cd(iHisto+1);
        hRatio->GetYaxis()->SetRangeUser(0.,3.);
        TLine line1(0,1.,0.45,1.);
        line1.SetLineStyle(2);
        hRatio->Draw();
        line1.DrawClone();
        statBoxPos = gPad->GetUxmax(); // Height of the stat box
        latex->DrawLatex(statBoxPos-0.35, 0.65, Form("%.0f < p_{T,jet} < %.0f GeV/c",jetptMin,jetptMax)); // Display jet pT cut applied

        //delete legend;
    }
    


    

    time(&end); // end instant of program execution
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    PlotRatios();
    return 0;
}



//
// Fitting parameter values
//
// 5-15 GeV/c jet pT initial fit parameters (deltaRbins = 100, massBins = 100)
/*parametersVectors.paramA = {5354.55, 4825.02, 3236.23, 1600.03, 926.955, 1684.15, 737.055, 572.857, 1.};
parametersVectors.paramB = {2.8335, 3.11021, 3.40724, 3.58429, 3.72211, 3.44699, 3.40815+0.5, 3.47049-4.5, 1.};
parametersVectors.paramC = {1.98347e+06, 2.20959e+06, 1.37818e+06, 920285, 1.18674e+06, 502167, 394194, 1.};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { 2*sigmaInitial, sigmaInitial, 2.5*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial, sigmaInitial };*/
// 15-30 GeV/c jet pT initial fit parameters (deltaRbins = 100, massBins = 100)
/*parametersVectors.paramA = {151.495, 160.988, 107.586, 95.7014, 82.6835, 221.485, 207.169, 891.023, 1448.5};
parametersVectors.paramB = {2.65389+0.5, 2.883, 2.96528, 3.19736, 3.23468, 3.3857, 3.29321, 3.21522};
parametersVectors.paramC = {21040.2, 30964.5, 29364.8, 29796.1, 35084.3, 98945.7, 124168, 464545, 673854};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { sigmaInitial, sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial };*/
// 5-30 GeV/c jet pT initial fit parameters (deltaRbins = 100, massBins = 100)
/*parametersVectors.paramA = {5497.6, 4963.91, 3334.83, 1681.7, 1006.26, 1904.17, 944.77, 1460.87, 1448.5};
parametersVectors.paramB = {3.10233, 3.39654, 3.56564, 3.69124, 3.42524, 3.40274, 3.36828, 3.21522, };
parametersVectors.paramC = {1.41624e+06, 2.01454e+06, 2.23897e+06, 1.40804e+06, 955357, 1.28567e+06, 626330, 858739, 673854};
parametersVectors.paramM0 = { m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter, m_0_parameter };
parametersVectors.paramSigma = { 3*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 1.5*sigmaInitial, 1.5*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial, 2*sigmaInitial };*/