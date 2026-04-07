/*
 *
 *
 * Macro for performing unfolding on prompt D0s measured 
 * data and applying correction to FeedDownSubtraction.C
 * resulting distributions.
 * 
 * 
 * 
 * 
**/

#include "../commonUtilities.h"
#include "sidebandClosure.h"
#include "efficiencyClosure.h"
#include "unfoldingClosure.h"

using namespace std;

struct ClosureTestData2 {

    // Input objects: pT,jet vs DeltaR vs pT,D0
    TH3D* hInputParticle = nullptr;
    TH3D* hInputDetector = nullptr;

    // Correction objects
    RooUnfoldResponse* response;                                        // response matrix
    std::vector<TH2D*> hKineEffParticle = {nullptr, nullptr, nullptr};  // particle level kinematic efficiency -> 0 = numerator, 1 = denominator, 2 = efficiency
    std::vector<TH2D*> hKineEffDetector = {nullptr, nullptr, nullptr};  // detector level kinematic efficiency -> 0 = numerator, 1 = denominator, 2 = efficiency
    
    // Unfolding objects
    std::vector<RooUnfoldBayes*> unfold;                                // unfolding objects, there are iterationNumber unfolding objects
    std::vector<TH2D*> hUnfolded;                                       // unfolded 2D histogram (jet pT vs DeltaR), there are iterationNumber unfolding objects
};

// 2 - Sideband subtraction + efficiency correction + unfolding closure test
void SecondClosureTest(){
    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution
    
    // Number of unfolding procedure iterations
    int iterationNumber = 4;

    // Load binning from reflections file
    TFile* fBinning = new TFile(Form("../1-SignalTreatment/Reflections/binningInfo.root"),"read");
    if (!fBinning || fBinning->IsZombie()) {
        std::cerr << "Error: Unable to open the first ROOT binning info file." << std::endl;
    }
    BinningStruct binning = retrieveBinningFromFile(fBinning);
    double jetptMin = binning.ptjetBinEdges_detector[0];
    double jetptMax = binning.ptjetBinEdges_detector[binning.ptjetBinEdges_detector.size() - 1];
    double minDeltaR = binning.deltaRBinEdges_detector[0];
    double maxDeltaR = binning.deltaRBinEdges_detector[binning.deltaRBinEdges_detector.size() - 1];
    double hfptMin = binning.ptHFBinEdges_detector[0]; //ptHFBinEdges[0] - should start from 0 or from the lowest pT,D value?
    double hfptMax = binning.ptHFBinEdges_detector[binning.ptHFBinEdges_detector.size() - 1];


    // Opening files
    TFile* fClosureInput = new TFile("mc_closure_input_data.root","read");
    if (!fClosureInput || fClosureInput->IsZombie()) {
        std::cerr << "Error: Unable to open 1st closure input ROOT file." << std::endl;
    }

    // ----1: Perform side-band subtraction
    // Example: selecting a model
    FitModelType modelToUse = FitModelType::SignalReflectionsOnly;
    TH3D* hBackgroundSubtracted = SidebandClosure(fClosureInput, binning, modelToUse);

    // ----2: Perform efficiency correction
    //EfficiencyData efficiencyDatacontainer = EfficiencyClosure(fClosureInput, hBackgroundSubtracted, binningStruct, bdtPtCuts);
    //std::vector<TH1D*> hSelEff_run3style = {efficiencyDatacontainer.hSelectionEfficiency.first,efficiencyDatacontainer.hSelectionEfficiency.second};
    //TH2D* hEfficiencyCorrected = efficiencyDatacontainer.hEfficiencyCorrected.second;

    // ----3: Perform unfolding
    //UnfoldData unfoldDataContainer = UnfoldingClosure(fClosureInput, hSelEff_run3style, hEfficiencyCorrected, binningStruct, bdtPtCuts);

    // ----4: Compare input (MC particle level) with output (background subtracted, efficiency corrected, unfolded) distributions
    //TH2D* inputMCP = CompareClosureTest(fClosureInput, unfoldDataContainer.hUnfoldedKinCorrected, binningStruct);

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    SecondClosureTest();
    return 0;
}
