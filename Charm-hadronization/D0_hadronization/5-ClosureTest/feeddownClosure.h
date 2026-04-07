/*
 * Macro for performing feed-down subtraction to third closure test
 * 
 * 
 * 
 * Author: Christian Reckziegel
**/


using namespace std;

struct FeedDownData {
    // Kinematic efficiencies: first = particle level, second = detector level
    std::pair<TH2D*, TH2D*> kinematicEfficiency;

    // B contribution distributions
    TH3D* hPowheg3DParticle;
    TH2D* hBPowheg2DParticle;
    TH2D* hBPowheg2DDetector

    // D0 feed-down subtractd distribution
    TH2D* hFeddown;
}
std::pair<TH2D*, TH2D*> calculateKinematicEfficiency(TFile* fClosureInputMatched, FeedDownData& dataContainer) {
    
    // Particle and detector level kinematic efficiency
    std::pair<TH2D*, TH2D*> kinematicEfficiency;

    // Create DeltaR vs pT,jet matched dstributions

    // Fill distributions

    // Calculate kinematic efficiencies

}

TH2D* createBContribution(TFile* fPowheg, FeedDownData& dataContainer, const std::vector<TH1D*>& hSelEff_run3style) {
    // Create 3D B distribution

    // Apply selection efficiency correction

    // Project and obtain 2D B distribution

    // Apply particle level kinematic efficiency correction

}

TH2D* smearBContribution(FeedDownData& dataContainer) {
    // Fold 2D B distribution

    // Apply detector level kinematic efficiency correction

    // Scale by POWHEG and data luminosities


}

TH2D* subtractBs(TH2D* hEfficiencyCorrected) {
    
    // Obtain efficiency corrected distribution

    // Subtract B contribution
}

TH2D* FeedDownClosure(TFile* fClosureInputNonMatched, TFile* fClosureInputMatched, TFile* fPowheg, const std::vector<TH1D*>& hSelEff_run3style, TH2D* hEfficiencyCorrected, const BinningStruct& binningStruct, const std::vector<std::pair<double, double>>& bdtPtCuts) {

    // Testing range of pT,jet bins
    //ptjetBinEdges_detector = {15., 30.};

    // Create data container to hold all distributions
    FeedDownData dataContainer;

    // Calculate kinematic efficiency
    std::pair<TH2D*, TH2D*> kinematicEfficiency = calculateKinematicEfficiency(fClosureInputMatched, dataContainer);

    // Fill and correct kinematically non-prompt contribution to detector level
    TH2D* hBPowheg2DParticle = createBContribution(fPowheg, dataContainer, hSelEff_run3style);

    // Smear, correct kinematically and scale B contribution
    TH2D* hBPowheg2DDetector = smearBContribution(dataContainer);

    // Subtract B contribution from measured D0 candidates data
    TH2D* hFeddown = subtractBs(hEfficiencyCorrected);


    TH2D* hFeddown = hEfficiencyCorrected;

    return hFeddown;
}
