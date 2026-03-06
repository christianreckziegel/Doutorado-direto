#include "fastjet/PseudoJet.hh"
using namespace fastjet;


class MyInfo: public PseudoJet::UserInfoBase {
    public:
    MyInfo(int id, int prompt, bool hasD0MesonConstituent, int eventNumber) : pdg(id), isPrompt(prompt), hasD0Constituent(hasD0MesonConstituent), pythiaEventNumber(eventNumber) {} // default value of isPrompt is -1, meaning not a D0 particle
    
    int GetPDG() const { return pdg; }
    int GetPromptness() const { return isPrompt; }
    bool hasD0() const { return hasD0Constituent; } // Check if the PseudoJet has a D0 or D0bar constituent
    int GetPythiaEventNumber() const { return pythiaEventNumber; } // Get the original Pythia event number for reference

    private:
    int pdg; // PDG code of the particle
    int isPrompt; // 1 if prompt, 0 if from b decay, -1 if not a D0 particle
    bool hasD0Constituent; // true if the PseudoJet has a D0 or D0bar constituent, false otherwise
    int pythiaEventNumber; // The event number from the original Pythia event, for reference
};