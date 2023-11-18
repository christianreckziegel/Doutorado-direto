# FastJet workflow stage

The main analysis and variations were implemented:
- **myAnalysis.cc**: main implementation for identification of heavy flavour decay candidates, checking of invariant mass, jet tagging and analysis.
- **myAnalysis2.cc**: same as main implementation with dedicated function for processing each event and making use of pointers for more efficient memory management.
- **myAnalysis3.cc**: initial analysis without user information associated to PseudoJets.

### Memory leak
All options seems to have memory leak issues which were investigated but not solved. It appears to be related to my Cling library but not specifically the analysis implemented.
The way it's implemented it can process 10k events but not 100k without consuming all my computer memory (32 GB).
For that reason, and as an imediate solution, the analysis was segmented in steps of 10k events and automatized in the workflow in /BigData_analysis.

