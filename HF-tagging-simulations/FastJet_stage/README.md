# FastJet workflow stage

The main analysis and variations were implemented:
- **myAnalysis.cc**: main implementation for identification of heavy flavour decay candidates, checking of invariant mass, jet tagging and analysis. ROOT running inside FastJet.
- **myAnalysis2.cc**: same as main implementation with dedicated function for processing each event and making use of pointers for more efficient memory management. ROOT running inside FastJet.
- **myAnalysis3.cc**: initial analysis without user information associated to PseudoJets. ROOT running inside FastJet.
- **myAnalysis_onRoot.C**: initial attempt of analysis running FastJet inside ROOT.

### Result files
Stored with a name patter as\
"AnalysisResults_<number_of_events>_events_option_<data_storing_option>.root"\
Obs.: more about data storing option of /Pythia_stage.

### Memory leak
All options seems to have memory leak issues which were investigated but not solved. It appears to be related to my Cling library but not specifically the analysis implemented. The output message after running with Valgrind can be found in errorOutAnalysis3.txt for the analysis number 3.\
The way it's implemented it can process 10k events but not 100k without consuming all my computer memory (32 GB).
For that reason, and as an imediate solution, the analysis was segmented in steps of 10k events and automatized in the workflow in /BigData_analysis.

