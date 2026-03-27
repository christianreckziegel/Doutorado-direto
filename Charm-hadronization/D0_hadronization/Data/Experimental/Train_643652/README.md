

This dataset has
- experimental data (not Monte Carlo simulation)
- D0 jets
- from derived data `JE_HF_LHC23_pass4_Thin_2P3PDstar_D0CJ_4_D0_1`
- was deployed at March 26th, 2026
- the working points optimal cuts (background score) from BDT for this dataset are:
```
std::vector<std::pair<double, double>> bdtPtCuts = {
        {0, 0.12}, {1, 0.12}, {2, 0.12}, {3, 0.16}, {4, 0.2}, {5, 0.25}, {6, 0.4}, {7, 0.6}, {8, 0.8}, {10, 0.8}, {12, 1}, {16, 1}, {50, -1}
    };
```
Got it from `https://alimonitor.cern.ch/hyperloop/view-wagon/42852/configuration?timestamp=1773169281439`. Concistency: 12 binsPtMl intervals and 12 bins for cutsMl.
Example, for pT,D⁰ \in [1;2], optimal score cut is 0.12, and so on.
Only consider the values until the last pT limit (do not consider the score for entries with pT,D⁰ > 50 GeV/c)