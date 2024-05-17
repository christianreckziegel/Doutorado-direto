#!/bin/bash
OPTION="-b --aod-file AO2D_LHC22o_pass4_minBias_small.root --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error"
o2-analysis-bc-converter {OPTION} | \
o2-analysis-timestamp {OPTION} | \
o2-analysis-event-selection {OPTION} | \
o2-analysis-ft0-corrected-table {OPTION} | \
o2-analysis-trackselection {OPTION} | \
o2-analysis-track-propagation {OPTION} | \
o2-analysis-tracks-extra-converter {OPTION} | \
o2-analysis-pid-tpc-base {OPTION} | \
o2-analysis-pid-tpc-full {OPTION} | \
o2-analysis-pid-tof-base {OPTION} | \
o2-analysis-pid-tof-full {OPTION} | \
o2-analysis-track-to-collision-associator {OPTION} | \
o2-analysis-hf-candidate-creator-2prong {OPTION} | \
o2-analysis-hf-candidate-selector-d0 {OPTION} | \
o2-analysis-hf-derived-data-creator-d0-to-k-pi {OPTION} | \
o2-analysis-multiplicity-table {OPTION} | \
o2-analysis-je-jet-deriveddata-producer {OPTION} | \
o2-analysis-je-jet-deriveddata-trigger-producer {OPTION} | \
o2-analysis-je-jet-finder-d0-data-charged {OPTION} | \
o2-analysis-je-jet-hf-fragmentation -b --aod-writer-keep AOD/JETDISTTABLE/0 --aod-file AO2D_LHC22o_pass4_minBias_small.root --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error
