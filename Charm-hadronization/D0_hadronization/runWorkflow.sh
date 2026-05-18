
#
# Script to run the entire analysis steps:
# 1 - MC fit templates + side-band subtraction
# 2 - Efficiency estimation: run 2 and run 3 style
# 3 - Feed-down subtraction
# 4 - Unfolding
#
# P.S.: closure tests are not included for now.
#
#!/bin/bash
# Make sure the script stops on the first error occurs
set -e

# Define working directory
WORKDIR=/home/christian/cernbox/Analyses/With_Nima/D0_hadronization
# Spawn a subshell to avoid polluting the current environment with environment variables
(
  # Enter O2 (ALICE software) environment
  #alienv enter O2Physics/latest RooUnfold/latest
  eval `alienv printenv O2Physics/latest RooUnfold/latest`
  
  # 0.0 - Find optimal BDT score thresholds
  cd $WORKDIR/1-SignalTreatment/BDTOptimization
  root -b -q BDTOptimization.C > log_bdt_optimization.txt 2>&1
  echo "0.0 - BDT optimization done!"
  cd $WORKDIR
  # 1.1 - MC fit templates
  cd $WORKDIR/1-SignalTreatment/Reflections
  root -b -q ReflectionsTreatment.C > log_reflections.txt 2>&1
  echo "1.1 - MC fit templates done!"
  cd $WORKDIR
  # 1.2 - Side-band subtraction
  cd $WORKDIR/1-SignalTreatment/SideBand
  root -b -q BackgroundSubtraction.C > log_sideband.txt 2>&1
  echo "1.2 - Background subtraction done!"
  cd $WORKDIR
  # 2.1 - Efficiency estimation: run 2 style
  cd $WORKDIR/2-Efficiency
  root -b -q Efficiency_run2_style.C > log_run2_style_efficiency.txt 2>&1
  echo "2.1 - Run 2 style efficiency done!"
  # 2.2 - Efficiency estimation: run 3 style
  root -b -q Efficiency_run3_style_detector_level.C > log_run3_style_efficiency.txt 2>&1
  echo "2.2 - Run 3 style efficiency done!"
  cd $WORKDIR
  # 3 - Feed-down subtraction
  cd $WORKDIR/3-Feed-Down
  root -b -q FeedDownSubtraction.C > log_feeddown.txt 2>&1
  echo "3 - Feed-down subtraction done!"
  cd $WORKDIR
  # 4 - Unfolding
  cd $WORKDIR/4-Unfolding
  root -b -q Unfolding.C > log_unfolding.txt 2>&1
  echo "4 - Unfolding procedure done!"
  cd $WORKDIR
  # 5.1 - First closure test: unfolding
  cd $WORKDIR/5-ClosureTest
  root -b -q FirstClosureTest.C > log_1st_closure_test.txt 2>&1
  echo "5.1 - First closure test done!"
  # 5.2 - Second closure test: background subtraction + efficiency correction + unfolding
  root -b -q SecondClosureTest.C > log_2nd_closure_test.txt 2>&1
  echo "5.2 - Second closure test done!"
  # 5.4 - Fourth closure test: my run 3 measurement vs. Emma Yeats' run 2 measurement
  root -b -q FourthClosureTest.C > log_4th_closure_test.txt 2>&1
  echo "5.4 - Fourth closure test done!"
  # End
  echo "All done!"
)
# End of subshell, back to original environment
echo "Back to original environment. Workflow finished."