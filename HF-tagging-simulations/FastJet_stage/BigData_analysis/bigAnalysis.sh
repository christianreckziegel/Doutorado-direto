#! /bin/bash

start_time=$(date +%s)

source /home/christian/Softwares/root/bin/thisroot.sh
g++ -g myAnalysisParameter.cc -o myAnalysisParameter `/home/christian/Softwares/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` `root-config --cflags --libs`
echo "Analysis compiled."

# total number of events to be analysed
numOfEvents=1000000
numberOfEventsInFileString="1m"

# number of events per run = 10k
eventsPerRun=10000 # it's safe until 10000
# running 10k events per analysis run
numOfRuns=$((numOfEvents/eventsPerRun))

echo "Starting analysis loop..."
# analyse 10k events in each run
for ((i = 0; i < numOfRuns; i++))
do
    #echo "Run number $i"
    # Add your analysis command here
    #./myAnalysisParameter <numberOfEventsPerRun> <startingEventNumber> <numberOfEventsInFile>
    ./myAnalysisParameter $eventsPerRun $(($i*$eventsPerRun)) "$numberOfEventsInFileString"
    #echo "$eventsPerRun events per run, starting event number: $(($i*$eventsPerRun))", input file string
    #example: ./myAnalysisParameter 100 0 1k

    # show processing progress
    if [ $(($i + 1)) -eq $(($numOfRuns*10/100)) ]
    then
        echo -en "\r10% done |==>                  |\n"
    elif [ $(($i + 1)) -eq $(($numOfRuns*20/100)) ]; then
        echo -en "\r20% done|====>                |\n"
    elif [ $(($i + 1)) -eq $(($numOfRuns*30/100)) ]; then
        echo -en "\r30% done|======>              |\n"
    elif [ $(($i + 1)) -eq $(($numOfRuns*40/100)) ]; then
        echo -en "\r40% done|========>            |\n"
    elif [ $(($i + 1)) -eq $(($numOfRuns*50/100)) ]; then
        echo -en "\r50% done|==========>          |\n"
    elif [ $(($i + 1)) -eq $(($numOfRuns*60/100)) ]; then
        echo -en "\r60% done|============>        |\n"
    elif [ $(($i + 1)) -eq $(($numOfRuns*70/100)) ]; then
        echo -en "\r70% done|==============>      |\n"
    elif [ $(($i + 1)) -eq $(($numOfRuns*80/100)) ];then
        echo -en "\r80% done|================>    |\n"
    elif [ $(($i + 1)) -eq $(($numOfRuns*90/100)) ]; then
        echo -en "\r90% done|==================>  |\n"
    elif [ $(($i + 1)) -eq $(($numOfRuns*100/100)) ]; then
        echo -en "\r100% done.|====================>|\n"
        echo -en "\rAnalysis loop done.\n"
        echo -en "" # jump line
    fi



    
done

echo "Merging files"
g++ -g bigAnalysisMerging.C -o bigAnalysisMerging `root-config --cflags --libs`
echo "Merging macro compiled."
./bigAnalysisMerging "$numberOfEventsInFileString" $numOfRuns
echo "Execution finalized."

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
elapsed_time_minutes=$((elapsed_time / 60))
# Print execution elapsed time
echo "Elapsed time: $elapsed_time_minutes minute(s)"