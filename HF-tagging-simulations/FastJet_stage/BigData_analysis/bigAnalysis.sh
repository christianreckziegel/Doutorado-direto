#! /bin/bash

source /home/christian/Softwares/root/bin/thisroot.sh
g++ -g myAnalysisParameter.cc -o myAnalysisParameter `/home/christian/Softwares/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` `root-config --cflags --libs`
echo "myAnalysisParameter.cc compiled."

# total number of events to be analysed
numOfEvents=1000
numberOfEventsInFileString="1k"

# number of events per run = 10k
eventsPerRun=100
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
    #echo "$eventsPerRun events per run, starting event number: $(($i*$eventsPerRun))"

    # show processing progress
    if [ $(($i + 1)) -eq $(($numOfRuns*10/100)) ]
    then
        echo "10% done"
        echo "|==>                  |"
    elif [ $(($i + 1)) -eq $(($numOfRuns*20/100)) ]; then
        echo "20% done"
        echo "|====>                |"
    elif [ $(($i + 1)) -eq $(($numOfRuns*30/100)) ]; then
        echo "30% done"
        echo "|======>              |"
    elif [ $(($i + 1)) -eq $(($numOfRuns*40/100)) ]; then
        echo "40% done"
        echo "|========>            |"
    elif [ $(($i + 1)) -eq $(($numOfRuns*50/100)) ]; then
        echo "50% done"
        echo "|==========>          |"
    elif [ $(($i + 1)) -eq $(($numOfRuns*60/100)) ]; then
        echo "60% done"
        echo "|============>        |"
    elif [ $(($i + 1)) -eq $(($numOfRuns*70/100)) ]; then
        echo "70% done"
        echo "|==============>      |"
    elif [ $(($i + 1)) -eq $(($numOfRuns*80/100)) ];then
        echo "80% done"
        echo "|================>    |"
    elif [ $(($i + 1)) -eq $(($numOfRuns*90/100)) ]; then
        echo "90% done"
        echo "|==================>  |"
    elif [ $(($i + 1)) -eq $(($numOfRuns*100/100)) ]; then
        echo "100% done."
        echo "|====================>|"
        echo "Analysis loop done."
        echo "" # jump line
    fi



    
done

echo "Merging files"
g++ -g bigAnalysisMerging.C -o bigAnalysisMerging `root-config --cflags --libs`
./bigAnalysisMerging "$numberOfEventsInFileString" $numOfRuns

