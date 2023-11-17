For compiling Pythia code with ROOT
$ g++ main91.cc -o main91 -w -I/home/christian/Softwares/pythia8310/include -O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread  -L/home/christian/Softwares/pythia8310/lib -Wl,-rpath,/home/christian/Softwares/pythia8310/lib -lpythia8 -ldl -L/home/christian/Softwares/root/lib -Wl,-rpath,/home/christian/Softwares/root/lib -lCore\
 -pthread -std=c++17 -m64 -I/home/christian/Softwares/root/include -L/home/christian/Softwares/root/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/home/christian/Softwares/root/lib -pthread -lm -ldl -rdynamic

For compiling Pythia code with FastJet
g++ main71.cc -o main71 -w -I/home/christian/Softwares/pythia8310/include -O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread  -L/home/christian/Softwares/pythia8310/lib -Wl,-rpath,/home/christian/Softwares/pythia8310/lib -lpythia8 -ldl -I/home/christian/Softwares/fastjet-install/include -L/home/christian/Softwares/fastjet-install/lib -Wl,-rpath,/home/christian/Softwares/fastjet-install/lib -lfastjet

For compiling Pythia code with ROOT and FastJet
g++ main95.cc -o main95 -w -I../include -O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread  -L../lib -Wl,-rpath,../lib -lpythia8 -ldl -I/home/christian/Softwares/fastjet-install/include\
     	-L/home/christian/Softwares/fastjet-install/lib -Wl,-rpath,/home/christian/Softwares/fastjet-install/lib -lfastjet -L/home/christian/Softwares/root/lib -Wl,-rpath,/home/christian/Softwares/root/lib -lCore -pthread -std=c++17 -m64 -I/home/christian/Softwares/root/include -L/home/christian/Softwares/root/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/home/christian/Softwares/root/lib -pthread -lm -ldl 

Or use the name of an example which uses both ROOT and FastJet for your macro, such as main95.cc and type
$ make main95

