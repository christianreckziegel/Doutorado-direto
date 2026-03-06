In order to compile PYTHIA code use
Pure PYTHIA code
```
$ g++ myprog.cc -o myprog `$PYTHIA8PATH/bin/pythia8-config --cxxflags --libs`
```
PYTHIA with ROOT code
```
$ g++ generateCharmEvents.cc -o generateCharmEvents `$PYTHIA8PATH/bin/pythia8-config --cxxflags --libs` `root-config --cflags --libs`
```
PYTHIA with FastJet code
```
$ g++
```
PYTHIA with FastJet and ROOT code
```
$ g++ reconstructCharmJets.cc -o reconstructCharmJets `$FASTJET_DIR/bin/fastjet-config --cxxflags --libs --plugins` `root-config --cflags --glibs`
```
Standalone ROOT
```
$ g++ myprog.cxx -o myprog `root-config --cflags --glibs`
```
Standalone FastJet
```
$ g++ myprog.cxx -o myprog `$FASTJET_DIR/bin/fastjet-config --cxxflags --libs --plugins`
```