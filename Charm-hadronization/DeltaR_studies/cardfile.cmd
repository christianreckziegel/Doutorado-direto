Main:numberOfEvents = 100000 ! Generate 10k or 100k events for each energy point (5.02 TeV and 13.6 TeV)
Beams:eCM = 5020. ! Set the center-of-mass energy for the collisions to 5.02 TeV and 13.6 TeV
Beams:idA = 2212 ! pp collisions
Beams:idB = 2212
HardQCD:all = on ! Minimum bias events with hard QCD processes, which include charm production
PhaseSpace:pTHatMin = 1.0 ! No minimum pT cut for hard processes
421:mayDecay = false ! Force D0 to be stable so we can analyze its kinematics in the tree
-421:mayDecay = false ! Force D0bar to be stable so we can analyze its kinematics in the tree