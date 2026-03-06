

bool eventHasD0(const Pythia8::Event& event) {
    //
    for (int iParticle = 0; iParticle < event.size(); iParticle++) {
        if (abs(event[iParticle].id()) == 421) { // Check for D0 or D0bar
            return true;
        }
    }
    return false;
}

// check if the particle is a B meson, baryon or an excited state of a B hadron
bool isBHadron(int pdg)
{
    int apdg = abs(pdg);

    if(apdg/1000 == 5) return true;
    if((apdg/100)%10 == 5) return true;

    return false;
}

// check if the particle with index idx in the event has a b hadron ancestor
bool hasBHadronAncestor(int idx, const Pythia8::Event& event)
{
    const Pythia8::Particle &particle = event[idx];

    int mother1 = particle.mother1();
    int mother2 = particle.mother2();

    if(mother1 == 0) return false;

    for(int iMother = mother1; iMother <= mother2; iMother++)
    {
        const Pythia8::Particle& mother = event[iMother];

        int pdg = mother.id();
        int status = mother.status();

        // Check if this is a B hadron
        if(isBHadron(pdg)) {
            return true;
        }

        // Stop recursion at hadronization stage: 81 - 89 : primary hadrons produced by hadronization process
        if(abs(status) >= 81 && abs(status) <= 89) {
            continue;
        }

        // Recursive search
        if(hasBHadronAncestor(iMother, event)) {
            return true;
        }
    }

    return false;
}