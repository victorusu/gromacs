#include <cstdio>

#include "gromacs/nblib/particletype.h"
#include "gromacs/nblib/molecules.h"
#include "gromacs/nblib/topology.h"
#include "gromacs/nblib/box.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nblib/nbkerneloptions.h"
#include "gromacs/nblib/forcecalculator.h"
#include "gromacs/nblib/integrator.h"

using namespace nblib;

//! \internal \brief Parameters from gromos43A1
struct ArAtom
{
    ParticleName name = "Ar";
    Mass         mass = 39.94800;
    C6           c6   = 0.0062647225;
    C12          c12  = 9.847044e-06;
};

ParticleTypesInteractions interactions;
Molecule                  argonMolecule("AR");
TopologyBuilder           topologyBuilder;
Box                       box(6.05449);
int                       numParticles = 12;
NBKernelOptions           options      = NBKernelOptions();

std::vector<gmx::RVec> coordinates = {
    { 0.794, 1.439, 0.610 }, { 1.397, 0.673, 1.916 }, { 0.659, 1.080, 0.573 },
    { 1.105, 0.090, 3.431 }, { 1.741, 1.291, 3.432 }, { 1.936, 1.441, 5.873 },
    { 0.960, 2.246, 1.659 }, { 0.382, 3.023, 2.793 }, { 0.053, 4.857, 4.242 },
    { 2.655, 5.057, 2.211 }, { 4.114, 0.737, 0.614 }, { 5.977, 5.104, 5.217 },
};

std::vector<gmx::RVec> velocities = {
    { 0.0055, -0.1400, 0.2127 },   { 0.0930, -0.0160, -0.0086 }, { 0.1678, 0.2476, -0.0660 },
    { 0.1591, -0.0934, -0.0835 },  { -0.0317, 0.0573, 0.1453 },  { 0.0597, 0.0013, -0.0462 },
    { 0.0484, -0.0357, 0.0168 },   { 0.0530, 0.0295, -0.2694 },  { -0.0550, -0.0896, 0.0494 },
    { -0.0799, -0.2534, -0.0079 }, { 0.0436, -0.1557, 0.1849 },  { -0.0214, 0.0446, 0.0758 },
};

std::vector<gmx::RVec> forces = {
    { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
    { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
    { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
    { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
};

int main()
{
    // ArAtom has all the parameters for an argon atom
    ArAtom arAtom;
    // Create a particle with argon
    ParticleType argonAtom(arAtom.name, arAtom.mass);
    // Add the argon particle to a molecule
    argonMolecule.addParticle(ParticleName("Ar"), argonAtom);
    // Add non-bonded interactions for argon
    interactions.add(argonAtom.name(), arAtom.c6, arAtom.c12);
    // Add the requested number of argon molecules to a topology
    topologyBuilder.addMolecule(argonMolecule, numParticles);
    // Add the argon interactions to the topology
    topologyBuilder.addParticleTypesInteractions(interactions);
    // Build the topology
    Topology topology = topologyBuilder.buildTopology();
    // A simulation state contains all the molecular information about the system
    SimulationState simState(coordinates, velocities, forces, box, topology);
    // Use a simple cutoff rule for Coulomb
    options.coulombType = BenchMarkCoulomb::Cutoff;
    // Some performance flags can be set a run time
    options.nbnxmSimd = BenchMarkKernels::Simd4XM;
    // The force calculator contains all the data needed to compute forces
    ForceCalculator forceCalculator(simState, options);
    // Print some diagnostic info
    printf("initial forces on particle 0: x %4f y %4f z %4f\n", forces[0][0], forces[0][1], forces[0][2]);
    // The forces are handed back to the user
    gmx::ArrayRef<gmx::RVec> userForces = forceCalculator.compute();
    // Print some diagnostic info
    printf("  final forces on particle 0: x %4f y %4f z %4f\n", userForces[0][0], userForces[0][1],
           userForces[0][2]);
    // The forces are not automatically updated in case the user wants to add their own
    std::copy(userForces.begin(), userForces.end(), begin(simState.forces()));
    // Integration requires masses, positions, and forces
    LeapFrog integrator(simState);
    // Print some diagnostic info
    printf("initial position of particle 0: x %4f y %4f z %4f\n", simState.coordinates()[0][0],
           simState.coordinates()[0][1], simState.coordinates()[0][2]);
    // Integrate with a time step of 1 fs
    integrator.integrate(1.0);
    // Print some diagnostic info
    printf("  final position of particle 0: x %4f y %4f z %4f\n", simState.coordinates()[0][0],
           simState.coordinates()[0][1], simState.coordinates()[0][2]);
    return 0;
}
