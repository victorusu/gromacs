//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_MDRUNNER_H
#define GROMACS_MDRUNNER_H

#include <tuple>
#include <unordered_map>
#include <string>
#include <vector>

#include "atoms.h"
#include "interactions.h"
#include "simulationstate.h"
#include "nbkerneloptions.h"
#include "forcecalculator.h"

#include "gromacs/math/vectypes.h"

class TopologyBuilder;

namespace nblib
{

    // BondsForce bondsf();
    // AnglesForce anglesf();
    // NonBondedForce nonBondedf;

    // ForceCalculator(bondsf, anglesf, nonBondedf);

    // MDRunner runner(...);

    // for (int i = 0; i < N; i++) {
    //     runner.run(10).updatePairlist();
    // }

    // for (int i = 0; i < N; i++) {
    //     runner.computeForces().updatePairlist().integrate();
    // }

    // for (int i = 0; i < N; i++) {
    //     runner.computeForces().computeStrangeForces().updatePairlist().integrate();
    // }

    // for (int i = 0; i < N; i++) {
    //     runner.computeBondedForces().myStrangeMethod().computeNonbondedForces().updatePairlist().integrate();
    // }
// ;

class MDRunner {
public:

    // Add an update class, but for now update is a function.
    MDRunner(SimulationState simulationState, NBKernelOptions nbKernelOptions);

    MDRunner& updatePairList();

    MDRunner& update();

    MDRunner& computeForces();

    // MDRunner& run(const int numSteps);

private:

    SimulationState simulationState_;
    NBKernelOptions nbKernelOptions_;
    std::unique_ptr<ForceCalculator> forceCalculator_;

};

} //namespace nblib

#endif //GROMACS_MDRUNNER_H
