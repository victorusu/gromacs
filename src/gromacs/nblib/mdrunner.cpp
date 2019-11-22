#include "mdrunner.h"

namespace nblib {

MDRunner::MDRunner(SimulationState simulationState, NBKernelOptions nbKernelOptions)
: simulationState_(simulationState), nbKernelOptions_(nbKernelOptions)
{
    // TODO
    // finalize this function
    forceCalculator_ = std::make_unique<ForceCalculator>(simulationState_, nbKernelOptions_);
    forceCalculator_->setupInstance(simulationState_, nbKernelOptions_);
}

MDRunner& MDRunner::updatePairList()
{
    // TODO
    // implement this function
    return *this;
}

MDRunner& MDRunner::update()
{
    // TODO
    // implement this function
    return *this;
}

MDRunner& MDRunner::computeForces()
{
    // TODO
    // implement this function
    return *this;
}

} // namespace nblib
