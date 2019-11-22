/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * This implements basic nblib tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 */

#ifndef GROMACS_FORCECALCULATOR_H
#define GROMACS_FORCECALCULATOR_H

#include "gromacs/timing/cyclecounter.h"
#include "gromacs/compat/optional.h"

#include "gromacs/compat/optional.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gridset.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"

// #include "nbkernelsystem.h"

#include "nbkerneloptions.h"
#include "nbkerneldef.h"
#include "simulationstate.h"

namespace nblib {

class ForceCalculator
{
public:

    ForceCalculator(SimulationState       &simulationState,
                    const NBKernelOptions &options);

    //! Sets up and runs the kernel calls
    //! TODO Refactor this function to return a handle to dispatchNonbondedKernel
    //!      that callers can manipulate directly.
    void compute(const int timestep, const bool printTimings = false);

    // Just passing SimulationState and NBKernelOptions to demonstrate dependency
    ForceCalculator& setupInstance(SimulationState& simulationState, const NBKernelOptions &nbKernelOptions);


private:

    // void printTimingsOutput(const NBKernelOptions &options,
    //                         const NBKernelSystem  &system,
    //                         const gmx::index      &numPairs,
    //                         gmx_cycles_t           cycles);

    SimulationState simulationState_;
    NBKernelOptions nbKernelOptions_;
    std::unique_ptr<nonbonded_verlet_t> nonbondedVerlet_;

    interaction_const_t setupInteractionConst(const NBKernelOptions &options);

    void expandSimdOptionAndPushBack(const NBKernelOptions        &options,
                            std::vector<NBKernelOptions> *optionsList);

    gmx::compat::optional<std::string> checkKernelSetup(const NBKernelOptions &options);

    Nbnxm::KernelSetup getKernelSetup(const NBKernelOptions &options);

    real ewaldCoeff(const real ewald_rtol, const real pairlistCutoff);
};

} // namespace nblib

#endif //GROMACS_FORCECALCULATOR_H
