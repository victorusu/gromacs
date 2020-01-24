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
 * Implements nblib ForceCalculator
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef GROMACS_FORCECALCULATOR_H
#define GROMACS_FORCECALCULATOR_H

#include "gromacs/timing/cyclecounter.h"

#include "nbkerneldef.h"
#include "nbkerneloptions.h"
#include "simulationstate.h"

namespace nblib {

enum class CombinationRule : int
{
    Geometric = 0,
    Count = 1
};

class ForceCalculator
{
public:

    // TODO: Depend on simulationState
    ForceCalculator(const SimulationState& system,
                    const NBKernelOptions& options);

    //! Sets up and runs the kernel calls
    //! TODO Refactor this function to return a handle to dispatchNonbondedKernel
    //!      that callers can manipulate directly.
    void compute(const bool printTimings = false);

private:

    void unpackTopologyToGmx();
    std::unique_ptr<nonbonded_verlet_t> setupNbnxmInstance();

    //void printTimingsOutput(const NBKernelOptions &options,
    //                        const SimulationState &system,
    //                        const gmx::index      &numPairs,
    //                        gmx_cycles_t           cycles);

    SimulationState system_;
    NBKernelOptions options_;

    //! Storage for parameters for short range interactions.
    std::vector<real> nonbondedParameters_;
    //! Atom masses
    std::vector<real> masses_;
    //! Atom info where all atoms are marked to have Van der Waals interactions
    std::vector<int> atomInfoAllVdw_;
    //! Legacy matrix for box
    matrix box_;
    //! Information about exclusions.
    t_blocka excls_;

};

} // namespace nblib

#endif //GROMACS_FORCECALCULATOR_H
