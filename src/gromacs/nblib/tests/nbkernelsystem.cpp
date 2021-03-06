/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * This implements topology setup tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "gmxpre.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/nblib/forcecalculator.h"
#include "gromacs/nblib/gmxsetup.h"
#include "gromacs/nblib/integrator.h"
#include "gromacs/nblib/topology.h"
#include "gromacs/topology/exclusionblocks.h"

#include "testhelpers.h"
#include "testsystems.h"

namespace nblib
{
namespace test
{
namespace
{

using ::testing::Eq;
using ::testing::Pointwise;

// This is defined in src/gromacs/mdtypes/forcerec.h but there is also a
// legacy C6 macro defined there that conflicts with the nblib C6 type.
// Todo: Once that C6 has been refactored into a regular function, this
//       file can just include forcerec.h
#define SET_CGINFO_HAS_VDW(cgi) (cgi) = ((cgi) | (1 << 23))

TEST(NBlibTest, SpcMethanolForcesAreCorrect)
{
    auto options        = NBKernelOptions();
    options.nbnxmSimd   = BenchMarkKernels::SimdNo;
    options.coulombType = BenchMarkCoulomb::Cutoff;

    SpcMethanolSimulationStateBuilder spcMethanolSystemBuilder;

    auto simState        = spcMethanolSystemBuilder.setupSimulationState();
    auto forceCalculator = ForceCalculator(simState, options);

    gmx::ArrayRef<gmx::RVec> forces;
    ASSERT_NO_THROW(forces = forceCalculator.compute());

    Vector3DTest forcesOutputTest;
    forcesOutputTest.testVectors(forces, "SPC-methanol forces");
}

TEST(NBlibTest, ExpectedNumberOfForces)
{
    auto options      = NBKernelOptions();
    options.nbnxmSimd = BenchMarkKernels::SimdNo;

    SpcMethanolSimulationStateBuilder spcMethanolSystemBuilder;

    auto simState        = spcMethanolSystemBuilder.setupSimulationState();
    auto forceCalculator = ForceCalculator(simState, options);

    gmx::ArrayRef<gmx::RVec> forces = forceCalculator.compute();
    EXPECT_EQ(simState.topology().numParticles(), forces.size());
}

TEST(NBlibTest, CanIntegrateSystem)
{
    auto options          = NBKernelOptions();
    options.nbnxmSimd     = BenchMarkKernels::SimdNo;
    options.numIterations = 1;

    SpcMethanolSimulationStateBuilder spcMethanolSystemBuilder;

    auto simState        = spcMethanolSystemBuilder.setupSimulationState();
    auto forceCalculator = ForceCalculator(simState, options);

    LeapFrog integrator(simState);

    for (int iter = 0; iter < options.numIterations; iter++)
    {
        gmx::ArrayRef<gmx::RVec> forces = forceCalculator.compute();
        // Todo: this makes it apparent that the SimState design is broken
        std::copy(forces.begin(), forces.end(), begin(simState.forces()));
        EXPECT_NO_THROW(integrator.integrate(1.0));
    }
}

TEST(NBlibTest, ArgonForcesAreCorrect)
{
    auto options        = NBKernelOptions();
    options.nbnxmSimd   = BenchMarkKernels::SimdNo;
    options.coulombType = BenchMarkCoulomb::Cutoff;

    ArgonSimulationStateBuilder argonSystemBuilder;

    auto simState        = argonSystemBuilder.setupSimulationState();
    auto forceCalculator = ForceCalculator(simState, options);

    gmx::ArrayRef<gmx::RVec> testForces;
    testForces = forceCalculator.compute();

    Vector3DTest forcesOutputTest;
    forcesOutputTest.testVectors(testForces, "Argon forces");
}

} // namespace
} // namespace test
} // namespace nblib
