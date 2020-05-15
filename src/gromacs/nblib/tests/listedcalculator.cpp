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
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include <gtest/gtest.h>

#include "gromacs/nblib/gmxsetup.h"
#include "gromacs/nblib/listedcalculator.h"

#include "testhelpers.h"
#include "testsystems.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NBlibTest, ListedForceCalculatorCanConstruct)
{
    ListedInteractionData interactions;
    EXPECT_NO_THROW(ListedForceCalculator listedForceCalculator(interactions));
}

TEST(NBlibTest, HarmonicScalarKernelCanCompute)
{
    real k = 1.1;
    real x0 = 1.0;
    real x = 1.2;

    real force, epot;
    std::tie(force, epot) = harmonicScalarForce(k, x0, x);

    EXPECT_REAL_EQ_TOL(-k* (x-x0), force, gmx::test::defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(0.5 * k* (x-x0)*(x-x0), epot, gmx::test::defaultFloatTolerance());
}

TEST(NBlibTest, CalcForces)
{
    HarmonicBondType bond1{376560, 0.136};
    HarmonicBondType bond2{313800, 0.1};
    std::vector<HarmonicBondType> bonds{bond1, bond2};
    std::vector<std::tuple<int, int, int>> indices{ {0,1,0}, {1,2,1} };

    std::vector<gmx::RVec> x{ {1.97, 1.46, 1.209}, {1.978, 1.415, 1.082}, {1.905, 1.46, 1.03}};
    std::vector<gmx::RVec> forces(3, gmx::RVec{0, 0, 0});

    real energy = calcForces(indices, bonds, x, &forces);

    EXPECT_REAL_EQ_TOL(energy, 0.2113433, gmx::test::defaultRealTolerance());
    std::vector<gmx::RVec> refForces{ {-22.8980637, 128.801575, 363.505951}, {-43.2698593, -88.0130997, -410.639252},
                                      {66.167923, -40.788475, 47.1333084}};

    for (int i = 0; i < refForces.size(); ++i)
    {
        for (int j = 0; j < DIM; ++j)
        {
            EXPECT_REAL_EQ_TOL(refForces[i][j], forces[i][j], gmx::test::defaultRealTolerance());
        }
    }
}

} // namespace
} // namespace test
} // namespace nblib
