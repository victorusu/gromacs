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
 * This implements molecule setup tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "gromacs/nblib/atomtype.h"
#include "gromacs/nblib/molecules.h"

#include <iostream>

#include "testutils/testasserts.h"

namespace nblib {
namespace test {
namespace {

struct ArAtom
{
    AtomName name = "Ar";
    Mass     mass = 1.0;
    C6       c6   = 1;
    C12      c12  = 1;
};

struct OwAtom
{
    AtomName name = "Ow";
    Mass     mass = 16;
    C6       c6   = 1.;
    C12      c12  = 1.;
};

struct HwAtom
{
    AtomName name = "Hw";
    Mass     mass = 1;
    C6       c6   = 1.;
    C12      c12  = 1.;
};

TEST(NBlibTest, CanConstructMoleculeWithoutChargeOrResidueName)
{
    ArAtom   arAtom;
    AtomType Ar(arAtom.name, arAtom.mass, arAtom.c6, arAtom.c12);
    Molecule argon("Ar");
    EXPECT_NO_THROW(argon.addAtom(AtomName("Ar"), Ar));
}

TEST(NBlibTest, CanConstructMoleculeWithChargeWithoutResidueName)
{
    ArAtom   arAtom;
    AtomType Ar(arAtom.name, arAtom.mass, arAtom.c6, arAtom.c12);
    Molecule argon("Ar");
    EXPECT_NO_THROW(argon.addAtom(AtomName("Ar"), Charge(0), Ar));
}

TEST(NBlibTest, CanConstructMoleculeWithoutChargeWithResidueName)
{
    ArAtom   arAtom;
    AtomType Ar(arAtom.name, arAtom.mass, arAtom.c6, arAtom.c12);
    Molecule argon("Ar");
    EXPECT_NO_THROW(argon.addAtom(AtomName("Ar"), ResidueName("ar2"), Ar));
}

TEST(NBlibTest, CanConstructMoleculeWithChargeWithResidueName)
{
    ArAtom   arAtom;
    AtomType Ar(arAtom.name, arAtom.mass, arAtom.c6, arAtom.c12);
    Molecule argon("Ar");
    EXPECT_NO_THROW(argon.addAtom(AtomName("Ar"), ResidueName("ar2"), Charge(0), Ar));
}

TEST(NBlibTest, CanGetNumAtomsInMolecule)
{
    //! Manually Create Molecule (Water)

    //! 1. Define Atom Type
    OwAtom   owAtom;
    AtomType Ow(owAtom.name, owAtom.mass, owAtom.c6, owAtom.c12);
    HwAtom   hwAtom;
    AtomType Hw(hwAtom.name, hwAtom.mass, hwAtom.c6, hwAtom.c12);

    //! 2. Define Molecule
    Molecule water("water");

    water.addAtom(AtomName("Oxygen"), Charge(-0.6), Ow);
    water.addAtom(AtomName("H1"), Charge(+0.3), Hw);
    water.addAtom(AtomName("H2"), Hw);

    auto numAtoms = water.numAtomsInMolecule();

    EXPECT_EQ(3, numAtoms);
}

TEST(NBlibTest, CanConstructExclusionListFromNames)
{
    //! Manually Create Molecule (Water)

    //! 1. Define Atom Type
    OwAtom   owAtom;
    AtomType Ow(owAtom.name, owAtom.mass, owAtom.c6, owAtom.c12);
    HwAtom   hwAtom;
    AtomType Hw(hwAtom.name, hwAtom.mass, hwAtom.c6, hwAtom.c12);

    //! 2. Define Molecule
    Molecule water("water");

    water.addAtom(AtomName("Oxygen"), Charge(-0.6), Ow);
    water.addAtom(AtomName("H1"), Charge(+0.3), Hw);
    water.addAtom(AtomName("H2"), Charge(+0.3), Hw);

    water.addExclusion("H1", "Oxygen");
    water.addExclusion("H1", "H2");
    water.addExclusion("H2", "Oxygen");

    std::vector<std::tuple<int, int>> exclusions = water.getExclusions();

    std::vector<std::tuple<int, int>> reference{ {0,0}, {0,1}, {0,2},
                                                 {1,0}, {1,1}, {1,2},
                                                 {2,0}, {2,1}, {2,2}};

    ASSERT_EQ(exclusions.size(), 9);
    for (std::size_t i = 0; i < exclusions.size(); ++i)
        EXPECT_EQ(exclusions[i], reference[i]);

}

TEST(NBlibTest, CanConstructExclusionListFromNamesAndIndicesMixed)
{
    //! Manually Create Molecule (Water)

    //! 1. Define Atom Type
    OwAtom   owAtom;
    AtomType Ow(owAtom.name, owAtom.mass, owAtom.c6, owAtom.c12);
    HwAtom   hwAtom;
    AtomType Hw(hwAtom.name, hwAtom.mass, hwAtom.c6, hwAtom.c12);

    //! 2. Define Molecule
    Molecule water("water");

    water.addAtom(AtomName("Oxygen"), Charge(-0.6), Ow);
    water.addAtom(AtomName("H1"), Charge(+0.3), Hw);
    water.addAtom(AtomName("H2"), Charge(+0.3), Hw);

    water.addExclusion("H1", "Oxygen");
    water.addExclusion("H1", "H2");
    water.addExclusion(2, 0);

    std::vector<std::tuple<int, int>> exclusions = water.getExclusions();

    std::vector<std::tuple<int, int>> reference{ {0,0}, {0,1}, {0,2},
                                                 {1,0}, {1,1}, {1,2},
                                                 {2,0}, {2,1}, {2,2}};

    ASSERT_EQ(exclusions.size(), 9);
    for (std::size_t i = 0; i < exclusions.size(); ++i)
        EXPECT_EQ(exclusions[i], reference[i]);

}

}  // namespace
}  // namespace test
}  // namespace nblib
