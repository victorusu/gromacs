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
 * This implements topology setup tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "gromacs/nblib/atomtype.h"
#include "gromacs/nblib/topology.h"

#include "testutils/testasserts.h"

namespace nblib {
namespace test {
namespace {

class TwoWaterMolecules
{
public:
    Topology buildTopology()
    {
        //! Manually Create Molecule (Water)

        //! Define Atom Type
        AtomType Ow("Ow", 16, 1., 1.);
        AtomType Hw("Hw", 1, 1., 1.);

        //! Define Molecule
        Molecule water("water");

        //! Add the atoms
        water.addAtom("Oxygen", Charge(-0.6), Ow);
        water.addAtom("H1", Charge(+0.3), Hw);
        water.addAtom("H2", Charge(+0.3), Hw);

        //! Add the exclusions
        water.addExclusion("Oxygen", "H1");
        water.addExclusion("Oxygen", "H2");
        water.addExclusion("H1", "H2");

        // Todo: Add bonds functionality so this can be used/tested
        //water.addHarmonicBond(HarmonicType{1, 2, "H1", "Oxygen"});
        //water.addHarmonicBond(HarmonicType{1, 2, "H2", "Oxygen"});

        //! Build the topology
        TopologyBuilder topologyBuilder;

        //! Add some molecules to the topology
        topologyBuilder.addMolecule(water, 2);
        Topology topology = topologyBuilder.buildTopology();
        return topology;
    }
};

TEST(NBlibTest, TopologyHasCharges)
{
    TwoWaterMolecules waters;
    Topology watersTopology = waters.buildTopology();
    std::vector<real> testCharges = watersTopology.getCharges();
    std::vector<real> refCharges = { -0.6, 0.3, 0.3};
    EXPECT_EQ(refCharges, testCharges);
}

TEST(NBlibTest, TopologyHasMasses)
{
    TwoWaterMolecules waters;
    Topology watersTopology = waters.buildTopology();
    std::vector<real> testMasses = watersTopology.getMasses();
    std::vector<real> refMasses = { 16., 1., 1.};
    EXPECT_EQ(refMasses, testMasses);
}

}  // namespace
}  // namespace test
}  // namespace nblib
