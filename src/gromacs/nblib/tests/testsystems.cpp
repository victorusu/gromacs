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
 * This implements nblib test systems
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "gmxpre.h"

#include "testsystems.h"

namespace nblib
{

//! User class to hold ParticleTypes
//! Note: this is not part of NBLIB, users should write their own
class ParticleLibrary
{
public:
    ParticleLibrary()
    {
        ParticleType Ow(ParticleName("Ow"), Mass(15.99940));
        ParticleType H(ParticleName("H"), Mass(1.008));
        ParticleType OMet(ParticleName("OMet"), Mass(15.999));
        ParticleType CMet(ParticleName("CMet"), Mass(15.035));
        ParticleType Ar(ParticleName("Ar"), Mass(39.94800));

        particles_.insert(std::make_pair(Ow.name(), Ow));
        particles_.insert(std::make_pair(H.name(), H));
        particles_.insert(std::make_pair(OMet.name(), OMet));
        particles_.insert(std::make_pair(CMet.name(), CMet));
        particles_.insert(std::make_pair(Ar.name(), Ar));

        c6_[Ow.name()]   = 0.0026173456;
        c6_[H.name()]    = 0;
        c6_[OMet.name()] = 0.0022619536;
        c6_[CMet.name()] = 0.0088755241;
        c6_[Ar.name()]   = 0.0062647225;

        c12_[Ow.name()]   = 2.634129e-06;
        c12_[H.name()]    = 0;
        c12_[OMet.name()] = 1.505529e-06;
        c12_[CMet.name()] = 2.0852922e-05;
        c12_[Ar.name()]   = 9.847044e-06;
    }

    ParticleType type(const ParticleName& particleName) const
    {
        return particles_.at(particleName);
    }
    C6  c6(const ParticleName& particleName) const { return c6_.at(particleName); }
    C12 c12(const ParticleName& particleName) const { return c12_.at(particleName); }

private:
    std::map<ParticleName, ParticleType> particles_;
    std::map<ParticleName, C6>           c6_;
    std::map<ParticleName, C12>          c12_;
};

std::unordered_map<std::string, Charge> Charges{ { "Ow", -0.82 },
                                                 { "Hw", +0.41 },
                                                 { "OMet", -0.574 },
                                                 { "CMet", +0.176 },
                                                 { "HMet", +0.398 } };

WaterMoleculeBuilder::WaterMoleculeBuilder() : water_("SOL")
{
    ParticleLibrary plib;

    //! Add the particles
    water_.addParticle(ParticleName("Oxygen"), Charges.at("Ow"), plib.type("Ow"));
    water_.addParticle(ParticleName("H1"), Charges.at("Hw"), plib.type("H"));
    water_.addParticle(ParticleName("H2"), Charges.at("Hw"), plib.type("H"));
}

Molecule WaterMoleculeBuilder::waterMolecule()
{
    addExclusionsFromNames();
    return water_;
}

Molecule WaterMoleculeBuilder::waterMoleculeWithoutExclusions()
{
    return water_;
}

void WaterMoleculeBuilder::addExclusionsFromNames()
{
    water_.addExclusion("H1", "Oxygen");
    water_.addExclusion("H2", "Oxygen");
    water_.addExclusion("H1", "H2");
}

MethanolMoleculeBuilder::MethanolMoleculeBuilder() : methanol_("MeOH")
{
    ParticleLibrary library;

    //! Add the particles
    methanol_.addParticle(ParticleName("Me1"), Charges.at("CMet"), library.type("CMet"));
    methanol_.addParticle(ParticleName("O2"), Charges.at("OMet"), library.type("OMet"));
    methanol_.addParticle(ParticleName("H3"), Charges.at("HMet"), library.type("H"));

    // Add the exclusions
    methanol_.addExclusion("Me1", "O2");
    methanol_.addExclusion("Me1", "H3");
    methanol_.addExclusion("H3", "O2");
}

Molecule MethanolMoleculeBuilder::methanolMolecule()
{
    return methanol_;
}


Topology WaterTopology::buildTopology(int numMolecules)
{
    ParticleLibrary library;

    ParticleTypesInteractions interactions;
    std::vector<std::string>  typeNames = { "Ow", "H" };
    for (const auto& name : typeNames)
    {
        interactions.add(name, library.c6(name), library.c12(name));
    }

    //! Add some molecules to the topology
    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(water(), numMolecules);

    //! Add non-bonded interaction information
    topologyBuilder.addParticleTypesInteractions(interactions);

    Topology topology = topologyBuilder.buildTopology();
    return topology;
}

Molecule WaterTopology::water()
{
    return waterMolecule_.waterMolecule();
}

Topology SpcMethanolTopologyBuilder::buildTopology(int numWater, int numMethanol)
{
    ParticleLibrary library;

    ParticleTypesInteractions interactions;
    std::vector<std::string>  typeNames = { "Ow", "H", "OMet", "CMet" };
    for (const auto& name : typeNames)
    {
        interactions.add(name, library.c6(name), library.c12(name));
    }

    //! Add some molecules to the topology
    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(methanol(), numMethanol);
    topologyBuilder.addMolecule(water(), numWater);

    //! Add non-bonded interaction information
    topologyBuilder.addParticleTypesInteractions(interactions);

    Topology topology = topologyBuilder.buildTopology();
    return topology;
}

Molecule SpcMethanolTopologyBuilder::methanol()
{
    return methanolMolecule_.methanolMolecule();
}

Molecule SpcMethanolTopologyBuilder::water()
{
    return waterMolecule_.waterMolecule();
}

ArgonTopologyBuilder::ArgonTopologyBuilder(const int& numParticles)
{
    ParticleLibrary library;

    ParticleTypesInteractions nbinteractions;
    nbinteractions.add("Ar", library.c6("Ar"), library.c12("Ar"));

    Molecule argonMolecule("AR");
    argonMolecule.addParticle(ParticleName("AR"), library.type("Ar"));

    topologyBuilder_.addMolecule(argonMolecule, numParticles);
    topologyBuilder_.addParticleTypesInteractions((nbinteractions));
}

Topology ArgonTopologyBuilder::argonTopology()
{
    return topologyBuilder_.buildTopology();
}

ArgonSimulationStateBuilder::ArgonSimulationStateBuilder() :
    box_(6.05449),
    topology_(ArgonTopologyBuilder(12).argonTopology())
{

    coordinates_ = {
        { 0.794, 1.439, 0.610 }, { 1.397, 0.673, 1.916 }, { 0.659, 1.080, 0.573 },
        { 1.105, 0.090, 3.431 }, { 1.741, 1.291, 3.432 }, { 1.936, 1.441, 5.873 },
        { 0.960, 2.246, 1.659 }, { 0.382, 3.023, 2.793 }, { 0.053, 4.857, 4.242 },
        { 2.655, 5.057, 2.211 }, { 4.114, 0.737, 0.614 }, { 5.977, 5.104, 5.217 },
    };

    velocities_ = {
        { 0.0055, -0.1400, 0.2127 },   { 0.0930, -0.0160, -0.0086 }, { 0.1678, 0.2476, -0.0660 },
        { 0.1591, -0.0934, -0.0835 },  { -0.0317, 0.0573, 0.1453 },  { 0.0597, 0.0013, -0.0462 },
        { 0.0484, -0.0357, 0.0168 },   { 0.0530, 0.0295, -0.2694 },  { -0.0550, -0.0896, 0.0494 },
        { -0.0799, -0.2534, -0.0079 }, { 0.0436, -0.1557, 0.1849 },  { -0.0214, 0.0446, 0.0758 },
    };
    forces_ = {
        { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
        { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
        { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
        { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
    };
}

void ArgonSimulationStateBuilder::setCoordinate(int particleNum, int dimension, real value)
{
    GMX_ASSERT((dimension >= 0 and dimension <= 2), "Must provide a valid dimension\n");
    coordinates_.at(particleNum)[dimension] = value;
}

void ArgonSimulationStateBuilder::setVelocity(int particleNum, int dimension, real value)
{
    GMX_ASSERT((dimension >= 0 and dimension <= 2), "Must provide a valid dimension\n");
    velocities_.at(particleNum)[dimension] = value;
}

SimulationState ArgonSimulationStateBuilder::setupSimulationState()
{
    return SimulationState(coordinates_, velocities_, forces_, box_, topology_);
}

const Topology& ArgonSimulationStateBuilder::topology() const
{
    return topology_;
}

Box& ArgonSimulationStateBuilder::box()
{
    return box_;
}

std::vector<gmx::RVec>& ArgonSimulationStateBuilder::coordinates()
{
    return coordinates_;
}

std::vector<gmx::RVec>& ArgonSimulationStateBuilder::velocities()
{
    return velocities_;
}

SpcMethanolSimulationStateBuilder::SpcMethanolSimulationStateBuilder() :
    box_(3.01000),
    topology_(SpcMethanolTopologyBuilder().buildTopology(1, 1))
{
    coordinates_ = {
        { 1.970, 1.460, 1.209 }, // Me1
        { 1.978, 1.415, 1.082 }, // O2
        { 1.905, 1.460, 1.030 }, // H3
        { 1.555, 1.511, 0.703 }, // Ow
        { 1.498, 1.495, 0.784 }, // Hw1
        { 1.496, 1.521, 0.623 }, // Hw2
    };

    velocities_ = {
        { -0.8587, -0.1344, -0.0643 }, { 0.0623, -0.1787, 0.0036 }, { -0.5020, -0.9564, 0.0997 },
        { 0.869, 1.245, 1.665 },       { 0.169, 0.275, 1.565 },     { 0.269, 2.275, 1.465 },
    };

    forces_ = {
        { 0.000, 0.000, 0.000 }, { 0.000, 0.000, 0.000 }, { 0.000, 0.000, 0.000 },
        { 0.000, 0.000, 0.000 }, { 0.000, 0.000, 0.000 }, { 0.000, 0.000, 0.000 },
    };
}

SimulationState SpcMethanolSimulationStateBuilder::setupSimulationState()
{
    return SimulationState(coordinates_, velocities_, forces_, box_, topology_);
}

std::vector<gmx::RVec>& SpcMethanolSimulationStateBuilder::coordinates()
{
    return coordinates_;
}

std::vector<gmx::RVec>& SpcMethanolSimulationStateBuilder::velocities()
{
    return velocities_;
}

} // namespace nblib
