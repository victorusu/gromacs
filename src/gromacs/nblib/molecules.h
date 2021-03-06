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
/*! \file
 * \brief
 * Implements nblib Molecule
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \inpublicapi
 * \ingroup nblib
 */
#ifndef GMX_NBLIB_MOLECULES_H
#define GMX_NBLIB_MOLECULES_H

#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "gromacs/nblib/bondtypes.h"
#include "gromacs/nblib/particletype.h"
#include "gromacs/nblib/util.h"

namespace nblib
{
class TopologyBuilder;

using ParticleName = std::string;
using Charge       = real;
using ResidueName  = std::string;

#define SUPPORTED_BOND_TYPES \
    HarmonicBondType, G96BondType, CubicBondType, FENEBondType, HalfAttractiveQuarticBondType

using SupportedBondTypes = TypeList<SUPPORTED_BOND_TYPES>;

class Molecule
{
    template<class Bond>
    struct BondData
    {
        using type = Bond;

        std::vector<Bond> interactionTypes_;
        std::vector<std::tuple<ParticleName, ResidueName, ParticleName, ResidueName>> interactions_;
    };

    // BondContainerTypes is TypeList<BondData<HarmonicBondType>, ...>
    using BondContainerTypes = Map<BondData, SupportedBondTypes>;
    // InteractionTuple is std::tuple<BondData<HarmonicBondType>, ...>
    using InteractionTuple = Reduce<std::tuple, BondContainerTypes>;

public:
    Molecule(std::string moleculeName);

    // Add a particle to the molecule with full specification of parameters.
    Molecule& addParticle(const ParticleName& particleName,
                          const ResidueName&  residueName,
                          const Charge&       charge,
                          ParticleType const& particleType);

    // Force explicit use of correct types
    template<typename T, typename U, typename V>
    Molecule& addParticle(const T&            particleName,
                          const U&            residueName,
                          const V&            charge,
                          ParticleType const& particleType) = delete;

    // Add a particle to the molecule with implicit charge of 0
    Molecule& addParticle(const ParticleName& particleName,
                          const ResidueName&  residueName,
                          ParticleType const& particleType);

    // Add a particle to the molecule with residueName set using particleName
    Molecule& addParticle(const ParticleName& particleName, const Charge& charge, ParticleType const& particleType);

    // Force explicit use of correct types, covers both implicit charge and residueName
    template<typename T, typename U>
    Molecule& addParticle(const T& particleName, const U& charge, ParticleType const& particleType) = delete;

    // Add a particle to the molecule with residueName set using particleName with implicit charge of 0
    Molecule& addParticle(const ParticleName& particleName, ParticleType const& particleType);

    // Force explicit use of correct types
    template<typename T>
    Molecule& addParticle(const T& particleName, ParticleType const& particleType) = delete;

    // TODO: add exclusions based on the unique ID given to the particle of the molecule
    void addExclusion(int particleIndex, int particleIndexToExclude);

    // Specify an exclusion with particle and residue names that have been added to molecule
    void addExclusion(std::tuple<std::string, std::string> particle,
                      std::tuple<std::string, std::string> particleToExclude);

    // Specify an exclusion with particle names that have been added to molecule
    void addExclusion(const std::string& particleName, const std::string& particleNameToExclude);

    //! add various types of interactions to the molecule
    //! Note: adding an interaction type not listed in SUPPORTED_BOND_TYPES results in a compilation error
    template<class Interaction>
    void addInteraction(const ParticleName& particleNameI,
                        const ResidueName&  residueNameI,
                        const ParticleName& particleNameJ,
                        const ResidueName&  residueNameJ,
                        const Interaction&  interaction);

    // add interactions with default residue name
    template<class Interaction>
    void addInteraction(const ParticleName& particleNameI,
                        const ParticleName& particleNameJ,
                        const Interaction&  interaction);

    // The number of molecules
    int numParticlesInMolecule() const;

    // Return the ParticleType data for a specific particle name that has been added to the molecule
    const ParticleType& at(const std::string& particlesTypeName) const;

    // convert exclusions given by name to indices and unify with exclusions given by indices
    // returns a sorted vector containing no duplicates of particles to exclude by indices
    std::vector<std::tuple<int, int>> getExclusions() const;

    // return all interactions stored in Molecule
    const InteractionTuple& interactionData() const;

    // return name of ith particle
    const ParticleName& particleName(int i) const;

    // return name of ith residue
    const ResidueName& residueName(int i) const;

    // The molecule name
    std::string name() const;

    friend class TopologyBuilder;

private:
    //! Name of the molecule
    std::string name_;

    struct ParticleData
    {
        std::string particleName_;
        std::string residueName_;
        std::string particleTypeName_;
        real        charge_;
    };

    //! one entry per particle in molecule
    std::vector<ParticleData> particles_;

    //! collection of distinct particle types in molecule
    std::unordered_map<std::string, ParticleType> particleTypes_;

    //! Used for calculated exclusions based on particle indices in molecule
    std::vector<std::tuple<int, int>> exclusions_;

    //! we cannot efficiently compute indices during the build-phase
    //! so we delay the conversion until TopologyBuilder requests it
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> exclusionsByName_;

    //! collection of data for all types of interactions
    InteractionTuple interactionData_;
};

#define ADD_INTERACTION_EXTERN_TEMPLATE(x)                                      \
    extern template void Molecule::addInteraction(                              \
            const ParticleName& particleNameI, const ResidueName& residueNameI, \
            const ParticleName& particleNameJ, const ResidueName& residueNameJ, const x& interaction);
MAP(ADD_INTERACTION_EXTERN_TEMPLATE, SUPPORTED_BOND_TYPES)
#undef ADD_INTERACTION_EXTERN_TEMPLATE

#define ADD_INTERACTION_EXTERN_TEMPLATE(x)         \
    extern template void Molecule::addInteraction( \
            const ParticleName& particleNameI, const ParticleName& particleNameJ, const x& interaction);
MAP(ADD_INTERACTION_EXTERN_TEMPLATE, SUPPORTED_BOND_TYPES)
#undef ADD_INTERACTION_EXTERN_TEMPLATE

} // namespace nblib
#endif // GMX_NBLIB_MOLECULES_H
