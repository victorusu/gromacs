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
/*! \inpublicapi \file
 * \brief
 * Implements nblib Topology and TopologyBuilder
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef GMX_NBLIB_TOPOLOGY_H
#define GMX_NBLIB_TOPOLOGY_H

#include <vector>

#include "gromacs/nblib/molecules.h"
#include "gromacs/utility/listoflists.h"

#include "interactions.h"

namespace gmx
{
struct ExclusionBlock;
}

namespace nblib
{

namespace detail
{

// Converts tuples of particle indices to exclude to the gmx::ExclusionBlock format
std::vector<gmx::ExclusionBlock> toGmxExclusionBlock(const std::vector<std::tuple<int, int>>& tupleList);

// Add offset to all indices in inBlock
std::vector<gmx::ExclusionBlock> offsetGmxBlock(std::vector<gmx::ExclusionBlock> inBlock, int offset);

// extract all bonds of type B from a vector of molecules. The second argument tuple element
// specifies multiples of the molecule given as first tuple element. Let (S, I) denote the return
// value tuple. Then J[i] = I[S[i]] for all i in 0...S.size() is the full sequence of BondType
// instances as they occur in the input tuple
template<class B>
std::tuple<std::vector<size_t>, std::vector<B>> collectBonds(const std::vector<std::tuple<Molecule, int>>&);

#define COLLECT_BONDS_EXTERN_TEMPLATE(x)                                          \
    extern template std::tuple<std::vector<size_t>, std::vector<x>> collectBonds( \
            const std::vector<std::tuple<Molecule, int>>&);
MAP(COLLECT_BONDS_EXTERN_TEMPLATE, SUPPORTED_BOND_TYPES)
#undef COLLECT_BONDS_EXTERN_TEMPLATE

// Return a list of unique BondType instances U and an index list S of size aggregatedBonds.size()
// such that the BondType instance at aggregatedBonds[i] is equal to U[S[i]]
// returns std::tuple(S, U)
template<class B>
std::tuple<std::vector<size_t>, std::vector<B>> eliminateDuplicateBonds(const std::vector<B>& collectedBonds);

#define ELIMINATE_DUPLICATE_EXTERN_TEMPLATE(x)                                               \
    extern template std::tuple<std::vector<size_t>, std::vector<x>> eliminateDuplicateBonds( \
            const std::vector<x>& collectedBonds);
MAP(ELIMINATE_DUPLICATE_EXTERN_TEMPLATE, SUPPORTED_BOND_TYPES)
#undef ELIMINATE_DUPLICATE_EXTERN_TEMPLATE

//! Helper class for Topology to keep track of particle IDs
class ParticleSequencer
{
    // store by (molecule name, molecule nr, residue name, particle name)
    using DataType =
            std::unordered_map<std::string, std::unordered_map<int, std::unordered_map<ResidueName, std::unordered_map<ParticleName, int>>>>;

public:
    //! build sequence from a list of molecules
    void build(const std::vector<std::tuple<Molecule, int>>& moleculesList);

    //! access ID by (molecule name, molecule nr, residue name, particle name)
    int operator()(const std::string&, int, const ResidueName&, const ParticleName&) const;

private:
    DataType data_;
};

template<class B>
std::vector<std::tuple<int, int>> sequencePairIDs(const std::vector<std::tuple<Molecule, int>>&,
                                                  const detail::ParticleSequencer&);

#define SEQUENCE_PAIR_ID_EXTERN_TEMPLATE(x)                               \
    extern template std::vector<std::tuple<int, int>> sequencePairIDs<x>( \
            const std::vector<std::tuple<Molecule, int>>&, const detail::ParticleSequencer&);
MAP(SEQUENCE_PAIR_ID_EXTERN_TEMPLATE, SUPPORTED_BOND_TYPES)
#undef SEQUENCE_PAIR_ID_EXTERN_TEMPLATE

} // namespace detail

/*! \inpublicapi
 * \ingroup nblib
 * \brief System Topology
 *
 * Contains all topology information meant to be used by the simulation
 * engine internally. Private constructor ensures that a Topology object
 * exists in a scope in a valid state after it has been built using a
 * Topology Builder.
 */
class Topology
{
    template<class B>
    struct BondData
    {
        typedef B type;

        // tuple format: <particleID i, particleID j, BondInstanceIndex>
        std::vector<std::tuple<int, int, int>> indices;
        // vector of unique BondType instances of type B
        std::vector<B> bondInstances;
    };

    // std::tuple<BondData<BondType1>, ...>
    using InteractionData = Reduce<std::tuple, Map<BondData, SupportedBondTypes>>;

public:
    //! Returns the total number of particles in the system
    const int& numParticles() const;

    //! Returns a vector of particle types
    const std::vector<ParticleType>& getParticleTypes() const;

    //! Return the ParticleType ID of all particles
    const std::vector<int>& getParticleTypeIdOfAllParticles() const;

    //! Returns a vector of particles partial charges
    const std::vector<real>& getCharges() const;

    //! Returns exclusions in proper, performant, GROMACS layout
    const gmx::ListOfLists<int>& getGmxExclusions() const { return exclusions_; }

    int sequenceID(std::string moleculeName, int moleculeNr, ResidueName residueName, ParticleName particleName);

    //! Returns a map of non-bonded force parameters indexed by ParticleType names
    const NonBondedInteractionMap& getNonBondedInteractionMap() const;

    //! Returns the interaction data
    const InteractionData& getInteractionData() const;

    //! Returns the combination rule used to generate the NonBondedInteractionMap
    CombinationRule getCombinationRule() const;

private:
    Topology() = default;

    friend class TopologyBuilder;

    //! Total number of particles in the system
    int numParticles_;
    //! unique collection of ParticleTypes
    std::vector<ParticleType> particleTypes_;
    //! store an ID of each particles's type
    std::vector<int> particleTypeIdOfAllParticles_;
    //! Storage for particles partial charges
    std::vector<real> charges_;
    //! Information about exclusions.
    gmx::ListOfLists<int> exclusions_;
    //! Associate molecule, residue and particle names with sequence numbers
    detail::ParticleSequencer particleSequencer_;
    //! Map that should hold all nonbonded interactions for all particle types
    NonBondedInteractionMap nonBondedInteractionMap_;
    //! data about bonds for all supported types
    InteractionData interactionData_;
    //! Combination Rule used to generate the nonbonded interactions
    CombinationRule combinationRule_;
};

/*! \brief Topology Builder
 *
 * \libinternal
 * \ingroup nblib
 *
 * A helper class to assist building of topologies. They also ensure that
 * topologies only exist in a valid state within the scope of the
 * simulation program.
 *
 */
class TopologyBuilder
{
public:
    //! Constructor
    TopologyBuilder();

    /*! \brief
     * Builds and Returns a valid Topology
     *
     * This function accounts for all the molecules added along with their
     * exclusions and returns a topology with a valid state that is usable
     * by the GROMACS back-end.
     */
    Topology buildTopology();

    // Adds a molecules of a certain type into the topology
    TopologyBuilder& addMolecule(const Molecule& moleculeType, int nMolecules);

    void addParticleTypesInteractions(const ParticleTypesInteractions& particleTypesInteractions);

private:
    //! Internally stored topology
    Topology topology_;

    //! Total number of particles in the system
    int numParticles_;

    //! List of molecule types and number of molecules
    std::vector<std::tuple<Molecule, int>> molecules_;

    // Builds a GROMACS-compliant performant exclusions list aggregating exclusions from all molecules
    gmx::ListOfLists<int> createExclusionsListOfLists() const;

    // Gather interaction data from molecules
    Topology::InteractionData createInteractionData(const detail::ParticleSequencer&);

    // Helper function to extract quantities like mass, charge, etc from the system
    template<typename T, class Extractor>
    std::vector<T> extractParticleTypeQuantity(Extractor&& extractor);

    //! distinct collection of ParticleTypes
    std::unordered_map<std::string, ParticleType> particleTypes_;

    //! ParticleType nonbonded parameters
    ParticleTypesInteractions particleTypesInteractions_;
};

//! utility function to extract Particle quantities and expand them to the full
//! array of length numParticles()
template<class F>
inline auto expandQuantity(const Topology& topology, F&& particleTypeExtractor)
{
    using ValueType = decltype((std::declval<ParticleType>().*std::declval<F>())());

    std::vector<ValueType> ret;
    ret.reserve(topology.numParticles());

    const std::vector<ParticleType>& particleTypes = topology.getParticleTypes();

    for (size_t id : topology.getParticleTypeIdOfAllParticles())
    {
        ret.push_back((particleTypes[id].*particleTypeExtractor)());
    }

    return ret;
}

} // namespace nblib

#endif // GMX_NBLIB_TOPOLOGY_H
