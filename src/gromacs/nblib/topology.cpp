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
 * Implements nblib Topology and TopologyBuilder
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "gmxpre.h"

#include "topology.h"

#include <numeric>
#include <sstream>

#include "gromacs/nblib/particletype.h"
#include "gromacs/nblib/util.h"
#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

namespace nblib
{

namespace detail
{

std::vector<gmx::ExclusionBlock> toGmxExclusionBlock(const std::vector<std::tuple<int, int>>& tupleList)
{
    std::vector<gmx::ExclusionBlock> ret;

    auto firstLowerThan = [](auto const& tup1, auto const& tup2) {
        return std::get<0>(tup1) < std::get<0>(tup2);
    };

    //! initialize pair of iterators delimiting the range of exclusions for
    //! the first particle in the list
    GMX_ASSERT(!tupleList.empty(), "tupleList must not be empty\n");
    auto range = std::equal_range(std::begin(tupleList), std::end(tupleList), tupleList[0], firstLowerThan);
    auto it1 = range.first;
    auto it2 = range.second;

    //! loop over all exclusions in molecule, linear in tupleList.size()
    while (it1 != std::end(tupleList))
    {
        gmx::ExclusionBlock localBlock;
        //! loop over all exclusions for current particle
        for (; it1 != it2; ++it1)
        {
            localBlock.atomNumber.push_back(std::get<1>(*it1));
        }

        ret.push_back(localBlock);

        //! update the upper bound of the range for the next particle
        if (it1 != end(tupleList))
        {
            it2 = std::upper_bound(it1, std::end(tupleList), *it1, firstLowerThan);
        }
    }

    return ret;
}

std::vector<gmx::ExclusionBlock> offsetGmxBlock(std::vector<gmx::ExclusionBlock> inBlock, int offset)
{
    //! shift particle numbers by offset
    for (auto& localBlock : inBlock)
    {
        std::transform(std::begin(localBlock.atomNumber), std::end(localBlock.atomNumber),
                       std::begin(localBlock.atomNumber), [offset](auto i) { return i + offset; });
    }

    return inBlock;
}

int ParticleSequencer::operator()(const std::string&  moleculeName,
                                  int                 moleculeNr,
                                  const ResidueName&  residueName,
                                  const ParticleName& particleName) const
{
    try
    {
        return data_.at(moleculeName).at(moleculeNr).at(residueName).at(particleName);
    }
    catch (const std::out_of_range& outOfRange)
    {
        // TODO: use string format function once we have it
        if (moleculeName != residueName)
        {
            printf("No particle %s in residue %s in molecule %s found\n", particleName.c_str(),
                   residueName.c_str(), moleculeName.c_str());
        }
        else
        {
            printf("No particle %s in molecule %s found\n", particleName.c_str(), moleculeName.c_str());
        }

        throw gmx::InvalidInputError(outOfRange.what());
    };
}

void ParticleSequencer::build(const std::vector<std::tuple<Molecule, int>>& moleculesList)
{
    int currentID = 0;
    for (auto& molNumberTuple : moleculesList)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        size_t          numMols  = std::get<1>(molNumberTuple);

        auto& moleculeMap = data_[molecule.name()];

        for (size_t i = 0; i < numMols; ++i)
        {
            auto& moleculeNrMap = moleculeMap[i];
            for (int j = 0; j < molecule.numParticlesInMolecule(); ++j)
            {
                moleculeNrMap[molecule.residueName(j)][molecule.particleName(j)] = currentID++;
            }
        }
    }
}

template<class B>
std::tuple<std::vector<size_t>, std::vector<B>> collectBonds(const std::vector<std::tuple<Molecule, int>>& molecules)
{
    std::vector<B>      collectedBonds;
    std::vector<size_t> expansionArray;
    for (auto& molNumberTuple : molecules)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        size_t          numMols  = std::get<1>(molNumberTuple);

        auto& interactions = pickType<B>(molecule.interactionData()).interactionTypes_;

        std::vector<size_t> moleculeExpansion(interactions.size());
        // assign indices to the bonds in the current molecule, continue counting from
        // the number of bonds seen so far (=collectedBonds.size())
        std::iota(begin(moleculeExpansion), end(moleculeExpansion), collectedBonds.size());

        std::copy(begin(interactions), end(interactions), std::back_inserter(collectedBonds));

        for (size_t i = 0; i < numMols; ++i)
        {
            std::copy(begin(moleculeExpansion), end(moleculeExpansion), std::back_inserter(expansionArray));
        }
    }
    return std::make_tuple(expansionArray, collectedBonds);
}

#define COLLECT_BONDS_INSTANTIATE_TEMPLATE(x)                              \
    template std::tuple<std::vector<size_t>, std::vector<x>> collectBonds( \
            const std::vector<std::tuple<Molecule, int>>&);
MAP(COLLECT_BONDS_INSTANTIATE_TEMPLATE, SUPPORTED_BOND_TYPES)

//! for each interaction, translate the (moleculeName, nr, residueName, particleName)-pair
//! to particle sequence ID pairs
template<class B>
std::vector<std::tuple<int, int>> sequencePairIDs(const std::vector<std::tuple<Molecule, int>>& molecules,
                                                  const detail::ParticleSequencer& particleSequencer)
{
    std::vector<std::tuple<int, int>> interactionPairIDs;

    // loop over all molecules
    for (auto& molNumberTuple : molecules)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        size_t          numMols  = std::get<1>(molNumberTuple);

        for (size_t i = 0; i < numMols; ++i)
        {
            auto& interactions = pickType<B>(molecule.interactionData()).interactions_;
            for (const auto& interactionPair : interactions)
            {
                // the particle sequence ID of the first particle
                int id1 = particleSequencer(molecule.name(), i, std::get<1>(interactionPair),
                                            std::get<0>(interactionPair));
                // the particle sequence ID of the second particle
                int id2 = particleSequencer(molecule.name(), i, std::get<3>(interactionPair),
                                            std::get<2>(interactionPair));

                // we choose to store the lower sequence ID first. this allows for better unit tests
                // that are agnostic to how the input was set up
                interactionPairIDs.emplace_back(std::min(id1, id2), std::max(id1, id2));
            }
        }
    }
    return interactionPairIDs;
}

#define SEQUENCE_PAIR_ID_INSTANTIATE_TEMPLATE(x)                   \
    template std::vector<std::tuple<int, int>> sequencePairIDs<x>( \
            const std::vector<std::tuple<Molecule, int>>&, const detail::ParticleSequencer&);
MAP(SEQUENCE_PAIR_ID_INSTANTIATE_TEMPLATE, SUPPORTED_BOND_TYPES)
#undef SEQUENCE_PAIR_ID_INSTANTIATE_TEMPLATE

template<class B>
std::tuple<std::vector<size_t>, std::vector<B>> eliminateDuplicateBonds(const std::vector<B>& aggregatedBonds)
{
    std::vector<size_t> uniqueIndices(aggregatedBonds.size());
    std::vector<B>      uniqueBondInstances;
    // if there are no interactions of type B we're done now
    if (aggregatedBonds.empty())
    {
        return std::make_tuple(uniqueIndices, uniqueBondInstances);
    }

    // create 0,1,2,... sequence
    std::iota(begin(uniqueIndices), end(uniqueIndices), 0);

    std::vector<std::tuple<B, size_t>> enumeratedBonds(aggregatedBonds.size());
    // append each interaction with its index
    std::transform(begin(aggregatedBonds), end(aggregatedBonds), begin(uniqueIndices),
                   begin(enumeratedBonds), [](B b, size_t i) { return std::make_tuple(b, i); });

    auto sortKey = [](const auto& t1, const auto& t2) { return std::get<0>(t1) < std::get<0>(t2); };
    // sort w.r.t bonds. the result will contain contiguous segments of identical bond instances
    // the associated int indicates the original index of each BondType instance in the input vector
    std::sort(begin(enumeratedBonds), end(enumeratedBonds), sortKey);

    // initialize it1 and it2 to delimit first range of equal BondType instances
    auto range = std::equal_range(begin(enumeratedBonds), end(enumeratedBonds), enumeratedBonds[0], sortKey);
    auto it1 = range.first;
    auto it2 = range.second;

    // number of unique instances of BondType B = number of contiguous segments in enumeratedBonds =
    //         number of iterations in the outer while loop below
    while (it1 != end(enumeratedBonds))
    {
        uniqueBondInstances.push_back(std::get<0>(*it1));

        // loop over all identical BondType instances;
        for (; it1 != it2; ++it1)
        {
            // we note down that the BondType instance at index <interactionIndex>
            // can be found in the uniqueBondInstances container at index <uniqueBondInstances.size()>
            int interactionIndex            = std::get<1>(*it1);
            uniqueIndices[interactionIndex] = uniqueBondInstances.size() - 1;
        }

        // Note it1 has been incremented and is now equal to it2
        if (it1 != end(enumeratedBonds))
        {
            it2 = std::upper_bound(it1, end(enumeratedBonds), *it1, sortKey);
        }
    }

    return make_tuple(uniqueIndices, uniqueBondInstances);
}

#define ELIMINATE_DUPLICATE_BONDS_INSTANTIATE_TEMPLATE(x)                             \
    template std::tuple<std::vector<size_t>, std::vector<x>> eliminateDuplicateBonds( \
            const std::vector<x>& aggregatedBonds);
MAP(ELIMINATE_DUPLICATE_BONDS_INSTANTIATE_TEMPLATE, SUPPORTED_BOND_TYPES)
#undef ELIMINATE_DUPLICATE_BONDS_INSTANTIATE_TEMPLATE

} // namespace detail

TopologyBuilder::TopologyBuilder() : numParticles_(0) {}

gmx::ListOfLists<int> TopologyBuilder::createExclusionsListOfLists() const
{
    const auto& moleculesList = molecules_;

    std::vector<gmx::ExclusionBlock> exclusionBlockGlobal;
    exclusionBlockGlobal.reserve(numParticles_);

    size_t particleNumberOffset = 0;
    for (const auto& molNumberTuple : moleculesList)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        size_t          numMols  = std::get<1>(molNumberTuple);

        std::vector<gmx::ExclusionBlock> exclusionBlockPerMolecule =
                detail::toGmxExclusionBlock(molecule.getExclusions());

        //! duplicate the exclusionBlockPerMolecule for the number of Molecules of (numMols)
        for (size_t i = 0; i < numMols; ++i)
        {
            auto offsetExclusions =
                    detail::offsetGmxBlock(exclusionBlockPerMolecule, particleNumberOffset);

            std::copy(std::begin(offsetExclusions), std::end(offsetExclusions),
                      std::back_inserter(exclusionBlockGlobal));

            particleNumberOffset += molecule.numParticlesInMolecule();
        }
    }

    gmx::ListOfLists<int> exclusionsListOfListsGlobal;
    for (const auto& block : exclusionBlockGlobal)
    {
        exclusionsListOfListsGlobal.pushBack(block.atomNumber);
    }

    return exclusionsListOfListsGlobal;
}

Topology::InteractionData TopologyBuilder::createInteractionData(const detail::ParticleSequencer& particleSequencer)
{
    Topology::InteractionData interactionData;

    // this code is doing the compile time equivalent of
    // for (int i = 0; i < interactionData.size(); ++i)
    //     create(get<i>(interactionData));

    auto create = [this, &particleSequencer](auto& interactionDataElement) {
        using BondType = typename std::decay_t<decltype(interactionDataElement)>::type;

        // first compression stage: each bond per molecule listed once,
        // eliminates duplicates from multiple identical molecules
        auto  compressedDataStage1 = detail::collectBonds<BondType>(this->molecules_);
        auto& expansionArrayStage1 = std::get<0>(compressedDataStage1);
        auto& moleculeBonds        = std::get<1>(compressedDataStage1);

        // second compression stage: recognize bond duplicates among bonds from all molecules put together
        auto  compressedDataStage2 = detail::eliminateDuplicateBonds(moleculeBonds);
        auto& expansionArrayStage2 = std::get<0>(compressedDataStage2);
        auto& uniqueBondInstances  = std::get<1>(compressedDataStage2);

        // combine stage 1 + 2 expansion arrays
        std::vector<size_t> expansionArray(expansionArrayStage1.size());
        std::transform(begin(expansionArrayStage1), end(expansionArrayStage1), begin(expansionArray),
                       [& S2 = expansionArrayStage2](size_t S1Element) { return S2[S1Element]; });

        // add data about BondType instances
        interactionDataElement.bondInstances = std::move(uniqueBondInstances);

        interactionDataElement.indices.resize(expansionArray.size());
        // pairIndices contains the particle sequence IDs (i,j) of all interaction pairs of type <BondType>
        auto pairIndices = detail::sequencePairIDs<BondType>(this->molecules_, particleSequencer);
        // zip pairIndices(i,j) + expansionArray(k) -> interactionDataElement.indices(i,j,k)
        std::transform(begin(pairIndices), end(pairIndices), begin(expansionArray),
                       begin(interactionDataElement.indices), [](auto pairIndex, auto bondIndex) {
                           return std::make_tuple(std::get<0>(pairIndex), std::get<1>(pairIndex), bondIndex);
                       });
    };

    for_each_tuple(create, interactionData);

    return interactionData;
}

template<typename T, class Extractor>
std::vector<T> TopologyBuilder::extractParticleTypeQuantity(Extractor&& extractor)
{
    auto& moleculesList = molecules_;

    //! returned object
    std::vector<T> ret;
    ret.reserve(numParticles_);

    for (auto& molNumberTuple : moleculesList)
    {
        Molecule& molecule = std::get<0>(molNumberTuple);
        size_t    numMols  = std::get<1>(molNumberTuple);

        for (size_t i = 0; i < numMols; ++i)
        {
            for (auto& particleData : molecule.particles_)
            {
                ret.push_back(extractor(particleData, molecule.particleTypes_));
            }
        }
    }

    return ret;
}

Topology TopologyBuilder::buildTopology()
{
    topology_.numParticles_ = numParticles_;

    topology_.exclusions_ = createExclusionsListOfLists();
    topology_.charges_    = extractParticleTypeQuantity<real>([](const auto& data, auto& map) {
        ignore_unused(map);
        return data.charge_;
    });

    // map unique ParticleTypes to IDs
    std::unordered_map<std::string, int> nameToId;
    for (auto& name_particleType_tuple : particleTypes_)
    {
        topology_.particleTypes_.push_back(name_particleType_tuple.second);
        nameToId[name_particleType_tuple.first] = nameToId.size();
    }

    topology_.particleTypeIdOfAllParticles_ =
            extractParticleTypeQuantity<int>([&nameToId](const auto& data, auto& map) {
                ignore_unused(map);
                return nameToId[data.particleTypeName_];
            });

    detail::ParticleSequencer particleSequencer;
    particleSequencer.build(molecules_);
    topology_.particleSequencer_ = std::move(particleSequencer);

    topology_.combinationRule_         = particleTypesInteractions_.getCombinationRule();
    topology_.nonBondedInteractionMap_ = particleTypesInteractions_.generateTable();

    topology_.interactionData_ = createInteractionData(topology_.particleSequencer_);

    // Check whether there is any missing term in the particleTypesInteractions compared to the
    // list of particletypes
    for (const auto& particleType1 : particleTypes_)
    {
        for (const auto& particleType2 : particleTypes_)
        {
            auto interactionKey = std::make_tuple(particleType1.first, particleType2.first);
            if (topology_.nonBondedInteractionMap_.count(interactionKey) == 0)
            {
                std::string message =
                        formatString("Missing nonbonded interaction parameters for pair {} {}",
                                     particleType1.first, particleType2.first);
                GMX_THROW(gmx::InvalidInputError(message));
            }
        }
    }

    return topology_;
}

TopologyBuilder& TopologyBuilder::addMolecule(const Molecule& molecule, const int nMolecules)
{
    /*!
     * 1. Push-back a tuple of molecule type and nMolecules
     * 2. Append exclusion list into the data structure
     */

    molecules_.emplace_back(molecule, nMolecules);
    numParticles_ += nMolecules * molecule.numParticlesInMolecule();

    for (const auto& name_type_tuple : molecule.particleTypes_)
    {
        //! If we already have the particleType, we need to make
        //! sure that the type's parameters are actually the same
        //! otherwise we would overwrite them
        if (particleTypes_.count(name_type_tuple.first) > 0)
        {
            if (!(particleTypes_.at(name_type_tuple.first) == name_type_tuple.second))
            {
                GMX_THROW(gmx::InvalidInputError(
                        "Differing ParticleTypes with identical names encountered"));
            }
        }
    }

    // Note: insert does nothing if the key already exists
    particleTypes_.insert(molecule.particleTypes_.begin(), molecule.particleTypes_.end());

    return *this;
}

void TopologyBuilder::addParticleTypesInteractions(const ParticleTypesInteractions& particleTypesInteractions)
{
    particleTypesInteractions_.merge(particleTypesInteractions);
}

const int& Topology::numParticles() const
{
    return numParticles_;
}

const std::vector<real>& Topology::getCharges() const
{
    return charges_;
}

const std::vector<ParticleType>& Topology::getParticleTypes() const
{
    return particleTypes_;
}

const std::vector<int>& Topology::getParticleTypeIdOfAllParticles() const
{
    return particleTypeIdOfAllParticles_;
}

int Topology::sequenceID(std::string moleculeName, int moleculeNr, ResidueName residueName, ParticleName particleName)
{
    return particleSequencer_(moleculeName, moleculeNr, residueName, particleName);
}

const NonBondedInteractionMap& Topology::getNonBondedInteractionMap() const
{
    return nonBondedInteractionMap_;
}

const Topology::InteractionData& Topology::getInteractionData() const
{
    return interactionData_;
}

CombinationRule Topology::getCombinationRule() const
{
    return combinationRule_;
}

} // namespace nblib
