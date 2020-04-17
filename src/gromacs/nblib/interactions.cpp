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
 * Implements nblib particle-types interactions
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#include "gmxpre.h"

#include "interactions.h"

#include <cmath>

#include "gromacs/nblib/util.h"
#include "gromacs/utility/exceptions.h"

namespace nblib
{

C6 NonBondedInteractionMap::getC6(const ParticleTypeName& name1, const ParticleTypeName& name2) const
{
    return std::get<0>(this->at(std::make_tuple(name1, name2)));
}

C12 NonBondedInteractionMap::getC12(const ParticleTypeName& name1, const ParticleTypeName& name2) const
{
    return std::get<1>(this->at(std::make_tuple(name1, name2)));
}

namespace detail
{
real combineNonbondedParameters(real v, real w, CombinationRule combinationRule)
{
    if (combinationRule == CombinationRule::Geometric)
    {
        return std::sqrt(v * w);
    }
    else
    {
        throw gmx::InvalidInputError("unknown LJ Combination rule specified\n");
    }
}

} // namespace detail

ParticleTypesInteractions::ParticleTypesInteractions(CombinationRule cr) : combinationRule_(cr) {}

void ParticleTypesInteractions::add(const ParticleTypeName& particleTypeName, C6 c6, C12 c12)
{
    auto insertLocation = singleParticleInteractionsMap_.insert(
            std::make_pair(particleTypeName, std::make_tuple(c6, c12)));

    if (!insertLocation.second) // if particleTypeName already existed
    {
        if (std::get<0>(insertLocation.first->second) != c6
            || std::get<1>(insertLocation.first->second) != c12)
        {
            std::string message = formatString(
                    "Attempting to add nonbonded interaction parameters for particle "
                    "type {} twice",
                    particleTypeName.c_str());
            GMX_THROW(gmx::InvalidInputError(message));
        }
    }
}

void ParticleTypesInteractions::add(const ParticleTypeName& particleTypeName1,
                                    const ParticleTypeName& particleTypeName2,
                                    C6                      c6,
                                    C12                     c12)
{
    auto interactionKey         = std::make_tuple(particleTypeName1, particleTypeName2);
    auto possibleInteractionKey = std::make_tuple(particleTypeName2, particleTypeName1);

    auto insertLocation = twoParticlesInteractionsMap_.insert(
            std::make_pair(interactionKey, std::make_tuple(c6, c12)));
    twoParticlesInteractionsMap_.insert(std::make_pair(possibleInteractionKey, std::make_tuple(c6, c12)));

    if (!insertLocation.second) // if particleTypeName already existed
    {
        if (std::get<0>(insertLocation.first->second) != c6
            || std::get<1>(insertLocation.first->second) != c12)
        {
            std::string message = formatString(
                    "Attempting to add nonbonded interaction parameters between the particle types "
                    "{} {} twice",
                    particleTypeName1.c_str(), particleTypeName2.c_str());
            GMX_THROW(gmx::InvalidInputError(message));
        }
    }
}

NonBondedInteractionMap ParticleTypesInteractions::generateTable()
{
    NonBondedInteractionMap nonbondedParameters_;

    // creating the combination rule based interaction matrix
    for (const auto& particleType1 : singleParticleInteractionsMap_)
    {
        real c6_1  = std::get<0>(particleType1.second);
        real c12_1 = std::get<1>(particleType1.second);

        for (const auto& particleType2 : singleParticleInteractionsMap_)
        {
            real c6_2  = std::get<0>(particleType2.second);
            real c12_2 = std::get<1>(particleType2.second);

            real c6_combo  = detail::combineNonbondedParameters(c6_1, c6_2, combinationRule_);
            real c12_combo = detail::combineNonbondedParameters(c12_1, c12_2, combinationRule_);

            auto interactionKey = std::make_tuple(particleType1.first, particleType2.first);
            nonbondedParameters_[interactionKey] = std::make_tuple(c6_combo, c12_combo);
        }
    }

    // updating the interaction matrix based on the user fine tuned parameters
    for (const auto& particleTypeTuple : twoParticlesInteractionsMap_)
    {
        real c6_combo  = std::get<0>(particleTypeTuple.second);
        real c12_combo = std::get<1>(particleTypeTuple.second);

        nonbondedParameters_[particleTypeTuple.first] = std::make_tuple(c6_combo, c12_combo);
    }

    std::set<ParticleTypeName> particleTypes;
    for (auto const& typeKey : nonbondedParameters_)
    { // we don't need to get<1> because the list is guaranteed to be symmetric
        particleTypes.insert(std::get<0>(typeKey.first));
    }

    // check whether there is any missing interaction
    for (const ParticleTypeName& particleTypeName1 : particleTypes)
    {
        for (const ParticleTypeName& particleTypeName2 : particleTypes)
        {
            auto interactionKey = std::make_tuple(particleTypeName1, particleTypeName2);
            if (nonbondedParameters_.count(interactionKey) == 0)
            {
                std::string message = formatString("Missing interaction between {} {}",
                                                   particleTypeName1.c_str(), particleTypeName2.c_str());
                GMX_THROW(gmx::InvalidInputError(message));
            }
        }
    }
    return nonbondedParameters_;
}

CombinationRule ParticleTypesInteractions::getCombinationRule() const
{
    return combinationRule_;
}

void ParticleTypesInteractions::merge(const ParticleTypesInteractions& other)
{
    for (const auto& keyval : other.singleParticleInteractionsMap_)
    {
        add(keyval.first, std::get<0>(keyval.second), std::get<1>(keyval.second));
    }

    for (const auto& keyval : other.twoParticlesInteractionsMap_)
    {
        add(std::get<0>(keyval.first), std::get<1>(keyval.first), std::get<0>(keyval.second),
            std::get<1>(keyval.second));
    }
}

} // namespace nblib
