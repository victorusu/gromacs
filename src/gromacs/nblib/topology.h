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
 * Implements nblib Topology and TopologyBuilder
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef GROMACS_TOPOLOGY_H
#define GROMACS_TOPOLOGY_H

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/topology/block.h"

#include "molecules.h"

struct t_blocka;

namespace gmx
{
struct ExclusionBlock;
}

namespace nblib
{

namespace detail
{
std::vector<gmx::ExclusionBlock> toGmxExclusionBlock(const std::vector<std::tuple<int, int>>& tupleList);
std::vector<gmx::ExclusionBlock> offsetGmxBlock(std::vector<gmx::ExclusionBlock> inBlock, int offset);
} // namespace detail

/*! \libinternal
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
public:
    //! Returns the total number of atoms in the system
    const int& numAtoms() const;

    //! Returns a vector of atom names
    const std::vector<AtomType>& getAtomTypes() const;

    //! Return the AtomType ID of all atoms
    const std::vector<int>& getAtomTypeIdOfallAtoms() const;

    //! Returns a vector of atom partial charges
    const std::vector<real>& getCharges() const;

    //! Returns a vector of atomic masses
    const std::vector<real>& getMasses() const;

    //! Returns full list of nonbondedParameters
    const std::vector<std::tuple<real, real>>& getNonbondedParameters() const;

    const std::vector<int>& getAtomInfoAllVdw() const;

    // TODO: This function is only needed for testing. Need
    //       another way for testing exclusion correctness
    const t_blocka& getGMXexclusions() const { return excls_; }

private:
    Topology() = default;

    friend class TopologyBuilder;

    //! Total number of atoms in the system
    int numAtoms_;
    //! Storage for parameters for short range interactions.
    std::vector<std::tuple<real, real>> nonbondedParameters_;
    //! unique collection of AtomTypes
    std::vector<AtomType> atomTypes_;
    //! store an ID of each atom's type
    std::vector<int> atomTypeIdOfAllAtoms_;
    //! Storage for atom partial charges.
    std::vector<real> charges_;
    //! Atom masses
    std::vector<real> masses_;
    //! Atom info where all atoms are marked to have Van der Waals interactions
    std::vector<int> atomInfoAllVdw_;
    //! Information about exclusions.
    t_blocka excls_;
};

/*! \libinternal
 * \ingroup nblib
 * \brief Topology Builder
 *
 * A helper class to assist building of topologies. They also ensure that
 * topologies only exist in a valid state within the scope of the
 * simulation program.
 */
class TopologyBuilder
{
public:
    //! Constructor
    TopologyBuilder();

    /*! \brief Builds and Returns a valid Topology
     *
     * This function accounts for all the molecules added along with their
     * exclusions and returns a topology with a valid state that is usable
     * by the GROMACS back-end.
     */
    Topology buildTopology();

    //! Adds a molecules of a certain type into the topology
    TopologyBuilder& addMolecule(const Molecule& moleculeType, int nMolecules);

private:
    //! Internally stored topology
    Topology topology_;

    //! Total number of atoms in the system
    int numAtoms_;

    //! List of molecule types and number of molecules
    std::vector<std::tuple<Molecule, int>> molecules_;

    //! Builds a GROMACS-compliant performant exclusions list aggregating exclusions from all molecules
    t_blocka createExclusionsList() const;

    //! Helper function to extract quantities like mass, charge, etc from the system
    template<typename T, class Extractor>
    std::vector<T> extractAtomTypeQuantity(Extractor extractor);

    //! distinct collection of AtomTypes
    std::unordered_map<std::string, AtomType> atomTypes_;
};

} // namespace nblib

#endif // GROMACS_TOPOLOGY_H
