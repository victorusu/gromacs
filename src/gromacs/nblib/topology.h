//
// Created by sebkelle on 19.11.19.
//

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/topology/block.h"
#include "molecules.h"

#ifndef GROMACS_TOPOLOGY_H
#define GROMACS_TOPOLOGY_H

struct t_blocka;

namespace nblib {

class Topology {
public:

    const std::vector<int>& atoms() const;

    const std::vector<real>& charges() const;

    const std::vector<real>& masses() const;

    const std::vector<real>& nonbondedParameters() const;

    const std::vector<int>& atomTypes() const;

    const std::vector<int>& atomInfoAllVdw() const;

    int numAtoms() const
    {
        // TODO
        // Revisit this return value
        return charges_.size();
    }

    // TODO: This function is only needed for testing. Need
    //       another way for testing exclusion correctness
    const t_blocka& getGMXexclusions() const
    {
        return excls_;
    }

private:
    Topology() = default;

    friend class TopologyBuilder;

    //! Storage for parameters for short range interactions.
    std::vector<real>      nonbondedParameters_;
    //! Storage for atom type parameters.
    std::vector<int>       atomTypes_;
    //! Storage for atom partial charges.
    std::vector<real>      charges_;
    //! Atom masses
    std::vector<real>      masses_;
    //! Atom info where all atoms are marked to have Van der Waals interactions
    std::vector<int>       atomInfoAllVdw_;
    //! Information about exclusions.
    t_blocka               excls_;


};

class TopologyBuilder {
public:
    TopologyBuilder() = default;

    Topology buildTopology();

    TopologyBuilder& addMolecule(Molecule moleculeType, int nMolecules);

private:
    Topology topology_;

    int numAtoms_;

    std::vector<std::tuple<Molecule, int>> molecules_;

    t_blocka createExclusionsList() const;

    template <class Extractor>
    std::vector<real> extractQuantity(Extractor extractor);

    void extractAllProperties();
};

} // namespace nblib

#endif //GROMACS_TOPOLOGY_H
