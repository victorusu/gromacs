//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_MOLECULES_H
#define GROMACS_MOLECULES_H

#include <tuple>
#include <unordered_map>
#include <string>
#include <vector>

#include "atoms.h"
#include "interactions.h"

#include "gromacs/math/vectypes.h"

class TopologyBuilder;

namespace nblib
{

class MoleculeType {
public:
    MoleculeType(std::string moleculeName);

    MoleculeType& addAtom(const std::string &atomName, const std::string &residueName, AtomType const &atomType);

    MoleculeType& addAtom(const std::string &atomName, AtomType const &atomType);

    void addHarmonicBond(HarmonicType harmonicBond);

    // TODO: add exclusions based on the unique ID given to the atom of the molecule
    void addExclusion(const int atomIndex, const int atomIndexToExclude);
    void addExclusion(std::string atomName, std::string atomNameToExclude);

    int numAtomsInMolecule() const;

    friend class TopologyBuilder;

private:
    std::string name_;

    //! one entry per atom in molecule
    std::vector<std::tuple<std::string, std::string>> atoms_;
    //! collection of distinct AtomTypes in molecule
    std::unordered_map<std::string, AtomType> atomTypes_;

    std::vector<std::tuple<int, int>> exclusions_;

    std::vector<HarmonicType> harmonicInteractions_;

    int atomNameToIndex(std::string atomName);

};

} //namespace nblib
#endif //GROMACS_MOLECULES_H
