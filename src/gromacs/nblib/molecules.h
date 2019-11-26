//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_MOLECULES_H
#define GROMACS_MOLECULES_H

#include <tuple>
#include <unordered_map>
#include <string>
#include <vector>

#include "interactions.h"

#include "gromacs/math/vectypes.h"

class AtomType;
class TopologyBuilder;

namespace nblib
{

using AtomName = std::string;
using Charge = real;
using ResidueName = std::string;

class Molecule {
public:
    Molecule(std::string moleculeName);

    Molecule& addAtom(const AtomName& atomName, const ResidueName& residueName, const Charge& charge, AtomType const &atomType);

    Molecule& addAtom(const AtomName& atomName, const ResidueName& residueName, AtomType const &atomType);

    Molecule& addAtom(const AtomName& atomName, const Charge& charge, AtomType const &atomType);

    Molecule& addAtom(const AtomName& atomName, AtomType const &atomType);

    void addHarmonicBond(HarmonicType harmonicBond);

    // TODO: add exclusions based on the unique ID given to the atom of the molecule
    void addExclusion(const int atomIndex, const int atomIndexToExclude);

    void addExclusion(std::tuple<std::string, std::string> atom, std::tuple<std::string, std::string> atomToExclude);

    void addExclusion(std::string atomName, std::string atomNameToExclude);

    int numAtomsInMolecule() const;

    friend class TopologyBuilder;

private:
    std::string name_;

    //! one entry per atom in molecule
    std::vector<std::tuple<std::string, std::string>> atoms_;

    //! collection of distinct Atoms in molecule
    std::unordered_map<std::string, std::tuple<AtomType, real>> atomTypes_;

    std::vector<std::tuple<int, int>> exclusions_;

    std::vector<HarmonicType> harmonicInteractions_;

    int atomNameAndResidueToIndex(std::tuple<std::string, std::string> atomResNameTuple);

    void addAtomSelfExclusion(std::string atomName, std::string resName);

};

} //namespace nblib
#endif //GROMACS_MOLECULES_H
