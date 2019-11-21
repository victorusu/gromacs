//
// Created by sebkelle on 20.11.19.
//

#include "molecules.h"

namespace nblib {

MoleculeType::MoleculeType(std::string moleculeName) : moleculeName_(std::move(moleculeName)) {}

MoleculeType& MoleculeType::addAtom(const std::string &moleculeAtomName, const std::string &residueName, AtomType const &atom)
{
    // check whether we already have the atom type
    if (!atomTypes_.count(moleculeAtomName))
    {
        atomTypes_[moleculeAtomName] = atom;
    }

    atoms_.emplace_back(std::make_tuple(moleculeAtomName, residueName));

    return *this;
}

MoleculeType& MoleculeType::addAtom(const std::string &name, AtomType const &atom)
{
    return this->addAtom(name, "", atom);
}

int MoleculeType::numAtomsInMolecule() const
{
    return atoms_.size();
}

void MoleculeType::addHarmonicBond(HarmonicType harmonicBond)
{
    harmonicInteractions_.push_back(harmonicBond);
}

void MoleculeType::addExclusion(int atomWithExclusion, int atomToExclude)
{
    exclusions_.emplace_back(std::make_tuple(atomWithExclusion, atomToExclude));
}

} // namespace nblib
