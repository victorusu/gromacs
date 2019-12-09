//
// Created by sebkelle on 20.11.19.
//

#include "atomtype.h"

namespace nblib {

AtomType::AtomType() noexcept :
  name_(""),
  mass_(0),
  c6_(0),
  c12_(0)
{}

AtomType::AtomType(AtomTypeName atomName, AtomicMass mass, C6Param c6, C12Param c12)
: name_(std::move(atomName.name_)),
  mass_(mass.mass_),
  c6_(c6.c6_),
  c12_(c12.c12_)
{}

std::string AtomType::name() const { return name_; }

real AtomType::mass() const { return mass_; }

real AtomType::c6() const { return c6_; }

real AtomType::c12() const { return c12_; }

} // namespace nblib
