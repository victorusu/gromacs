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
 * Implements nblib supported bondtypes
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef GMX_NBLIB_BONDTYPES_H
#define GMX_NBLIB_BONDTYPES_H

#include "gromacs/nblib/particletype.h"

namespace nblib
{
using Name          = std::string;
using ForceConstant = real;
using EquilDistance = real;
using Exponent      = real;

/*! \brief Harmonic bond type
 *
 *  It represents the interaction of the form
 * V(r; forceConstant, equilDistance) = 0.5 * forceConstant * (r - equilDistance)^2
 */
class HarmonicBondType
{
public:
    HarmonicBondType() = default;

    HarmonicBondType(ForceConstant forceConstant, EquilDistance equilDistance);

    ForceConstant forceConstant() const { return forceConstant_; }

    EquilDistance equilDistance() const { return equilDistance_; }

    //! returns a tuple containing const references to the data members
    decltype(auto) properties() const;

    bool operator==(const HarmonicBondType& other) const;
    bool operator<(const HarmonicBondType& other) const;

private:
    ForceConstant forceConstant_;
    EquilDistance equilDistance_;
};


/*! \brief GROMOS bond type
 *
 * It represents the interaction of the form
 * V(r; forceConstant, equilDistance) = 0.25 * forceConstant * (r^2 - equilDistance^2)^2
 */
class G96BondType
{
public:
    G96BondType() = default;

    G96BondType(ForceConstant forceConstant, EquilDistance equilDistance);

    ForceConstant forceConstant() const { return forceConstant_; }

    EquilDistance equilDistance() const { return equilDistance_; }

    //! returns a tuple containing const references to the data members
    decltype(auto) properties() const;

    bool operator==(const G96BondType& other) const;
    bool operator<(const G96BondType& other) const;

private:
    ForceConstant forceConstant_;
    EquilDistance equilDistance_;
};


/*! \brief Cubic bond type
 *
 * It represents the interaction of the form
 * V(r; quadraticForceConstant, cubicForceConstant, equilDistance) = quadraticForceConstant * (r -
 * equilDistance)^2 + quadraticForceConstant * cubicForceConstant * (r - equilDistance)
 */
class CubicBondType
{
public:
    CubicBondType() = default;

    CubicBondType(ForceConstant quadraticForceConstant,
                  ForceConstant cubicForceConstant,
                  EquilDistance equilDistance);

    ForceConstant quadraticForceConstant() const { return quadraticForceConstant_; }

    ForceConstant cubicForceConstant() const { return cubicForceConstant_; }

    EquilDistance equilDistance() const { return equilDistance_; }

    //! returns a tuple containing const references to the data members
    decltype(auto) properties() const;

    bool operator==(const CubicBondType& other) const;
    bool operator<(const CubicBondType& other) const;

private:
    ForceConstant quadraticForceConstant_;
    ForceConstant cubicForceConstant_;
    EquilDistance equilDistance_;
};


/*! \brief FENE bond type
 *
 * It represents the interaction of the form
 * V(r; forceConstant, equilDistance) = - 0.5 * forceConstant * equilDistance^2 * log( 1 - (r / equilDistance)^2)
 */
class FENEBondType
{
public:
    FENEBondType() = default;

    FENEBondType(ForceConstant forceConstant, EquilDistance equilDistance);

    ForceConstant forceConstant() const { return forceConstant_; }

    EquilDistance equilDistance() const { return equilDistance_; }

    //! returns a tuple containing const references to the data members
    decltype(auto) properties() const;

    bool operator==(const FENEBondType& other) const;
    bool operator<(const FENEBondType& other) const;

private:
    ForceConstant forceConstant_;
    EquilDistance equilDistance_;
};


/*! \brief Morse bond type
 *
 * It represents the interaction of the form
 * V(r; forceConstant, exponent, equilDistance) = forceConstant * ( 1 - exp( -exponent * (r - equilDistance))
 */
class MorseBondType
{
public:
    MorseBondType() = default;

    MorseBondType(ForceConstant forceConstant, Exponent exponent, EquilDistance equilDistance);

    ForceConstant forceConstant() const { return forceConstant_; }

    Exponent exponent() const { return exponent_; }

    EquilDistance equilDistance() const { return equilDistance_; }

    //! returns a tuple containing const references to the data members
    decltype(auto) properties() const;

    bool operator==(const MorseBondType& other) const;
    bool operator<(const MorseBondType& other) const;

private:
    ForceConstant forceConstant_;
    Exponent      exponent_;
    EquilDistance equilDistance_;
};


/*! \brief Half-attractive quartic bond type
 *
 * It represents the interaction of the form
 * V(r; forceConstant, equilDistance) = 0.5 * forceConstant * (r - equilDistance)^4
 */
class HalfAttractiveQuarticBondType
{
public:
    HalfAttractiveQuarticBondType() = default;

    HalfAttractiveQuarticBondType(ForceConstant forceConstant, EquilDistance equilDistance);

    ForceConstant forceConstant() const { return forceConstant_; }

    EquilDistance equilDistance() const { return equilDistance_; }

    //! returns a tuple containing const references to the data members
    decltype(auto) properties() const;

    bool operator==(const HalfAttractiveQuarticBondType& other) const;
    bool operator<(const HalfAttractiveQuarticBondType& other) const;

private:
    ForceConstant forceConstant_;
    EquilDistance equilDistance_;
};

} // namespace nblib
#endif // GMX_NBLIB_BONDTYPES_H
