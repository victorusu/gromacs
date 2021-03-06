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
#include "gmxpre.h"

#include "bondtypes.h"

#include "gromacs/nblib/util.h"


namespace nblib
{

HarmonicBondType::HarmonicBondType(ForceConstant forceConstant, EquilDistance equilDistance) :
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

decltype(auto) HarmonicBondType::properties() const
{
    return std::tie(forceConstant_, equilDistance_);
}

bool HarmonicBondType::operator==(const HarmonicBondType& other) const
{
    return properties() == other.properties();
}

bool HarmonicBondType::operator<(const HarmonicBondType& other) const
{
    return properties() < other.properties();
}

G96BondType::G96BondType(ForceConstant forceConstant, EquilDistance equilDistance) :
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

decltype(auto) G96BondType::properties() const
{
    return std::tie(forceConstant_, equilDistance_);
}

bool G96BondType::operator==(const G96BondType& other) const
{
    return properties() == other.properties();
}

bool G96BondType::operator<(const G96BondType& other) const
{
    return properties() < other.properties();
}

CubicBondType::CubicBondType(ForceConstant quadraticForceConstant,
                             ForceConstant cubicForceConstant,
                             EquilDistance equilDistance) :
    quadraticForceConstant_(quadraticForceConstant),
    cubicForceConstant_(cubicForceConstant),
    equilDistance_(equilDistance)
{
}

decltype(auto) CubicBondType::properties() const
{
    return std::tie(quadraticForceConstant_, cubicForceConstant_, equilDistance_);
}

bool CubicBondType::operator==(const CubicBondType& other) const
{
    return properties() == other.properties();
}

bool CubicBondType::operator<(const CubicBondType& other) const
{
    return properties() < other.properties();
}


FENEBondType::FENEBondType(ForceConstant forceConstant, EquilDistance equilDistance) :
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

decltype(auto) FENEBondType::properties() const
{
    return std::tie(forceConstant_, equilDistance_);
}

bool FENEBondType::operator==(const FENEBondType& other) const
{
    return properties() == other.properties();
}

bool FENEBondType::operator<(const FENEBondType& other) const
{
    return properties() < other.properties();
}

MorseBondType::MorseBondType(ForceConstant forceConstant, Exponent exponent, EquilDistance equilDistance) :
    forceConstant_(forceConstant),
    exponent_(exponent),
    equilDistance_(equilDistance)
{
}

decltype(auto) MorseBondType::properties() const
{
    return std::tie(forceConstant_, exponent_, equilDistance_);
}

bool MorseBondType::operator==(const MorseBondType& other) const
{
    return properties() == other.properties();
}

bool MorseBondType::operator<(const MorseBondType& other) const
{
    return properties() < other.properties();
}

HalfAttractiveQuarticBondType::HalfAttractiveQuarticBondType(ForceConstant forceConstant,
                                                             EquilDistance equilDistance) :
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

decltype(auto) HalfAttractiveQuarticBondType::properties() const
{
    return std::tie(forceConstant_, equilDistance_);
}

bool HalfAttractiveQuarticBondType::operator==(const HalfAttractiveQuarticBondType& other) const
{
    return properties() == other.properties();
}

bool HalfAttractiveQuarticBondType::operator<(const HalfAttractiveQuarticBondType& other) const
{
    return properties() < other.properties();
}

} // namespace nblib
