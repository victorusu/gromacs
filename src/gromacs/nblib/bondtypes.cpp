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

namespace detail
{

template<class B, class... Ps>
bool commonCompare(const B& a, const B& b, const std::tuple<Ps...>& properties)
{
    auto elementComparison = transform_tuple([&a, &b](const auto& f) { return f(a) == f(b); }, properties);
    auto logicalAnd        = [](auto&& a, auto&& b) { return a && b; };
    return std17::apply([f = logicalAnd](auto&&... args) { return binary_fold(f, args...); },
                        elementComparison);
}

template<class B, class Op, class... Ps>
bool commonRelational(const B& a, const B& b, Op&& op, const std::tuple<Ps...>& properties)
{
    auto compare = [](const auto& a, const auto& b, auto&& op, const auto&... properties) {
        auto start = [](auto&& recurse, const auto& a, const auto& b, auto&& op,
                        const auto& property, const auto&... properties) {
            if (op(property(a), property(b)))
            {
                return true;
            }
            else if (op(property(b), property(a)))
            {
                return false;
            }
            else
            {
                return std::get<(sizeof...(properties) == 0)>(recurse)(
                        std::forward<decltype(recurse)>(recurse), a, b,
                        std::forward<decltype(op)>(op), properties...);
            }
        };

        auto stop = [](auto&&... args) {
            ignore_unused(args...);
            return false;
        };

        return start(std::make_tuple(start, stop), a, b, std::forward<decltype(op)>(op), properties...);
    };

    return std17::apply(
            [&a, &b, &op, f = compare](const auto&... args) {
                return f(a, b, std::forward<Op>(op), args...);
            },
            properties);
}

} // namespace detail

HarmonicBondType::HarmonicBondType(Name name, ForceConstant forceConstant, EquilDistance equilDistance) :
    name_(std::move(name)),
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

auto HarmonicBondType::properties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.equilDistance(); },
                           [](const auto& x) { return x.forceConstant(); });
}

bool operator==(const HarmonicBondType& a, const HarmonicBondType& b)
{
    return detail::commonCompare(a, b, HarmonicBondType::properties());
}

bool operator<(const HarmonicBondType& a, const HarmonicBondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan), HarmonicBondType::properties());
}

G96BondType::G96BondType(Name name, ForceConstant forceConstant, EquilDistance equilDistance) :
    name_(std::move(name)),
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

auto G96BondType::properties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.equilDistance(); },
                           [](const auto& x) { return x.forceConstant(); });
}

bool operator==(const G96BondType& a, const G96BondType& b)
{
    return detail::commonCompare(a, b, G96BondType::properties());
}

bool operator<(const G96BondType& a, const G96BondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan), G96BondType::properties());
}

CubicBondType::CubicBondType(Name          name,
                             ForceConstant quadraticForceConstant,
                             ForceConstant cubicForceConstant,
                             EquilDistance equilDistance) :
    name_(std::move(name)),
    quadraticForceConstant_(quadraticForceConstant),
    cubicForceConstant_(cubicForceConstant),
    equilDistance_(equilDistance)
{
}

auto CubicBondType::properties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.quadraticForceConstant(); },
                           [](const auto& x) { return x.cubicForceConstant(); },
                           [](const auto& x) { return x.equilDistance(); });
}

bool operator==(const CubicBondType& a, const CubicBondType& b)
{
    return detail::commonCompare(a, b, CubicBondType::properties());
}

bool operator<(const CubicBondType& a, const CubicBondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan), CubicBondType::properties());
}

FENEBondType::FENEBondType(Name name, ForceConstant forceConstant, EquilDistance equilDistance) :
    name_(std::move(name)),
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

auto FENEBondType::properties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.equilDistance(); },
                           [](const auto& x) { return x.forceConstant(); });
}

bool operator==(const FENEBondType& a, const FENEBondType& b)
{
    return detail::commonCompare(a, b, FENEBondType::properties());
}

bool operator<(const FENEBondType& a, const FENEBondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan), FENEBondType::properties());
}

MorseBondType::MorseBondType(Name name, ForceConstant forceConstant, Exponent exponent, EquilDistance equilDistance) :
    name_(std::move(name)),
    forceConstant_(forceConstant),
    exponent_(exponent),
    equilDistance_(equilDistance)
{
}

auto MorseBondType::properties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.forceConstant(); },
                           [](const auto& x) { return x.exponent(); },
                           [](const auto& x) { return x.equilDistance(); });
}

bool operator==(const MorseBondType& a, const MorseBondType& b)
{
    return detail::commonCompare(a, b, MorseBondType::properties());
}

bool operator<(const MorseBondType& a, const MorseBondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan), MorseBondType::properties());
}

HalfAttractiveQuarticBondType::HalfAttractiveQuarticBondType(Name          name,
                                                             ForceConstant forceConstant,
                                                             EquilDistance equilDistance) :
    name_(std::move(name)),
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

auto HalfAttractiveQuarticBondType::properties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.equilDistance(); },
                           [](const auto& x) { return x.forceConstant(); });
}

bool operator==(const HalfAttractiveQuarticBondType& a, const HalfAttractiveQuarticBondType& b)
{
    return detail::commonCompare(a, b, HalfAttractiveQuarticBondType::properties());
}

bool operator<(const HalfAttractiveQuarticBondType& a, const HalfAttractiveQuarticBondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan),
                                    HalfAttractiveQuarticBondType::properties());
}

} // namespace nblib
