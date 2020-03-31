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

#include "util.h"
#include "bondtypes.h"

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
            else if (op(property(a), property(b)))
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
    name_(name),
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

auto HarmonicBondType::getProperties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.equilDistance(); },
                           [](const auto& x) { return x.forceConstant(); });
}

bool operator==(const HarmonicBondType& a, const HarmonicBondType& b)
{
    static_assert(sizeof(HarmonicBondType)
                          == sizeof(std::string) + sizeof(ForceConstant) + sizeof(EquilDistance),
                  "HarmonicBondType operator == incorrect");

    return detail::commonCompare(a, b, HarmonicBondType::getProperties());
}

bool operator<(const HarmonicBondType& a, const HarmonicBondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan), HarmonicBondType::getProperties());
}

G96BondType::G96BondType(Name name, ForceConstant forceConstant, EquilDistance equilDistance) :
    name_(name),
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

auto G96BondType::getProperties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.equilDistance(); },
                           [](const auto& x) { return x.forceConstant(); });
}

bool operator==(const G96BondType& a, const G96BondType& b)
{
    static_assert(sizeof(G96BondType) == sizeof(std::string) + sizeof(ForceConstant) + sizeof(EquilDistance),
                  "G96BondType operator == incorrect");

    return detail::commonCompare(a, b, G96BondType::getProperties());
}

bool operator<(const G96BondType& a, const G96BondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan), G96BondType::getProperties());
}

CubicBondType::CubicBondType(Name          name,
                             ForceConstant quadraticForceConstant,
                             ForceConstant cubicForceConstant,
                             EquilDistance equilDistance) :
    name_(name),
    quadraticForceConstant_(quadraticForceConstant),
    cubicForceConstant_(cubicForceConstant),
    equilDistance_(equilDistance)
{
}

auto CubicBondType::getProperties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.quadraticForceConstant(); },
                           [](const auto& x) { return x.cubicForceConstant(); },
                           [](const auto& x) { return x.equilDistance(); });
}

bool operator==(const CubicBondType& a, const CubicBondType& b)
{
    static_assert(sizeof(CubicBondType)
                          == sizeof(std::string) + 3 * sizeof(ForceConstant) + sizeof(EquilDistance),
                  "CubicBondType operator == incorrect");

    return detail::commonCompare(a, b, CubicBondType::getProperties());
}

bool operator<(const CubicBondType& a, const CubicBondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan), CubicBondType::getProperties());
}

FENEBondType::FENEBondType(Name name, ForceConstant forceConstant, EquilDistance equilDistance) :
    name_(name),
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

auto FENEBondType::getProperties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.equilDistance(); },
                           [](const auto& x) { return x.forceConstant(); });
}

bool operator==(const FENEBondType& a, const FENEBondType& b)
{
    static_assert(sizeof(FENEBondType) == sizeof(std::string) + sizeof(ForceConstant) + sizeof(EquilDistance),
                  "G96BondType operator == incorrect");

    return detail::commonCompare(a, b, FENEBondType::getProperties());
}

bool operator<(const FENEBondType& a, const FENEBondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan), FENEBondType::getProperties());
}

MorseBondType::MorseBondType(Name name, ForceConstant forceConstant, Exponent exponent, EquilDistance equilDistance) :
    name_(name),
    forceConstant_(forceConstant),
    exponent_(exponent),
    equilDistance_(equilDistance)
{
}

auto MorseBondType::getProperties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.forceConstant(); },
                           [](const auto& x) { return x.exponent(); },
                           [](const auto& x) { return x.equilDistance(); });
}

bool operator==(const MorseBondType& a, const MorseBondType& b)
{
    static_assert(sizeof(MorseBondType)
                          == sizeof(std::string) + 3 * sizeof(ForceConstant) + sizeof(EquilDistance),
                  "CubicBondType operator == incorrect");

    return detail::commonCompare(a, b, MorseBondType::getProperties());
}

bool operator<(const MorseBondType& a, const MorseBondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan), MorseBondType::getProperties());
}

HalfAttractiveQuarticBondType::HalfAttractiveQuarticBondType(Name          name,
                                                             ForceConstant forceConstant,
                                                             EquilDistance equilDistance) :
    name_(name),
    forceConstant_(forceConstant),
    equilDistance_(equilDistance)
{
}

auto HalfAttractiveQuarticBondType::getProperties()
{
    return std::make_tuple([](const auto& x) { return x.name(); },
                           [](const auto& x) { return x.equilDistance(); },
                           [](const auto& x) { return x.forceConstant(); });
}

bool operator==(const HalfAttractiveQuarticBondType& a, const HalfAttractiveQuarticBondType& b)
{
    static_assert(sizeof(HalfAttractiveQuarticBondType)
                          == sizeof(std::string) + sizeof(ForceConstant) + sizeof(EquilDistance),
                  "G96BondType operator == incorrect");

    return detail::commonCompare(a, b, HalfAttractiveQuarticBondType::getProperties());
}

bool operator<(const HalfAttractiveQuarticBondType& a, const HalfAttractiveQuarticBondType& b)
{
    auto lessThan = [](const auto& a, const auto& b) { return a < b; };
    return detail::commonRelational(a, b, std::move(lessThan),
                                    HalfAttractiveQuarticBondType::getProperties());
}

} // namespace nblib
