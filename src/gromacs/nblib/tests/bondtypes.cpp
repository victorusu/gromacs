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
/*! \internal \file
 * \brief
 * This implements basic nblib box tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "gromacs/nblib/bondtypes.h"

#include "gromacs/nblib/util.h"

#include "testutils/testasserts.h"

namespace nblib
{

namespace test_detail
{

template<class B>
void testThreeParameterBondEquality(const B& deduceType)
{
    ignore_unused(deduceType);
    B a("someName", 1, 2);
    B b("someName", 1, 2);
    EXPECT_TRUE(a == b);

    B c("someName", 1, 3);
    EXPECT_FALSE(a == c);
}

template<class B>
void testFourParameterBondEquality(const B& deduceType)
{
    ignore_unused(deduceType);
    B a("someName", 1, 2, 3);
    B b("someName", 1, 2, 3);
    EXPECT_TRUE(a == b);

    B c("someName", 2, 3, 4);
    EXPECT_FALSE(a == c);
}

template<class B>
void testThreeParameterBondLessThan(const B& deduceType)
{
    ignore_unused(deduceType);
    B a("h1", 1, 2);
    B b("h1", 1, 3);
    EXPECT_TRUE(a < b);
    EXPECT_FALSE(b < a);

    B c("h1", 1, 2);
    B d("h1", 1, 2);
    EXPECT_FALSE(c < d);

    B e("a", 1, 3);
    B f("b", 1, 2);
    EXPECT_TRUE(e < f);
    EXPECT_FALSE(f < e);
}

template<class B>
void testFourParameterBondLessThan(const B& deduceType)
{
    ignore_unused(deduceType);
    B a("h1", 1, 2, 1);
    B b("h1", 1, 3, 1);
    EXPECT_TRUE(a < b);
    EXPECT_FALSE(b < a);

    B c("h1", 1, 2, 3);
    B d("h1", 1, 2, 3);
    EXPECT_FALSE(c < d);
}

} // namespace test_detail

TEST(NBlibTest, BondTypesOperatorEqualWorks)
{
    auto bondList3 = std::make_tuple(HarmonicBondType(), G96BondType(), FENEBondType(),
                                     HalfAttractiveQuarticBondType());
    for_each_tuple([](const auto& b) { test_detail::testThreeParameterBondEquality(b); }, bondList3);

    auto bondList4 = std::make_tuple(CubicBondType(), MorseBondType());
    for_each_tuple([](const auto& b) { test_detail::testFourParameterBondEquality(b); }, bondList4);
}

TEST(NBlibTest, BondTypesLessThanWorks)
{
    auto bondList3 = std::make_tuple(HarmonicBondType(), G96BondType(), FENEBondType(),
                                     HalfAttractiveQuarticBondType());
    for_each_tuple([](const auto& b) { test_detail::testThreeParameterBondLessThan(b); }, bondList3);

    auto bondList4 = std::make_tuple(CubicBondType(), MorseBondType());
    for_each_tuple([](const auto& b) { test_detail::testFourParameterBondLessThan(b); }, bondList4);
}


} // namespace nblib
