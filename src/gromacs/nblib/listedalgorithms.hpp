/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Implements a force calculator based on GROMACS data structures.
 *
 * Intended for internal use inside the ForceCalculator.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#ifndef GROMACS_LISTED_ALGORITHMS_HPP
#define GROMACS_LISTED_ALGORITHMS_HPP

#include <tuple>
#include <vector>

#include "gromacs/nblib/listedinteractions.h"
#include "gromacs/nblib/listedkernels.hpp"

namespace nblib
{


/*! implement a loop over bonds for a given BondType and Kernel
 *  corresponds to e.g. the "bonds" function at Gromacs:bonded.cpp@450
 *
 * \tparam BondType
 * \tparam Kernel
 * \param indices
 * \param bondInstances
 * \param x
 * \param kernel
 * \return
 */
template <class BondType, class Kernel>
real calcForces(const std::vector<std::tuple<int, int, int>>& indices,
                const std::vector<BondType>& bondInstances,
                const std::vector<gmx::RVec>& x,
                std::vector<gmx::RVec>& force,
                Kernel&& kernel)
{
    real Epot = 0.0;

    for (const auto& index : indices)
    {
        int i = std::get<0>(index);
        int j = std::get<1>(index);
        const gmx::RVec& x1 = x[i];
        const gmx::RVec& x2 = x[j];
        const BondType& bond = bondInstances[std::get<2>(index)];

        real dr2 = dot(x1, x2);
        real dr  = std::sqrt(dr2);

        real v;
        std::tie(force[i], v) = kernel(x1, x2, bond);
        Epot += v;
        //force[j] = -1 * force[i];
    }

   return Epot;
}

} // namespace nblib

#endif // GROMACS_GMXCALCULATOR_H
