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
 * Implements kernels for nblib supported bondtypes
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef GMX_NBLIB_LISTEDKERNELS_HPP
#define GMX_NBLIB_LISTEDKERNELS_HPP

#include <tuple>

#include "gromacs/nblib/basicdefinitions.h"

namespace nblib
{

/*! \brief kernel to calculate the scalar part for simple harmonic bond forces
 *         for lambda = 0
 *
 * \param k spring constant
 * \param x0 equilibrium distance
 * \param x  input bond length
 *
 * \return tuple<force, potential energy>
 */
template <class T>
std::tuple<T, T> harmonicScalarForce(T k, T x0, T x)
{
    real dx  = x - x0;
    real dx2 = dx * dx;

    real force = -k * dx;
    real epot = 0.5 * k * dx2;

    return std::make_tuple(force, epot);

    /* That was 6 flops */
}


/*! \brief kernel to calculate the scalar part for simple harmonic bond forces
 *         for non-zero lambda to interpolate between A and B states
 *
 * \param kA spring constant state A
 * \param kB spring constant state B
 * \param xA equilibrium distance state A
 * \param xB equilibrium distance state B
 * \param x  input bond length
 * \param lambda interpolation factor between A and B state
 *
 * \return tuple<force, potential energy, lambda-interpolated energy>
 */
template <class T>
std::tuple<T, T, T> harmonicScalarForce(T kA, T kB, T xA, T xB, T x, T lambda)
{
    // code unchanged relative to Gromacs

    real L1 = 1.0 - lambda;
    real kk = L1 * kA + lambda * kB;
    real x0 = L1 * xA + lambda * xB;

    real dx  = x - x0;
    real dx2 = dx * dx;

    real force     = -kk * dx;
    real epot      = 0.5 * kk * dx2;
    real dvdlambda = 0.5 * (kB - kA) * dx2 + (xA - xB) * kk * dx;

    return std::make_tuple(force, epot, dvdlambda);

    /* That was 19 flops */
}

//! WIP: abstraction layer for different 2-center bonds
template <class T>
auto bondKernel(T dr, const HarmonicBondType& bond)
{
    return harmonicScalarForce(bond.forceConstant(), bond.equilDistance(), dr);
}

template <class T>
auto bondKernel(T dr, const G96BondType& bond)
{
    // Todo: implement me
    return std::make_tuple(real(0.0), real(0.0));
}

template <class T>
auto bondKernel(T dr, const MorseBondType& bond)
{
    // Todo: implement me
    return std::make_tuple(real(0.0), real(0.0));
}

template <class T>
auto bondKernel(T dr, const FENEBondType& bond)
{
    // Todo: implement me
    return std::make_tuple(real(0.0), real(0.0));
}

template <class T>
auto bondKernel(T dr, const CubicBondType& bond)
{
    // Todo: implement me
    return std::make_tuple(real(0.0), real(0.0));
}
// ... more bond types

} // namespace nblib
#endif // GMX_NBLIB_LISTEDINTERACTIONS_H
