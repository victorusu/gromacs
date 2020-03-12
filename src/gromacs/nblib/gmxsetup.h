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
 * \brief Translation layer to GROMACS data structures for force calculations.
 *
 * Implements the translation layer between the user scope and
 * GROMACS data structures for force calculations. Sets up the
 * non-bonded verlet.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef GROMACS_GMXSETUP_H
#define GROMACS_GMXSETUP_H

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nbnxm/nbnxm.h"

#include "gmxcalculator.h"
#include "interactions.h"
#include "nbkerneloptions.h"

namespace nblib
{

class NbvSetupUtil
{
public:
    NbvSetupUtil(SimulationState system, const NBKernelOptions& options);

    //! Sets up and returns a GmxForceCalculator
    std::unique_ptr<GmxForceCalculator> setupGmxForceCalculator();

private:
    //! Sets non-bonded parameters to be used to build GMX data structures
    void setNonBondedParameters(const Topology& topology);

    //! Marks particles to have Van der Waals interactions
    void setParticleInfoAllVdv(size_t numParticles);

    //! Returns the kernel setup
    Nbnxm::KernelSetup getKernelSetup(const NBKernelOptions& options);

    //! Sets Particle Types and Charges and VdW params
    void setAtomProperties(std::unique_ptr<nonbonded_verlet_t>& nbv, t_mdatoms& mdatoms);

    //! Sets up non-bonded verlet on the GmxForceCalculator
    std::unique_ptr<nonbonded_verlet_t> setupNbnxmInstance(const Topology&        topology,
                                                           const NBKernelOptions& options);

    SimulationState                  system_;
    std::shared_ptr<NBKernelOptions> options_;

    //! Storage for parameters for short range interactions.
    std::vector<real> nonbondedParameters_;

    //! Particle info where all particles are marked to have Van der Waals interactions
    std::vector<int> particleInfoAllVdw_;
};

//! Set up StepWorkload data
gmx::StepWorkload setupStepWorkload(std::shared_ptr<NBKernelOptions> options);

//! Return an interaction constants struct with members set appropriately
interaction_const_t setupInteractionConst(std::shared_ptr<NBKernelOptions> options);


} // namespace nblib
#endif // GROMACS_GMXSETUP_H
