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
 * This implements basic nblib tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 */

#include "forcecalculator.h"
#include "integrator.h"


#include "gromacs/compat/optional.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gridset.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"

#include "gromacs/math/vec.h"


namespace nblib {

static Nbnxm::KernelType translateBenchmarkEnum(const BenchMarkKernels &kernel)
{
    int kernelInt = static_cast<int>(kernel);
    return static_cast<Nbnxm::KernelType>(kernelInt);
}


ForceCalculator::ForceCalculator(SimulationState         &simulationState,
                                 const NBKernelOptions   &options)
                                 : simulationState_(simulationState), nbKernelOptions_(options)
{}

// //! Print timings outputs
// void ForceCalculator::printTimingsOutput(const NBKernelOptions &options,
//                                          const NBKernelSystem  &simulationState,
//                                          const gmx::index      &numPairs,
//                                          gmx_cycles_t           cycles)
// {
//     const gmx::EnumerationArray<BenchMarkKernels, std::string>  kernelNames   = { "auto", "no", "4xM", "2xMM" };
//     const gmx::EnumerationArray<BenchMarkCombRule, std::string> combruleNames = { "geom.", "LB", "none" };

//     // Generate an, accurate, estimate of the number of non-zero pair interactions
//     const real                          atomDensity          = simulationState.coordinates.size()/det(simulationState.box);
//     const real                          numPairsWithinCutoff = atomDensity*4.0/3.0*M_PI*std::pow(options.pairlistCutoff, 3);
//     const real                          numUsefulPairs       = simulationState.coordinates.size()*0.5*(numPairsWithinCutoff + 1);
// #if GMX_SIMD
//     if (options.nbnxmSimd != BenchMarkKernels::SimdNo)
//     {
//         fprintf(stdout, "SIMD width:           %d\n", GMX_SIMD_REAL_WIDTH);
//     }
// #endif
//     fprintf(stdout, "System size:          %zu atoms\n", simulationState.coordinates.size());
//     fprintf(stdout, "Cut-off radius:       %g nm\n", options.pairlistCutoff);
//     fprintf(stdout, "Number of threads:    %d\n", options.numThreads);
//     fprintf(stdout, "Number of iterations: %d\n", options.numIterations);
//     fprintf(stdout, "Compute energies:     %s\n",
//             options.computeVirialAndEnergy ? "yes" : "no");
//     if (options.coulombType != BenchMarkCoulomb::ReactionField)
//     {
//         fprintf(stdout, "Ewald excl. corr.:    %s\n",
//                 options.nbnxmSimd == BenchMarkKernels::SimdNo || options.useTabulatedEwaldCorr ? "table" : "analytical");
//     }
//     printf("\n");

//     fprintf(stdout, "Coulomb LJ   comb. SIMD    Mcycles  Mcycles/it.   %s\n",
//             options.cyclesPerPair ? "cycles/pair" : "pairs/cycle");
//     fprintf(stdout, "                                                total    useful\n");

//     fprintf(stdout, "%-7s %-4s %-5s %-4s ",
//             options.coulombType == BenchMarkCoulomb::Pme ? "Ewald" : "RF",
//             options.useHalfLJOptimization ? "half" : "all",
//             combruleNames[options.ljCombinationRule].c_str(),
//             kernelNames[options.nbnxmSimd].c_str());

//     const double dCycles = static_cast<double>(cycles);
//     if (options.cyclesPerPair)
//     {
//         fprintf(stdout, "%10.3f %10.4f %8.4f %8.4f\n",
//                 cycles*1e-6,
//                 dCycles/options.numIterations*1e-6,
//                 dCycles/(options.numIterations*numPairs),
//                 dCycles/(options.numIterations*numUsefulPairs));
//     }
//     else
//     {
//         fprintf(stdout, "%10.3f %10.4f %8.4f %8.4f\n",
//                 dCycles*1e-6,
//                 dCycles/options.numIterations*1e-6,
//                 options.numIterations*numPairs/dCycles,
//                 options.numIterations*numUsefulPairs/dCycles);
//     }
// }


//! Sets up and runs the kernel calls
//! TODO Refactor this function to return a handle to dispatchNonbondedKernel
//!      that callers can manipulate directly.
void ForceCalculator::compute(const int timestep, const bool printTimings)
{
    // gmx_omp_nthreads_set(emntPairsearch, options.numThreads);
    // gmx_omp_nthreads_set(emntNonbonded, options.numThreads);

    std::vector<NBKernelOptions> optionsList;
    expandSimdOptionAndPushBack(nbKernelOptions_, &optionsList);
    GMX_RELEASE_ASSERT(!optionsList.empty(), "Expect at least one benchmark setup");

    // We set the interaction cut-off to the pairlist cut-off
    interaction_const_t   ic   = setupInteractionConst(nbKernelOptions_);
    t_nrnb                nrnb = { 0 };
    gmx_enerdata_t        enerd(1, 0);

    gmx::StepWorkload     stepWork;
    stepWork.computeForces = true;
    if (nbKernelOptions_.computeVirialAndEnergy)
    {
        stepWork.computeVirial = true;
        stepWork.computeEnergy = true;
    }

    setupInstance(simulationState_, nbKernelOptions_);

    const PairlistSet                  &pairlistSet   = nonbondedVerlet_->pairlistSets().pairlistSet(gmx::InteractionLocality::Local);
    // const gmx::index                    num  Pairs      = pairlistSet.natpair_ljq_ + pairlistSet.natpair_lj_ + pairlistSet.natpair_q_;
    gmx_cycles_t                        cycles        = gmx_cycles_read();

    matrix boxMatrix;
    convertBoxToGMXFormat(simulationState_.box(), boxMatrix);


    t_forcerec forceRec;
    forceRec.ntype = simulationState_.topology().numAtoms();
    forceRec.nbfp  = simulationState_.topology().nonbondedParameters();
    snew(forceRec.shift_vec, SHIFTS);
    calc_shifts(boxMatrix, forceRec.shift_vec);

    std::vector<gmx::RVec> &currentCoords = simulationState_.coordinates();
    for (int iter = 0; iter < nbKernelOptions_.numIterations; iter++)
    {
        // Run the kernel without force clearing
        nonbondedVerlet_->dispatchNonbondedKernel(gmx::InteractionLocality::Local,
                                     ic, stepWork, enbvClearFNo, forceRec,
                                     &enerd,
                                     &nrnb);
        // There is one output data structure per thread
        std::vector<nbnxn_atomdata_output_t> nbvAtomsOut = nonbondedVerlet_->nbat.get()->out;
        integrateCoordinates(nbvAtomsOut, nbKernelOptions_, boxMatrix, currentCoords);
    }
    // Not copying back because currentCoords is a ref to simulationState_.coordinates();
    // simulationState_.coordinates = currentCoords;

    cycles = gmx_cycles_read() - cycles;
    // TODO
    // Enable printTimings
    // if (printTimings)
    // {
    //     printTimingsOutput(nbKernelOptions_, simulationState_, numPairs, cycles);
    // }
}

ForceCalculator& ForceCalculator::setupInstance(SimulationState &simulationState, const NBKernelOptions &nbKernelOptions)
{
    const auto         pinPolicy       = (nbKernelOptions.useGpu ? gmx::PinningPolicy::PinnedIfSupported : gmx::PinningPolicy::CannotBePinned);
    const int          numThreads      = nbKernelOptions.numThreads;
    // Note: the nbKernelOptions and Nbnxm combination rule enums values should match
    const int          combinationRule = static_cast<int>(nbKernelOptions.ljCombinationRule);

    auto               messageWhenInvalid = checkKernelSetup(nbKernelOptions);
    if (messageWhenInvalid)
    {
        gmx_fatal(FARGS, "Requested kernel is unavailable because %s.",
                  messageWhenInvalid->c_str());
    }
    Nbnxm::KernelSetup        kernelSetup = getKernelSetup(nbKernelOptions);

    PairlistParams            pairlistParams(kernelSetup.kernelType, false, nbKernelOptions.pairlistCutoff, false);

    Nbnxm::GridSet            gridSet(epbcXYZ, false, nullptr, nullptr, pairlistParams.pairlistType, false, numThreads, pinPolicy);

    auto                      pairlistSets = std::make_unique<PairlistSets>(pairlistParams, false, 0);

    auto                      pairSearch   = std::make_unique<PairSearch>(epbcXYZ, false, nullptr, nullptr,
                                                                          pairlistParams.pairlistType,
                                                                          false, numThreads, pinPolicy);

    auto atomData     = std::make_unique<nbnxn_atomdata_t>(pinPolicy);

    // Put everything together
    nonbondedVerlet_ = std::make_unique<nonbonded_verlet_t>(std::move(pairlistSets),
                                                    std::move(pairSearch),
                                                    std::move(atomData),
                                                    kernelSetup,
                                                    nullptr,
                                                    nullptr);

    auto numAtoms = simulationState.coordinates().size();

    nbnxn_atomdata_init(gmx::MDLogger(),
                        nonbondedVerlet_->nbat.get(), kernelSetup.kernelType,
                        combinationRule, numAtoms, simulationState.topology().nonbondedParameters(),
                        1, numThreads);

    t_nrnb nrnb;

    matrix boxMatrix;
    convertBoxToGMXFormat(simulationState.box(), boxMatrix);

    GMX_RELEASE_ASSERT(!TRICLINIC(boxMatrix), "Only rectangular unit-cells are supported here");
    const rvec               lowerCorner = { 0, 0, 0 };
    const rvec               upperCorner = {
        boxMatrix[XX][XX],
        boxMatrix[YY][YY],
        boxMatrix[ZZ][ZZ]
    };


    gmx::ArrayRef<const int> atomInfo;
    atomInfo = simulationState.topology().atomInfoAllVdw();

    const real atomDensity = simulationState.coordinates().size()/det(boxMatrix);

    nbnxn_put_on_grid(nonbondedVerlet_.get(),
                      boxMatrix, 0, lowerCorner, upperCorner,
                      nullptr, {0, int(simulationState.coordinates().size())}, atomDensity,
                      atomInfo, simulationState.coordinates(),
                      0, nullptr);

    nonbondedVerlet_->constructPairlist(gmx::InteractionLocality::Local,
                           &simulationState.topology().getGMXexclusions(), 0, &nrnb);

    t_mdatoms mdatoms;
    // We only use (read) the atom type and charge from mdatoms
    mdatoms.typeA   = const_cast<int *>(simulationState.topology().atomTypes().data());
    mdatoms.chargeA = const_cast<real *>(simulationState.topology().charges().data());
    nonbondedVerlet_->setAtomProperties(mdatoms, atomInfo);

    return *this;
}

interaction_const_t ForceCalculator::setupInteractionConst(const NBKernelOptions &options)
{
    interaction_const_t ic;

    ic.vdwtype          = evdwCUT;
    ic.vdw_modifier     = eintmodPOTSHIFT;
    ic.rvdw             = options.pairlistCutoff;

    ic.eeltype          = (options.coulombType == BenchMarkCoulomb::Pme ? eelPME : eelRF);
    ic.coulomb_modifier = eintmodPOTSHIFT;
    ic.rcoulomb         = options.pairlistCutoff;

    // Reaction-field with epsilon_rf=inf
    // TODO: Replace by calc_rffac() after refactoring that
    ic.k_rf             = 0.5*std::pow(ic.rcoulomb, -3);
    ic.c_rf             = 1/ic.rcoulomb + ic.k_rf*ic.rcoulomb*ic.rcoulomb;

    if (EEL_PME_EWALD(ic.eeltype))
    {
        // Ewald coefficients, we ignore the potential shift
        ic.ewaldcoeff_q = ewaldCoeff(1e-5, options.pairlistCutoff);
        GMX_RELEASE_ASSERT(ic.ewaldcoeff_q > 0, "Ewald coefficient should be > 0");
        ic.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
        init_interaction_const_tables(nullptr, &ic);
    }

    return ic;
}

void ForceCalculator::expandSimdOptionAndPushBack(const NBKernelOptions        &options,
                            std::vector<NBKernelOptions> *optionsList)
{
    if (options.nbnxmSimd == BenchMarkKernels::SimdAuto)
    {
        bool addedInstance = false;
#ifdef GMX_NBNXN_SIMD_4XN
        optionsList->push_back(options);
        optionsList->back().nbnxmSimd = BenchMarkKernels::Simd4XM;
        addedInstance = true;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
        optionsList->push_back(options);
        optionsList->back().nbnxmSimd = BenchMarkKernels::Simd2XMM;
        addedInstance = true;
#endif
        if (!addedInstance)
        {
            optionsList->push_back(options);
            optionsList->back().nbnxmSimd = BenchMarkKernels::SimdNo;
        }
    }
    else
    {
        optionsList->push_back(options);
    }
}

gmx::compat::optional<std::string> ForceCalculator::checkKernelSetup(const NBKernelOptions &options)
{
    GMX_RELEASE_ASSERT(options.nbnxmSimd < BenchMarkKernels::Count &&
                       options.nbnxmSimd != BenchMarkKernels::SimdAuto, "Need a valid kernel SIMD type");

    // Check SIMD support
    if ((options.nbnxmSimd != BenchMarkKernels::SimdNo && !GMX_SIMD)
#ifndef GMX_NBNXN_SIMD_4XN
        || options.nbnxmSimd == BenchMarkKernels::Simd4XM
#endif
#ifndef GMX_NBNXN_SIMD_2XNN
        || options.nbnxmSimd == BenchMarkKernels::Simd2XMM
#endif
        )
    {
        return "the requested SIMD kernel was not set up at configuration time";
    }

    return {};
}

Nbnxm::KernelSetup ForceCalculator::getKernelSetup(const NBKernelOptions &options)
{
    auto messageWhenInvalid = checkKernelSetup(options);
    GMX_RELEASE_ASSERT(!messageWhenInvalid, "Need valid options");

    Nbnxm::KernelSetup kernelSetup;

    //The int enum options.nbnxnSimd is set up to match Nbnxm::KernelType + 1
    kernelSetup.kernelType         = translateBenchmarkEnum(options.nbnxmSimd);
    // The plain-C kernel does not support analytical ewald correction
    if (kernelSetup.kernelType == Nbnxm::KernelType::Cpu4x4_PlainC)
    {
        kernelSetup.ewaldExclusionType = Nbnxm::EwaldExclusionType::Table;
    }
    else
    {
        kernelSetup.ewaldExclusionType = options.useTabulatedEwaldCorr ? Nbnxm::EwaldExclusionType::Table : Nbnxm::EwaldExclusionType::Analytical;
    }

    return kernelSetup;
}

real ForceCalculator::ewaldCoeff(const real ewald_rtol, const real pairlistCutoff)
{
    return calc_ewaldcoeff_q(pairlistCutoff, ewald_rtol);
}

} // namespace nblib
