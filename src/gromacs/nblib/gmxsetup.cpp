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
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "gmxsetup.h"

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/nblib/particletype.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"

namespace nblib
{

//! Helper to translate between the different enumeration values.
static Nbnxm::KernelType translateBenchmarkEnum(const BenchMarkKernels& kernel)
{
    int kernelInt = static_cast<int>(kernel);
    return static_cast<Nbnxm::KernelType>(kernelInt);
}

/*! \brief Checks the kernel setup
 *
 * Returns an error string when the kernel is not available.
 */
static std::optional<std::string> checkKernelSetup(const NBKernelOptions& options)
{
    GMX_RELEASE_ASSERT(options.nbnxmSimd < BenchMarkKernels::Count
                               && options.nbnxmSimd != BenchMarkKernels::SimdAuto,
                       "Need a valid kernel SIMD type");

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

void NbvSetupUtil::setExecutionContext(const NBKernelOptions& options)
{
    // Todo: find a more general way to initialize hardware
    gmx_omp_nthreads_set(emntPairsearch, options.numThreads);
    gmx_omp_nthreads_set(emntNonbonded, options.numThreads);
}

Nbnxm::KernelSetup NbvSetupUtil::getKernelSetup(const NBKernelOptions& options)
{
    auto messageWhenInvalid = checkKernelSetup(options);
    GMX_RELEASE_ASSERT(!messageWhenInvalid, "Need valid options");

    Nbnxm::KernelSetup kernelSetup;

    // The int enum options.nbnxnSimd is set up to match Nbnxm::KernelType + 1
    kernelSetup.kernelType = translateBenchmarkEnum(options.nbnxmSimd);
    // The plain-C kernel does not support analytical ewald correction
    if (kernelSetup.kernelType == Nbnxm::KernelType::Cpu4x4_PlainC)
    {
        kernelSetup.ewaldExclusionType = Nbnxm::EwaldExclusionType::Table;
    }
    else
    {
        kernelSetup.ewaldExclusionType = options.useTabulatedEwaldCorr
                                                 ? Nbnxm::EwaldExclusionType::Table
                                                 : Nbnxm::EwaldExclusionType::Analytical;
    }

    return kernelSetup;
}

void NbvSetupUtil::setParticleInfoAllVdv(const size_t numParticles)

{
    particleInfoAllVdw_.resize(numParticles);
    for (size_t particleI = 0; particleI < numParticles; particleI++)
    {
        SET_CGINFO_HAS_VDW(particleInfoAllVdw_[particleI]);
        SET_CGINFO_HAS_Q(particleInfoAllVdw_[particleI]);
    }
}

void NbvSetupUtil::setNonBondedParameters(const std::vector<ParticleType>& particleTypes,
                                          const NonBondedInteractionMap&   nonBondedInteractionMap)
{
    /* Todo: Refactor nbnxm to take nonbondedParameters_ directly
     *
     * initial self-handling of combination rules
     * size: 2*(numParticleTypes^2)
     */
    nonbondedParameters_.reserve(2 * particleTypes.size() * particleTypes.size());

    constexpr real c6factor  = 6.0;
    constexpr real c12factor = 12.0;

    for (const ParticleType& particleType1 : particleTypes)
    {
        for (const ParticleType& particleType2 : particleTypes)
        {
            nonbondedParameters_.push_back(
                    nonBondedInteractionMap.getC6(particleType1.name(), particleType2.name()) * c6factor);
            nonbondedParameters_.push_back(
                    nonBondedInteractionMap.getC12(particleType1.name(), particleType2.name()) * c12factor);
        }
    }
}

void NbvSetupUtil::setAtomProperties(const std::vector<int>&  particleTypeIdOfAllParticles,
                                     const std::vector<real>& charges)
{
    // We only use (read) the atom type and charge from mdatoms
    gmxForceCalculator_->mdatoms_.typeA   = const_cast<int*>(particleTypeIdOfAllParticles.data());
    gmxForceCalculator_->mdatoms_.chargeA = const_cast<real*>(charges.data());

    gmxForceCalculator_->nbv_->setAtomProperties(gmxForceCalculator_->mdatoms_, particleInfoAllVdw_);
}

//! Sets up and returns a Nbnxm object for the given options and system
void NbvSetupUtil::setupNbnxmInstance(const size_t numParticleTypes, const NBKernelOptions& options)
{
    const auto pinPolicy  = (options.useGpu ? gmx::PinningPolicy::PinnedIfSupported
                                           : gmx::PinningPolicy::CannotBePinned);
    const int  numThreads = options.numThreads;
    // Note: the options and Nbnxm combination rule enums values should match
    const int combinationRule = static_cast<int>(options.ljCombinationRule);

    auto messageWhenInvalid = checkKernelSetup(options);
    if (messageWhenInvalid)
    {
        gmx_fatal(FARGS, "Requested kernel is unavailable because %s.", messageWhenInvalid->c_str());
    }

    Nbnxm::KernelSetup kernelSetup = getKernelSetup(options);

    PairlistParams pairlistParams(kernelSetup.kernelType, false, options.pairlistCutoff, false);
    Nbnxm::GridSet gridSet(PbcType::Xyz, false, nullptr, nullptr, pairlistParams.pairlistType,
                           false, numThreads, pinPolicy);
    auto           pairlistSets = std::make_unique<PairlistSets>(pairlistParams, false, 0);
    auto           pairSearch =
            std::make_unique<PairSearch>(PbcType::Xyz, false, nullptr, nullptr,
                                         pairlistParams.pairlistType, false, numThreads, pinPolicy);

    auto atomData = std::make_unique<nbnxn_atomdata_t>(pinPolicy);

    // Put everything together
    auto nbv = std::make_unique<nonbonded_verlet_t>(std::move(pairlistSets), std::move(pairSearch),
                                                    std::move(atomData), kernelSetup, nullptr, nullptr);

    // Needs to be called with the number of unique ParticleTypes
    nbnxn_atomdata_init(gmx::MDLogger(), nbv->nbat.get(), kernelSetup.kernelType, combinationRule,
                        numParticleTypes, nonbondedParameters_, 1, numThreads);

    gmxForceCalculator_->nbv_ = std::move(nbv);
}

static real ewaldCoeff(const real ewald_rtol, const real pairlistCutoff)
{
    return calc_ewaldcoeff_q(pairlistCutoff, ewald_rtol);
}

void NbvSetupUtil::setupStepWorkload(const NBKernelOptions& options)
{
    gmx::StepWorkload stepWork;
    stepWork.computeForces          = true;
    stepWork.computeNonbondedForces = true;

    if (options.computeVirialAndEnergy)
    {
        stepWork.computeVirial = true;
        stepWork.computeEnergy = true;
    }

    gmxForceCalculator_->stepWork_ = stepWork;
}

void NbvSetupUtil::setupInteractionConst(const NBKernelOptions& options)
{
    interaction_const_t interactionConst;
    interactionConst.vdwtype      = evdwCUT;
    interactionConst.vdw_modifier = eintmodPOTSHIFT;
    interactionConst.rvdw         = options.pairlistCutoff;

    switch (options.coulombType)
    {
        case BenchMarkCoulomb::Pme: interactionConst.eeltype = eelPME; break;
        case BenchMarkCoulomb::Cutoff: interactionConst.eeltype = eelCUT; break;
        case BenchMarkCoulomb::ReactionField: interactionConst.eeltype = eelRF; break;
        case BenchMarkCoulomb::Count:
            GMX_THROW(gmx::InvalidInputError("Unsupported electrostatic interaction"));
    }
    interactionConst.coulomb_modifier = eintmodPOTSHIFT;
    interactionConst.rcoulomb         = options.pairlistCutoff;
    // Note: values correspond to ic.coulomb_modifier = eintmodPOTSHIFT
    interactionConst.dispersion_shift.cpot = -1.0 / gmx::power6(interactionConst.rvdw);
    interactionConst.repulsion_shift.cpot  = -1.0 / gmx::power12(interactionConst.rvdw);

    // These are the initialized values but we leave them here so that later
    // these can become options.
    interactionConst.epsilon_r  = 1.0;
    interactionConst.epsilon_rf = 1.0;

    /* Set the Coulomb energy conversion factor */
    if (interactionConst.epsilon_r != 0)
    {
        interactionConst.epsfac = ONE_4PI_EPS0 / interactionConst.epsilon_r;
    }
    else
    {
        /* eps = 0 is infinite dieletric: no Coulomb interactions */
        interactionConst.epsfac = 0;
    }

    calc_rffac(nullptr, interactionConst.epsilon_r, interactionConst.epsilon_rf,
               interactionConst.rcoulomb, &interactionConst.k_rf, &interactionConst.c_rf);

    if (EEL_PME_EWALD(interactionConst.eeltype))
    {
        // Ewald coefficients, we ignore the potential shift
        interactionConst.ewaldcoeff_q = ewaldCoeff(1e-5, options.pairlistCutoff);
        GMX_RELEASE_ASSERT(interactionConst.ewaldcoeff_q > 0, "Ewald coefficient should be > 0");
        interactionConst.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
        init_interaction_const_tables(nullptr, &interactionConst, 0);
    }

    gmxForceCalculator_->interactionConst_ = std::move(interactionConst);
}

void NbvSetupUtil::setupForceRec(const matrix& box)
{
    gmxForceCalculator_->forcerec_.nbfp = nonbondedParameters_;
    snew(gmxForceCalculator_->forcerec_.shift_vec, SHIFTS);
    calc_shifts(box, gmxForceCalculator_->forcerec_.shift_vec);
}

void NbvSetupUtil::setParticlesOnGrid(const std::vector<gmx::RVec>& coordinates, const Box& box)
{
    gmxForceCalculator_->setParticlesOnGrid(particleInfoAllVdw_, coordinates, box);
}

void NbvSetupUtil::constructPairList(const gmx::ListOfLists<int>& exclusions)
{
    t_nrnb nrnb;
    gmxForceCalculator_->nbv_->constructPairlist(gmx::InteractionLocality::Local, exclusions, 0, &nrnb);
}

void NbvSetupUtil::setForcesToZero(size_t numParticles)
{
    gmxForceCalculator_->verletForces_ = std::vector<gmx::RVec>(numParticles, gmx::RVec(0, 0, 0));
}

std::unique_ptr<GmxForceCalculator> GmxSetupDirector::setupGmxForceCalculator(const SimulationState& system,
                                                                              const NBKernelOptions& options)
{
    NbvSetupUtil nbvSetupUtil;
    nbvSetupUtil.setExecutionContext(options);
    nbvSetupUtil.setNonBondedParameters(system.topology().getParticleTypes(),
                                        system.topology().getNonBondedInteractionMap());
    nbvSetupUtil.setParticleInfoAllVdv(system.topology().numParticles());

    nbvSetupUtil.setupInteractionConst(options);
    nbvSetupUtil.setupStepWorkload(options);
    nbvSetupUtil.setupNbnxmInstance(system.topology().getParticleTypes().size(), options);
    nbvSetupUtil.setParticlesOnGrid(system.coordinates(), system.box());
    nbvSetupUtil.constructPairList(system.topology().getGmxExclusions());
    nbvSetupUtil.setAtomProperties(system.topology().getParticleTypeIdOfAllParticles(),
                                   system.topology().getCharges());

    const matrix& box = system.box().legacyMatrix();
    nbvSetupUtil.setupForceRec(box);
    nbvSetupUtil.setForcesToZero(system.topology().numParticles());

    return nbvSetupUtil.getGmxForceCalculator();
}

} // namespace nblib
