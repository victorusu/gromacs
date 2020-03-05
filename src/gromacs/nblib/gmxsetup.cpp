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

#include "gromacs/compat/optional.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/nblib/particletype.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/pbcutil/ishift.h"
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
static gmx::compat::optional<std::string> checkKernelSetup(const NBKernelOptions& options)
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

NbvSetupUtil::NbvSetupUtil(SimulationState system, const NBKernelOptions& options) :
        system_(std::move(system))
{
    options_ = std::make_shared<NBKernelOptions>(options);

    //! Todo: find a more general way to initialize hardware
    gmx_omp_nthreads_set(emntPairsearch, options.numThreads);
    gmx_omp_nthreads_set(emntNonbonded, options.numThreads);

    unpackTopologyToGmx(system_.topology());
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

void NbvSetupUtil::unpackTopologyToGmx(const Topology& topology)

{
    const std::vector<ParticleType>& particleTypes = topology.getParticleTypes();

    size_t numParticles = topology.numParticles();

    //! Todo: Refactor nbnxm to take this (nonbondedParameters_) directly
    //!
    //! initial self-handling of combination rules
    //! size: 2*(numParticleTypes^2)
    nonbondedParameters_.reserve(2 * particleTypes.size() * particleTypes.size());

    constexpr real c6factor  = 6.0;
    constexpr real c12factor = 12.0;

    for (const ParticleType& particleType1 : particleTypes)
    {
        real c6_1  = particleType1.c6() * c6factor;
        real c12_1 = particleType1.c12() * c12factor;
        for (const ParticleType& particleType2 : particleTypes)
        {
            real c6_2  = particleType2.c6() * c6factor;
            real c12_2 = particleType2.c12() * c12factor;

            real c6_combo = detail::combineNonbondedParameters(c6_1, c6_2, CombinationRule::Geometric);
            real c12_combo = detail::combineNonbondedParameters(c12_1, c12_2, CombinationRule::Geometric);
            nonbondedParameters_.push_back(c6_combo);
            nonbondedParameters_.push_back(c12_combo);
        }
    }

    particleInfoAllVdw_.resize(numParticles);
    for (size_t particleI = 0; particleI < numParticles; particleI++)
    {
        SET_CGINFO_HAS_VDW(particleInfoAllVdw_[particleI]);
        SET_CGINFO_HAS_Q(particleInfoAllVdw_[particleI]);
    }
}

static t_mdatoms setAtomData(SimulationState simulationState)
{
    t_mdatoms mdatoms;
    // We only use (read) the atom type and charge from mdatoms
    mdatoms.typeA = const_cast<int*>(simulationState.topology().getParticleTypeIdOfAllParticles().data());
    mdatoms.chargeA = const_cast<real*>(simulationState.topology().getCharges().data());
    return mdatoms;
}

void NbvSetupUtil::setAtomProperties(std::unique_ptr<nonbonded_verlet_t>& nbv, t_mdatoms& mdatoms)
{
    nbv->setAtomProperties(mdatoms, particleInfoAllVdw_);;
}

//! Sets up and returns a Nbnxm object for the given options and system
std::unique_ptr<nonbonded_verlet_t> NbvSetupUtil::setupNbnxmInstance(const Topology& topology, const NBKernelOptions& options)
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
                        topology.getParticleTypes().size(), nonbondedParameters_, 1, numThreads);

    setParticlesOnGrid(nbv, particleInfoAllVdw_, system_.coordinates(), system_.box());

    t_nrnb nrnb;
    nbv->constructPairlist(gmx::InteractionLocality::Local, topology.getGmxExclusions(), 0, &nrnb);

    return nbv;
}

std::unique_ptr<GmxForceCalculator> NbvSetupUtil::setupGmxForceCalculator()
{
    auto gmxForceCalculator_p = std::make_unique<GmxForceCalculator>(system_, options_);

    gmxForceCalculator_p->nbv_ = setupNbnxmInstance(system_.topology(), *options_);

    gmxForceCalculator_p->mdatoms_ = setAtomData(system_);
    setAtomProperties(gmxForceCalculator_p->nbv_, gmxForceCalculator_p->mdatoms_);

    // const PairlistSet& pairlistSet = nbv->pairlistSets().pairlistSet(gmx::InteractionLocality::Local);
    // const gmx::index numPairs = pairlistSet.natpair_ljq_ + pairlistSet.natpair_lj_ + pairlistSet.natpair_q_;
    // gmx_cycles_t cycles = gmx_cycles_read();

    const matrix& box_ = system_.box().legacyMatrix;

    gmxForceCalculator_p->forcerec_.nbfp = nonbondedParameters_;
    snew(gmxForceCalculator_p->forcerec_.shift_vec, SHIFTS);
    calc_shifts(box_, gmxForceCalculator_p->forcerec_.shift_vec);

    put_atoms_in_box(PbcType::Xyz, box_, system_.coordinates());

    gmxForceCalculator_p->verletForces_ =
            gmx::PaddedHostVector<gmx::RVec>(system_.topology().numParticles(), gmx::RVec(0, 0, 0));

    return gmxForceCalculator_p;
}

} // namespace nblib
