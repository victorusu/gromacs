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
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef GROMACS_UTIL_H
#define GROMACS_UTIL_H

#include <functional>
#include <tuple>
#include <type_traits>
#include <vector>

#include "gromacs/math/vectypes.h"

namespace nblib
{

//! generate Velocites from a Maxwell Boltzmann distro, masses should be the
//! same as the ones specified for the Topology object
std::vector<gmx::RVec> generateVelocity(real Temperature, unsigned int seed, std::vector<real> const& masses);

bool checkNumericValues(const std::vector<gmx::RVec>& values);

inline void ignore_unused()
{
}

template<class T, class... Ts>
inline void ignore_unused(T& x, Ts&... xs)
{
    static_cast<void>(x);
    ignore_unused(xs...);
}

namespace detail
{

template<class T, class U>
struct disjunction : std::integral_constant<bool, T::value || U::value>
{
};

template<class T>
struct void_t
{
    typedef void type;
};

template<class T, class Enable = void>
struct HasTypeMember
{
    typedef T type;
};

template<class T>
struct HasTypeMember<T, typename void_t<typename T::type>::type>
{
    typedef typename T::type type;
};

template<int N, typename T, typename Tuple>
struct CompareField :
    disjunction<std::is_same<T, typename std::tuple_element<N, Tuple>::type>,
                std::is_same<T, typename HasTypeMember<typename std::tuple_element<N, Tuple>::type>::type>>
{
};

template<int N, class T, class Tuple, bool Match = false>
struct MatchingField
{
    static decltype(auto) get(Tuple& tp)
    {
        // check next element
        return MatchingField<N + 1, T, Tuple, CompareField<N + 1, T, Tuple>::value>::get(tp);
    }
};

template<int N, class T, class Tuple>
struct MatchingField<N, T, Tuple, true>
{
    static decltype(auto) get(Tuple& tp) { return std::get<N>(tp); }
};

} // namespace detail

//! Function to return the element in Tuple whose type matches T
//! Note: if there are more than one, the first occurrence will be returned
template<typename T, typename Tuple>
decltype(auto) pickType(Tuple& tup)
{
    return detail::MatchingField<0, T, Tuple, detail::CompareField<0, T, Tuple>::value>::get(tup);
}

template<class... Ts>
struct TypeList
{
};

template<template<class...> class P, class L>
struct Map_
{
};

template<template<class...> class P, template<class...> class L, class... Ts>
struct Map_<P, L<Ts...>>
{
    typedef TypeList<P<Ts>...> type;
};

template<template<class...> class P, class L>
using Map = typename Map_<P, L>::type;


template<template<class...> class P, class L>
struct Reduce_
{
};

template<template<class...> class P, template<class...> class L, class... Ts>
struct Reduce_<P, L<Ts...>>
{
    typedef P<Ts...> type;
};

template<template<class...> class P, class L>
using Reduce = typename Reduce_<P, L>::type;

} // namespace nblib

namespace std17
{

// std::apply and std::invoke are part of C++17
// this code can be removed once the project moves to C++17

namespace detail
{

// stripped down version of std::invoke, does not handle class member functions
template<class F, class... Args>
constexpr decltype(auto) invoke(F&& f, Args&&... args)
{
    // call f with all arguments
    return std::forward<F>(f)(std::forward<Args>(args)...);
}

template<class F, class Tuple, size_t... Is>
constexpr decltype(auto) apply_impl(F&& f, Tuple&& t, std::index_sequence<Is...>)
{
    // unpack the tuple t into a function parameter pack and forward to invoke
    return invoke(std::forward<F>(f), std::get<Is>(std::forward<Tuple>(t))...);
}

} // namespace detail

// call f with all tuple elements as arguments
template<class F, class Tuple>
constexpr decltype(auto) apply(F&& f, Tuple&& t)
{
    // generate an index sequence to access all tuple elements
    return detail::apply_impl(std::forward<F>(f), std::forward<Tuple>(t),
                              std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>{}>{});
}

} // namespace std17

namespace nblib
{

// calls func with each element in tuple_
template<class F, class... Ts>
decltype(auto) for_each_tuple(F&& func, std::tuple<Ts...>& tuple_)
{
    // TODO: change return type to void and add C++17 [[ maybe_unused ]] attribute
    std17::apply([f = func](auto&... args) { return std::initializer_list<int>{ (f(args), 0)... }; }, tuple_);
}

// applies func to each element in tuple_ and returns result
template<class F, class... Ts>
decltype(auto) transform_tuple(F&& func, const std::tuple<Ts...>& tuple_)
{
    return std17::apply([f = func](auto&... args) { return std::make_tuple(f(args)...); }, tuple_);
}

// binary fold operations are natively supported in C++17
// and this helpers will become obsolete
template<class Op, class T1, class T2>
auto binary_fold(Op&& op, T1&& a, T2&& b)
{
    return op(a, b);
}

template<class Op, class First, class... Ts>
auto binary_fold(Op&& op, First&& first, Ts&&... args)
{
    return op(first, binary_fold(op, args...));
}

} // namespace nblib

#endif // GROMACS_UTIL_H
