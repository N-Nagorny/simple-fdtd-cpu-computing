#pragma once
#include <cstdint>
#include <tuple>
#include <boost/optional.hpp>
#include "rvlm/core/HalfOpenRange.hh"
#include "rvlm/core/SolidArray3d.hh"

#define RVLM_FDTD_DASSERT(expr)

template <typename T>
using Optional = boost::optional<T>;

enum class AxialDirection {
    positive,
    negative,
};

using Index = std::ptrdiff_t;

template <typename T>
using Triple = std::tuple<T, T, T>;

template <typename T>
Triple<T> make_triple(T const& val1, T const& val2, T const& val3) {
    return std::tuple<T, T, T>(val1, val2, val3);
}

template <typename T>
T get(Triple<T> const& triple, Index idx) {
    return idx == 0 ? std::get<0>(triple)
         : idx == 1 ? std::get<1>(triple)
         : idx == 2 ? std::get<2>(triple)
         : throw std::invalid_argument("");
}



using HalfOpenIndexRange = rvlm::core::HalfOpenRange<Index>;

using HalfOpenIndexRanges = Triple<HalfOpenIndexRange>;

inline
Triple<Index> getStarts(HalfOpenIndexRanges const& ranges) {
    return make_triple(std::get<0>(ranges).start,
                       std::get<1>(ranges).start,
                       std::get<0>(ranges).start);
}

template <typename T>
using Array = rvlm::core::SolidArray3d<T, Index>;
