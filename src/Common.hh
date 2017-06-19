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



template <typename T>
class XorMatrix3 {
    std::tuple<T&, T&, T&, T&, T&, T&> data;
public:

    XorMatrix3(std::nullptr_t, T& x01, T& x02, T& x10, std::nullptr_t, T& x12, T& x20, T& x21, std::nullptr_t)
        :data(x01, x02, x10, x12, x20, x21)
    {}

    template <int row, int col>
    T& get() {
        static_assert(row != col, "");
        static_assert(0 <= row && row <= 2, "");
        static_assert(0 <= col && col <= 2, "");

        return row == 0 && col == 1 ? std::ref(std::get<0>(data))
             : row == 0 && col == 2 ? std::ref(std::get<1>(data))
             : row == 1 && col == 0 ? std::ref(std::get<2>(data))
             : row == 1 && col == 2 ? std::ref(std::get<3>(data))
             : row == 2 && col == 0 ? std::ref(std::get<4>(data))
             : row == 2 && col == 1 ? std::ref(std::get<5>(data))
             : throw std::invalid_argument("");
    }

    template <int row, int col>
    T const& get() const {
        using ThisType = XorMatrix3<T>;
        return const_cast<T const&>(
               const_cast<ThisType const*>(this)->template get<row, col>());
    }


};

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
