#pragma once
#include <cmath>
#include <boost/units/get_system.hpp>
#include <boost/units/get_dimension.hpp>
#include <boost/units/static_rational.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/dimension.hpp>

/*
 * This file provides fixed versions for boost::units::pow<N, M>()
 * and boost::units::sqrt().
 */

template<typename Unit, int N, int M = 1>
using static_power_helper = boost::units::unit<
        typename boost::units::static_power<
                typename boost::units::get_dimension<Unit>::type,
                typename boost::units::static_rational<N, M>::type>::type,
        typename boost::units::get_system<Unit>::type>;

template<int N, typename Unit, typename Y>
auto ratpow(boost::units::quantity<Unit, Y> const &x)
        -> boost::units::quantity<static_power_helper<Unit, N>, Y> {

    using std::pow;

    #pragma ide diagnostic ignored "OCDFAInspection"
    return decltype(ratpow<N>(x))::from_value(pow(x.value(), N));
}

template<int N, int M, typename Unit, typename Y>
auto ratpow(boost::units::quantity<Unit, Y> const &x)
        -> boost::units::quantity<static_power_helper<Unit, N, M>, Y> {

    using std::pow;

    #pragma ide diagnostic ignored "OCDFAInspection"
    return decltype(ratpow<N, M>(x))::from_value(pow(x.value(), (Y) N / M));
}

template<typename Unit, typename Y>
auto usqrt(boost::units::quantity<Unit, Y> const &x)
        -> decltype(ratpow<1, 2>(x)) {

    return ratpow<1, 2>(x);
}

template<typename Unit, typename Y>
auto reciprocal(boost::units::quantity<Unit, Y> const &x)
        -> decltype(ratpow<-1>(x)) {

    return ratpow<-1>(x);
}

