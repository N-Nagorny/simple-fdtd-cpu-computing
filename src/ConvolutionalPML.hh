#pragma once
#include <iostream>
#include "Constants.hh"
#include "BoostUnitsHelpers.hh"
#include "Common.hh"
#include "Dimensions.hh"
#include "yee_grid.hh"
#include <boost/units/cmath.hpp>

template <typename Y>
class ConvolutionalPML {
public:

    using Const = Constants<Y>;

    template <typename YY>
    using Array = rvlm::core::SolidArray3d<YY, Index>;

    ConvolutionalPML(YeeGrid<Y> *grid,
                     HalfOpenIndexRanges         const& presentCells,
                     Triple<Optional<AxialDirection>> const& gradientDirection,
                     Triple<Dimensionless<Y>>    const& reflectionCoeffs,
                     Triple<Dimensionless<Y>>    const& parametersG)
        : mYeeGrid(grid)
        , mPresentCells(presentCells)
        , mGradientDirection(gradientDirection)
        , mReflectionCoeffs(reflectionCoeffs)
        , mParametersG(parametersG)
        , mFieldExy(presentCells, ElectricIntensity<Y>::from_value(0))
        , mFieldExz(presentCells, ElectricIntensity<Y>::from_value(0))
        , mFieldEyx(presentCells, ElectricIntensity<Y>::from_value(0))
        , mFieldEyz(presentCells, ElectricIntensity<Y>::from_value(0))
        , mFieldEzx(presentCells, ElectricIntensity<Y>::from_value(0))
        , mFieldEzy(presentCells, ElectricIntensity<Y>::from_value(0))
        , mFieldHxy(presentCells, MagneticIntensity<Y>::from_value(0))
        , mFieldHxz(presentCells, MagneticIntensity<Y>::from_value(0))
        , mFieldHyx(presentCells, MagneticIntensity<Y>::from_value(0))
        , mFieldHyz(presentCells, MagneticIntensity<Y>::from_value(0))
        , mFieldHzx(presentCells, MagneticIntensity<Y>::from_value(0))
        , mFieldHzy(presentCells, MagneticIntensity<Y>::from_value(0))
        {}

    void setup() {
        Optional<AxialDirection> gradX, gradY, gradZ;
        std::tie(gradX, gradY, gradZ) = mGradientDirection;

        if (gradX) {
            precalculateCoefficients<0>();
            fillSigmas<0>(*gradX);
        }

        if (gradY) {
            precalculateCoefficients<1>();
            fillSigmas<1>(*gradY);
        }

        if (gradZ) {
            precalculateCoefficients<2>();
            fillSigmas<2>(*gradZ);
        }
    }

    void calcH() {

    }

    void calcE() {

    }

private:

    template <int AxisNum>
    ElectricConductivity<Y> conductivityProfile(Index itwice) {
        auto const& condt0 = std::get<AxisNum>(mCachedConductivity0);
        auto const& condt1 = std::get<AxisNum>(mCachedConductivity1);
        auto const&      g = std::get<AxisNum>(mParametersG);

        if (itwice <= 0)
            return condt0;

        using boost::units::pow;

        Dimensionless<Y> exponent = Dimensionless<Y>(itwice)/Dimensionless<Y>(2);
        Dimensionless<Y> powVal   = pow(g, exponent);
        return condt1 * powVal;
    }

    template <int AxisNum>
    void precalculateCoefficients() {
        auto const& g = std::get<AxisNum>(mParametersG);
        auto const& r = std::get<AxisNum>(mReflectionCoeffs);
        auto const& d = std::get<AxisNum>(mYeeGrid->spatialSteps);
        auto const& N = std::get<AxisNum>(mPresentCells).size();
        auto & condt0 = std::get<AxisNum>(mCachedConductivity0);
        auto & condt1 = std::get<AxisNum>(mCachedConductivity1);


        using boost::units::log;
        using boost::units::sqrt;
        using boost::units::pow;

        Dimensionless<Y> logR   (( log(r) ));
        Dimensionless<Y> logG   (( log(g) ));
        Dimensionless<Y> sqrtG  (( usqrt(g) ));
        Dimensionless<Y> powG   (( pow(g, Dimensionless<Y>(N)) ));

        ElectricConductivity<Y> sigma0 {
            - Const::C() * Const::EPS_0() / (Dimensionless<Y>(2) * d)
                    * logG / (powG - Dimensionless<Y>(1.0)) * logR
        };

        condt0 = sigma0 * (sqrtG - Dimensionless<Y>(1)) / logG;
        condt1 = sigma0 * (g - Dimensionless<Y>(1)) / (sqrtG * logG);
    }

    template <int Axis>
    void fillSigmas(AxialDirection dir) {
        constexpr int A0 = Axis;
        constexpr int A1 = (A0 + 1) % 3;
        constexpr int A2 = (A0 + 2) % 3;

        Index const s2 = std::get<A2>(mPresentCells).start;
        auto const& arrEps    = std::get<A0>(mYeeGrid->epsilons);
        auto const& arrMu     = std::get<A0>(mYeeGrid->mus);
        auto & arrSigmaE = std::get<A0>(mYeeGrid->sigmasE);
        auto & arrSigmaH = std::get<A0>(mYeeGrid->sigmasH);

        Index itwiceE;
        Index itwiceH;
        Index itwiceStep;

        if (dir == AxialDirection::positive) {
            itwiceE = 0;
            itwiceH = 1;
            itwiceStep = 2;
        } else {
            Index n = std::get<A0>(mPresentCells).size();
            itwiceH = 2*(n - 1);
            itwiceE = itwiceH - 1;
            itwiceStep = -2;
        }

        for (auto i0: std::get<A0>(mPresentCells)) {
            (void)i0;
            RVLM_FDTD_DASSERT(itwiceE >= 0);
            RVLM_FDTD_DASSERT(itwiceH >= 0);

            ElectricConductivity<Y> cdty0 = conductivityProfile<A0>(itwiceE);
            ElectricConductivity<Y> cdty1 = conductivityProfile<A0>(itwiceH);
            std::cerr << "PROFILE for axis " << A0 << " [" << itwiceE << "] " << cdty0 << '\n';

            for (auto i1: std::get<A1>(mPresentCells)) {
                (void)i1;

                for (auto i2: std::get<A2>(mPresentCells)) {
                    (void)i2;

                    Index cx, cy, cz;
                    auto cur1 = arrEps.template getCursorX<A0,A1,A2>(i0, i1, i2);
                    arrEps.cursorCoordinates(cur1, cx, cy, cz);
                    auto cur2 = arrEps.getCursor(cx, cy, cz);
                    if (cur1 != cur2)
                        std::cerr << "SHIT\n";

                    auto const& eps = arrEps   .template at<A0,A1,A2>(i0, i1, i2);
                    auto const& mu  = arrMu    .template at<A0,A1,A2>(i0, i1, i2);
                    auto & sigmaE   = arrSigmaE.template at<A0,A1,A2>(i0, i1, i2);
                    auto & sigmaH   = arrSigmaH.template at<A0,A1,A2>(i0, i1, i2);

                    sigmaE = cdty0;
                    sigmaH = cdty1 * mu / eps;
                }
            }

            itwiceE += itwiceStep;
            itwiceH += itwiceStep;
        }
    }

    YeeGrid<Y> *mYeeGrid;

    HalfOpenIndexRanges
        mPresentCells;

    Triple<Optional<AxialDirection>>
        mGradientDirection;

    Triple<Dimensionless<Y>>
        mReflectionCoeffs,
        mParametersG;

    Triple<ElectricConductivity<Y>>
        mCachedConductivity0,
        mCachedConductivity1;

    Array<ElectricIntensity<Y>>
        mFieldExy, mFieldExz,
        mFieldEyx, mFieldEyz,
        mFieldEzx, mFieldEzy;

    Array<MagneticIntensity<Y>>
        mFieldHxy, mFieldHxz,
        mFieldHyx, mFieldHyz,
        mFieldHzx, mFieldHzy;
};
