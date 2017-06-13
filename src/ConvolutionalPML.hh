#pragma once
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
                     HalfOpenIndexRanges      const& presentCells,
                     Triple<int>              const& gradientDirection,
                     Triple<Dimensionless<Y>> const& reflectionCoeffs,
                     Triple<Dimensionless<Y>> const& parametersG)
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
        int gradX, gradY, gradZ;
        std::tie(gradX, gradY, gradZ) = mGradientDirection;

        //if (gradX > 0)
        //    fillSigmasTowardXmax();

        //if (gradX < 0)
        //    fillSigmasTowardXmin();

        // ...

        precalculateCoefficients<0>();
        fillSigmasForward<0>();
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
    void fillSigmasForward() {
        constexpr int A0 = Axis;
        constexpr int A1 = (A0 + 1) % 3;
        constexpr int A2 = (A0 + 2) % 3;

        Index const s0 = std::get<A0>(mPresentCells).start;
        Index const s1 = std::get<A1>(mPresentCells).start;
        Index const s2 = std::get<A2>(mPresentCells).start;
        auto const& arrEps    = std::get<A0>(mYeeGrid->epsilons);
        auto const& arrMu     = std::get<A0>(mYeeGrid->mus);
        auto const& arrSigmaE = std::get<A0>(mYeeGrid->sigmasE);
        auto const& arrSigmaH = std::get<A0>(mYeeGrid->sigmasH);
        auto curEps    = arrEps   .template getCursor<A0,A1,A2>(s0, s1, s2);
        auto curMu     = arrMu    .template getCursor<A0,A1,A2>(s0, s1, s2);
        auto curSigmaE = arrSigmaE.template getCursor<A0,A1,A2>(s0, s1, s2);
        auto curSigmaH = arrSigmaH.template getCursor<A0,A1,A2>(s0, s1, s2);

        Index itwice = 0;
        for (auto i0: std::get<A0>(mPresentCells)) {
            (void)i0;

            ElectricConductivity<Y> cdty0 = conductivityProfile<A0>(itwice);
            ElectricConductivity<Y> cdty1 = conductivityProfile<A0>(itwice+1);

            for (auto i1: std::get<A1>(mPresentCells)) {
                (void)i1;

                for (auto i2: std::get<A2>(mPresentCells)) {
                    (void)i2;

                    auto const& eps = arrEps   .at(curEps);
                    auto const& mu  = arrMu    .at(curMu);
                    auto & sigmaE   = arrSigmaE.at(curSigmaE);
                    auto & sigmaH   = arrSigmaH.at(curSigmaH);

                    sigmaE = cdty0;
                    sigmaH = cdty1 * mu / eps;

                    arrEps   .template cursorMoveToNext<A2>(curEps);
                    arrMu    .template cursorMoveToNext<A2>(curMu);
                    arrSigmaE.template cursorMoveToNext<A2>(curSigmaE);
                    arrSigmaH.template cursorMoveToNext<A2>(curSigmaH);
                }

                arrEps   .template cursorMoveToNext<A1>(curEps);
                arrMu    .template cursorMoveToNext<A1>(curMu);
                arrSigmaE.template cursorMoveToNext<A1>(curSigmaE);
                arrSigmaH.template cursorMoveToNext<A1>(curSigmaH);
            }

            itwice += 2;

            arrEps   .template cursorMoveToNext<A0>(curEps);
            arrMu    .template cursorMoveToNext<A0>(curMu);
            arrSigmaE.template cursorMoveToNext<A0>(curSigmaE);
            arrSigmaH.template cursorMoveToNext<A0>(curSigmaH);
        }
    }

    YeeGrid<Y> *mYeeGrid;

    HalfOpenIndexRanges
        mPresentCells;

    Triple<int>
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
