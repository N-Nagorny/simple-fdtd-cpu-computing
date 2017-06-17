#pragma once
#include <iostream>
#include <boost/units/cmath.hpp>
#include "Constants.hh"
#include "BoostUnitsHelpers.hh"
#include "Common.hh"
#include "Dimensions.hh"
#include "yee_grid.hh"

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
        , mFieldCax(presentCells, Dimensionless<Y>::from_value(0))
        , mFieldCay(presentCells, Dimensionless<Y>::from_value(0))
        , mFieldCaz(presentCells, Dimensionless<Y>::from_value(0))
        , mFieldDax(presentCells, Dimensionless<Y>::from_value(0))
        , mFieldDay(presentCells, Dimensionless<Y>::from_value(0))
        , mFieldDaz(presentCells, Dimensionless<Y>::from_value(0))
        , mFieldCbx(presentCells, ElectricCurlCoefficient<Y>::from_value(0))
        , mFieldCby(presentCells, ElectricCurlCoefficient<Y>::from_value(0))
        , mFieldCbz(presentCells, ElectricCurlCoefficient<Y>::from_value(0))
        , mFieldDbx(presentCells, MagneticCurlCoefficient<Y>::from_value(0))
        , mFieldDby(presentCells, MagneticCurlCoefficient<Y>::from_value(0))
        , mFieldDbz(presentCells, MagneticCurlCoefficient<Y>::from_value(0))
        , mFieldsCa(mFieldCax, mFieldCay, mFieldCaz)
        , mFieldsCb(mFieldCbx, mFieldCby, mFieldCbz)
        , mFieldsDa(mFieldDax, mFieldDay, mFieldDaz)
        , mFieldsDb(mFieldDbx, mFieldDby, mFieldDbz)
        , mFieldsPmlE( {mFieldExy,   mFieldExy, mFieldExz},
                       {mFieldEyx, mFieldExy,   mFieldEyz},
                       {mFieldEzx, mFieldEzy, mFieldExy})
        , mFieldsPmlH( {mFieldHxy,   mFieldHxy, mFieldHxz},
                       {mFieldHyx, mFieldHxy,   mFieldHyz},
                       {mFieldHzx, mFieldHzy, mFieldHxy})
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

    void precalculatePmlCoefficients() {
        precalculatePmlCoefficients<0>();
        precalculatePmlCoefficients<1>();
        precalculatePmlCoefficients<2>();
    }

    template <int Axis>
    void precalculatePmlCoefficients() {
        auto const& arrEps    = std::get<Axis>(mYeeGrid->epsilons);
        auto const& arrMu     = std::get<Axis>(mYeeGrid->mus);
        auto const& arrSigmaE = std::get<Axis>(mYeeGrid->sigmasE);
        auto const& arrSigmaH = std::get<Axis>(mYeeGrid->sigmasH);
        auto      & arrCa     = std::get<Axis>(mFieldsCa);
        auto      & arrCb     = std::get<Axis>(mFieldsCb);
        auto      & arrDa     = std::get<Axis>(mFieldsDa);
        auto      & arrDb     = std::get<Axis>(mFieldsDb);

        Time<Y> const& deltaT = mYeeGrid->delta_t;
        Dimensionless<Y> one = Y(1);
        Dimensionless<Y> two = Y(2);
        for (Index ix: std::get<0>(mPresentCells))
        for (Index iy: std::get<1>(mPresentCells))
        for (Index iz: std::get<2>(mPresentCells)) {
            auto const& eps    = arrEps    .at(ix, iy, iz);
            auto const& mu     = arrMu     .at(ix, iy, iz);
            auto const& sigmaE = arrSigmaE .at(ix, iy, iz);
            auto const& sigmaH = arrSigmaH .at(ix, iy, iz);
            auto      & Ca     = arrCa     .at(ix, iy, iz);
            auto      & Cb     = arrCb     .at(ix, iy, iz);
            auto      & Da     = arrDa     .at(ix, iy, iz);
            auto      & Db     = arrDb     .at(ix, iy, iz);

            auto epsSubexpr = one + sigmaE * deltaT / (two * eps);
            auto muSubexpr  = one + sigmaH * deltaT / (two * mu);
            Ca = (one - sigmaE * deltaT / (two * eps) ) / epsSubexpr;
            Da = (one - sigmaH * deltaT / (two * mu) )  / muSubexpr;
            Cb = deltaT / eps / epsSubexpr;
            Db = deltaT / mu  / muSubexpr;
        }
    }

    template <int Axis0, int Axis1>
    void calcPmlH() {
        static_assert(Axis0 != Axis1, "");
        constexpr int A0 = Axis0;
        constexpr int A1 = Axis1;
        constexpr int A2 = 3 - (Axis1 + Axis0);

        auto const& arrDa = std::get<A0>(mFieldsDa);
        auto const& arrDb = std::get<A1>(mFieldsDb);
        auto const& arrE  = std::get<A2>(mYeeGrid->fieldsE);
        auto const& delta = std::get<A1>(mYeeGrid->spatialSteps);
        auto      & arrH  = std::get<A0>(std::get<A1>(mFieldsPmlH));

        for (Index i0: std::get<A0>(mPresentCells))
        for (Index i1: std::get<A1>(mPresentCells))
        for (Index i2: std::get<A2>(mPresentCells)) {
            auto const& Da = arrDa.template at<A0,A1,A2>(i0, i1, i2);
            auto const& Db = arrDb.template at<A0,A1,A2>(i0, i1, i2);
            auto const& E  = arrE .template at<A0,A1,A2>(i0, i1, i2);
            auto const& E1 = arrE .template at<A0,A1,A2>(i0, i1+1, i2);
            auto      & H  = arrH .template at<A0,A1,A2>(i0, i1, i2);

            H = H * Da - Db / delta * (E1 - E);
        }
    }

    void calcH() {
        calcPmlH<0, 1>();
        calcPmlH<0, 2>();
        sumPmlFieldsH<0>();
        calcPmlH<1, 0>();
        calcPmlH<1, 2>();
        sumPmlFieldsH<1>();
        calcPmlH<2, 0>();
        calcPmlH<2, 1>();
        sumPmlFieldsH<1>();
    }

    template <int Axis0, int Axis1>
    void calcPmlE() {
        static_assert(Axis0 != Axis1, "");
        constexpr int A0 = Axis0;
        constexpr int A1 = Axis1;
        constexpr int A2 = 3 - (Axis1 + Axis0);

        auto const& arrCa = std::get<A0>(mFieldsCa);
        auto const& arrCb = std::get<A1>(mFieldsCb);
        auto const& arrH  = std::get<A2>(mYeeGrid->fieldsH);
        auto const& delta = std::get<A1>(mYeeGrid->spatialSteps);
        auto      & arrE  = std::get<A0>(std::get<A1>(mFieldsPmlE));

        for (Index i0: std::get<A0>(mPresentCells))
        for (Index i1: std::get<A1>(mPresentCells))
        for (Index i2: std::get<A2>(mPresentCells)) {
            auto const& Ca = arrCa.template at<A0,A1,A2>(i0, i1, i2);
            auto const& Cb = arrCb.template at<A0,A1,A2>(i0, i1, i2);
            auto const& H  = arrH .template at<A0,A1,A2>(i0, i1, i2);
            auto const& H1 = arrH .template at<A0,A1,A2>(i0, i1-1, i2);
            auto      & E  = arrE .template at<A0,A1,A2>(i0, i1, i2);

            E = E * Ca - Cb / delta * (H - H1);
        }
    }

    template <int Axis>
    void sumPmlFieldsE() {
        constexpr int A0 = Axis;
        constexpr int A1 = (Axis + 1) % 3;
        constexpr int A2 = (Axis + 2) % 3;

        auto const& arrE01 = std::get<A0>(std::get<A1>(mFieldsPmlE));
        auto const& arrE02 = std::get<A0>(std::get<A2>(mFieldsPmlE));
        auto      & arrE0  = std::get<A0>(mYeeGrid->fieldsE);

        for (Index i0: std::get<A0>(mPresentCells))
        for (Index i1: std::get<A1>(mPresentCells))
        for (Index i2: std::get<A2>(mPresentCells)) {
            auto const& E01 = arrE01.template at<A0,A1,A2>(i0,i1,i2);
            auto const& E02 = arrE02.template at<A0,A1,A2>(i0,i1,i2);
            auto      & E0  = arrE0 .template at<A0,A1,A2>(i0,i1,i2);

            E0 = E01 + E02;
        }
    }

    template <int Axis>
    void sumPmlFieldsH() {
        constexpr int A0 = Axis;
        constexpr int A1 = (Axis + 1) % 3;
        constexpr int A2 = (Axis + 2) % 3;

        auto const& arrH01 = std::get<A0>(std::get<A1>(mFieldsPmlH));
        auto const& arrH02 = std::get<A0>(std::get<A2>(mFieldsPmlH));
        auto      & arrH0  = std::get<A0>(mYeeGrid->fieldsH);

        for (Index i0: std::get<A0>(mPresentCells))
        for (Index i1: std::get<A1>(mPresentCells))
        for (Index i2: std::get<A2>(mPresentCells)) {
            auto const& H01 = arrH01.template at<A0,A1,A2>(i0,i1,i2);
            auto const& H02 = arrH02.template at<A0,A1,A2>(i0,i1,i2);
            auto      & H0  = arrH0 .template at<A0,A1,A2>(i0,i1,i2);

            H0 = H01 + H02;
        }
    }

    void calcE() {
        calcPmlE<0, 1>();
        calcPmlE<0, 2>();
        sumPmlFieldsE<0>();
        calcPmlE<1, 0>();
        calcPmlE<1, 2>();
        sumPmlFieldsE<1>();
        calcPmlE<2, 0>();
        calcPmlE<2, 1>();
        sumPmlFieldsE<2>();
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

    Array<Dimensionless<Y>>
        mFieldDax, mFieldDay, mFieldDaz;

    Array<MagneticCurlCoefficient<Y>>
        mFieldDbx, mFieldDby, mFieldDbz;

    Array<Dimensionless<Y>>
        mFieldCax, mFieldCay, mFieldCaz;

    Array<ElectricCurlCoefficient<Y>>
        mFieldCbx, mFieldCby, mFieldCbz;

    Triple<Triple<Array<MagneticIntensity<Y>>&>> mFieldsPmlH;
    Triple<Triple<Array<ElectricIntensity<Y>>&>> mFieldsPmlE;

    Triple<Array<Dimensionless<Y>>&> mFieldsCa;
    Triple<Array<ElectricCurlCoefficient<Y>>&> mFieldsCb;
    Triple<Array<Dimensionless<Y>>&> mFieldsDa;
    Triple<Array<MagneticCurlCoefficient<Y>>&> mFieldsDb;
};
