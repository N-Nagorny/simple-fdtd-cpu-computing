#pragma once
#include "yee_grid.hh"
#include "Constants.hh"
#include "Dimensions.hh"
#include "rvlm/core/SolidArray3d.hh"

template <typename Y>
class CurrentSource {
public:

    using Const     = Constants<Y>;
    using ValueType = Y;
    using ArrayType = rvlm::core::SolidArray3d<Y>;

    CurrentSource(int ix, int iy, int iz) {
        mPosX = ix;
        mPosY = iy;
        mPosZ = iz;
    }

    void calcCoefs(YeeGrid<Y> const& grid) {
        auto const delta_t = grid.delta_t;
        auto const delta_x = grid.delta_x;
        auto const delta_y = grid.delta_y;
        auto const epsilon = grid.epsilon_Ez.at(mPosX, mPosY, mPosZ);

        mCoeff = delta_t / (epsilon * delta_x * delta_y);
    }

    void updateFields(YeeGrid<Y>& grid, Current<Y> const& current) const {
        int ix  = mPosX;
        int iy  = mPosY;
        int iz  = mPosZ;

        grid.Ez.at(ix, iy, iz) -= mCoeff * current;
    }

private:

    using CoefficientUnit = boost::units::divide_typeof_helper<
            boost::units::si::resistance,
            boost::units::si::length>::type;

    int mPosX, mPosY, mPosZ;
    boost::units::quantity<CoefficientUnit, Y> mCoeff;
};
