#pragma once
#include "yee_grid.hh"
#include "Constants.hh"
#include "Dimensions.hh"

template <typename T>
class ResistiveSource {
public:
    using Const = rvlm::core::Constants<T>;

    typedef T                           ValueType;
    typedef rvlm::core::SolidArray3d<T> ArrayType;
    typedef typename ArrayType::IndexType       IndexType;

public:

    ResistiveSource(IndexType ix, IndexType iy, IndexType iz, Resistance<T> r) {
        mPosX = ix;
        mPosY = iy;
        mPosZ = iz;
        mResistance  = r;
    }

    void calcCoefs(YeeGrid<T> const& grid) {
        auto const delta_t = grid.delta_t;
        auto const delta_x = grid.delta_x;
        auto const delta_y = grid.delta_y;
        auto const delta_z = grid.delta_z;
        auto const epsilon = grid.epsilon_Ez.at(mPosX, mPosY, mPosZ);
        auto const sigmaE  = grid.sigma_Ez.at(mPosX, mPosY, mPosZ);
        auto const subexpr = delta_t * delta_z * sigmaE /
                            ((T)2 * mResistance * epsilon * delta_x * delta_y);

        mCE = ((T)1 - subexpr) / ((T)1 + subexpr);

        mDE = ElectricCurlCoefficient<T>();
        // TODO: !!!
        // mDE = (delta_t / epsilon) / ((T)1 + subexpr);
    }

    void resqueFields(YeeGrid<T> const& grid) {
        mPrevE = grid.Ez.at(mPosX, mPosY, mPosZ);
    }

    void updateFields(YeeGrid<T>& grid, ElectricPotential<T> voltage) const {
        IndexType ix  = mPosX;
        IndexType iy  = mPosY;
        IndexType iz  = mPosZ;
        ValueType delta_x = grid.delta_x;
        ValueType delta_y = grid.delta_y;

        ValueType& curEz      = grid.Ez          .at(ix,   iy,   iz);
        ValueType  curHy0     = grid.Hy          .at(ix,   iy,   iz);
        ValueType  curHy1     = grid.Hy          .at(ix-1, iy,   iz);
        ValueType  curHx0     = grid.Hx          .at(ix,   iy,   iz);
        ValueType  curHx1     = grid.Hx          .at(ix,   iy-1, iz);

        curEz = mCE * mPrevE + mDE * ((curHy0 - curHy1) / delta_x -
                                  (curHx0 - curHx1) / delta_y +
                                  voltage / (mResistance * delta_x * delta_y));
    }

    ValueType getC() const { return mCE; }
    ValueType getD() const { return mDE; }
    Resistance<T> getResistance() const { return mResistance; }

private:
    IndexType mPosX, mPosY, mPosZ;
    T mCE;
    ElectricCurlCoefficient<T> mDE;

    Resistance<T> mResistance;
    ElectricIntensity<ValueType> mPrevE;
};
