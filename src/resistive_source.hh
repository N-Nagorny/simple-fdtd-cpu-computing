#pragma once
#include "yee_grid.hh"
#include "rvlm/core/Constants.hh"

template <typename valueType>
class ResistiveSource {
public:
    using Const = rvlm::core::Constants<valueType>;

    typedef valueType                           ValueType;
    typedef rvlm::core::SolidArray3d<valueType> ArrayType;
    typedef typename ArrayType::IndexType       IndexType;

public:

    ResistiveSource(IndexType ix, IndexType iy, IndexType iz, ValueType r) {
        mPosX = ix;
        mPosY = iy;
        mPosZ = iz;
        mResistance  = r;
    }

    void calcCoefs(YeeGrid<ValueType> const& grid) {
        ValueType delta_t = grid.delta_t;
        ValueType delta_x = grid.delta_x;
        ValueType delta_y = grid.delta_y;
        ValueType delta_z = grid.delta_z;
        ValueType epsilon = grid.epsilon_Ez.at(mPosX, mPosY, mPosZ);
        ValueType subexpr = delta_t * delta_z /
                            (2 * mResistance * delta_x * delta_y);

        mCE = (1 - subexpr) / (1 + subexpr);
        mDE = (delta_t / epsilon) / (1 + subexpr);
    }

    void resqueFields(YeeGrid<ValueType> const& grid) {
        mPrevE = grid.Ez.at(mPosX, mPosY, mPosZ);
    }

    void updateFields(YeeGrid<ValueType>& grid, ValueType voltage) const {
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
    ValueType getResistance() const { return mResistance; }

private:
    IndexType mPosX, mPosY, mPosZ;
    ValueType mCE, mDE;
    ValueType mResistance;
    ValueType mPrevE;
};
