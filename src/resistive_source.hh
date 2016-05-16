#pragma once
#include "yee_grid.hh"

class ResistiveSource {
public:

    ResistiveSource(int ix, int iy, int iz, float resistance) {
        mPosX = ix;
        mPosY = iy;
        mPosZ = iz;
        mResistance  = resistance;
    }

    void calcCoefs(YeeGrid const& grid) {
        float delta_t = grid.delta_t;
        float delta_x = grid.delta_x;
        float delta_y = grid.delta_y;
        float delta_z = grid.delta_z;
        float epsilon = grid.epsilon_Ez.at(mPosX, mPosY, mPosZ);
        float subexpr = delta_t * delta_z /
                        (2 * mResistance * delta_x * delta_y);

        mCE = (1 - subexpr) / (1 + subexpr);
        mDE = (delta_t / epsilon) / (1 + subexpr);
    }

    void resqueFields(YeeGrid const& grid) {
        mPrevE = grid.Ez.at(mPosX, mPosY, mPosZ);
    }

    void updateFields(YeeGrid& grid, float voltage) const {
        int ix = mPosX;
        int iy = mPosY;
        int iz = mPosZ;
        float delta_x = grid.delta_x;
        float delta_y = grid.delta_y;

        float& curEz      = grid.Ez          .at(ix,   iy,   iz);
        float  curHy0     = grid.Hy          .at(ix,   iy,   iz);
        float  curHy1     = grid.Hy          .at(ix-1, iy,   iz);
        float  curHx0     = grid.Hx          .at(ix,   iy,   iz);
        float  curHx1     = grid.Hx          .at(ix,   iy-1, iz);

        curEz = mCE * mPrevE + mDE * ((curHy0 - curHy1) / delta_x -
                                     (curHx0 - curHx1) / delta_y +
                                     voltage / (mResistance * delta_x * delta_y));
    }

    float getC() const { return mCE; }
    float getD() const { return mDE; }

private:
    int mPosX, mPosY, mPosZ;
    float mCE, mDE;
    float mResistance;
    float mPrevE;
};
