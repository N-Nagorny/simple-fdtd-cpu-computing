#pragma once
#include <initializer_list>
#include "yee_grid.hh"
#include "ConvolutionalPML.hh"

template <typename Y>
class PmlShell {
public:

    PmlShell(YeeGrid<Y> *grid, int h,
             Dimensionless<Y> const& r,
             Dimensionless<Y> const& g)
        : mYeeGrid     (grid)
        , mParametersR (r, r, r)
        , mParametersG (g, g, g)
    {
        HalfOpenIndexRanges rangesX = splitRange(grid->numX, h);
        HalfOpenIndexRanges rangesY = splitRange(grid->numY, h);
        HalfOpenIndexRanges rangesZ = splitRange(grid->numZ, h);

        for (int nX: {0, 1, 2})
        for (int nY: {0, 1, 2})
        for (int nZ: {0, 1, 2}) {
            if (nX == 1 && nY == 1 && nZ == 1)
                continue;

            auto pml = new ConvolutionalPML<Y>(mYeeGrid,
                    make_triple(get(rangesX, nX),
                                get(rangesY, nY),
                                get(rangesZ, nZ)),
                    make_triple(numToDirection(nX),
                                numToDirection(nY),
                                numToDirection(nZ)),
                    mParametersR,
                    mParametersG);

            mPmls.push_back(std::unique_ptr<ConvolutionalPML<Y>>(pml));
        }
    }

    void setup() {
        for (auto& pml: mPmls)
            pml->setup();
    }

    void calcE() {
        for (auto& pml: mPmls)
            pml->calcE();
    }

    void calcH() {
        for (auto& pml: mPmls)
            pml->calcH();
    }

private:

    static Triple<HalfOpenIndexRange> splitRange(int n, int h) {
        return { {1, 1 + h}, {1+h, n-h-1}, {n-h-1, n-1} };
    }

    static constexpr Optional<AxialDirection> numToDirection(int nRange) {
        return nRange == 0 ? AxialDirection::negative
             : nRange == 1 ? Optional<AxialDirection>()
             : nRange == 2 ? AxialDirection::positive
             : throw std::invalid_argument("");
    }

    YeeGrid<Y> *mYeeGrid;
    Triple<Dimensionless<Y>> mParametersR, mParametersG;
    std::vector<std::unique_ptr<ConvolutionalPML<Y>>> mPmls;
};
