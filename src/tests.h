#pragma once
#include "focus.h"
#include <memory>

class GaussianMean {
public:
    std::unique_ptr<Info> info;

    GaussianMean() {
        // newP for Gaussian
        auto newP = [](double St, int tau, double m0){
            std::unique_ptr<Piece> p = std::make_unique<PieceGau>();
            p->St = St;
            p->tau = tau;
            p->m0 = m0;
            return p;
        };
        info = std::make_unique<Info>(newP, NAN);
    }

    void update(double y) {
        info->update(y);
    }

    double statistic() const {
        return info->statistic();
    }
};