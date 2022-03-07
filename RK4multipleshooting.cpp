//
// Created by ohmy on 2022-03-07.
//

#include "RK4multipleshooting.h"
using namespace casadi;

MX RK4multipleshooting::rk4(const Function &f, const MX &dt, const MX &x, const MX &u)
{
    auto k1 = f({{x}, {u}});
    auto k2 = f({{x + ((double)1 / 2) * dt * k1[0]}, {u}});
    auto k3 = f({{x + ((double)1 / 2) * dt * k2[0]}, {u}});
    auto k4 = f({{x + dt * k3[0]}, {u}});
    return  ((double)1 / 6) * dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
}

void RK4multipleshooting::generate_constraints() {
    for (int k = 0; k < N; ++k) {
        const auto &x_next
            =  X(all, k)
                    + rk4(
                            f,
                            (1 / (double) N) * T,
                            X(all, k),
                            U(all, k));
        g.push_back(X(all, k + 1) - x_next);
    }
}

MX RK4multipleshooting::integrated_cost(MX t0, MX tf, int N) {
    //Integrated cost
    MX _J = 0;
    for (int k = 0; k < N; ++k) {
        _J = _J + rk4(J, (tf - t0) / (double) N, X(all, k),U(all, k));
    }
    MX params = MX::vertcat({t0, tf});
    return _J;
}

