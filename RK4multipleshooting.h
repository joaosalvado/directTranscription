//
// Created by ohmy on 2022-03-07.
//

#ifndef DIRECTTRANSCRIPTION_RK4MULTIPLESHOOTING_H
#define DIRECTTRANSCRIPTION_RK4MULTIPLESHOOTING_H

#include<casadi/casadi.hpp>
using namespace casadi;

class RK4multipleshooting {
public:
    int N;
    double T;
    casadi::MX &X;
    casadi::MX &U;
    casadi::Function f;
    casadi::Function J;
    Slice all;
    std::vector<MX> g;
    RK4multipleshooting(
            casadi::MX &X_, casadi::MX &U_, int N, int T,
            Function &f_, Function J_)
        : X(X_), U(U_), f(f_), J(J_) {
        this->N = N; this->T = T;
        generate_constraints();
    }
    void generate_constraints();
    MX integrated_cost(MX t0, MX tf, int N);
    MX rk4(const Function &f, const MX &dt, const MX &x, const MX &u);

    virtual ~RK4multipleshooting(){}

};


#endif //DIRECTTRANSCRIPTION_RK4MULTIPLESHOOTING_H
