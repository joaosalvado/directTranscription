//
// Created by ohmy on 2022-03-09.
//

#ifndef DIRECTTRANSCRIPTION_CGLMS_H
#define DIRECTTRANSCRIPTION_CGLMS_H

#include <casadi/casadi.hpp>
using namespace casadi;

class CGLms {
public:
    int N;
    double T;
    // For a Gauss polynomial of degree n, and N shooting intervals ones has:
    // X (nx * N)
    // U (nu * N)
    // Tau ( nx * n+1 * N )
    // Let nx be the number of states and nu of controls, e.g. se2 nx = 3 (x,y,o)
    casadi::MX &X; // State Shooting vars
    casadi::MX &U; // Control vars
    casadi::MX W; // States + collocation states
    casadi::MX Tau; // Collocation nodes vars
    casadi::Function f;
    casadi::Function J;
    Slice all;
    std::vector<MX> g;
    Opti &ocp;


    // Legendre-Gauss-Lobato parameters
    DM tau; // Collocation nodes in [-1, 1]
    DM D;   // Differentiation matrix
    DM w;   // Weight quadrature
    int n;  // Polynomial degree
    MX J_;

    CGLms(casadi::MX &X_, casadi::MX &U_, int N, int T, int n,
          Function &f_, Function J_, Opti & ocp_)
    : X(X_), U(U_), f(f_), J(J_), ocp(ocp_) {
        this->N = N; this->T = T; this->n = n;
        create_GGL_params(n);
        generate_constraints();
    }

    void create_GGL_params(int degree);
    MX integrated_cost(MX t0, MX tf, int N);
    void generate_constraints();
    MX getStates(){return W;};

    virtual ~CGLms(){}
};


#endif //DIRECTTRANSCRIPTION_CGLMS_H
