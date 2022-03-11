//
// Created by ohmy on 2022-03-09.
//

#include "CGLms.h"


void CGLms::create_GGL_params(int degree) {
    if (degree == 1) {
        this->tau = DM( {-1,0,1});
        this->D = DM ({
                              {-0.500000000000000, 0.500000000000000, 0},
                              {-0.500000000000000, 0.500000000000000, 0}   });
        this->D = this->D(Slice(0,2));
        this->w = DM({1,1});
    }
    if (degree == 2) {
        this->tau = DM( {-1,0,1});
        this->D = DM ({
            {-1.50000000000000,	2.00000000000000,	-0.500000000000000},
            {-0.500000000000000,	0,	0.500000000000000},
            {0.500000000000000,	-2,	1.50000000000000} });
        this->w = DM({0.333333333333333,1.333333333333333,0.333333333333333});
    }
    if (degree == 3) {
        this->tau = DM( {-1,-0.500000000000000,0.500000000000000,1});
        this->D = DM ({
            {-3.16666666666667,	4.00000000000000,	-1.33333333333333,	0.500000000000000},
            {-1.00000000000000,	0.333333333333333,	1.00000000000000,	-0.333333333333333},
            {0.333333333333333,	-1.00000000000000,	-0.333333333333333,	1.00000000000000},
            {-0.500000000000000,	1.33333333333333,	-4.00000000000000,	3.16666666666667} });
        this->w = DM({0.111111111111111, 0.888888888888889, 0.888888888888889, 0.111111111111111});
    }
}


void CGLms::generate_constraints() {
    int nx = X.size1();
    double h = T/(N-1);
    J_ = 0;

    W = MX::vertcat({W, X(all,0)});
    for (int k = 0; k < N; ++k) {
        // Create collocation states
        auto o = ocp.variable(nx, n);


        auto xo = MX::horzcat({ X(all,k), o});
        W = MX::horzcat({W, o});
        MX u = MX::repmat(U(all, k), 1, n+1);
        // Defect constraints
        auto x_dot_ = f({{o}, {U(all, k)}});
        auto x_dot = MX::vertcat(x_dot_);
        auto F = 0.5*h  * x_dot ;
        auto g_d = mtimes(D(Slice(1,n+1), all), transpose(xo)) - transpose(F);
        // Shooting gap
        auto g_s = X(all,k+1) - o(all,n-1);

        g.push_back(g_s);
        g.push_back(MX::vec(g_d));



        // Cost
        auto l_k_ = J({{xo}, {U(all,k)}});
        auto l_k = MX::vertcat(l_k_);
        J_ += 0.5* h * mtimes(l_k , w); // quadrature
    }

}

MX CGLms::integrated_cost(MX t0, MX tf, int N) {
    return J_;
}

