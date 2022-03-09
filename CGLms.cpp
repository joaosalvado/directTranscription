//
// Created by ohmy on 2022-03-09.
//

#include "CGLms.h"


void CGLms::create_GGL_params(int degree) {
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

    for (int k = 0; k < N; ++k) {
        // Create collocation states
        auto o = ocp.variable(nx, n);


        auto xo = MX::horzcat({ X(all,k), o});
        MX u = MX::repmat(U(all, k), 1, n+1);
        // Defect constraints
        auto x_dot_ = f({{o}, {U(all, k)}});
        auto x_dot = MX::vertcat(x_dot_);
        std::cout << x_dot.size1() << std::endl;
        std::cout << x_dot.size2() << std::endl;
        auto F = 0.5*h  * x_dot ;
        auto g_d = mtimes(D(Slice(1,n+1), all), transpose(xo)) - transpose(F);
        std::cout << g_d.size1() << std::endl;
        std::cout << g_d.size2() << std::endl;
        auto h = MX::vec(g_d);
        std::cout << h.size1() << std::endl;
        std::cout << h.size2() << std::endl;
        // Shooting gap
        auto g_s = X(all,k+1) - o(all,n-1);

        g.push_back(g_s);
        g.push_back(MX::vec(g_d));



        // Cost
        auto l_k_ = J({{xo}, {U(all,k)}});
        auto l_k = MX::vertcat(l_k_);
        std::cout << l_k.size1() << std::endl;
        std::cout << l_k.size2() << std::endl;
        std::cout << w.size1() << std::endl;
        std::cout << w.size2() << std::endl;
        J_ += MX::sum2(l_k) ;//0.5* h * mtimes(l_k , w); // quadrature
    }

}

MX CGLms::integrated_cost(MX t0, MX tf, int N) {
    return J_;
}

