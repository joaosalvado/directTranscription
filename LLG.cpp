//
// Created by ohmy on 2022-03-11.
//

#include "LLG.h"
//
// Created by ohmy on 2022-03-08.
//


void LLG::generate_constraints() {
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
        auto F = h  * x_dot ;
        auto g_d = mtimes(transpose(D), transpose(xo)) - transpose(F);
        // Shooting gap
        auto X_end = mtimes(xo,E);
        auto g_s = X(all,k+1) - X_end;

        g.push_back(g_s);
        g.push_back(MX::vec(g_d));

        // Cost
        auto l_k_ = J({{o}, {U(all,k)}});
        auto l_k = MX::vertcat(l_k_);
        J_ += h * mtimes(l_k , w); // quadrature
    }

}

MX LLG::integrated_cost(MX t0, MX tf, int N) {
    return J_;
}

void LLG::create_LLG_params(int degree){
    tau = casadi::collocation_points(degree, "legendre");
    casadi::collocation_coeff(tau.get_elements(), D, E, w);
}



