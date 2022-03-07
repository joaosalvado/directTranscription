#include <iostream>
#include <casadi/casadi.hpp>
using namespace casadi;

#include "RK4multipleshooting.h"

int main() {
    // 1 - Problem setup

    // 1.1 - Params
    double T = 10;
    int N = 20;
    double L = 0.5;
    casadi::Opti ocp;
    DM x0 { 1, 1, 0 };
    DM xf { 2, 5, M_PI/2};

    // 1.2 - State: SE2
    struct se2{
        casadi::MX X; //Discrete states
        casadi::SX x_ode{SX::sym("x")};
        casadi::SX y_ode{SX::sym("y")};
        casadi::SX o_ode{SX::sym("o")};
        casadi::SX X_ode() { return casadi::SX::vertcat({x_ode, y_ode, o_ode}); };
    };
    auto x = std::make_shared<se2>();
    x->X =  ocp.variable(3, N+1) ;

    // 1.3 - Control: Speed and Angular velocity
    struct vw{
        casadi::MX U; // Discrete controls
        casadi::SX v_ode{SX::sym("v")};
        casadi::SX w_ode{SX::sym("w")};
        casadi::SX U_ode()  { return casadi::SX::vertcat({v_ode, w_ode}); };
    };
    auto u = std::make_shared<vw>();
    u->U = ocp.variable(2, N) ;

    // 1.4 - Dynamics Model
    SX X_dot = SX::vertcat(
            {u->v_ode * cos(x->o_ode),
             u->v_ode * sin(x->o_ode),
             (1 / L) * u->w_ode});
    auto f = Function("f", {x->X_ode(), u->U_ode()}, {X_dot});

    // 1.5 - Cost
    SX l = u->v_ode*u->v_ode + u->w_ode*u->w_ode;
    auto J = Function("l", {x->X_ode(), u->U_ode()}, {l});

    // 2 - Transcription Methods

    // 2.0 - Boundary Constraints
    Slice all;
    auto a = x->X.size1();
    auto b = x->X.size2();
    ocp.subject_to( x->X(all, 0 ) - x0 == 0);
    ocp.subject_to( x->X(all, N ) - xf == 0 );

    // 2.1 - Direct Local Collocation Multiple-shooting RK4
    RK4multipleshooting rk4ms = RK4multipleshooting(x->X, u->U, N, T, f, J);
    MX cost = rk4ms.integrated_cost(0, T, N);

    ocp.subject_to(rk4ms.g);
    ocp.minimize(cost);

    // 3 - Solve
    Dict p_opts, s_opts;
    p_opts["expand"] = true;
    std::string solver{"ipopt"};
    s_opts["print_level"] = 0;
    s_opts["linear_solver"] = "ma27";

    ocp.solver(solver, p_opts, s_opts);

    ocp.solve();

    return 0;
}
