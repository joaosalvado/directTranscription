#include <iostream>
#include <casadi/casadi.hpp>
using namespace casadi;

#include "RK4multipleshooting.h"
#include "LGLms.h"
#include "CGLms.h"
#include  "LLG.h"
#include "Plotter_dt.h"

int main() {
    // 1 - Problem setup

    // 1.1 - Params
    int n = 3;
    double T = 10;
    int N = 20; //3*T for RK4
    double L = 0.2;
    casadi::Opti ocp;
    DM x0 = DM::vertcat({ 2, 2, -M_PI + 0.25*M_PI });
    DM xf = DM::vertcat({ 1, 1, -M_PI/2});
    Slice all;

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
    //SX l = u->v_ode*u->v_ode + u->w_ode*u->w_ode;
    SX l = u->v_ode*u->v_ode + u->w_ode*u->w_ode;
    auto J = Function("l", {x->X_ode(), u->U_ode()}, {l});

    // 1.6 - Bounds on controls
    double v_bound = 1;
    double w_bound = 2;
    DM u_bound = DM::vertcat({v_bound, w_bound});
    for(int k = 0; k < N; ++k) {
        ocp.subject_to(u->U(all, k) <= u_bound);
        ocp.subject_to(u->U(all, k) >= -u_bound);
    }

    // 1.7 - Boundary Constraints
    ocp.subject_to( x->X(all, 0 ) - x0 == 0);
    //ocp.subject_to( x->X(all, N ) - xf == 0 );

    // Cost
    MX cost;
    for(int k = 0; k < N-1; k++){
        //cost = cost + mtimes(transpose(u->U(all, k+1) -  u->U(all, k)), u->U(all, k+1) -  u->U(all, k)) ;
        cost = cost + sum1(sum2(u->U(all, k+1) -  u->U(all, k)));
    }
    // 2 - Transcription Methods



    // 2.1 - Direct Local Collocation Multiple-shooting RK4
//    RK4multipleshooting rk4ms = RK4multipleshooting(x->X, u->U, N, T, f, J);
//    MX cost = rk4ms.integrated_cost(0, T, N);
//
//    for(auto g_i : rk4ms.g){
//        ocp.subject_to(g_i == 0 );
//    }
//    cost = cost + mtimes(transpose(x->X(all, N ) - xf),(x->X(all, N ) - xf));
//
//    ocp.minimize(cost);
    // 2.2 - Direct Global Collocation Multiple-shooting LGL
    LGLms lgl_ms = LGLms(x->X, u->U, N, T, n, f, J, ocp);
    cost = lgl_ms.integrated_cost(0, T, N);

    for(auto g_i : lgl_ms.g){
        ocp.subject_to(g_i == 0 );
    }
    cost = cost + mtimes(transpose(x->X(all, N ) - xf),(x->X(all, N ) - xf));
    ocp.minimize(cost);

    // 2.3 - Direct Global Collocation Multiple-shooting CGL
/*    auto cgl_ms = CGLms(x->X, u->U, N, T, n, f, J, ocp);
    MX cost = cgl_ms.integrated_cost(0, T, N);

    for(auto g_i : cgl_ms.g){
        ocp.subject_to(g_i == 0 );
    }

    //cost = cost + mtimes(transpose(x->X(all, N ) - xf),(x->X(all, N ) - xf));
    ocp.minimize(cost);*/

    // 2.4 - Direct Global Collocation Multiple-shooting LLG
/*    LLG llg_ms = LLG(x->X, u->U, N, T, n, f, J, ocp);
    MX cost = llg_ms.integrated_cost(0, T, N);

    for(auto g_i : llg_ms.g){
        ocp.subject_to(g_i == 0 );
    }
    ocp.minimize(cost);*/


    // 3 - Solve
    Dict p_opts, s_opts;
    p_opts["expand"] = true;
    std::string solver{"ipopt"};
    s_opts["print_level"] = 0;
    s_opts["linear_solver"] = "ma27";
    ocp.solver(solver, p_opts, s_opts);

    auto solution = ocp.solve();
    DM Xsol = solution.value(x->X);
    DM Usol = solution.value(u->U);
    std::cout << Xsol << std::endl;
    std::cout << Usol << std::endl;
    std::cout << solution.value(cost) << std::endl;

    // 4 - Plot Solution
    Plotter_dt plotter;
    plotter.plot_path(Xsol(0,all), Xsol(1,all));
    x->X = lgl_ms.getStates();
    Xsol = solution.value(x->X);
    plotter.plot_path_heading(Xsol(0,all), Xsol(1,all), Xsol(2, all));
    return 0;
}
