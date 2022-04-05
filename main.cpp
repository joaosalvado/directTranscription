#include <iostream>
#include <casadi/casadi.hpp>
using namespace casadi;

#include "RK4multipleshooting.h"
#include "LGLms.h"
#include "CGLms.h"
#include  "LLG.h"
#include "Plotter_dt.h"


double estimate_distance( std::vector<double> x0, std::vector<double> xf){
    double alpha = 0.5;
    double beta = 1.0;
    auto euclidean_dist =
            std::sqrt(std::pow(x0[0]-xf[0],2.0) + std::pow(x0[1]-xf[1], 2.0));
    auto angle_dist = std::abs(x0[2] -xf[2]);
    if(angle_dist > M_PI) angle_dist = 2*M_PI - angle_dist;
    return beta * (euclidean_dist + alpha  * angle_dist);
}

double estimate_T(double d, double v_std){
    double beta = 1.0;
    return beta * d / v_std;
}

double estimate_T1(int  N, double L, double v_max){
    double beta = 1.0;
    return N * L / (beta * v_max);
}

// 0.2 - State: SE2
struct se1{
    casadi::MX X; //Discrete states
    casadi::SX x_ode{SX::sym("x")};
    casadi::SX y_ode{SX::sym("y")};
    casadi::SX o_ode{SX::sym("o")};
    casadi::SX X_ode() { return casadi::SX::vertcat({x_ode, y_ode, o_ode}); };
};
// 0.3 - Control: Speed and Angular velocity
struct vw{
    casadi::MX U; // Discrete controls
    casadi::SX v_ode{SX::sym("v")};
    casadi::SX w_ode{SX::sym("w")};
    casadi::SX U_ode()  { return casadi::SX::vertcat({v_ode, w_ode}); };
};

void rh(std::vector<double> x0, std::vector <double> u0, std::vector<double> xf, double Tf,
std::vector<double> &xf_new, std::vector<double> &u0_new ){
    // 1.1 - Params
    int n = 3;


    int N = 10; //3*T for RK4
    double L = 0.5;
    double v_std = 0.5;
    double v_max = 1.0;
    casadi::Opti ocp;
    Slice all;
    double T = estimate_T1(3*N, L, v_max);

    auto x = std::make_shared<se1>();
    x->X =  ocp.variable(3, N+1) ;
    auto u = std::make_shared<vw>();
    u->U = ocp.variable(2, N) ;

    // 1.4 - Dynamics Model
    SX X_dot = SX::vertcat(
            {u->v_ode *(1 - x->o_ode*x->o_ode),
             u->v_ode * (2 * x->o_ode),
             u->w_ode});
    auto f = Function("f", {x->X_ode(), u->U_ode()}, {X_dot});

    // 1.5 - Cost
    SX l = (u->v_ode-v_std)*(u->v_ode-v_std) + u->w_ode*u->w_ode;
    auto J = Function("l", {x->X_ode(), u->U_ode()}, {l});

    // 1.6 - Bounds on controls
    double v_bound = v_max;
    double w_bound = 0.5;
    DM u_bound = DM::vertcat({v_bound, w_bound});
    for(int k = 0; k < N; ++k) {
        ocp.subject_to(u->U(all, k) <= w_bound);
        ocp.subject_to(u->U(all, k) >= -w_bound);
        ocp.subject_to(u->U(0,k) * ( 1 + x->X(2, k)*x->X(2,k)) <= v_bound);
        ocp.subject_to(u->U(0,k) * ( 1 + x->X(2, k)*x->X(2,k)) >= -v_bound);
    }
    // 1.7 - Boundary Constraints
    ocp.subject_to( x->X(all, 0 ) - x0 == 0);
    //ocp.subject_to( u->U(all, 0 ) - u0 == 0);

    // 2 - Transcription Methods
    // 2.2 - Direct Global Collocation Multiple-shooting LGL
    LGLms lgl_ms = LGLms(x->X, u->U, N, T, n, f, J, ocp);
    MX cost = lgl_ms.integrated_cost(0, T, N);

    for(auto g_i : lgl_ms.g){
        ocp.subject_to(g_i == 0 );
    }
    // cost = cost + mtimes(transpose(x->X(all, N ) - xf),(x->X(all, N ) - xf));

    // Second Part
    int n2 = 4;
    double T2 = (Tf-T);
    int N2 = 5;

    auto xend = ocp.variable(3,N2+1);
    auto uend = ocp.variable(2,N2);
    for(int k = 0; k < N2; ++k) {
        ocp.subject_to(uend(all, k) <= u_bound);
        ocp.subject_to(uend(all, k) >= -u_bound);
    }
    for(int k = 0; k < N2; ++k) {
        ocp.subject_to(xend(all, k) >= 0);
    }

    SX l2 = (u->v_ode-v_std)*(u->v_ode-v_std) + u->w_ode*u->w_ode;
/*    SX l2 = (x->x_ode-xf(0).scalar())*(x->x_ode-xf(0).scalar())
            + (x->y_ode-xf(1).scalar())*(x->y_ode-xf(1).scalar());*/
    auto J2 = Function("l", {x->X_ode(), u->U_ode()}, {l2});
    LGLms lgl_ms2 = LGLms(xend, uend, N2, T2, n2, f, J2, ocp);
    MX cost2 = lgl_ms2.integrated_cost(0, T2, N2);
    for(auto g_i : lgl_ms2.g){
        ocp.subject_to(g_i == 0 );
    }

    // Patch condition
    ocp.subject_to(xend(all, 0) - x->X(all,N) ==0);
    // End goal
    //ocp.subject_to( xend(all, N2) - xf == 0 );

    auto x_plot1 = lgl_ms.getStates();
    auto x_plot2 = lgl_ms2.getStates();
    auto x_plot = MX::horzcat({x_plot1, x_plot2});

    for(int k = 0; k < N-1; k++){
        cost = cost + mtimes(transpose(u->U(all, k+1) -  u->U(all, k)), u->U(all, k+1) -  u->U(all, k)) ;
        //cost = cost + sum1(sum2(u->U(all, k+1) -  u->U(all, k)));
        //cost = cost + mtimes(transpose(x->X(all, k ) - xf),(x->X(all, k ) - xf));
    }
    cost = cost + (1/10)*mtimes(transpose(x->X(all, N ) - xf),(x->X(all, N ) - xf));
    for(int k = 0; k < N2-1; k++){
        cost2 = cost2 + mtimes(transpose(uend(all, k+1) -  uend(all, k)), uend(all, k+1) -  uend(all, k)) ;
        //cost2 = cost2 + sum1(sum2(uend(all, k+1) -  uend(all, k)));
        cost2 = cost2 + mtimes(transpose(xend(all, k ) - xf),(xend(all, k ) - xf));
    }
    auto slicexy = Slice(0,2);
    cost2 = cost2 + mtimes(transpose(xend(all, N2 ) - xf),(xend(all, N2 ) - xf));
    //ocp.minimize((T*T)*cost+(T2*T2)*cost2);
    ocp.minimize((T)*cost+(T2)*cost2);



    // 3 - Solve
    Dict p_opts, s_opts;
    p_opts["expand"] = true;
    std::string solver{"ipopt"};
    s_opts["print_level"] = 0;
    s_opts["linear_solver"] = "ma27";
    ocp.solver(solver, p_opts, s_opts);

    auto solution = ocp.solve();

    DM Usol = solution.value(u->U);
    auto Nu = Usol.size2();
    u0_new.push_back(Usol(0, Nu-1).scalar());
    u0_new.push_back(Usol(1, Nu-1).scalar());

    x->X = MX::horzcat({x->X, xend});
    u->U = MX::horzcat({u->U, uend});

    Plotter_dt plotter;
    DM Xsol = solution.value(x->X);
    Usol = solution.value(u->U);
    std::cout << Xsol << std::endl;
    std::cout << Usol << std::endl;
    std::cout << solution.value(cost) << std::endl;

    // 4 - Plot Solution
    plotter.plot_path(Xsol(0,all), Xsol(1,all));
    Xsol = solution.value(x_plot1);
    plotter.plot_more_points_path(Xsol(0,all), Xsol(1,all), Xsol(2,all), T,n, lgl_ms.tau.get_elements());

    //plotter.plot_path_heading(Xsol(0,all), Xsol(1,all), Xsol(2, all));

    auto N_ = Xsol.size2();
    xf_new.push_back(Xsol(0, N_-1).scalar());
    xf_new.push_back(Xsol(1, N_-1).scalar());
    xf_new.push_back(Xsol(2, N_-1).scalar());

}


int main() {
    double v_std = 0.5;
    DM xf = DM::vertcat({ 1, 20, tan(0.9*M_PI/2)});
    DM x0 = DM::vertcat({ 1, 1, tan(0.9*M_PI/2)});

    // Iter 1
    std::vector<double> xf_new, uf_new;
    auto d = estimate_distance(x0.get_elements(), xf.get_elements());
    rh(x0.get_elements(), {0.0,0.0},xf.get_elements(), estimate_T(d,v_std), xf_new, uf_new );

    // Iter 2
    std::vector<double> xf_new2, uf_new2;
    d = estimate_distance(xf_new, xf.get_elements());
    rh(xf_new, uf_new,xf.get_elements(), estimate_T(d, v_std), xf_new2, uf_new2 );

    // Iter 3
    std::vector<double> xf_new3, uf_new3;
    d = estimate_distance(xf_new2, xf.get_elements());
    rh(xf_new2, uf_new2,xf.get_elements(), estimate_T(d, v_std), xf_new2, uf_new2 );
    return 0;
}
