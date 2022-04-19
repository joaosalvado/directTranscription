#include <iostream>
#include <casadi/casadi.hpp>
using namespace casadi;

#include "RK4multipleshooting.h"
#include "LGLms.h"
#include "CGLms.h"
#include  "LLG.h"
#include "Plotter_dt.h"


int N_init = 20; // 3*T for RK4
double L = 0.5;
double v_std = 0.5;
double v_max = 1.0;


double estimate_T(double d, double v_std){
    double beta = 1;
    return beta * d / v_std;
}

double estimate_T1(int  N, double L, double v_max){
    double beta = 1.0;
    return N * L / (beta * v_max);
}

double estimate_N( double T1, double L, double v_max){
    double beta = 1.0;
    return beta * T1 * v_max / L;
}


double T1 = estimate_T1(3*N_init, L, v_max);

double estimate_distance( std::vector<double> x0, std::vector<double> xf){
    double alpha = 0;
    double beta = 1.0;
    auto euclidean_dist =
            std::sqrt(std::pow(x0[0]-xf[0],2.0) + std::pow(x0[1]-xf[1], 2.0));
    auto angle_dist = std::abs(x0[2] -xf[2]);
    /*double o_std = atan2(xf[1]-x0[1], xf[0] - x0[0]);
    double d = std::sqrt((xf[1]-x0[1])*(xf[1]-x0[1]) +  (xf[0]-x0[0]) *(xf[0]-x0[0])) - v_std*T1;
    double tg2_std = sin(o_std)/(1+cos(o_std));
    auto angle_dist = x0[2] - tg2_std;
    alpha = d;*/
    //if(angle_dist > M_PI) angle_dist = 2*M_PI - angle_dist;
    return beta * (euclidean_dist + alpha  * angle_dist);
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
std::vector<double> &xf_new, std::vector<double> &u0_new , int final = 0){
    // 1.1 - Params
    int n =  3;

    casadi::Opti ocp;
    double T;
    int N;
    if(final){
        T = T1; //std::max(2*Tf, 2.0);
        N = N_init; //std::max(estimate_N(T,L, v_max), 3.0);
    } else{
        T = T1;
        N = N_init;
    }
    Slice all;

    auto x = std::make_shared<se1>();
    x->X =  ocp.variable(3, N+1) ;
    auto u = std::make_shared<vw>();
    u->U = ocp.variable(2, N) ;

    // 1.4 - Dynamics Model
    SX X_dot = SX::vertcat(
            {u->v_ode *(1 - x->o_ode*x->o_ode)/(1 + x->o_ode*x->o_ode),
             u->v_ode * (2 * x->o_ode)/(1 + x->o_ode*x->o_ode),
             u->w_ode});
    auto f = Function("f", {x->X_ode(), u->U_ode()}, {X_dot});

    // 1.5 - Cost
    SX l;
    Function J;
    if(final == 1 || final == 2){
        l = (u->v_ode)*(u->v_ode) + u->w_ode*u->w_ode;
        J = Function("l", {x->X_ode(), u->U_ode()}, {l});
    } else{
        l = (u->v_ode-v_std)*(u->v_ode-v_std) + u->w_ode*u->w_ode;
        J = Function("l", {x->X_ode(), u->U_ode()}, {l});
    }

    // 1.6 - Bounds on controls
    double v_bound = v_max;
    double w_bound = 0.5;
    DM u_bound = DM::vertcat({v_bound, w_bound});
    for(int k = 0; k < u->U.size2(); ++k) {
        ocp.subject_to(u->U(1, k) <= w_bound);
        ocp.subject_to(u->U(1, k) >= -w_bound);
        //ocp.subject_to(2*u->U(1,k) - (1 + x->X(2,k)*x->X(2,k))*(T)*w_bound < 0);
        //ocp.subject_to(-2*u->U(1,k) + (1 + x->X(2,k)*x->X(2,k))*(T)*(-w_bound) < 0);
        ocp.subject_to(u->U(0,k) <= v_bound);
        ocp.subject_to(u->U(0,k) >= -0.5*v_bound);
        //ocp.subject_to(u->U(0,k) * ( 1 + x->X(2, k)*x->X(2,k)) <= v_bound);
        //ocp.subject_to(u->U(0,k) * ( 1 + x->X(2, k)*x->X(2,k)) >= -0.5*v_bound);
    }
    // 1.7 - Boundary Constraints
    ocp.subject_to( x->X(all, 0 ) - x0 == 0);
    //ocp.subject_to( u->U(1, 0 ) - u0[0] == 0);
    //ocp.subject_to(u->U(1,u->U.size2()-1) == 0);

    // 2 - Transcription Methods
    // 2.2 - Direct Global Collocation Multiple-shooting LGL
    LGLms lgl_ms = LGLms(x->X, u->U, N, T, n, f, J, ocp);
    MX cost = lgl_ms.integrated_cost(0, T, N);

    for(auto g_i : lgl_ms.g){
        ocp.subject_to(g_i == 0 );
    }
    // cost = cost + mtimes(transpose(x->X(all, N ) - xf),(x->X(all, N ) - xf));

    // Second Part
/*    int n2 = 3;
    double T2 = Tf-T;
    int N2 = N;*/

    MX xend, uend, x_plot1;
    x_plot1 = lgl_ms.getStates();
    if(final == 0) { // patching polygon
/*        xend = ocp.variable(3, N2 + 1);
        uend = ocp.variable(2, N2);
        for (int k = 0; k < N2; ++k) {
            ocp.subject_to(uend(all, k) <= u_bound);
            ocp.subject_to(uend(all, k) >= -u_bound);
        }*/
//        for (int k = 0; k < N2; ++k) {
//            ocp.subject_to(xend(all, k) >= 0);
//        }

/*        SX l2 = (u->v_ode - v_std) * (u->v_ode - v_std) + u->w_ode * u->w_ode;
*//*    SX l2 = (x->x_ode-xf(0).scalar())*(x->x_ode-xf(0).scalar())
            + (x->y_ode-xf(1).scalar())*(x->y_ode-xf(1).scalar());*//*
        auto J2 = Function("l", {x->X_ode(), u->U_ode()}, {T2*l2});
        LGLms lgl_ms2 = LGLms(xend, uend, N2, T2, n2, f, J2, ocp);
        MX cost2 = lgl_ms2.integrated_cost(0, T2, N2);*/
/*        for (auto g_i: lgl_ms2.g) {
            ocp.subject_to(g_i == 0);
        }*/

        // Patch condition
/*        ocp.subject_to(xend(all, 0) - x->X(all, N) == 0);*/
        // End goal
        //ocp.subject_to( xend(all, N2) - xf == 0 );
/*
        auto x_plot2 = lgl_ms2.getStates();*/
        auto x_plot = MX::horzcat({x_plot1});
        for( int k = 0; k < u->U.size2()-1; ++k){
            //cost = cost + mtimes(transpose(u->U(all, k + 1) - u->U(all, k)), u->U(all, k + 1) - u->U(all, k));
        }
        for (int k = 0; k < x->X.size2() - 1; k++) {
            //cost = cost + sum1(sum2(u->U(all, k+1) -  u->U(all, k)));
            //cost = cost + mtimes(transpose(x->X(all, k ) - xf),(x->X(all, k ) - xf));
            //cost = cost + mtimes(transpose(x->X(Slice(0,2), k) - DM({xf[0],xf[1]}) ), (x->X(Slice(0,2), k) - DM({xf[0],xf[1]})));
        }
        cost = cost + (1/Tf)*mtimes(transpose(x->X(Slice(0,2), x->X.size2()-1) - DM({xf[0],xf[1]})), (x->X(Slice(0,2), x->X.size2()-1) - DM({xf[0],xf[1]}) ));
        double o_std = atan2(xf[1]-x0[1], xf[0] - x0[0]);
        double d = std::sqrt((xf[1]-x0[1])*(xf[1]-x0[1]) +  (xf[0]-x0[0]) *(xf[0]-x0[0])) - v_max*T;
        double tg2_std = sin(o_std)/(1+cos(o_std));
        if(d > 0) {
            cost = cost + (1/Tf)*d*mtimes(transpose(x->X(2, x->X.size2() - 1) - tg2_std), (x->X(2, x->X.size2() - 1) - tg2_std));
        }
        /*for (int k = 0; k < xend.size2() - 1; k++) {
            //cost2 = cost2 + mtimes(transpose(uend(all, k + 1) - uend(all, k)), uend(all, k + 1) - uend(all, k));
            //cost2 = cost2 + sum1(sum2(uend(all, k+1) -  uend(all, k)));
            cost2 = cost2 + mtimes(transpose(xend(Slice(0,2), k) - DM({xf[0],xf[1]}) ), (xend(Slice(0,2), k) - DM({xf[0],xf[1]})));
            //cost2 = cost2 + mtimes(transpose(xend(all, k) - xf ), (xend( all, k) - xf));
        }
        auto slicexy = Slice(0, 2);
        cost2 = cost2 + mtimes(transpose(xend(Slice(0,2),xend.size2()-1) - DM({xf[0],xf[1]}) ), (xend(Slice(0,2), xend.size2()-1) - DM({xf[0],xf[1]})));*/
        //cost2 = cost2 + mtimes(transpose(xend(all, xend.size2()-1) - xf), (xend(all, xend.size2()-1) - xf));
        //ocp.minimize((T*T)*cost+(T2*T2)*cost2);
        ocp.minimize( cost );
    }

    if(final == 1){ // Final Trajectory
        ocp.subject_to(u->U(0,u->U.size2()-1) == 0);
        //ocp.subject_to(x->X(0, x->X.size2()-1) - xf[0]== 0);
        //ocp.subject_to(x->X(2, x->X.size2()-1) - xf[2] <=  0.1);
        //ocp.subject_to(x->X(2, x->X.size2()-1) - xf[2] >=  0.1);
        //ocp.subject_to(x->X(1, x->X.size2()-1) - xf[1] <=  0.1);
        //ocp.subject_to(x->X(1, x->X.size2()-1) - xf[1] >=  0.1);
        //ocp.subject_to(mtimes(transpose(x->X(all, x->X.size2()-1) - xf), (x->X(all, x->X.size2()-1) - xf)) <= 0.5);
        //ocp.subject_to( mtimes(transpose(x->X(Slice(0,2), x->X.size2()-1) - DM({xf[0],xf[1]})), (x->X(Slice(0,2), x->X.size2()-1) - DM({xf[0],xf[1]}) )) <= 0.2 );
        for (int k = 0; k < x->X.size2() - 1; k++) {
            //cost = cost + mtimes(transpose(u->U(all, k + 1) - u->U(all, k)), u->U(all, k + 1) - u->U(all, k));
            //cost = cost + sum1(sum2(u->U(all, k+1) -  u->U(all, k)));
            //cost = cost + mtimes(transpose(x->X(all, k ) - xf),(x->X(all, k ) - xf));
            //cost = cost + mtimes(transpose(x->X(Slice(0,2), k ) - DM({xf[0],xf[1]})),(x->X(Slice(0,2), k ) - DM({xf[0],xf[1]})));
        }
        //cost = cost + mtimes(transpose(x->X(Slice(0,2), x->X.size2()-1) - DM({xf[0],xf[1]})), (x->X(Slice(0,2), x->X.size2()-1) - DM({xf[0],xf[1]}) ));
        cost = cost + mtimes(transpose(x->X(all, x->X.size2()-1) - xf), (x->X(all, x->X.size2()-1) - xf));
        double o_std = xf[2];
        double tg2_std = xf[2];//sin(o_std)/(1+cos(o_std));
        //cost = cost + mtimes(transpose(x->X(2, x->X.size2()-1) - tg2_std), (x->X(2, x->X.size2()-1) - tg2_std ));
        ocp.minimize(cost );
    }

    if(final ==2){

        ocp.subject_to(u->U(0,u->U.size2()-1) == 0);
        ocp.subject_to( mtimes(transpose(x->X(Slice(0,2), x->X.size2()-1) - DM({xf[0],xf[1]})), (x->X(Slice(0,2), x->X.size2()-1) - DM({xf[0],xf[1]}) )) <= 0.2 );
        //ocp.subject_to(mtimes(transpose(x->X(all, x->X.size2()-1) - xf), (x->X(all, x->X.size2()-1) - xf)) <= 0.1);
        for (int k = 0; k < x->X.size2() - 1; k++) {
            //cost = cost + mtimes(transpose(u->U(all, k + 1) - u->U(all, k)), u->U(all, k + 1) - u->U(all, k));
            //cost = cost + sum1(sum2(u->U(all, k+1) -  u->U(all, k)));
            //cost = cost + mtimes(transpose(x->X(all, k ) - xf),(x->X(all, k ) - xf));
            //cost = cost + mtimes(transpose(x->X(Slice(0,2), k ) - DM({xf[0],xf[1]})),(x->X(Slice(0,2), k ) - DM({xf[0],xf[1]})));
        }
        //cost = cost + mtimes(transpose(x->X(Slice(0,2), x->X.size2()-1) - DM({xf[0],xf[1]})), (x->X(Slice(0,2), x->X.size2()-1) - DM({xf[0],xf[1]}) ));
        //cost = cost + mtimes(transpose(x->X(all, x->X.size2()-1) - xf), (x->X(all, x->X.size2()-1) - xf));
        //double o_std = atan2(xf[1]-x0[1], xf[0] - x0[0]);
        //double tg2_std = sin(o_std)/(1+cos(o_std));
        double tg2_std = xf[2];
        cost = cost + T*mtimes(transpose(x->X(2, x->X.size2()-1) - tg2_std), (x->X(2, x->X.size2()-1) - tg2_std ));
        ocp.minimize(cost );
    }


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

    if(!final){
        x->X = MX::horzcat({x->X, xend});
        u->U = MX::horzcat({u->U, uend});
    };

    Plotter_dt plotter;
    DM Xsol = solution.value(x->X);
    Usol = solution.value(u->U);

    // Trajectory solution
    std::cout << Xsol << std::endl;
    std::cout << Usol << std::endl;
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

void rh_full(std::vector<double> x0, std::vector<double> xf){
    double T_prev, T_curr, d_prev, d_curr;
    double v_prev = v_std, v_curr, v_avg;
    std::vector<double> xf_new, uf_new;
    std::vector<double> xf_prev = x0, uf_prev = { 0.0, 0.0 };
    bool first = true, final = false;
    while(true) {
        d_curr = estimate_distance(xf_prev, xf);
        if(first){
            first = false;
            T_curr = estimate_T(d_curr, v_std);
        } else{
            v_prev = d_prev/T_prev;
            v_curr = (d_prev-d_curr)/T1;
            v_avg = 0.5*v_prev + 0.5*v_curr;
            T_curr = d_curr / v_avg;
        }

        xf_new.clear();
        uf_new.clear();
        if(T_curr < T1 || final) {
            if(!final){
                rh(xf_prev, uf_prev, xf, T1, xf_new, uf_new, 1);
            } else{
                if(d_curr < 0.5 && std::abs(xf_new[2]-xf[2]) < 0.1) break;
                rh(xf_prev, uf_prev, xf, T1, xf_new, uf_new, 2);
                break;
            }
            final = true;
            //break; // final trajectory
        } else {
            rh(xf_prev, uf_prev, xf, T_curr, xf_new, uf_new);
        }
        d_prev = d_curr;
        T_prev = T_curr;
        xf_prev = xf_new;
        uf_prev = uf_new;
    }
    // Final trajectory
    //rh(xf_prev, uf_prev, xf, T1, xf_new, uf_new, true);
}


int main() {
    DM x0 = DM::vertcat({ 10, 10, tan(0*M_PI/2)});
    DM xf = DM::vertcat({ 30, 30, tan(-0.5*M_PI/2)});

    rh_full( x0.get_elements(), xf.get_elements());

    // Iter 1
    //   std::vector<double> xf_new, uf_new;
//    auto d1 = estimate_distance(x0.get_elements(), xf.get_elements());
//    auto T_total_1 = estimate_T(d1, v_std);
//    rh(x0.get_elements(), {0.0,0.0},xf.get_elements(), T_total_1, xf_new, uf_new );
//
//    // Iter 2
//    std::vector<double> xf_new2, uf_new2;
//    auto d2 = estimate_distance(xf_new, xf.get_elements());
//    auto v_1 = d1/T_total_1;
//    auto v_2 = (d1 - d2)/T1;
//    auto v_avg = 0.5*v_1 + 0.5*v_2;
//    auto T_total_2 = d2 / v_avg;
//    rh(xf_new, uf_new,xf.get_elements(), T_total_2, xf_new2, uf_new2 );
//
//    // Iter 3
//    std::vector<double> xf_new3, uf_new3;
//    auto d3  = estimate_distance(xf_new2, xf.get_elements());
//    auto v_3 = (d2 - d3)/T1;
//    v_avg = 0.5*v_2 + 0.5*v_3;
//    auto T_total_3 = d3 / v_avg;
//    rh(xf_new2, uf_new2,xf.get_elements(), T_total_3, xf_new3, uf_new3 , true );
//    return 0;
}


