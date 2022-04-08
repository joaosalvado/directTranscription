//
// Created by ohmy on 2022-03-07.
//

#include "Plotter_dt.h"
#include "matplot/matplot.h"
using namespace matplot;

// note dt > t
std::vector<double> Plotter_dt::basis(
        double t, double dt, int n, std::vector<double> tau_v){
    auto phi = std::vector<double>(n+1,1);
    double tau = 2*t / dt - 1;
    for(int n_i = 0; n_i < n+1; ++n_i){
        for(int tau_i = 0; tau_i < n+1; ++tau_i) {
            if(tau_i == n_i) continue;
            phi[n_i] = phi[n_i] * (tau - tau_v[tau_i]) / (tau_v[n_i]- tau_v[tau_i]);
        }
    }
    return phi;
}

void Plotter_dt::plot_more_points_path(DM x, DM y, DM o, double Tf, int n, std::vector<double> tau_v){
    auto N = x.size2();
    auto dt = n*Tf / (N-1); // polynomial period

    // Linear Spaced time
    float space = Tf/(2*N);
    std::vector<float> t(2*N);
    std::generate(t.begin(), t.end(), [n = 0, &space]() mutable { return n++ * space; });

    std::vector<double> x_f, y_f, o_f;
    int pol_i = 0;
    for(auto t_ : t){
        if( t_ - pol_i*dt > dt){ pol_i++;}

        auto phi = basis(t_ - pol_i*dt, dt, n, tau_v);
        double x_t_ = 0, y_t_ = 0, o_t_ = 0;
        for( int n_i = 0; n_i < n+1; ++n_i){
            x_t_ = x_t_ + x(n_i + n*pol_i).scalar() * phi[n_i];
            y_t_ = y_t_ + y(n_i + n*pol_i).scalar() * phi[n_i];
            o_t_ = o_t_ + o(n_i + n*pol_i).scalar() * phi[n_i];
        }
        x_f.push_back(x_t_); y_f.push_back(y_t_); o_f.push_back(o_t_);

    }


    auto p = plot(x_f, y_f);
    p->line_width(1);
    p->marker(line_spec::marker_style::diamond);
    p->color("r");

    show();
}

void Plotter_dt::plot_path(DM x, DM y){
    auto p = plot(x.get_elements(), y.get_elements());
    p->line_width(2);
    p->marker(line_spec::marker_style::asterisk);
    matplot::hold(true);
    matplot::axis(matplot::equal);
    //show();
}

void Plotter_dt::plot_path_heading(DM x, DM y, DM o){
    int N = x.size2();
    double L= 0.025;
    auto x_head = std::vector<double>(N);
    auto y_head = std::vector<double>(N);
    for(int k = 0; k < N; ++k){
        x_head[k] = x(k).scalar() +  L*cos(o(k).scalar());
        y_head[k] = y(k).scalar() + L*sin(o(k).scalar());
/*        x_head[k] = x(k+1).scalar() - x(k).scalar();
        y_head[k] = y(k+1).scalar() - y(k).scalar();*/
        auto a = arrow(x(k).scalar(), y(k).scalar(), x_head[k], y_head[k]);
        a->line_width(2);
    }
    double x0 = x(0).scalar();
    double x1 = x(1).scalar();
    double y0 = y(0).scalar();
    double y1 = y(1).scalar();
    double a = L*cos(o(0).scalar());
    double b = L*sin(o(0).scalar());
    double c = L*cos(o(1).scalar());
    double d = L*sin(o(1).scalar());
/*    matplot::quiver(
            std::vector<double>{x(0).scalar(), x(1).scalar()},
            std::vector<double>{y(0).scalar(), y(1).scalar()},
            std::vector<double>{ 0.1, 0.2},
            std::vector<double>{0.3,  0.4});*/
/*    auto q = quiver(x.get_elements(), y.get_elements(), x_head, y_head);
    q->line_width(2);*/
    matplot::axis(matplot::equal);
    show();
}

void Plotter_dt::plot_var(DM var){

}