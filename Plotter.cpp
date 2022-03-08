//
// Created by ohmy on 2022-03-07.
//

#include "Plotter.h"
#include "matplot/matplot.h"
using namespace matplot;
void Plotter::plot_path(DM x, DM y){
    auto p = plot(x.get_elements(), y.get_elements());
    p->line_width(2);
    p->marker(line_spec::marker_style::asterisk);
    matplot::axis(matplot::equal);
    show();
}

void Plotter::plot_path_heading(DM x, DM y, DM o){
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

void Plotter::plot_var(DM var){

}