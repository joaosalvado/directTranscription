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
    show();
}

void Plotter::plot_path_heading(DM x, DM y, DM o){
    int N = x.size2();
    double L= 10;
    auto x_head = std::vector<double>(N);
    auto y_head = std::vector<double>(N);
    for(int k = 0; k < N; ++k){
        x_head[k] =  L*cos(o(k).scalar());
        y_head[k] =  L*sin(o(k).scalar());
    }
    matplot::quiver(
            std::vector<double>{x(0).scalar(), x(1).scalar()},
            std::vector<double>{y(0).scalar(), y(1).scalar()},
            std::vector<double>{ L*cos(o(0).scalar()),  L*cos(o(1).scalar())},
            std::vector<double>{L*sin(o(0).scalar()),  L*sin(o(1).scalar())});
   /* auto q = quiver(x.get_elements(), y.get_elements(), x_head, y_head);
    q->line_width(2);*/
    show();
}

void Plotter::plot_var(DM var){

}