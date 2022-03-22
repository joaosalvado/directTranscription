//
// Created by ohmy on 2022-03-07.
//

#ifndef DIRECTTRANSCRIPTION_PLOTTER_DT_H
#define DIRECTTRANSCRIPTION_PLOTTER_DT_H

#include <casadi/casadi.hpp>
using namespace casadi;

class Plotter_dt {
public:
    std::vector<std::function<void( std::vector<double> )>> plot_robot;
    void plot_path(DM x, DM y);
    void plot_path_heading(DM x, DM y, DM o);
    void plot_var(DM var);
};


#endif //DIRECTTRANSCRIPTION_PLOTTER_DT_H
