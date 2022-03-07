//
// Created by ohmy on 2022-03-07.
//

#ifndef DIRECTTRANSCRIPTION_PLOTTER_H
#define DIRECTTRANSCRIPTION_PLOTTER_H

#include <casadi/casadi.hpp>
using namespace casadi;

class Plotter {
public:
    void plot_path(DM x, DM y);
    void plot_path_heading(DM x, DM y, DM o);
    void plot_var(DM var);
};


#endif //DIRECTTRANSCRIPTION_PLOTTER_H
