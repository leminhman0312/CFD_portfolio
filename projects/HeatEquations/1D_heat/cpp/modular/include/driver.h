#ifndef DRIVER_H
#define DRIVER_H

#include "headers.h"

void sim(double dx, double dt, int imax,
         double t0, double tb, double F);

void plot(double dt);

void convergence_study(double dx, int imax, double alpha,
                       double t0, double tb, double t_target);

void plot_error_schemes(double dt, double t_target);
void plot_convergence(double t_target);


#endif
