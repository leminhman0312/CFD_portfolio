#ifndef METHODS_H
#define METHODS_H

#include "headers.h"

void thomasTriDiagonal(int N,
                       const double a[], const double b[],
                       const double c[], const double d[],
                       double u[]);

void exact_solution_profile(int imax, const double x[], double t,
                            double alpha, double L,
                            double tb, double t0,
                            int nterms, double Tout[]);

void FTCS_explicit(int nmax, double F, double tb, double t0,
                   int imax, double Tout[]);

void implicit_interior_one_step(double F, double tb, double t0,
                                int imax, double u_full[]);

void DufortFrankel(int nmax, double F, double tb, double t0,
                   int imax, double Tout[]);

void FTCS_implicit(int nmax, double F, double tb, double t0,
                   int imax, double Tout[]);

void CrankNicolson(int nmax, double F, double tb, double t0,
                   int imax, double Tout[]);

#endif
