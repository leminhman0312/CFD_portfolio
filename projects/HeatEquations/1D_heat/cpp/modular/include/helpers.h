#ifndef HELPERS_H
#define HELPERS_H

#include "headers.h"

void create_dir();
void make_dt_tag(double dt, char tag[8]);
void build_grid(int imax, double dx, double x[]);
void write_block_T(FILE* f, double t, int imax,
                   const double x[], const double T[]);

double err_L2(int imax, const double a[], const double b[]);
double err_Linf(int imax, const double a[], const double b[]);


#endif
