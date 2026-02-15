#include "helpers.h"

void create_dir() {
  struct stat st;
  if (stat("data", &st) != 0) mkdir("data", 0755);
  if (stat("plot", &st) != 0) mkdir("plot", 0755);
}

void make_dt_tag(double dt, char tag[8]) {
  int code = (int)lround(dt * 100.0);
  std::snprintf(tag, 8, "%03d", code);
}

void build_grid(int imax, double dx, double x[]) {
  x[1] = 0.0;
  for (int i = 1; i <= imax - 1; i++) x[i + 1] = x[i] + dx;
}

void write_block_T(FILE* f, double t, int imax,
                   const double x[], const double T[]) {
  std::fprintf(f, "# t = %.2f hr\n# x\tT\n", t);
  for (int i = 1; i <= imax; i++) {
    std::fprintf(f, "%.6f\t%.6f\n", x[i], T[i]);
  }
  std::fprintf(f, "\n\n");
}

double err_L2(int imax, const double a[], const double b[]) {
  double s = 0.0;
  for (int i = 1; i <= imax; i++) {
    double e = a[i] - b[i];
    s += e * e;
  }
  return std::sqrt(s / (double)imax);
}

double err_Linf(int imax, const double a[], const double b[]) {
  double m = 0.0;
  for (int i = 1; i <= imax; i++) {
    double e = std::fabs(a[i] - b[i]);
    if (e > m) m = e;
  }
  return m;
}
