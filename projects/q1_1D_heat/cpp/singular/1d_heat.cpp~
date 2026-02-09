// ============================================================================
//  Heat equation 1D (multiple schemes)  -  Organized layout
// ============================================================================

// ---- Standard / STL includes (keep what you actually use) -------------------
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#include <sys/stat.h>

// ============================================================================
//  Forward declarations
// ============================================================================

// ---------------------------------------------------------------------------
// Driver routines
// ---------------------------------------------------------------------------
void sim(double delta_x, double delta_t, int imax,
         double t0, double tboundary, double F_num);

void plot(double delta_t);

void compare_schemes(double delta_t, double t_target = -1.0);

void compare_dt(const char* scheme, double dt1, double dt2, double t_target);
void compare_dt_all_times(const char* scheme, double dt1, double dt2);

void compare_error_schemes(double delta_t);
void compare_error_schemes(double delta_t, double t_target);

void convergence_study(double delta_x, int imax, double alpha,
                       double t0, double tboundary, double t_target,
                       const double dt_list[], int ndt);

// ---------------------------------------------------------------------------
// Numerical schemes
// ---------------------------------------------------------------------------
void FTCS_explicit(int nmax, double F_num, double tboundary,
                   double t0, int imax, double out_un[]);

void DufortFrankel(int nmax, double F_num, double tboundary,
                   double t0, int imax, double out_un[]);

void FTCS_implicit(int nmax, double F_num, double tboundary,
                   double t0, int imax, double out_un[]);

void CrankNicolson(int nmax, double F_num, double tboundary,
                   double t0, int imax, double out_un[]);

// Linear solver
void thomasTriDiagonal(int imax, double a[], double b[],
                       double c[], double d[], double u[]);

// Helper used by Dufort-Frankel starter step
void FTCS_return_implicit_interior(int nmax, double F_num, double tboundary,
                                   double t0, int imax, double u_full[]);

// Exact solution (series)
void exact_solution_profile(int imax, const double x[], double t,
                            double alpha, double L,
                            double tboundary, double t0,
                            int nterms, double out_T[]);

// ---------------------------------------------------------------------------
// File / plotting helpers
// ---------------------------------------------------------------------------
void ensure_plot_dir();
void ensure_compare_dir();
void ensure_data_dir();

void make_dt_tag(double dt, char tag[8]);

bool file_exists(const char* path);
bool scheme_data_exists(const char* scheme, double dt);

int last_block_index_and_time(const char* path, double* last_time_out);
double latest_time_in_file(const char* path);

// Error norms
double err_L2(int imax, const double a[], const double b[]);
double err_Linf(int imax, const double a[], const double b[]);

// ============================================================================
//  Main
// ============================================================================
int main() {
  int imax = 21;             // number of grid points
  double delta_x = 0.05;     // dx
  double delta_t = 0.05;     // dt
  double alpha = 0.1;        // thermal diffusivity
  double t0 = 100.0;         // initial temperature
  double tboundary = 300.0;  // boundary temperature

  double Fourier_num = (alpha * delta_t) / (std::pow(delta_x, 2.0));

  sim(delta_x, delta_t, imax, t0, tboundary, Fourier_num);

  std::printf("Plotting all schemes at 1 dt\n\n");
  plot(delta_t);

  std::printf("Comparing different schemes\n\n");
  compare_error_schemes(delta_t, 0.4);

  double dt_list[] = {0.10, 0.05, 0.02, 0.01, 0.005};
  std::printf("Convergence study\n\n");
  convergence_study(delta_x, imax, alpha, t0, tboundary, 0.4, dt_list, 5);

  return 0;
}

// ============================================================================
//  Numerical schemes
// ============================================================================

// ---- FTCS explicit ----------------------------------------------------------
void FTCS_explicit(int nmax, double F_num, double tboundary, double t0,
                   int imax, double out_un[]) {
  double u[imax + 1];
  double un[imax + 1];

  for (int i = 2; i <= imax - 1; i++) {
    u[i]  = t0;
    un[i] = t0;
  }

  u[1] = tboundary;
  u[imax] = tboundary;
  un[1] = tboundary;
  un[imax] = tboundary;

  for (int n = 1; n <= nmax; n++) {
    for (int i = 1; i <= imax; i++) u[i] = un[i];

    for (int i = 2; i <= imax - 1; i++) {
      un[i] = u[i] + F_num * (u[i + 1] - 2.0 * u[i] + u[i - 1]);
    }

    un[1] = tboundary;
    un[imax] = tboundary;
  }

  for (int i = 1; i <= imax; i++) out_un[i] = un[i];
}

// ---- Dufort–Frankel ---------------------------------------------------------
void DufortFrankel(int nmax, double F_num, double tboundary, double t0,
                   int imax, double out_un[]) {
  double d = 2.0 * F_num;

  double un_m1[imax + 1] = {0.0};
  double un_p1[imax + 1] = {0.0};
  double un[imax + 1];

  un_m1[1] = tboundary;
  un_m1[imax] = tboundary;
  for (int i = 2; i <= imax - 1; i++) un_m1[i] = t0;

  // starter step using implicit interior helper
  double u1[imax + 1];
  FTCS_return_implicit_interior(1, F_num, tboundary, t0, imax, u1);
  for (int i = 1; i <= imax; i++) un[i] = u1[i];

  for (int n = 1; n <= nmax - 1; n++) {
    for (int i = 2; i <= imax - 1; i++) {
      un_p1[i] =
          ((1.0 - d) * un_m1[i] + d * (un[i + 1] + un[i - 1])) / (1.0 + d);
    }

    un_p1[1] = tboundary;
    un_p1[imax] = tboundary;

    for (int i = 1; i <= imax; i++) {
      un_m1[i] = un[i];
      un[i] = un_p1[i];
    }
  }

  for (int i = 1; i <= imax; i++) out_un[i] = un[i];
}

// ---- Implicit interior helper (for starter step) ----------------------------
void FTCS_return_implicit_interior(int nmax, double F_num, double tboundary,
                                   double t0, int imax, double u_full[]) {
  int N = imax - 2;

  double a[N + 1], b[N + 1], c[N + 1], d[N + 1], u_int[N + 1];

  for (int i = 0; i <= N; i++) {
    a[i] = b[i] = c[i] = d[i] = u_int[i] = 0.0;
  }

  u_full[1] = tboundary;
  u_full[imax] = tboundary;
  for (int i = 2; i <= imax - 1; i++) u_full[i] = t0;

  for (int j = 1; j <= N; j++) {
    d[j] = 1.0 + 2.0 * F_num;
    a[j] = -F_num;
    b[j] = -F_num;
  }
  b[1] = 0.0;
  a[N] = 0.0;

  for (int t = 1; t <= nmax; t++) {
    for (int j = 1; j <= N; j++) {
      int i = j + 1;
      c[j] = u_full[i];
    }

    c[1] += F_num * tboundary;
    c[N] += F_num * tboundary;

    thomasTriDiagonal(N, a, b, c, d, u_int);

    u_full[1] = tboundary;
    u_full[imax] = tboundary;
    for (int j = 1; j <= N; j++) {
      int i = j + 1;
      u_full[i] = u_int[j];
    }
  }
}

// ---- FTCS implicit ----------------------------------------------------------
void FTCS_implicit(int nmax, double F_num, double tboundary, double t0,
                   int imax, double out_un[]) {
  double a[imax + 1], b[imax + 1], c[imax + 1], d[imax + 1], un[imax + 1];

  un[1] = tboundary;
  un[imax] = tboundary;
  for (int i = 2; i <= imax - 1; i++) un[i] = t0;

  d[1] = 1.0;
  a[1] = 0.0;
  b[1] = 0.0;

  d[imax] = 1.0;
  a[imax] = 0.0;
  b[imax] = 0.0;

  for (int i = 2; i <= imax - 1; i++) {
    d[i] = 1.0 + 2.0 * F_num;
    a[i] = -F_num;
    b[i] = -F_num;
  }

  for (int t = 1; t <= nmax; t++) {
    c[1] = tboundary;
    c[imax] = tboundary;
    for (int i = 2; i <= imax - 1; i++) c[i] = un[i];

    thomasTriDiagonal(imax, a, b, c, d, un);
  }

  for (int i = 1; i <= imax; i++) out_un[i] = un[i];
}

// ---- Crank–Nicolson ---------------------------------------------------------
void CrankNicolson(int nmax, double F_num, double tboundary, double t0,
                   int imax, double out_un[]) {
  double d_NC = F_num / 2.0;

  double a[imax + 1], b[imax + 1], c[imax + 1], d[imax + 1];
  double u0[imax + 1];
  double u_half[imax + 1];

  u0[1] = tboundary;
  u0[imax] = tboundary;
  for (int i = 2; i <= imax - 1; i++) u0[i] = t0;

  u_half[1] = tboundary;
  u_half[imax] = tboundary;

  d[1] = 1.0;
  a[1] = 0.0;
  b[1] = 0.0;

  d[imax] = 1.0;
  a[imax] = 0.0;
  b[imax] = 0.0;

  for (int i = 2; i <= imax - 1; i++) {
    d[i] = 1.0 + 2.0 * d_NC;
    a[i] = -d_NC;
    b[i] = -d_NC;
  }

  for (int n = 1; n <= nmax; n++) {
    u_half[1] = tboundary;
    u_half[imax] = tboundary;

    for (int i = 2; i <= imax - 1; i++) {
      u_half[i] = u0[i] + d_NC * (u0[i + 1] - 2.0 * u0[i] + u0[i - 1]);
    }

    c[1] = tboundary;
    c[imax] = tboundary;
    for (int i = 2; i <= imax - 1; i++) c[i] = u_half[i];

    thomasTriDiagonal(imax, a, b, c, d, u0);

    u0[1] = tboundary;
    u0[imax] = tboundary;
  }

  for (int i = 1; i <= imax; i++) out_un[i] = u0[i];
}

// ---- Thomas tridiagonal solver ---------------------------------------------
void thomasTriDiagonal(int imax, double a[], double b[], double c[], double d[],
                       double u[]) {
  double dprime[imax + 1];
  double cprime[imax + 1];

  dprime[1] = d[1];
  cprime[1] = c[1];

  for (int i = 2; i <= imax; i++) {
    dprime[i] = d[i] - (b[i] * a[i - 1]) / dprime[i - 1];
    cprime[i] = c[i] - (cprime[i - 1] * b[i]) / dprime[i - 1];
  }

  u[imax] = cprime[imax] / dprime[imax];

  for (int i = imax - 1; i >= 1; i--) {
    u[i] = (cprime[i] - a[i] * u[i + 1]) / dprime[i];
  }
}

// ---- Exact solution ---------------------------------------------------------
void exact_solution_profile(int imax, const double x[], double t, double alpha,
                            double L, double tboundary, double t0, int nterms,
                            double out_T[]) {
  const double pi = 4.0 * std::atan(1.0);
  const double A = (t0 - tboundary);

  for (int i = 1; i <= imax; i++) {
    double sum = 0.0;

    for (int n = 1; n <= nterms; n++) {
      double coeff = (2.0 * A / (n * pi)) * (1.0 - std::pow(-1.0, n));
      if (coeff == 0.0) continue;

      double k = n * pi / L;
      sum += coeff * std::sin(k * x[i]) * std::exp(-alpha * k * k * t);
    }

    out_T[i] = tboundary + sum;
  }

  out_T[1] = tboundary;
  out_T[imax] = tboundary;
}

// ============================================================================
//  Helpers: directories, filenames, file parsing, error norms
// ============================================================================

void ensure_plot_dir() {
  struct stat st;
  if (stat("plot", &st) != 0) mkdir("plot", 0755);
}

void ensure_compare_dir() {
  struct stat st;
  if (stat("compare", &st) != 0) mkdir("compare", 0755);
}

void ensure_data_dir() {
  struct stat st;
  if (stat("data", &st) != 0) mkdir("data", 0755);
}

void make_dt_tag(double dt, char tag[8]) {
  int code = (int)lround(dt * 100.0);
  std::snprintf(tag, 8, "%03d", code);
}

bool file_exists(const char* path) {
  struct stat st;
  return (stat(path, &st) == 0);
}

bool scheme_data_exists(const char* scheme, double dt) {
  char tag[8];
  make_dt_tag(dt, tag);

  char fname[128];
  std::snprintf(fname, sizeof(fname), "data/%s_%s.txt", scheme, tag);
  return file_exists(fname);
}

int last_block_index_and_time(const char* path, double* last_time_out) {
  std::ifstream in(path);
  if (!in) {
    if (last_time_out) *last_time_out = 0.0;
    return 0;
  }

  std::string line;
  int idx = -1;
  double last_t = 0.0;

  while (std::getline(in, line)) {
    if (line.rfind("# t =", 0) == 0) {
      idx += 1;
      const char* s = line.c_str();
      const char* eq = std::strchr(s, '=');
      if (eq) last_t = std::atof(eq + 1);
    }
  }

  if (idx < 0) idx = 0;
  if (last_time_out) *last_time_out = last_t;
  return idx;
}

double latest_time_in_file(const char* path) {
  FILE* f = std::fopen(path, "r");
  if (!f) return 0.0;

  char line[256];
  double last_t = 0.0;

  while (std::fgets(line, sizeof(line), f)) {
    if (std::strncmp(line, "# t =", 5) == 0) {
      double t = 0.0;
      if (std::sscanf(line, "# t = %lf", &t) == 1) last_t = t;
    }
  }

  std::fclose(f);
  return last_t;
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

// ============================================================================
//  Plot / compare drivers (gnuplot wrappers)
// ============================================================================

// ---- compare dt at a single time ------------------------------------------
void compare_dt(const char* scheme, double dt1, double dt2, double t_target) {
  ensure_compare_dir();

  if (!scheme_data_exists(scheme, dt1) || !scheme_data_exists(scheme, dt2))
    return;

  char tag1[8], tag2[8];
  make_dt_tag(dt1, tag1);
  make_dt_tag(dt2, tag2);

  int idx = (int)lround(t_target / 0.1);

  char cmd[512];
  std::snprintf(
      cmd, sizeof(cmd),
      "gnuplot -e \"scheme='%s'; tag1='%s'; dt1='%.2f'; tag2='%s'; dt2='%.2f'; "
      "idx=%d; tlabel='%.1f'\" gnuplot_scripts/compare_dt.gp",
      scheme, tag1, dt1, tag2, dt2, idx, t_target);

  std::system(cmd);
}

void compare_dt_all_times(const char* scheme, double dt1, double dt2) {
  double times[4] = {0.1, 0.2, 0.3, 0.4};
  for (int k = 0; k < 4; k++) compare_dt(scheme, dt1, dt2, times[k]);
}

// ---- compare schemes at one time -------------------------------------------
void compare_schemes(double delta_t, double t_target) {
  ensure_compare_dir();

  char tag[8];
  make_dt_tag(delta_t, tag);

  char exact_path[128];
  std::snprintf(exact_path, sizeof(exact_path), "data/exact_%s.txt", tag);

  int idx = 0;
  double tlabel = 0.0;

  if (t_target < 0.0) {
    idx = last_block_index_and_time(exact_path, &tlabel);
  } else {
    tlabel = t_target;
    idx = (int)lround(t_target / 0.1);
  }

  char cmd[256];
  std::snprintf(cmd, sizeof(cmd),
                "gnuplot -e \"tag='%s'; dt='%.2f'; idx=%d; tlabel='%.1f'\" "
                "gnuplot_scripts/compare_schemes.gp",
                tag, delta_t, idx, tlabel);

  std::system(cmd);
}

// ---- error scheme compare wrappers -----------------------------------------
void compare_error_schemes(double delta_t) {
  char tag[8];
  make_dt_tag(delta_t, tag);

  char exact_path[128];
  std::snprintf(exact_path, sizeof(exact_path), "data/exact_%s.txt", tag);

  double t_latest = latest_time_in_file(exact_path);
  if (t_latest <= 0.0) return;

  compare_error_schemes(delta_t, t_latest);
}

void compare_error_schemes(double delta_t, double t_target) {
  ensure_compare_dir();

  char tag[8];
  make_dt_tag(delta_t, tag);

  int idx = (int)lround(t_target / 0.1);

  char cmd[512];
  std::snprintf(cmd, sizeof(cmd),
                "gnuplot -e \"tag='%s'; dt=%.2f; idx=%d; tlabel='%.1f'; "
                "outpng='plot/error_schemes_%s_t%.1f.png'\" "
                "gnuplot_scripts/compare_error_schemes.gp",
                tag, delta_t, idx, t_target, tag, t_target);

  std::system(cmd);
}

// ---- plot all schemes for one dt -------------------------------------------
void plot(double delta_t) {
  ensure_plot_dir();

  char tag[8];
  make_dt_tag(delta_t, tag);

  char cmd[256];

  std::snprintf(
      cmd, sizeof(cmd),
      "gnuplot -e \"scheme='ftcs_explicit'; tag='%s'; dt='%.2f'; "
      "outpng='plot/ftcs_explicit_%s.png'\" gnuplot_scripts/plot_scheme.gp",
      tag, delta_t, tag);
  std::system(cmd);

  std::snprintf(
      cmd, sizeof(cmd),
      "gnuplot -e \"scheme='dufort'; tag='%s'; dt='%.2f'; "
      "outpng='plot/dufort_%s.png'\" gnuplot_scripts/plot_scheme.gp",
      tag, delta_t, tag);
  std::system(cmd);

  std::snprintf(
      cmd, sizeof(cmd),
      "gnuplot -e \"scheme='ftcs_implicit'; tag='%s'; dt='%.2f'; "
      "outpng='plot/ftcs_implicit_%s.png'\" gnuplot_scripts/plot_scheme.gp",
      tag, delta_t, tag);
  std::system(cmd);

  std::snprintf(
      cmd, sizeof(cmd),
      "gnuplot -e \"scheme='cn'; tag='%s'; dt='%.2f'; "
      "outpng='plot/cn_%s.png'\" gnuplot_scripts/plot_scheme.gp",
      tag, delta_t, tag);
  std::system(cmd);
}

// ============================================================================
//  Simulation driver and convergence
// ============================================================================

// ---- run one dt and write all scheme data files -----------------------------
void sim(double delta_x, double delta_t, int imax, double t0, double tboundary,
         double F_num) {
  ensure_data_dir();

  char tag[8];
  make_dt_tag(delta_t, tag);

  char f_explicit[64], f_dufort[64], f_implicit[64], f_cn[64], f_exact[64];
  std::snprintf(f_explicit, 64, "data/ftcs_explicit_%s.txt", tag);
  std::snprintf(f_dufort,   64, "data/dufort_%s.txt", tag);
  std::snprintf(f_implicit, 64, "data/ftcs_implicit_%s.txt", tag);
  std::snprintf(f_cn,       64, "data/cn_%s.txt", tag);
  std::snprintf(f_exact,    64, "data/exact_%s.txt", tag);

  FILE* outfile_explicit = std::fopen(f_explicit, "w");
  FILE* outfile_dufort   = std::fopen(f_dufort, "w");
  FILE* outfile_implicit = std::fopen(f_implicit, "w");
  FILE* outfile_cn       = std::fopen(f_cn, "w");
  FILE* outfile_exact    = std::fopen(f_exact, "w");

  char f_err_exp[64], f_err_imp[64], f_err_df[64], f_err_cn[64];
  std::snprintf(f_err_exp, 64, "data/error_ftcs_explicit_%s.txt", tag);
  std::snprintf(f_err_imp, 64, "data/error_ftcs_implicit_%s.txt", tag);
  std::snprintf(f_err_df,  64, "data/error_dufort_%s.txt", tag);
  std::snprintf(f_err_cn,  64, "data/error_cn_%s.txt", tag);

  FILE* outfile_err_exp = std::fopen(f_err_exp, "w");
  FILE* outfile_err_imp = std::fopen(f_err_imp, "w");
  FILE* outfile_err_df  = std::fopen(f_err_df, "w");
  FILE* outfile_err_cn  = std::fopen(f_err_cn, "w");

  double x_vector[imax + 1];
  x_vector[1] = 0.0;
  for (int i = 1; i <= imax - 1; i++) x_vector[i + 1] = x_vector[i] + delta_x;

  double L = (imax - 1) * delta_x;
  double alpha = F_num * (delta_x * delta_x) / delta_t;

  double T_exact[imax + 1];
  double T_explicit[imax + 1], T_dufort[imax + 1],
         T_implicit[imax + 1], T_cn[imax + 1];

  for (int k = 0; k < 5; k++) {
    double t_target = 0.1 * k;
    int nmax = (int)lround(t_target / delta_t);

    FTCS_explicit(nmax, F_num, tboundary, t0, imax, T_explicit);
    DufortFrankel(nmax, F_num, tboundary, t0, imax, T_dufort);
    FTCS_implicit(nmax, F_num, tboundary, t0, imax, T_implicit);
    CrankNicolson(nmax, F_num, tboundary, t0, imax, T_cn);

    exact_solution_profile(imax, x_vector, t_target, alpha, L,
                           tboundary, t0, 200, T_exact);

    std::fprintf(outfile_explicit, "# t = %.2f hr\n# x\tT\n", t_target);
    std::fprintf(outfile_dufort,   "# t = %.2f hr\n# x\tT\n", t_target);
    std::fprintf(outfile_implicit, "# t = %.2f hr\n# x\tT\n", t_target);
    std::fprintf(outfile_cn,       "# t = %.2f hr\n# x\tT\n", t_target);
    std::fprintf(outfile_exact,    "# t = %.2f hr\n# x\tT\n", t_target);

    std::fprintf(outfile_err_exp, "# t = %.2f hr\n# x\terr\n", t_target);
    std::fprintf(outfile_err_imp, "# t = %.2f hr\n# x\terr\n", t_target);
    std::fprintf(outfile_err_df,  "# t = %.2f hr\n# x\terr\n", t_target);
    std::fprintf(outfile_err_cn,  "# t = %.2f hr\n# x\terr\n", t_target);

    for (int i = 1; i <= imax; i++) {
      std::fprintf(outfile_explicit, "%.6f\t%.6f\n", x_vector[i], T_explicit[i]);
      std::fprintf(outfile_dufort,   "%.6f\t%.6f\n", x_vector[i], T_dufort[i]);
      std::fprintf(outfile_implicit, "%.6f\t%.6f\n", x_vector[i], T_implicit[i]);
      std::fprintf(outfile_cn,       "%.6f\t%.6f\n", x_vector[i], T_cn[i]);
      std::fprintf(outfile_exact,    "%.6f\t%.6f\n", x_vector[i], T_exact[i]);

      std::fprintf(outfile_err_exp, "%.6f\t%.6f\n", x_vector[i],
                   T_explicit[i] - T_exact[i]);
      std::fprintf(outfile_err_imp, "%.6f\t%.6f\n", x_vector[i],
                   T_implicit[i] - T_exact[i]);
      std::fprintf(outfile_err_df,  "%.6f\t%.6f\n", x_vector[i],
                   T_dufort[i] - T_exact[i]);
      std::fprintf(outfile_err_cn,  "%.6f\t%.6f\n", x_vector[i],
                   T_cn[i] - T_exact[i]);
    }

    std::fprintf(outfile_explicit, "\n\n");
    std::fprintf(outfile_dufort,   "\n\n");
    std::fprintf(outfile_implicit, "\n\n");
    std::fprintf(outfile_cn,       "\n\n");
    std::fprintf(outfile_exact,    "\n\n");

    std::fprintf(outfile_err_exp, "\n\n");
    std::fprintf(outfile_err_imp, "\n\n");
    std::fprintf(outfile_err_df,  "\n\n");
    std::fprintf(outfile_err_cn,  "\n\n");
  }

  std::fclose(outfile_explicit);
  std::fclose(outfile_dufort);
  std::fclose(outfile_implicit);
  std::fclose(outfile_cn);
  std::fclose(outfile_exact);

  std::fclose(outfile_err_exp);
  std::fclose(outfile_err_imp);
  std::fclose(outfile_err_df);
  std::fclose(outfile_err_cn);
}

// ---- Convergence study ------------------------------------------------------
void convergence_study(double delta_x, int imax, double alpha, double t0,
                       double tboundary, double t_target,
                       const double dt_list[], int ndt) {
  ensure_data_dir();
  ensure_compare_dir();

  char outdata[128];
  std::snprintf(outdata, sizeof(outdata), "data/convergence_t%03d.txt",
                (int)lround(t_target * 100.0));

  FILE* f = std::fopen(outdata, "w");
  if (!f) return;

  std::fprintf(f, "# convergence at t = %.2f hr\n", t_target);
  std::fprintf(f,
               "# dt\tL2_exp\tL2_imp\tL2_df\tL2_cn\t"
               "Linf_exp\tLinf_imp\tLinf_df\tLinf_cn\n");

  double x_vector[imax + 1];
  x_vector[1] = 0.0;
  for (int i = 1; i <= imax - 1; i++) x_vector[i + 1] = x_vector[i] + delta_x;

  double L = (imax - 1) * delta_x;

  for (int k = 0; k < ndt; k++) {
    double dt = dt_list[k];
    double F_num = (alpha * dt) / (delta_x * delta_x);
    int nmax = (int)lround(t_target / dt);

    double T_exp[imax + 1], T_imp[imax + 1], T_df[imax + 1],
           T_cn[imax + 1],  T_ex[imax + 1];

    FTCS_explicit(nmax, F_num, tboundary, t0, imax, T_exp);
    FTCS_implicit(nmax, F_num, tboundary, t0, imax, T_imp);
    DufortFrankel(nmax, F_num, tboundary, t0, imax, T_df);
    CrankNicolson(nmax, F_num, tboundary, t0, imax, T_cn);

    exact_solution_profile(imax, x_vector, t_target, alpha, L,
                           tboundary, t0, 200, T_ex);

    double L2_exp = err_L2(imax, T_exp, T_ex);
    double L2_imp = err_L2(imax, T_imp, T_ex);
    double L2_df  = err_L2(imax, T_df,  T_ex);
    double L2_cn  = err_L2(imax, T_cn,  T_ex);

    double Li_exp = err_Linf(imax, T_exp, T_ex);
    double Li_imp = err_Linf(imax, T_imp, T_ex);
    double Li_df  = err_Linf(imax, T_df,  T_ex);
    double Li_cn  = err_Linf(imax, T_cn,  T_ex);

    std::fprintf(
        f, "%.6f\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
        dt, L2_exp, L2_imp, L2_df, L2_cn, Li_exp, Li_imp, Li_df, Li_cn);
  }

  std::fclose(f);

  char cmd[512];
  std::snprintf(
      cmd, sizeof(cmd),
      "gnuplot -e \"infile='%s'; tlabel='%.2f'; "
      "outpng='compare/convergence_t%.2f.png'\" gnuplot_scripts/convergence.gp",
      outdata, t_target, t_target);
  std::system(cmd);
}
