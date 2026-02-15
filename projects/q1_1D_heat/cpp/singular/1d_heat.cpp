// ============================================================================
//  heat_1d_clean.cpp
//
//  Learning focused 1D heat equation code
//
//  Structure
//    main
//    helpers
//    tridiagonal solver
//    exact solution
//    numerical schemes
//    sim (writes data files)
//    convergence_study (writes one table file)
//    plot (calls gnuplot)
//
//  Build
//    g++ -O2 -std=c++17 heat_1d_clean.cpp -o heat1d
//
//  Run
//    ./heat1d
// ============================================================================

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <sys/stat.h>

// ============================================================================
//  Helpers
// ============================================================================

static void create_dir() {
  struct stat st;

  if (stat("data", &st) != 0) {
    mkdir("data", 0755);
  }

  if (stat("plot", &st) != 0) {
    mkdir("plot", 0755);
  }
}

static void make_dt_tag(double dt, char tag[8]) {
  int code = (int)lround(dt * 100.0);
  std::snprintf(tag, 8, "%03d", code);
}

static void make_t_tag(double t, char tag[8]) {
  int code = (int)lround(t * 100.0);
  std::snprintf(tag, 8, "%03d", code);
}

static void build_grid(int imax, double dx, double x[]) {
  x[1] = 0.0;
  for (int i = 1; i <= imax - 1; i++) {
    x[i + 1] = x[i] + dx;
  }
}

static void write_block_T(FILE* f, double t, int imax,
                          const double x[], const double T[]) {
  std::fprintf(f, "# t = %.2f hr\n# x\tT\n", t);
  for (int i = 1; i <= imax; i++) {
    std::fprintf(f, "%.6f\t%.6f\n", x[i], T[i]);
  }
  std::fprintf(f, "\n\n");
}

static void write_block_error(FILE* f, double t, int imax,
                              const double x[],
                              const double T[],
                              const double Texact[]) {
  std::fprintf(f, "# t = %.2f hr\n# x\terr\n", t);
  for (int i = 1; i <= imax; i++) {
    std::fprintf(f, "%.6f\t%.6f\n", x[i], (T[i] - Texact[i]));
  }
  std::fprintf(f, "\n\n");
}

static double error_L2(int imax, const double a[], const double b[]) {
  double s = 0.0;
  for (int i = 1; i <= imax; i++) {
    double e = a[i] - b[i];
    s += e * e;
  }
  return std::sqrt(s / (double)imax);
}

static double error_Linf(int imax, const double a[], const double b[]) {
  double m = 0.0;
  for (int i = 1; i <= imax; i++) {
    double e = std::fabs(a[i] - b[i]);
    if (e > m) m = e;
  }
  return m;
}

// ============================================================================
//  Tridiagonal solver
//  Solves for u[1..N] given a,b,c,d arrays indexed 1..N
// ============================================================================

static void thomasTriDiagonal(int N,
                             const double a[], const double b[],
                             const double c[], const double d[],
                             double u[]) {
  std::vector<double> dprime(N + 1, 0.0);
  std::vector<double> cprime(N + 1, 0.0);

  dprime[1] = d[1];
  cprime[1] = c[1];

  for (int i = 2; i <= N; i++) {
    dprime[i] = d[i] - (b[i] * a[i - 1]) / dprime[i - 1];
    cprime[i] = c[i] - (cprime[i - 1] * b[i]) / dprime[i - 1];
  }

  u[N] = cprime[N] / dprime[N];
  for (int i = N - 1; i >= 1; i--) {
    u[i] = (cprime[i] - a[i] * u[i + 1]) / dprime[i];
  }
}

// ============================================================================
//  Exact solution
// ============================================================================

static void exact_solution_profile(int imax, const double x[], double t,
                                   double alpha, double L,
                                   double tb, double t0,
                                   int nterms, double Tout[]) {
  const double pi = 4.0 * std::atan(1.0);
  const double A  = t0 - tb;

  for (int i = 1; i <= imax; i++) {
    double sum = 0.0;

    for (int n = 1; n <= nterms; n++) {
      double coeff = (2.0 * A / (n * pi)) * (1.0 - std::pow(-1.0, n));
      if (coeff == 0.0) continue;

      double k = n * pi / L;
      sum += coeff * std::sin(k * x[i]) * std::exp(-alpha * k * k * t);
    }

    Tout[i] = tb + sum;
  }

  Tout[1] = tb;
  Tout[imax] = tb;
}

// ============================================================================
//  Numerical schemes
// ============================================================================

static void FTCS_explicit(int nmax, double F, double tb, double t0,
                          int imax, double Tout[]) {
  std::vector<double> u(imax + 1, 0.0);
  std::vector<double> un(imax + 1, 0.0);

  un[1] = tb;
  un[imax] = tb;
  for (int i = 2; i <= imax - 1; i++) un[i] = t0;

  for (int n = 1; n <= nmax; n++) {
    u = un;
    for (int i = 2; i <= imax - 1; i++) {
      un[i] = u[i] + F * (u[i + 1] - 2.0 * u[i] + u[i - 1]);
    }
    un[1] = tb;
    un[imax] = tb;
  }

  for (int i = 1; i <= imax; i++) Tout[i] = un[i];
}

static void implicit_interior_one_step(double F, double tb, double t0,
                                       int imax, double u_full[]) {
  int N = imax - 2;

  std::vector<double> a(N + 1), b(N + 1), c(N + 1), d(N + 1), u_int(N + 1);

  u_full[1] = tb;
  u_full[imax] = tb;
  for (int j = 1; j <= N; j++) u_full[j + 1] = t0;

  for (int j = 1; j <= N; j++) {
    d[j] = 1.0 + 2.0 * F;
    a[j] = -F;
    b[j] = -F;
  }
  b[1] = 0.0;
  a[N] = 0.0;

  for (int j = 1; j <= N; j++) c[j] = u_full[j + 1];
  c[1] += F * tb;
  c[N] += F * tb;

  thomasTriDiagonal(N, a.data(), b.data(), c.data(), d.data(), u_int.data());

  u_full[1] = tb;
  u_full[imax] = tb;
  for (int j = 1; j <= N; j++) u_full[j + 1] = u_int[j];
}

static void DufortFrankel(int nmax, double F, double tb, double t0,
                          int imax, double Tout[]) {
  const double dval = 2.0 * F;

  std::vector<double> un_m1(imax + 1), un(imax + 1),
                      un_p1(imax + 1), u1(imax + 1);

  un_m1[1] = tb;
  un_m1[imax] = tb;
  for (int i = 2; i <= imax - 1; i++) un_m1[i] = t0;

  implicit_interior_one_step(F, tb, t0, imax, u1.data());
  un = u1;

  for (int step = 1; step <= nmax - 1; step++) {
    for (int i = 2; i <= imax - 1; i++) {
      un_p1[i] =
        ((1.0 - dval) * un_m1[i] + dval * (un[i + 1] + un[i - 1]))
        / (1.0 + dval);
    }
    un_p1[1] = tb;
    un_p1[imax] = tb;

    un_m1 = un;
    un    = un_p1;
  }

  for (int i = 1; i <= imax; i++) Tout[i] = un[i];
}

static void FTCS_implicit(int nmax, double F, double tb, double t0,
                          int imax, double Tout[]) {
  std::vector<double> a(imax + 1, 0.0), b(imax + 1, 0.0),
                      c(imax + 1, 0.0), d(imax + 1, 0.0),
                      un(imax + 1, 0.0);

  un[1] = tb;
  un[imax] = tb;
  for (int i = 2; i <= imax - 1; i++) un[i] = t0;

  d[1] = 1.0;
  d[imax] = 1.0;

  for (int i = 2; i <= imax - 1; i++) {
    d[i] = 1.0 + 2.0 * F;
    a[i] = -F;
    b[i] = -F;
  }

  for (int step = 1; step <= nmax; step++) {
    c[1] = tb;
    c[imax] = tb;
    for (int i = 2; i <= imax - 1; i++) c[i] = un[i];
    thomasTriDiagonal(imax, a.data(), b.data(), c.data(), d.data(), un.data());

    un[1] = tb;
    un[imax] = tb;
  }

  for (int i = 1; i <= imax; i++) Tout[i] = un[i];
}

static void CrankNicolson(int nmax, double F, double tb, double t0,
                          int imax, double Tout[]) {
  const double dnc = F / 2.0;

  std::vector<double> a(imax + 1, 0.0), b(imax + 1, 0.0),
                      c(imax + 1, 0.0), d(imax + 1, 0.0),
                      u0(imax + 1, 0.0), uhalf(imax + 1, 0.0);

  u0[1] = tb;
  u0[imax] = tb;
  for (int i = 2; i <= imax - 1; i++) u0[i] = t0;

  d[1] = 1.0;
  d[imax] = 1.0;

  for (int i = 2; i <= imax - 1; i++) {
    d[i] = 1.0 + 2.0 * dnc;
    a[i] = -dnc;
    b[i] = -dnc;
  }

  for (int step = 1; step <= nmax; step++) {
    uhalf[1] = tb;
    uhalf[imax] = tb;
    for (int i = 2; i <= imax - 1; i++) {
      uhalf[i] = u0[i] + dnc *
        (u0[i + 1] - 2.0 * u0[i] + u0[i - 1]);
    }

    c[1] = tb;
    c[imax] = tb;
    for (int i = 2; i <= imax - 1; i++) c[i] = uhalf[i];

    thomasTriDiagonal(imax, a.data(), b.data(), c.data(), d.data(), u0.data());

    u0[1] = tb;
    u0[imax] = tb;
  }

  for (int i = 1; i <= imax; i++) Tout[i] = u0[i];
}

// ============================================================================
//  sim
//  Writes solution files and error files for one dt
// ============================================================================

static void sim(double dx, double dt, int imax,
                double t0, double tb, double F) {
  char tag[8];
  make_dt_tag(dt, tag);

  char fexp[128], fdf[128], fimp[128], fcn[128], fex[128];
  std::snprintf(fexp, sizeof(fexp), "data/ftcs_explicit_%s.txt", tag);
  std::snprintf(fdf,  sizeof(fdf),  "data/dufort_%s.txt", tag);
  std::snprintf(fimp, sizeof(fimp), "data/ftcs_implicit_%s.txt", tag);
  std::snprintf(fcn,  sizeof(fcn),  "data/cn_%s.txt", tag);
  std::snprintf(fex,  sizeof(fex),  "data/exact_%s.txt", tag);

  char feexp[128], fedf[128], feimp[128], fecn[128];
  std::snprintf(feexp, sizeof(feexp), "data/error_ftcs_explicit_%s.txt", tag);
  std::snprintf(fedf,  sizeof(fedf),  "data/error_dufort_%s.txt", tag);
  std::snprintf(feimp, sizeof(feimp), "data/error_ftcs_implicit_%s.txt", tag);
  std::snprintf(fecn,  sizeof(fecn),  "data/error_cn_%s.txt", tag);

  FILE* uexp  = std::fopen(fexp,  "w");
  FILE* udf   = std::fopen(fdf,   "w");
  FILE* uimp  = std::fopen(fimp,  "w");
  FILE* ucn   = std::fopen(fcn,   "w");
  FILE* uex   = std::fopen(fex,   "w");

  FILE* ueexp = std::fopen(feexp, "w");
  FILE* uedf  = std::fopen(fedf,  "w");
  FILE* ueimp = std::fopen(feimp, "w");
  FILE* uecn  = std::fopen(fecn,  "w");

  if (!uexp || !udf || !uimp || !ucn || !uex || !ueexp || !uedf || !ueimp || !uecn) {
    std::perror("fopen");
    std::exit(1);
  }

  std::vector<double> x(imax + 1);
  build_grid(imax, dx, x.data());

  double L = (imax - 1) * dx;
  double alpha_loc = F * (dx * dx) / dt;

  std::vector<double> Tex(imax + 1),
                      Texp(imax + 1), Tdf(imax + 1),
                      Timp(imax + 1), Tcn(imax + 1);

  for (int k = 0; k <= 4; k++) {
    double t = 0.1 * k;
    int nmax = (int)lround(t / dt);

    FTCS_explicit(nmax, F, tb, t0, imax, Texp.data());
    DufortFrankel(nmax, F, tb, t0, imax, Tdf.data());
    FTCS_implicit(nmax, F, tb, t0, imax, Timp.data());
    CrankNicolson(nmax, F, tb, t0, imax, Tcn.data());

    exact_solution_profile(imax, x.data(), t, alpha_loc, L, tb, t0, 200, Tex.data());

    write_block_T(uexp, t, imax, x.data(), Texp.data());
    write_block_T(udf,  t, imax, x.data(), Tdf.data());
    write_block_T(uimp, t, imax, x.data(), Timp.data());
    write_block_T(ucn,  t, imax, x.data(), Tcn.data());
    write_block_T(uex,  t, imax, x.data(), Tex.data());

    write_block_error(ueexp, t, imax, x.data(), Texp.data(), Tex.data());
    write_block_error(uedf,  t, imax, x.data(), Tdf.data(),  Tex.data());
    write_block_error(ueimp, t, imax, x.data(), Timp.data(), Tex.data());
    write_block_error(uecn,  t, imax, x.data(), Tcn.data(),  Tex.data());
  }

  std::fclose(uexp);
  std::fclose(udf);
  std::fclose(uimp);
  std::fclose(ucn);
  std::fclose(uex);

  std::fclose(ueexp);
  std::fclose(uedf);
  std::fclose(ueimp);
  std::fclose(uecn);
}

// ============================================================================
//  convergence study
//  Writes data/convergence_tXXX.txt where XXX = t_target*100
// ============================================================================

static void convergence_study(double dx, int imax, double alpha,
                              double t0, double tb, double t_target) {
  const int ndt = 5;
  const double dt_list[ndt] = {0.10, 0.05, 0.02, 0.01, 0.005};

  char ttag[8];
  make_t_tag(t_target, ttag);

  char outdata[128];
  std::snprintf(outdata, sizeof(outdata), "data/convergence_t%s.txt", ttag);

  FILE* f = std::fopen(outdata, "w");
  if (!f) {
    std::perror("fopen convergence");
    std::exit(1);
  }

  std::fprintf(f, "# convergence at t = %.2f hr\n", t_target);
  std::fprintf(f,
    "# dt\tL2_exp\tL2_imp\tL2_df\tL2_cn\tLinf_exp\tLinf_imp\tLinf_df\tLinf_cn\n");

  std::vector<double> x(imax + 1);
  build_grid(imax, dx, x.data());
  double L = (imax - 1) * dx;

  std::vector<double> Tex(imax + 1), Texp(imax + 1), Timp(imax + 1),
                      Tdf(imax + 1), Tcn(imax + 1);

  for (int k = 0; k < ndt; k++) {
    double dt = dt_list[k];
    double F  = alpha * dt / (dx * dx);
    int nmax  = (int)lround(t_target / dt);

    FTCS_explicit(nmax, F, tb, t0, imax, Texp.data());
    FTCS_implicit(nmax, F, tb, t0, imax, Timp.data());
    DufortFrankel(nmax, F, tb, t0, imax, Tdf.data());
    CrankNicolson(nmax, F, tb, t0, imax, Tcn.data());

    exact_solution_profile(imax, x.data(), t_target, alpha, L, tb, t0, 200, Tex.data());

    double L2e = error_L2(imax, Texp.data(), Tex.data());
    double L2i = error_L2(imax, Timp.data(), Tex.data());
    double L2d = error_L2(imax, Tdf.data(),  Tex.data());
    double L2c = error_L2(imax, Tcn.data(),  Tex.data());

    double Lie = error_Linf(imax, Texp.data(), Tex.data());
    double Lii = error_Linf(imax, Timp.data(), Tex.data());
    double Lid = error_Linf(imax, Tdf.data(),  Tex.data());
    double Lic = error_Linf(imax, Tcn.data(),  Tex.data());

    std::fprintf(f,
      "%.6f\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
      dt, L2e, L2i, L2d, L2c, Lie, Lii, Lid, Lic);
  }

  std::fclose(f);
}

// ============================================================================
//  plot
//  Calls gnuplot scripts to generate png files into plot/
// ============================================================================

static void plot(double dt, double t_target) {
  char tag[8];
  make_dt_tag(dt, tag);

  char ttag[8];
  make_t_tag(t_target, ttag);

  char cmd[512];

  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"scheme='ftcs_explicit'; tag='%s'; dt='%.6f'; "
    "outpng='plot/ftcs_explicit_%s.png'\" gnuplot_scripts/plot_scheme.gp",
    tag, dt, tag);
  std::system(cmd);

  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"scheme='dufort'; tag='%s'; dt='%.6f'; "
    "outpng='plot/dufort_%s.png'\" gnuplot_scripts/plot_scheme.gp",
    tag, dt, tag);
  std::system(cmd);

  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"scheme='ftcs_implicit'; tag='%s'; dt='%.6f'; "
    "outpng='plot/ftcs_implicit_%s.png'\" gnuplot_scripts/plot_scheme.gp",
    tag, dt, tag);
  std::system(cmd);

  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"scheme='cn'; tag='%s'; dt='%.6f'; "
    "outpng='plot/cn_%s.png'\" gnuplot_scripts/plot_scheme.gp",
    tag, dt, tag);
  std::system(cmd);

  int idx = (int)lround(t_target / 0.1);

  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"tag='%s'; dt=%.6f; idx=%d; tlabel='%.1f'; "
    "outpng='plot/error_schemes_%s_t%.1f.png'\" "
    "gnuplot_scripts/compare_error_schemes.gp",
    tag, dt, idx, t_target, tag, t_target);
  std::system(cmd);

  char infile[256];
  std::snprintf(infile, sizeof(infile), "data/convergence_t%s.txt", ttag);

  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"infile='%s'; tlabel='%.2f'; "
    "outpng='plot/convergence_t%.2f.png'\" "
    "gnuplot_scripts/convergence.gp",
    infile, t_target, t_target);
  std::system(cmd);
}

// ============================================================================
//  main
// ============================================================================

int main() {
  int imax = 21;
  double dx = 0.05;
  double dt = 0.01;
  double alpha = 0.1;
  double t0 = 100.0;
  double tb = 300.0;

  double F = alpha * dt / (dx * dx);

  double t_target = 0.4;

  create_dir();
  sim(dx, dt, imax, t0, tb, F);
  convergence_study(dx, imax, alpha, t0, tb, t_target);
  plot(dt, t_target);

  std::printf("Done\n");
  return 0;
}
