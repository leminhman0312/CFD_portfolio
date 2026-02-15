#include "driver.h"
#include "helpers.h"
#include "methods.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

// -----------------------------------------------------------------------------
// sim: write scheme profiles and per scheme error files
// -----------------------------------------------------------------------------
void sim(double dx, double dt, int imax, double t0, double tb, double F) {
  char tag[8];
  make_dt_tag(dt, tag);

  char fexp[128], fdf[128], fimp[128], fcn[128], fex[128];
  std::snprintf(fexp, sizeof(fexp), "data/ftcs_explicit_%s.txt", tag);
  std::snprintf(fdf,  sizeof(fdf),  "data/dufort_%s.txt", tag);
  std::snprintf(fimp, sizeof(fimp), "data/ftcs_implicit_%s.txt", tag);
  std::snprintf(fcn,  sizeof(fcn),  "data/cn_%s.txt", tag);
  std::snprintf(fex,  sizeof(fex),  "data/exact_%s.txt", tag);

  FILE* uexp = std::fopen(fexp, "w");
  FILE* udf  = std::fopen(fdf,  "w");
  FILE* uimp = std::fopen(fimp, "w");
  FILE* ucn  = std::fopen(fcn,  "w");
  FILE* uex  = std::fopen(fex,  "w");

  if (!uexp || !udf || !uimp || !ucn || !uex) {
    std::perror("fopen");
    std::exit(1);
  }

  char ferr_exp[128], ferr_imp[128], ferr_df[128], ferr_cn[128];
  std::snprintf(ferr_exp, sizeof(ferr_exp), "data/error_ftcs_explicit_%s.txt", tag);
  std::snprintf(ferr_imp, sizeof(ferr_imp), "data/error_ftcs_implicit_%s.txt", tag);
  std::snprintf(ferr_df,  sizeof(ferr_df),  "data/error_dufort_%s.txt", tag);
  std::snprintf(ferr_cn,  sizeof(ferr_cn),  "data/error_cn_%s.txt", tag);

  FILE* uerr_exp = std::fopen(ferr_exp, "w");
  FILE* uerr_imp = std::fopen(ferr_imp, "w");
  FILE* uerr_df  = std::fopen(ferr_df,  "w");
  FILE* uerr_cn  = std::fopen(ferr_cn,  "w");

  if (!uerr_exp || !uerr_imp || !uerr_df || !uerr_cn) {
    std::perror("fopen");
    std::exit(1);
  }

  std::vector<double> x(imax + 1, 0.0);
  build_grid(imax, dx, x.data());

  double L = (imax - 1) * dx;
  double alpha = F * (dx * dx) / dt;

  std::vector<double> Tex(imax + 1, 0.0);
  std::vector<double> Texp(imax + 1, 0.0);
  std::vector<double> Tdf(imax + 1, 0.0);
  std::vector<double> Timp(imax + 1, 0.0);
  std::vector<double> Tcn(imax + 1, 0.0);

  for (int k = 0; k <= 4; k++) {
    double t = 0.1 * (double)k;
    int nmax = (int)lround(t / dt);

    FTCS_explicit(nmax, F, tb, t0, imax, Texp.data());
    DufortFrankel(nmax, F, tb, t0, imax, Tdf.data());
    FTCS_implicit(nmax, F, tb, t0, imax, Timp.data());
    CrankNicolson(nmax, F, tb, t0, imax, Tcn.data());

    exact_solution_profile(imax, x.data(), t, alpha, L, tb, t0, 200, Tex.data());

    write_block_T(uexp, t, imax, x.data(), Texp.data());
    write_block_T(udf,  t, imax, x.data(), Tdf.data());
    write_block_T(uimp, t, imax, x.data(), Timp.data());
    write_block_T(ucn,  t, imax, x.data(), Tcn.data());
    write_block_T(uex,  t, imax, x.data(), Tex.data());

    std::fprintf(uerr_exp, "# t = %.2f hr\n# x\terr\n", t);
    std::fprintf(uerr_imp, "# t = %.2f hr\n# x\terr\n", t);
    std::fprintf(uerr_df,  "# t = %.2f hr\n# x\terr\n", t);
    std::fprintf(uerr_cn,  "# t = %.2f hr\n# x\terr\n", t);

    for (int i = 1; i <= imax; i++) {
      std::fprintf(uerr_exp, "%.6f\t%.6f\n", x[i], Texp[i] - Tex[i]);
      std::fprintf(uerr_imp, "%.6f\t%.6f\n", x[i], Timp[i] - Tex[i]);
      std::fprintf(uerr_df,  "%.6f\t%.6f\n", x[i], Tdf[i]  - Tex[i]);
      std::fprintf(uerr_cn,  "%.6f\t%.6f\n", x[i], Tcn[i]  - Tex[i]);
    }

    std::fprintf(uerr_exp, "\n\n");
    std::fprintf(uerr_imp, "\n\n");
    std::fprintf(uerr_df,  "\n\n");
    std::fprintf(uerr_cn,  "\n\n");
  }

  std::fclose(uexp);
  std::fclose(udf);
  std::fclose(uimp);
  std::fclose(ucn);
  std::fclose(uex);

  std::fclose(uerr_exp);
  std::fclose(uerr_imp);
  std::fclose(uerr_df);
  std::fclose(uerr_cn);
}

// -----------------------------------------------------------------------------
// plot: profile plots for one dt (scheme by scheme)
// -----------------------------------------------------------------------------
void plot(double dt) {
  char tag[8];
  make_dt_tag(dt, tag);

  char cmd[512];

  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"scheme='ftcs_explicit'; tag='%s'; dt='%.6f'; outpng='plot/ftcs_explicit_%s.png'\" gnuplot_scripts/plot_scheme.gp",
    tag, dt, tag);
  std::system(cmd);

  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"scheme='dufort'; tag='%s'; dt='%.6f'; outpng='plot/dufort_%s.png'\" gnuplot_scripts/plot_scheme.gp",
    tag, dt, tag);
  std::system(cmd);

  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"scheme='ftcs_implicit'; tag='%s'; dt='%.6f'; outpng='plot/ftcs_implicit_%s.png'\" gnuplot_scripts/plot_scheme.gp",
    tag, dt, tag);
  std::system(cmd);

  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"scheme='cn'; tag='%s'; dt='%.6f'; outpng='plot/cn_%s.png'\" gnuplot_scripts/plot_scheme.gp",
    tag, dt, tag);
  std::system(cmd);
}

// -----------------------------------------------------------------------------
// plot error vs exact at one t_target (writes png using gnuplot script)
// -----------------------------------------------------------------------------
void plot_error_schemes(double dt, double t_target) {
  char tag[8];
  make_dt_tag(dt, tag);

  int idx = (int)lround(t_target / 0.1);

  char cmd[512];
  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"tag='%s'; dt=%.6f; idx=%d; tlabel='%.1f'; outpng='plot/error_schemes_%s_t%.1f.png'\" gnuplot_scripts/compare_error_schemes.gp",
    tag, dt, idx, t_target, tag, t_target);
  std::system(cmd);
}

// -----------------------------------------------------------------------------
// convergence table + plot
// -----------------------------------------------------------------------------
void convergence_study(double dx, int imax, double alpha,
                       double t0, double tb, double t_target) {
  const int ndt = 5;
  const double dt_list[ndt] = {0.10, 0.05, 0.02, 0.01, 0.005};

  char outdata[128];
  std::snprintf(outdata, sizeof(outdata), "data/convergence_t%03d.txt",
                (int)lround(t_target * 100.0));

  FILE* f = std::fopen(outdata, "w");
  if (!f) {
    std::perror("fopen");
    std::exit(1);
  }

  std::fprintf(f, "# convergence at t = %.2f hr\n", t_target);
  std::fprintf(f,
               "# dt\tL2_exp\tL2_imp\tL2_df\tL2_cn\t"
               "Linf_exp\tLinf_imp\tLinf_df\tLinf_cn\n");

  std::vector<double> x(imax + 1, 0.0);
  build_grid(imax, dx, x.data());

  double L = (imax - 1) * dx;

  std::vector<double> Tex(imax + 1, 0.0);
  std::vector<double> Texp(imax + 1, 0.0);
  std::vector<double> Timp(imax + 1, 0.0);
  std::vector<double> Tdf(imax + 1, 0.0);
  std::vector<double> Tcn(imax + 1, 0.0);

  for (int k = 0; k < ndt; k++) {
    double dt = dt_list[k];
    double F  = alpha * dt / (dx * dx);
    int nmax  = (int)lround(t_target / dt);

    FTCS_explicit(nmax, F, tb, t0, imax, Texp.data());
    FTCS_implicit(nmax, F, tb, t0, imax, Timp.data());
    DufortFrankel(nmax, F, tb, t0, imax, Tdf.data());
    CrankNicolson(nmax, F, tb, t0, imax, Tcn.data());

    exact_solution_profile(imax, x.data(), t_target, alpha, L, tb, t0, 200, Tex.data());

    double L2_exp = err_L2(imax, Texp.data(), Tex.data());
    double L2_imp = err_L2(imax, Timp.data(), Tex.data());
    double L2_df  = err_L2(imax, Tdf.data(),  Tex.data());
    double L2_cn  = err_L2(imax, Tcn.data(),  Tex.data());

    double Li_exp = err_Linf(imax, Texp.data(), Tex.data());
    double Li_imp = err_Linf(imax, Timp.data(), Tex.data());
    double Li_df  = err_Linf(imax, Tdf.data(),  Tex.data());
    double Li_cn  = err_Linf(imax, Tcn.data(),  Tex.data());

    std::fprintf(f,
      "%.6f\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
      dt, L2_exp, L2_imp, L2_df, L2_cn, Li_exp, Li_imp, Li_df, Li_cn);
  }

  std::fclose(f);
}

void plot_convergence(double t_target) {
  char infile[128];
  std::snprintf(infile, sizeof(infile), "data/convergence_t%03d.txt",
                (int)lround(t_target * 100.0));

  char cmd[512];
  std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"infile='%s'; tlabel='%.2f'; outpng='plot/convergence_t%.2f.png'\" gnuplot_scripts/convergence.gp",
    infile, t_target, t_target);
  std::system(cmd);
}
