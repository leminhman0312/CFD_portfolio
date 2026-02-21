// ============================================================================
// 2D Heat Conduction Solver C++
//
// Solves the transient two dimensional heat equation
//
//     ∂T/∂t = α ( ∂²T/∂x² + ∂²T/∂y² )
//
// Numerical methods
//   FTCS explicit scheme
//   Implicit ADI scheme using Thomas tridiagonal solvers
//
// Run modes
//   ./heat2d sim
//   ./heat2d convergence
//   ./heat2d all
//
// Directories created automatically
//   data
//   plot
//   gnuplot_scripts
// ============================================================================

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>

using std::vector;

// ---------------------------------------------------------------------------
// Prototypes
// ---------------------------------------------------------------------------

// Directory helpers
void ensure_dir(const char* name);
void ensure_project_dirs();

// Error norms
double error_L2(const vector<vector<double>>& u,
                const vector<vector<double>>& uref);
double error_Linf(const vector<vector<double>>& u,
                  const vector<vector<double>>& uref);

// High-level drivers
void sim(const vector<vector<double>>& u0,
         double t_end,
         double deltax, double deltay,
         double alpha,
         double dt_implicit,
         double dt_explicit_given,
         double t1, double t2, double t3, double t4);

void convergence(const vector<vector<double>>& u0,
                 double t_end,
                 double deltax, double deltay,
                 double alpha,
                 double t1, double t2, double t3, double t4);

// Solvers
vector<vector<double>> initializeField(int imax, int jmax,
                                       double t0, double t1, double t2,
                                       double t3, double t4);

vector<vector<double>> FTCS_Explicit(const vector<vector<double>>& u0,
                                     double t_end,
                                     double deltax, double deltay,
                                     double dt,
                                     double alpha,
                                     double t1, double t2, double t3, double t4);

vector<vector<double>> FTCS_implicit_ADI(const vector<vector<double>>& u0,
                                        int nmax,
                                        double deltax, double deltay,
                                        double dt,
                                        double alpha,
                                        double t1, double t2, double t3, double t4);

// Linear solver
void thomasTriDiagonal(int n,
                       const vector<double>& a,
                       const vector<double>& b,
                       vector<double>& c,
                       const vector<double>& d,
                       vector<double>& u);

// IO and plotting
void write_field_xyz(const char* filename,
                     const vector<vector<double>>& u,
                     double deltax,
                     double deltay);

void plotContourMatlabLike(const char* datafile,
                           const char* outpng,
                           double time_hr,
                           const char* scheme);

void plot_time_convergence(const char* infile,
                           const char* outpng);

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {

  ensure_project_dirs();

  const bool do_sim  = (argc == 1) || (argc > 1 && std::strcmp(argv[1], "sim") == 0);
  const bool do_conv = (argc > 1 && std::strcmp(argv[1], "convergence") == 0);
  const bool do_all  = (argc > 1 && std::strcmp(argv[1], "all") == 0);

  const double t_end = 0.5;

  const double deltax = 0.1;
  const double deltay = 0.1;
  const double alpha  = 0.645;

  const double dt_implicit = 0.01;
  const double dt_explicit_given = 0.01;

  const double xmin = 0.0, xmax = 3.5;
  const double ymin = 0.0, ymax = 3.5;

  const int imax = (int)std::ceil((xmax - xmin) / deltax + 1.0);
  const int jmax = (int)std::ceil((ymax - ymin) / deltay + 1.0);

  const double t0 = 0.0;
  const double t1 = 200.0;
  const double t2 = 200.0;
  const double t3 = 0.0;
  const double t4 = 0.0;

  auto u0 = initializeField(imax, jmax, t0, t1, t2, t3, t4);

  if (do_sim || do_all) {
    sim(u0, t_end, deltax, deltay,
        alpha, dt_implicit, dt_explicit_given,
        t1, t2, t3, t4);
  }

  if (do_conv || do_all) {
    convergence(u0, t_end, deltax, deltay,
                alpha, t1, t2, t3, t4);
  }

  std::printf("\nDone\n");
  return 0;
}

// ---------------------------------------------------------------------------
// Directory helpers
// ---------------------------------------------------------------------------

void ensure_dir(const char* name) {
  struct stat st;
  if (stat(name, &st) != 0) {
    mkdir(name, 0755);
  }
}

void ensure_project_dirs() {
  ensure_dir("data");
  ensure_dir("plot");
  ensure_dir("gnuplot_scripts");
}

// ---------------------------------------------------------------------------
// Error norms
// ---------------------------------------------------------------------------

double error_L2(const vector<vector<double>>& u,
                const vector<vector<double>>& uref) {

  const int imax = (int)u.size() - 1;
  const int jmax = (int)u[1].size() - 1;

  double sum = 0.0;
  int n = 0;

  for (int i = 2; i <= imax - 1; i++) {
    for (int j = 2; j <= jmax - 1; j++) {
      const double e = u[i][j] - uref[i][j];
      sum += e * e;
      n++;
    }
  }

  return std::sqrt(sum / (double)n);
}

double error_Linf(const vector<vector<double>>& u,
                  const vector<vector<double>>& uref) {

  const int imax = (int)u.size() - 1;
  const int jmax = (int)u[1].size() - 1;

  double max_err = 0.0;

  for (int i = 2; i <= imax - 1; i++) {
    for (int j = 2; j <= jmax - 1; j++) {
      const double e = std::fabs(u[i][j] - uref[i][j]);
      if (e > max_err) max_err = e;
    }
  }

  return max_err;
}

// ---------------------------------------------------------------------------
// High-level simulation driver
// ---------------------------------------------------------------------------

void sim(const vector<vector<double>>& u0,
         double t_end,
         double deltax, double deltay,
         double alpha,
         double dt_implicit,
         double dt_explicit_given,
         double t1, double t2, double t3, double t4) {

  std::printf("\nSIMULATING\n");

  write_field_xyz("data/initial.dat", u0, deltax, deltay);
  plotContourMatlabLike("data/initial.dat",
                        "plot/contour_initial.png",
                        0.0, "Initial conditions");

  const int nmax_implicit = (int)std::lround(t_end / dt_implicit);

  auto u_implicit =
    FTCS_implicit_ADI(u0, nmax_implicit,
                      deltax, deltay, dt_implicit,
                      alpha, t1, t2, t3, t4);

  write_field_xyz("data/implicit.dat", u_implicit, deltax, deltay);
  plotContourMatlabLike("data/implicit.dat",
                        "plot/contour_implicit.png",
                        t_end, "Implicit ADI");

  const double fx = alpha * dt_explicit_given / (deltax * deltax);
  const double fy = alpha * dt_explicit_given / (deltay * deltay);
  const double sum = fx + fy;

  if (sum <= 0.5) {

    auto u_explicit =
      FTCS_Explicit(u0, t_end, deltax, deltay,
                    dt_explicit_given, alpha,
                    t1, t2, t3, t4);

    write_field_xyz("data/explicit.dat", u_explicit, deltax, deltay);
    plotContourMatlabLike("data/explicit.dat",
                          "plot/contour_explicit.png",
                          t_end, "Explicit FTCS");

  } else {

    const double dt_safe =
      0.5 / (alpha * (1.0/(deltax*deltax) + 1.0/(deltay*deltay)));

    auto u_fail =
      FTCS_Explicit(u0, t_end, deltax, deltay,
                    dt_explicit_given, alpha,
                    t1, t2, t3, t4);

    auto u_safe =
      FTCS_Explicit(u0, t_end, deltax, deltay,
                    dt_safe, alpha,
                    t1, t2, t3, t4);

    write_field_xyz("data/explicit_failed.dat", u_fail, deltax, deltay);
    write_field_xyz("data/explicit_safe.dat",   u_safe, deltax, deltay);

    plotContourMatlabLike("data/explicit_failed.dat",
                          "plot/contour_explicit_failed.png",
                          t_end, "Explicit unstable dt");

    plotContourMatlabLike("data/explicit_safe.dat",
                          "plot/contour_explicit_safe.png",
                          t_end, "Explicit safe dt");
  }
}

// ---------------------------------------------------------------------------
// Time convergence study (implicit ADI)
// ---------------------------------------------------------------------------

void convergence(const vector<vector<double>>& u0,
                 double t_end,
                 double deltax, double deltay,
                 double alpha,
                 double t1, double t2, double t3, double t4) {

  std::printf("\nTIME CONVERGENCE STUDY (implicit ADI)\n");
  std::printf("dt        L2 error      Linf error\n");
  std::printf("----------------------------------\n");

  const double dt_start = 0.04;
  const double dt_ratio = 0.5;
  const double dt_min   = 1e-6;

  vector<double> dt_list;
  for (double dt = dt_start; dt >= dt_min; dt *= dt_ratio) {
    dt_list.push_back(dt);
  }

  const double dt_ref = dt_list.back();
  const int nmax_ref = (int)std::lround(t_end / dt_ref);

  auto u_ref =
    FTCS_implicit_ADI(u0, nmax_ref,
                      deltax, deltay, dt_ref,
                      alpha, t1, t2, t3, t4);

  FILE* f = std::fopen("data/time_convergence.dat", "w");
  if (!f) {
    std::perror("data/time_convergence.dat");
    std::exit(EXIT_FAILURE);
  }

  for (double dt : dt_list) {

    const int nmax = (int)std::lround(t_end / dt);

    auto u =
      FTCS_implicit_ADI(u0, nmax,
                        deltax, deltay, dt,
                        alpha, t1, t2, t3, t4);

    const double L2   = error_L2(u, u_ref);
    const double Linf = error_Linf(u, u_ref);

    std::printf("%-8.5f  %.6e  %.6e\n", dt, L2, Linf);
    std::fprintf(f, "%.10f %.10e %.10e\n", dt, L2, Linf);
  }

  std::fclose(f);

  std::printf("\nWrote data/time_convergence.dat\n");
  plot_time_convergence("data/time_convergence.dat",
                        "plot/time_convergence.png");
}

// ---------------------------------------------------------------------------
// Initialize field
// ---------------------------------------------------------------------------

vector<vector<double>> initializeField(int imax, int jmax,
                                       double t0, double t1, double t2,
                                       double t3, double t4) {

  vector<vector<double>> u(imax + 1, vector<double>(jmax + 1, t0));

  for (int j = 1; j <= jmax; j++) {
    u[1][j]    = t2;
    u[imax][j] = t4;
  }
  for (int i = 1; i <= imax; i++) {
    u[i][1]    = t1;
    u[i][jmax] = t3;
  }

  return u;
}

// ---------------------------------------------------------------------------
// Explicit FTCS solver
// ---------------------------------------------------------------------------

vector<vector<double>> FTCS_Explicit(const vector<vector<double>>& u0,
                                     double t_end,
                                     double deltax, double deltay,
                                     double dt,
                                     double alpha,
                                     double t1, double t2, double t3, double t4) {

  vector<vector<double>> u = u0;

  const int imax = (int)u.size() - 1;
  const int jmax = (int)u[1].size() - 1;

  vector<vector<double>> u_new(imax + 1, vector<double>(jmax + 1, 0.0));

  const int nfull = (int)std::floor(t_end / dt);
  const double dt_last = t_end - nfull * dt;
  const int nsteps = (dt_last > 0.0) ? nfull + 1 : nfull;

  for (int n = 1; n <= nsteps; n++) {

    const double dt_n = (n == nsteps && dt_last > 0.0) ? dt_last : dt;

    const double fx = alpha * dt_n / (deltax * deltax);
    const double fy = alpha * dt_n / (deltay * deltay);

    for (int i = 1; i <= imax; i++) {
      for (int j = 1; j <= jmax; j++) {
        u_new[i][j] = u[i][j];
      }
    }

    for (int i = 2; i <= imax - 1; i++) {
      for (int j = 2; j <= jmax - 1; j++) {
        u_new[i][j] =
          (1.0 - 2.0 * fx - 2.0 * fy) * u[i][j]
          + fx * (u[i + 1][j] + u[i - 1][j])
          + fy * (u[i][j + 1] + u[i][j - 1]);
      }
    }

    for (int j = 1; j <= jmax; j++) {
      u_new[1][j]    = t2;
      u_new[imax][j] = t4;
    }
    for (int i = 1; i <= imax; i++) {
      u_new[i][1]    = t1;
      u_new[i][jmax] = t3;
    }

    u.swap(u_new);
  }

  return u;
}

// ---------------------------------------------------------------------------
// Implicit ADI solver
// ---------------------------------------------------------------------------

vector<vector<double>> FTCS_implicit_ADI(const vector<vector<double>>& u0,
                                        int nmax,
                                        double deltax, double deltay,
                                        double dt,
                                        double alpha,
                                        double t1, double t2, double t3, double t4) {

  vector<vector<double>> u = u0;

  const int imax = (int)u.size() - 1;
  const int jmax = (int)u[1].size() - 1;

  vector<vector<double>> u_dummy(imax + 1, vector<double>(jmax + 1, 0.0));

  const double fx = alpha * dt / (deltax * deltax);
  const double fy = alpha * dt / (deltay * deltay);

  vector<double> ax(imax + 1, 0.0), bx(imax + 1, 0.0);
  vector<double> cx(imax + 1, 0.0), dx(imax + 1, 0.0);

  for (int i = 2; i <= imax - 1; i++) {
    ax[i] = -fx / 2.0;
    bx[i] = -fx / 2.0;
    dx[i] =  1.0 + fx;
  }
  dx[1] = dx[imax] = 1.0;

  vector<double> ay(jmax + 1, 0.0), by(jmax + 1, 0.0);
  vector<double> cy(jmax + 1, 0.0), dy(jmax + 1, 0.0);

  for (int j = 2; j <= jmax - 1; j++) {
    ay[j] = -fy / 2.0;
    by[j] = -fy / 2.0;
    dy[j] =  1.0 + fy;
  }
  dy[1] = dy[jmax] = 1.0;

  for (int n = 1; n <= nmax; n++) {

    for (int j = 1; j <= jmax; j++) {
      u[1][j]    = t2;
      u[imax][j] = t4;
    }
    for (int i = 1; i <= imax; i++) {
      u[i][1]    = t1;
      u[i][jmax] = t3;
    }

    // X sweep: solve for u_dummy
    for (int j = 2; j <= jmax - 1; j++) {

      for (int i = 2; i <= imax - 1; i++) {
        cx[i] = (1.0 - fy) * u[i][j]
              + (fy / 2.0) * (u[i][j + 1] + u[i][j - 1]);
      }

      cx[2]       += (fx / 2.0) * u[1][j];
      cx[imax - 1] += (fx / 2.0) * u[imax][j];

      vector<double> solx(imax + 1, 0.0);
      solx[1] = u[1][j];
      solx[imax] = u[imax][j];

      cx[1] = u[1][j];
      cx[imax] = u[imax][j];

      thomasTriDiagonal(imax, ax, bx, cx, dx, solx);

      for (int i = 1; i <= imax; i++) {
        u_dummy[i][j] = solx[i];
      }
    }

    // Enforce BCs on u_dummy
    for (int j = 1; j <= jmax; j++) {
      u_dummy[1][j]    = t2;
      u_dummy[imax][j] = t4;
    }
    for (int i = 1; i <= imax; i++) {
      u_dummy[i][1]    = t1;
      u_dummy[i][jmax] = t3;
    }

    // Y sweep: solve for u
    for (int i = 2; i <= imax - 1; i++) {

      for (int j = 2; j <= jmax - 1; j++) {
        cy[j] = (1.0 - fx) * u_dummy[i][j]
              + (fx / 2.0) * (u_dummy[i + 1][j] + u_dummy[i - 1][j]);
      }

      cy[2]       += (fy / 2.0) * u_dummy[i][1];
      cy[jmax - 1] += (fy / 2.0) * u_dummy[i][jmax];

      vector<double> soly(jmax + 1, 0.0);
      soly[1] = u_dummy[i][1];
      soly[jmax] = u_dummy[i][jmax];

      thomasTriDiagonal(jmax, ay, by, cy, dy, soly);

      for (int j = 1; j <= jmax; j++) {
        u[i][j] = soly[j];
      }
    }

    // Enforce BCs on u
    for (int j = 1; j <= jmax; j++) {
      u[1][j]    = t2;
      u[imax][j] = t4;
    }
    for (int i = 1; i <= imax; i++) {
      u[i][1]    = t1;
      u[i][jmax] = t3;
    }
  }

  return u;
}

// ---------------------------------------------------------------------------
// Thomas tridiagonal solver
// ---------------------------------------------------------------------------

void thomasTriDiagonal(int n,
                       const vector<double>& a,
                       const vector<double>& b,
                       vector<double>& c,
                       const vector<double>& d,
                       vector<double>& u) {

  vector<double> dprime(n + 1, 0.0);
  vector<double> cprime(n + 1, 0.0);

  dprime[1] = d[1];
  cprime[1] = c[1];

  for (int i = 2; i <= n; i++) {
    dprime[i] = d[i] - (b[i] * a[i - 1]) / dprime[i - 1];
    cprime[i] = c[i] - (cprime[i - 1] * b[i]) / dprime[i - 1];
  }

  u[n] = cprime[n] / dprime[n];

  for (int i = n - 1; i >= 1; i--) {
    u[i] = (cprime[i] - a[i] * u[i + 1]) / dprime[i];
  }
}

// ---------------------------------------------------------------------------
// Write field to file
// ---------------------------------------------------------------------------

void write_field_xyz(const char* filename,
                     const vector<vector<double>>& u,
                     double deltax,
                     double deltay) {

  const int imax = (int)u.size() - 1;
  const int jmax = (int)u[1].size() - 1;

  FILE* f = std::fopen(filename, "w");
  if (!f) return;

  for (int j = 1; j <= jmax; j++) {
    const double y = (j - 1) * deltay;
    for (int i = 1; i <= imax; i++) {
      const double x = (i - 1) * deltax;
      std::fprintf(f, "%.8f %.8f %.10f\n", x, y, u[i][j]);
    }
    std::fprintf(f, "\n");
  }

  std::fclose(f);
}

// ---------------------------------------------------------------------------
// Plot helpers
// ---------------------------------------------------------------------------

void plotContourMatlabLike(const char* datafile,
                           const char* outpng,
                           double time_hr,
                           const char* scheme) {

  char cmd[1024];
  std::snprintf(cmd, sizeof(cmd),
                "gnuplot -e \""
                "datafile=\\\"%s\\\"; "
                "outpng=\\\"%s\\\"; "
                "tlabel=%.3f; "
                "scheme=\\\"%s\\\"\" "
                "gnuplot_scripts/plot_contour_2d_matlab_like.gp",
                datafile, outpng, time_hr, scheme);

  std::printf("\n%s\n", scheme);
  std::system(cmd);
}

void plot_time_convergence(const char* datafile,
                           const char* outpng) {

  char cmd[1024];
  std::snprintf(cmd, sizeof(cmd),
                "gnuplot -e \""
                "datafile=\\\"%s\\\"; "
                "outpng=\\\"%s\\\"\" "
                "gnuplot_scripts/plot_convergence.gp",
                datafile, outpng);

  std::printf("\nTime convergence plot\n");
  std::system(cmd);
}
