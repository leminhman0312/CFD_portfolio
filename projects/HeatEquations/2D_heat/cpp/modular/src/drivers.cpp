#include "heat_2d.hpp"


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
