#include "heat_2d.hpp"


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
