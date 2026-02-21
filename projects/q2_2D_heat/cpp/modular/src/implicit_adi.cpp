#include "heat_2d.hpp"


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
