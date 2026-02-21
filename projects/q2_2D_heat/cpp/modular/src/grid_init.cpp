#include "heat_2d.hpp"


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
