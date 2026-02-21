#include "heat_2d.hpp"


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
