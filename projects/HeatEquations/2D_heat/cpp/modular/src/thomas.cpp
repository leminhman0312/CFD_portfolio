#include "heat_2d.hpp"


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
