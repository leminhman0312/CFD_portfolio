#include "methods.h"

void thomasTriDiagonal(int N,
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

void exact_solution_profile(int imax, const double x[], double t,
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

void FTCS_explicit(int nmax, double F, double tb, double t0,
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

void implicit_interior_one_step(double F, double tb, double t0,
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

void DufortFrankel(int nmax, double F, double tb, double t0,
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

void FTCS_implicit(int nmax, double F, double tb, double t0,
                   int imax, double Tout[]) {
  std::vector<double> a(imax + 1), b(imax + 1),
                      c(imax + 1), d(imax + 1), un(imax + 1);

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
  }

  for (int i = 1; i <= imax; i++) Tout[i] = un[i];
}

void CrankNicolson(int nmax, double F, double tb, double t0,
                   int imax, double Tout[]) {
  const double dnc = F / 2.0;

  std::vector<double> a(imax + 1), b(imax + 1),
                      c(imax + 1), d(imax + 1),
                      u0(imax + 1), uhalf(imax + 1);

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
      uhalf[i] = u0[i] + dnc * (u0[i + 1] - 2.0 * u0[i] + u0[i - 1]);
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
