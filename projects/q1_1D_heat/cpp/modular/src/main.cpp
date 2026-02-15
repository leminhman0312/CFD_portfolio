#include "headers.h"
#include "helpers.h"
#include "driver.h"

int main() {
  int imax = 21;
  double dx = 0.05;
  double dt = 0.01;
  double alpha = 0.1;
  double t0 = 100.0;
  double tb = 300.0;

  double F = alpha * dt / (dx * dx);


  create_dir();
  sim(dx, dt, imax, t0, tb, F);
  plot(dt);
  plot_error_schemes(dt, 0.4);
  convergence_study(dx, imax, alpha, t0, tb, 0.4);
  plot_convergence(0.4);


  std::printf("Done\n");
  return 0;
}
