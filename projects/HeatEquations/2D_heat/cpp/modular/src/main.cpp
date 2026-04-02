#include "heat_2d.hpp"



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
