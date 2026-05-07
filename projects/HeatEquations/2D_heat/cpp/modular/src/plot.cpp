#include "heat_2d.hpp"

void plotContourMatlabLike(const char* datafile, const char* outpng,
                           double time_hr, const char* scheme) {
  char cmd[1024];
  std::snprintf(cmd, sizeof(cmd),
                "python3 plot_heat_2d.py contour \"%s\" \"%s\" %.6f \"%s\"",
                datafile, outpng, time_hr, scheme);

  std::printf("\n%s\n", scheme);

  int status = std::system(cmd);
  if (status != 0) {
    std::fprintf(stderr, "Python contour plotting failed.\n");
    std::exit(EXIT_FAILURE);
  }
}

void plot_time_convergence(const char* datafile, const char* outpng) {
  char cmd[1024];
  std::snprintf(cmd, sizeof(cmd),
                "python3 plot_heat_2d.py convergence \"%s\" \"%s\"", datafile,
                outpng);

  std::printf("\nTime convergence plot\n");

  int status = std::system(cmd);
  if (status != 0) {
    std::fprintf(stderr, "Python convergence plotting failed.\n");
    std::exit(EXIT_FAILURE);
  }
}
