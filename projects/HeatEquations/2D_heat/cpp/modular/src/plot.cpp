#include "heat_2d.hpp"


void plotContourMatlabLike(const char* datafile,
                           const char* outpng,
                           double time_hr,
                           const char* scheme) {

  char cmd[1024];
  std::snprintf(cmd, sizeof(cmd),
                "gnuplot -e \""
                "datafile=\\\"%s\\\"; "
                "outpng=\\\"%s\\\"; "
                "tlabel=%.3f; "
                "scheme=\\\"%s\\\"\" "
                "gnuplot_scripts/plot_contour_2d_matlab_like.gp",
                datafile, outpng, time_hr, scheme);

  std::printf("\n%s\n", scheme);
  std::system(cmd);
}

void plot_time_convergence(const char* datafile,
                           const char* outpng) {

  char cmd[1024];
  std::snprintf(cmd, sizeof(cmd),
                "gnuplot -e \""
                "datafile=\\\"%s\\\"; "
                "outpng=\\\"%s\\\"\" "
                "gnuplot_scripts/plot_convergence.gp",
                datafile, outpng);

  std::printf("\nTime convergence plot\n");
  std::system(cmd);
}
