#include "heat_2d.hpp"


void write_field_xyz(const char* filename,
                     const vector<vector<double>>& u,
                     double deltax,
                     double deltay) {

  const int imax = (int)u.size() - 1;
  const int jmax = (int)u[1].size() - 1;

  FILE* f = std::fopen(filename, "w");
  if (!f) return;

  for (int j = 1; j <= jmax; j++) {
    const double y = (j - 1) * deltay;
    for (int i = 1; i <= imax; i++) {
      const double x = (i - 1) * deltax;
      std::fprintf(f, "%.8f %.8f %.10f\n", x, y, u[i][j]);
    }
    std::fprintf(f, "\n");
  }

  std::fclose(f);
}
