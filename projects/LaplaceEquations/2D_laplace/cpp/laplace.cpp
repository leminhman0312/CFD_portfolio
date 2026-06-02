#include <cmath>
#include <cstdio>

const int MAX = 128;

// prototypes
void grid_refinement_study(double error_max, int iter_max);
void init(double psi[][MAX], double psi_old[][MAX], int imax, int jmax, double dx);
void print_array(double psi[][MAX], int imax, int jmax);
void writeFile(const char filename[], double psi[][MAX], int imax, int jmax,
    double dx, double dy);


static FILE *openHistoryFile(const char filename[]);


static void writeHistory(FILE *fp, int iter, double error);


static bool diverged(double error);

int Point_GS(double psi[][MAX], double psi_old[][MAX], int imax, int jmax,
             double dx, double dy, double error_max, int iter_max,
             const char history_file[]);

int Point_SOR(double w, double psi[][MAX], double psi_old[][MAX], int imax,
              int jmax, double dx, double dy, double error_max, int iter_max,
              const char history_file[]);

int Line_GS(double psi[][MAX], double psi_old[][MAX], int imax, int jmax,
            double dx, double dy, double error_max, int iter_max,
            const char history_file[]);

int Line_SOR(double w, double psi[][MAX], double psi_old[][MAX], int imax,
             int jmax, double dx, double dy, double error_max, int iter_max,
             const char history_file[]);

void sweep_w_LineSOR(int imax, int jmax, double dx, double dy,
                     double error_max, int iter_max);

static void thomasTriDiagonal(int N, const double a[], const double b[],
                              const double rhs[], const double diag[],
                              double u[], double diag_prime[],
                              double rhs_prime[]);


void writeCenterline(const char filename[],
                     double psi[][MAX],
                     int imax, int jmax,
                     double dx, double dy);


// functions
void grid_refinement_study(double error_max, int iter_max)
{
    double grids[] = {0.4, 0.2, 0.1};

    for(int k = 0; k < 3; k++){

        // back out dx and dy from the grid arrays
        double dx = grids[k];
        double dy = grids[k]; // assume dx = dy

        //calculate imax and jmax from dx and dy
        int imax = int(6.0 / dx) + 1;
        int jmax = int(4.0 / dy) + 1;

        double psi[MAX][MAX];
        double psi_old[MAX][MAX];

        double w = 1.35;

        printf("\nGrid refinement study: dx = dy = %f\n", dx);
        printf("imax = %d, jmax = %d\n", imax, jmax);

        Line_SOR(w, psi, psi_old,
                 imax, jmax,
                 dx, dy,
                 error_max, iter_max,
                 NULL);

        if(dx == 0.4){
            writeCenterline("data/centerline_dx04.dat",
                            psi, imax, jmax, dx, dy);
        }
        else if(dx == 0.2){
            writeCenterline("data/centerline_dx02.dat",
                            psi, imax, jmax, dx, dy);
        }
        else if(dx == 0.1){
            writeCenterline("data/centerline_dx01.dat",
                            psi, imax, jmax, dx, dy);
        }
    }
}

void writeCenterline(const char filename[],
                     double psi[][MAX],
                     int imax, int jmax,
                     double dx, double dy)
{
    FILE *fp = fopen(filename, "w");

    if(fp == NULL){
            printf("Failed to open file: %s\n", filename);
            return;
    }

    int jmid = int(2.0 / dy) + 1;

    fprintf(fp, "# dx %f\n", dx);
    fprintf(fp, "# dy %f\n", dy);
    fprintf(fp, "# imax %d\n", imax);
    fprintf(fp, "# jmax %d\n", jmax);
    fprintf(fp, "# x psi\n");

    for(int i = 1; i <= imax; i++){

        double x = (i - 1) * dx;

        fprintf(fp, "%f %f\n", x, psi[i][jmid]);
    }

    fclose(fp);
}





static FILE *openHistoryFile(const char filename[])
{
    if(filename == NULL){
        return NULL;
    }

    FILE *fp = fopen(filename, "w");

    if(fp == NULL){
        printf("Failed to open history file: %s\n", filename);
        return NULL;
    }

    fprintf(fp, "# iter max_error\n");
    return fp;
}

static void writeHistory(FILE *fp, int iter, double error)
{
    if(fp != NULL){
        fprintf(fp, "%d %e\n", iter, error);
    }
}

static bool diverged(double error)
{
    return !std::isfinite(error);
}

// relaxation sweep for Line SOR
void sweep_w_LineSOR(int imax, int jmax, double dx, double dy,
                     double error_max, int iter_max)
{
    FILE *fp = fopen("data/line_sor_w_sweep.dat", "w");

    if(fp == NULL){
        printf("Failed to open data/line_sor_w_sweep.dat\n");
        return;
    }

    fprintf(fp, "# w iterations\n");
    printf("\nRelaxation iteration study\n");

    for(double w = 0.0; w <= 2.0001; w += 0.05){

        double psi[MAX][MAX];
        double psi_old[MAX][MAX];

        int iter = Line_SOR(w, psi, psi_old,
                            imax, jmax,
                            dx, dy,
                            error_max, iter_max,
                            NULL);

        fprintf(fp, "%6.2f %d\n", w, iter);
        printf("w = %6.2f   iter = %d\n", w, iter);
    }

    fclose(fp);
}

// init
void init(double psi[][MAX], double psi_old[][MAX], int imax, int jmax,
          double dx)
{
    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax; j++){
            psi[i][j] = 0.0;
            psi_old[i][j] = 0.0;
        }
    }

    int left_inlet = int(1.0 / dx) + 1;
    int right_inlet = int(1.2 / dx) + 1;

    // left boundary
    for(int j = 1; j <= jmax; j++){
        psi[1][j] = 0.0;
    }

    // top boundary
    for(int i = 1; i <= imax; i++){
        psi[i][jmax] = 0.0;
    }

    // bottom boundary
    for(int i = 1; i <= left_inlet; i++){
        psi[i][1] = 0.0;
    }

    for(int i = right_inlet; i <= imax; i++){
        psi[i][1] = 100.0;
    }
}

// print array
void print_array(double psi[][MAX], int imax, int jmax)
{
    for(int j = jmax; j >= 1; j--){
        for(int i = 1; i <= imax; i++){
            printf("%6.1f ", psi[i][j]);
        }
        printf("\n");
    }
}

// Point GS
int Point_GS(double psi[][MAX], double psi_old[][MAX], int imax, int jmax,
             double dx, double dy, double error_max, int iter_max,
             const char history_file[])
{
    init(psi, psi_old, imax, jmax, dx);

    FILE *fp = openHistoryFile(history_file);

    double error = 10.0;
    double diff = 0.0;
    int iter = 0;

    double beta = dx / dy;
    double beta2 = beta * beta;
    double coeff = 1.0 / (2.0 * (1.0 + beta2));

    printf("\nSolving Point GS\n");

    while(error > error_max && iter < iter_max){

        iter++;
        error = 0.0;

        for(int i = 2; i <= imax - 1; i++){
            for(int j = 2; j <= jmax - 1; j++){

                psi_old[i][j] = psi[i][j];

                psi[i][j] = coeff * (psi[i + 1][j] + psi[i - 1][j]
                                    + beta2 * (psi[i][j + 1] + psi[i][j - 1]));

                diff = fabs(psi[i][j] - psi_old[i][j]);

                if(diff > error){
                    error = diff;
                }
            }
        }

        // right boundary: dpsi/dx = 0
        for(int j = 2; j <= jmax - 1; j++){
            psi[imax][j] = psi[imax - 1][j];
        }

        // divergence check
        if(diverged(error)){
            printf("Point GS diverged at iter = %d, error = %e\n", iter, error);
            if(fp != NULL) fclose(fp);
            return iter_max;
        }

        writeHistory(fp, iter, error);
    }

    if(fp != NULL) fclose(fp);

    printf("Point GS iter = %d\n", iter);
    printf("Point GS error = %f\n", error);

    return iter;
}

// Point SOR
int Point_SOR(double w, double psi[][MAX], double psi_old[][MAX], int imax,
              int jmax, double dx, double dy, double error_max, int iter_max,
              const char history_file[])
{
    init(psi, psi_old, imax, jmax, dx);

    FILE *fp = openHistoryFile(history_file);

    double error = 10.0;
    double diff = 0.0;
    int iter = 0;

    double beta = dx / dy;
    double beta2 = beta * beta;
    double coeff = 1.0 / (2.0 * (1.0 + beta2));

    printf("\nSolving Point SOR\n");
    printf("w = %f\n", w);

    while(error > error_max && iter < iter_max){

        iter++;
        error = 0.0;

        for(int i = 2; i <= imax - 1; i++){
            for(int j = 2; j <= jmax - 1; j++){

                psi_old[i][j] = psi[i][j];

                double residual = psi[i + 1][j] + psi[i - 1][j]
                                + beta2 * (psi[i][j + 1] + psi[i][j - 1])
                                - 2.0 * (1.0 + beta2) * psi[i][j];

                psi[i][j] += w * coeff * residual;

                diff = fabs(psi[i][j] - psi_old[i][j]);

                if(diff > error){
                    error = diff;
                }
            }
        }

        // right boundary: dpsi/dx = 0
        for(int j = 2; j <= jmax - 1; j++){
            psi[imax][j] = psi[imax - 1][j];
        }

        // divergence check
        if(diverged(error)){
            printf("Point SOR diverged at iter = %d, error = %e\n", iter, error);
            if(fp != NULL) fclose(fp);
            return iter_max;
        }

        writeHistory(fp, iter, error);
    }

    if(fp != NULL) fclose(fp);

    printf("Point SOR iter = %d\n", iter);
    printf("Point SOR error = %f\n", error);

    return iter;
}

// Thomas tridiagonal solver
// Unknowns are indexed from 2 to N.
// b = lower diagonal, diag = main diagonal, a = upper diagonal, rhs = RHS.
static void thomasTriDiagonal(int N, const double a[], const double b[],
                              const double rhs[], const double diag[],
                              double u[], double diag_prime[],
                              double rhs_prime[])
{
    diag_prime[2] = diag[2];
    rhs_prime[2] = rhs[2];

    for(int i = 3; i <= N; i++){
        diag_prime[i] = diag[i] - (b[i] * a[i - 1]) / diag_prime[i - 1];
        rhs_prime[i] = rhs[i] - (rhs_prime[i - 1] * b[i]) / diag_prime[i - 1];
    }

    u[N] = rhs_prime[N] / diag_prime[N];

    for(int i = N - 1; i >= 2; i--){
        u[i] = (rhs_prime[i] - a[i] * u[i + 1]) / diag_prime[i];
    }
}

// Line GS
int Line_GS(double psi[][MAX], double psi_old[][MAX], int imax, int jmax,
            double dx, double dy, double error_max, int iter_max,
            const char history_file[])
{
    init(psi, psi_old, imax, jmax, dx);

    FILE *fp = openHistoryFile(history_file);

    double beta = dx / dy;
    double beta2 = beta * beta;

    int N = imax - 1;

    double a[MAX], b[MAX], rhs[MAX], diag[MAX], u[MAX];
    double diag_prime[MAX], rhs_prime[MAX];

    double error = 10.0;
    double diff = 0.0;
    int iter = 0;

    printf("\nSolving Line GS\n");

    while(error > error_max && iter < iter_max){

        iter++;
        error = 0.0;

        for(int j = 2; j <= jmax - 1; j++){

            for(int i = 2; i <= imax - 1; i++){
                psi_old[i][j] = psi[i][j];
            }

            for(int i = 2; i <= imax - 1; i++){
                b[i] = -1.0;
                diag[i] = 2.0 * (1.0 + beta2);
                a[i] = -1.0;

                rhs[i] = beta2 * (psi[i][j - 1] + psi[i][j + 1]);
            }

            b[2] = 0.0;
            rhs[2] += psi[1][j];

            a[imax - 1] = 0.0;
            rhs[imax - 1] += psi[imax][j];

            thomasTriDiagonal(N, a, b, rhs, diag, u, diag_prime, rhs_prime);

            for(int i = 2; i <= imax - 1; i++){
                psi[i][j] = u[i];

                diff = fabs(psi[i][j] - psi_old[i][j]);

                if(diff > error){
                    error = diff;
                }
            }
        }

        // right boundary: dpsi/dx = 0
        for(int j = 2; j <= jmax - 1; j++){
            psi[imax][j] = psi[imax - 1][j];
        }

        // divergence check
        if(diverged(error)){
            printf("Line GS diverged at iter = %d, error = %e\n", iter, error);
            if(fp != NULL) fclose(fp);
            return iter_max;
        }

        writeHistory(fp, iter, error);
    }

    if(fp != NULL) fclose(fp);

    printf("Line GS iter = %d\n", iter);
    printf("Line GS error = %f\n", error);

    return iter;
}

// Line SOR
int Line_SOR(double w, double psi[][MAX], double psi_old[][MAX], int imax,
             int jmax, double dx, double dy, double error_max, int iter_max,
             const char history_file[])
{
    init(psi, psi_old, imax, jmax, dx);

    FILE *fp = openHistoryFile(history_file);

    double beta = dx / dy;
    double beta2 = beta * beta;

    int N = imax - 1;

    double a[MAX], b[MAX], rhs[MAX], diag[MAX], u[MAX];
    double diag_prime[MAX], rhs_prime[MAX];

    double error = 10.0;
    double diff = 0.0;
    int iter = 0;

    printf("\nSolving Line SOR\n");
    printf("w = %f\n", w);

    while(error > error_max && iter < iter_max){

        iter++;
        error = 0.0;

        for(int j = 2; j <= jmax - 1; j++){

            for(int i = 2; i <= imax - 1; i++){
                psi_old[i][j] = psi[i][j];
            }

            // Build the Line GS tridiagonal system first.
            // SOR relaxation is applied after the Thomas solve.
            for(int i = 2; i <= imax - 1; i++){
                b[i] = -1.0;
                diag[i] = 2.0 * (1.0 + beta2);
                a[i] = -1.0;

                rhs[i] = beta2 * (psi[i][j - 1] + psi[i][j + 1]);
            }

            b[2] = 0.0;
            rhs[2] += psi[1][j];

            a[imax - 1] = 0.0;
            rhs[imax - 1] += psi[imax][j];

            thomasTriDiagonal(N, a, b, rhs, diag, u, diag_prime, rhs_prime);

            for(int i = 2; i <= imax - 1; i++){
                psi[i][j] = psi_old[i][j] + w * (u[i] - psi_old[i][j]);

                diff = fabs(psi[i][j] - psi_old[i][j]);

                if(diff > error){
                    error = diff;
                }
            }
        }

        // right boundary: dpsi/dx = 0
        for(int j = 2; j <= jmax - 1; j++){
            psi[imax][j] = psi[imax - 1][j];
        }

        // divergence check
        if(diverged(error)){
            printf("Line SOR diverged at iter = %d, error = %e\n", iter, error);
            if(fp != NULL) fclose(fp);
            return iter_max;
        }

        writeHistory(fp, iter, error);
    }

    if(fp != NULL) fclose(fp);

    printf("Line SOR iter = %d\n", iter);
    printf("Line SOR error = %f\n", error);

    return iter;
}

// write file
void writeFile(const char filename[], double psi[][MAX], int imax, int jmax,
               double dx, double dy)
{
    FILE *fp = fopen(filename, "w");

    if(fp == NULL){
        printf("Failed to open file: %s\n", filename);
        return;
    }

    for(int j = 1; j <= jmax; j++){
        for(int i = 1; i <= imax; i++){

            double xcoord = (i - 1) * dx;
            double ycoord = (j - 1) * dy;

            fprintf(fp, "%f %f %f\n", xcoord, ycoord, psi[i][j]);
        }

        fprintf(fp, "\n");
    }

    fclose(fp);

    printf("File written: %s\n", filename);
}

// main
int main()
{
    double dx = 0.2;
    double dy = 0.2;

    const int imax = 31;
    const int jmax = 21;

    double psi[MAX][MAX];
    double psi_old[MAX][MAX];

    double error_max = 0.01;
    int iter_max = 10000;

    double w_point = 1.35;
    double w_line = 1.35;

    // Run all solvers to generate solution files and convergence histories.
    Point_GS(psi, psi_old, imax, jmax, dx, dy, error_max, iter_max,
             "data/point_gs_history.dat");
    writeFile("data/point_gs.dat", psi, imax, jmax, dx, dy);

    Point_SOR(w_point, psi, psi_old, imax, jmax, dx, dy, error_max, iter_max,
              "data/point_sor_history.dat");
    writeFile("data/point_sor.dat", psi, imax, jmax, dx, dy);

    Line_GS(psi, psi_old, imax, jmax, dx, dy, error_max, iter_max,
            "data/line_gs_history.dat");
    writeFile("data/line_gs.dat", psi, imax, jmax, dx, dy);

    Line_SOR(w_line, psi, psi_old, imax, jmax, dx, dy, error_max, iter_max,
             "data/line_sor_history.dat");
    writeFile("data/line_sor.dat", psi, imax, jmax, dx, dy);

    // Iteration vs relaxation-factor study.
    sweep_w_LineSOR(imax, jmax, dx, dy, error_max, iter_max);

    // Grid refinement study for LineSOR
    grid_refinement_study(error_max, iter_max);


    return 0;
}
