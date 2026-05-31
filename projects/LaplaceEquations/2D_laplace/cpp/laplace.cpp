#include <cstdio>
#include <cmath>

// prototypes
void init(double psi[][22], double psi_old[][22],
          int imax, int jmax, double dx);

void print_array(double psi[][22], int imax, int jmax);

void writeFile(double psi[][22],
               int imax, int jmax,
               double dx, double dy);

void Point_GS(double psi[][22], double psi_old[][22],
              int imax, int jmax,
              double dx, double dy,
              double error_max, int iter_max);

void Point_SOR(double w,
               double psi[][22], double psi_old[][22],
               int imax, int jmax,
               double dx, double dy,
               double error_max, int iter_max);

void Line_GS(double psi[][22], double psi_old[][22],
             int imax, int jmax,
             double dx, double dy,
             double error_max, int iter_max);

static void thomasTriDiagonal(int N,
                              const double a[],
                              const double b[],
                              const double c[],
                              const double d[],
                              double u[],
                              double dprime[],
                              double cprime[]);


// init
void init(double psi[][22], double psi_old[][22],
          int imax, int jmax,
          double dx){

    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax; j++){
            psi[i][j] = 0.0;
            psi_old[i][j] = 0.0;
        }
    }

    int left_inlet  = int(1.0/dx) + 1;
    int right_inlet = int(1.2/dx) + 1;

    for(int j = 1; j <= jmax; j++){
        psi[1][j] = 0.0;
    }

    for(int i = 1; i <= imax; i++){
        psi[i][jmax] = 0.0;
    }

    for(int i = 1; i <= left_inlet; i++){
        psi[i][1] = 0.0;
    }

    for(int i = right_inlet; i <= imax; i++){
        psi[i][1] = 100.0;
    }
}


// print array
void print_array(double psi[][22], int imax, int jmax){

    for(int j = jmax; j >= 1; j--){
        for(int i = 1; i <= imax; i++){
            printf("%6.1f ", psi[i][j]);
        }
        printf("\n");
    }
}


// Point GS
void Point_GS(double psi[][22], double psi_old[][22],
              int imax, int jmax,
              double dx, double dy,
              double error_max, int iter_max){

    init(psi, psi_old, imax, jmax, dx);

    double error = 10.0;
    int iter = 0;

    double beta = dx/dy;
    double beta2 = beta*beta;
    double coeff = 1.0/(2.0*(1.0 + beta2));

    printf("\nSolving Point GS\n");

    while(error > error_max && iter < iter_max){

        iter++;
        error = 0.0;

        for(int i = 2; i <= imax-1; i++){
            for(int j = 2; j <= jmax-1; j++){

                psi_old[i][j] = psi[i][j];

                psi[i][j] = coeff*(psi[i+1][j]
                                + psi[i-1][j]
                                + beta2*(psi[i][j+1]
                                + psi[i][j-1]));

                error += fabs(psi[i][j] - psi_old[i][j]);
            }
        }

        for(int j = 2; j <= jmax-1; j++){
            psi[imax][j] = psi[imax-1][j];
        }
    }

    printf("Point GS iter = %d\n", iter);
    printf("Point GS error = %f\n", error);
}


// Point SOR
void Point_SOR(double w,
               double psi[][22], double psi_old[][22],
               int imax, int jmax,
               double dx, double dy,
               double error_max, int iter_max){

    init(psi, psi_old, imax, jmax, dx);

    double error = 10.0;
    int iter = 0;

    double beta = dx/dy;
    double beta2 = beta*beta;
    double coeff = 1.0/(2.0*(1.0 + beta2));

    printf("\nSolving Point SOR\n");
    printf("w = %f\n", w);

    while(error > error_max && iter < iter_max){

        iter++;
        error = 0.0;

        for(int i = 2; i <= imax-1; i++){
            for(int j = 2; j <= jmax-1; j++){

                psi_old[i][j] = psi[i][j];

                double residual =
                      psi[i+1][j]
                    + psi[i-1][j]
                    + beta2*(psi[i][j+1] + psi[i][j-1])
                    - 2.0*(1.0 + beta2)*psi[i][j];

                psi[i][j] += w*coeff*residual;

                error += fabs(psi[i][j] - psi_old[i][j]);
            }
        }

        for(int j = 2; j <= jmax-1; j++){
            psi[imax][j] = psi[imax-1][j];
        }
    }

    printf("Point SOR iter = %d\n", iter);
    printf("Point SOR error = %f\n", error);
}


// Thomas tridiagonal solver
static void thomasTriDiagonal(int N,
                              const double a[],
                              const double b[],
                              const double c[],
                              const double d[],
                              double u[],
                              double dprime[],
                              double cprime[]){

    dprime[2] = d[2];
    cprime[2] = c[2];

    for(int i = 3; i <= N; i++){
        dprime[i] = d[i] - (b[i]*a[i-1])/dprime[i-1];
        cprime[i] = c[i] - (cprime[i-1]*b[i])/dprime[i-1];
    }

    u[N] = cprime[N]/dprime[N];

    for(int i = N-1; i >= 2; i--){
        u[i] = (cprime[i] - a[i]*u[i+1])/dprime[i];
    }
}


// Line GS
void Line_GS(double psi[][22], double psi_old[][22],
             int imax, int jmax,
             double dx, double dy,
             double error_max, int iter_max){

    init(psi, psi_old, imax, jmax, dx);

    double beta = dx/dy;
    double beta2 = beta*beta;

    int N = imax - 1;

    double a[32], b[32], c[32], d[32], u[32];
    double dprime[32], cprime[32];

    double error = 10.0;
    int iter = 0;

    printf("\nSolving Line GS\n");

    while(error > error_max && iter < iter_max){

        iter++;
        error = 0.0;

        for(int j = 2; j <= jmax-1; j++){

            for(int i = 2; i <= imax-1; i++){
                psi_old[i][j] = psi[i][j];
            }

            for(int i = 2; i <= imax-1; i++){

                b[i] = -1.0;
                d[i] =  2.0*(1.0 + beta2);
                a[i] = -1.0;

                c[i] = beta2*(psi[i][j-1] + psi[i][j+1]);
            }

            b[2] = 0.0;
            c[2] += psi[1][j];

            a[imax-1] = 0.0;
            c[imax-1] += psi[imax][j];

            thomasTriDiagonal(N, a, b, c, d, u, dprime, cprime);

            for(int i = 2; i <= imax-1; i++){
                psi[i][j] = u[i];
                error += fabs(psi[i][j] - psi_old[i][j]);
            }
        }

        for(int j = 2; j <= jmax-1; j++){
            psi[imax][j] = psi[imax-1][j];
        }
    }

    printf("Line GS iter = %d\n", iter);
    printf("Line GS error = %f\n", error);
}


// write file
void writeFile(double psi[][22],
               int imax, int jmax,
               double dx, double dy){

    FILE *fp = fopen("data/psi.dat", "w");

    if(fp == NULL){
        printf("File failed to open\n");
        return;
    }

    for(int j = 1; j <= jmax; j++){
        for(int i = 1; i <= imax; i++){

            double xcoord = (i-1)*dx;
            double ycoord = (j-1)*dy;

            fprintf(fp, "%f %f %f\n",
                    xcoord, ycoord, psi[i][j]);
        }

        fprintf(fp, "\n");
    }

    fclose(fp);

    printf("File written: data/psi.dat\n");
}


// main
int main(){

    double dx = 0.2;
    double dy = 0.2;

    const int imax = 31;
    const int jmax = 21;

    double psi[imax+1][jmax+1];
    double psi_old[imax+1][jmax+1];

    double error_max = 0.01;
    int iter_max = 10000;

    double w = 1.00;

    // choose one solver
    Point_GS(psi, psi_old, imax, jmax, dx, dy, error_max, iter_max);

    Point_SOR(w, psi, psi_old, imax, jmax, dx, dy, error_max, iter_max);

    Line_GS(psi, psi_old, imax, jmax, dx, dy, error_max, iter_max);

    // writeFile(psi, imax, jmax, dx, dy);

    return 0;
}
