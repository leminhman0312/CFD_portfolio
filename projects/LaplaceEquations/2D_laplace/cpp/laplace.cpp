#include <cstdio>
#include <iostream>
#include <cmath>

void print_array(double psi[][22], int imax, int jmax){
    for(int j = jmax; j >= 1; j--){
        for(int i = 1; i <= imax; i++){
            printf("%6.1f ", psi[i][j]);
        }
        printf("\n");
    }
}

void Point_GS(double psi[][22], double psi_old[][22],
              int imax, int jmax,
              double error_max, int iter_max){

    double error = 10.0;
    int iter = 0;

    while (error > error_max && iter < iter_max){

        iter++;
        error = 0.0;

        for (int i = 2; i <= imax-1; i++) {
            for (int j = 2; j <= jmax-1; j++) {

                psi_old[i][j] = psi[i][j];

                psi[i][j] = 0.25*(psi[i+1][j]
                                + psi[i-1][j]
                                + psi[i][j+1]
                                + psi[i][j-1]);

                error = error + fabs(psi[i][j] - psi_old[i][j]);
            }
        }

        // right exit boundary: dpsi/dx = 0
        for (int j = 2; j <= jmax-1; j++) {
            psi[imax][j] = psi[imax-1][j];
        }
    }

    printf("iter = %d\n", iter);
    printf("error = %f\n", error);
}

void writeFile(double psi[][22],
               int imax, int jmax,
               double dx, double dy){

    FILE *fp = fopen("data/PointGS.dat", "w");

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

    printf("File written: data/PointGS.dat\n");
}


void init(double psi[][22], double psi_old[][22],
          int imax, int jmax,
          double dx){

    // initialize grid
    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax; j++){

            psi[i][j] = 0.0;
            psi_old[i][j] = 0.0;

        }
    }

    // inlet indices
    int left_inlet  = int(1.0/dx) + 1;
    int right_inlet = int(1.2/dx) + 1;

    // left boundary
    for (int j = 1; j <= jmax; j++) {
        psi[1][j] = 0.0;
    }

    // top boundary
    for (int i = 1; i <= imax; i++) {
        psi[i][jmax] = 0.0;
    }

    // bottom left section
    for (int i = 1; i <= left_inlet; i++) {
        psi[i][1] = 100.0;
    }

    // bottom right section
    for (int i = right_inlet; i <= imax; i++) {
        psi[i][1] = 100.0;
    }
}


int main(){

    // grid spacing
    double dx = 0.2;
    double dy = 0.2;

    // number of grid points
    const int imax = 31;
    const int jmax = 21;

    // 2D array
    double psi[imax+1][jmax+1]; // ignore index 0
    double psi_old[imax+1][jmax+1]; // ignore index 0

    // error terms
    double error_max = 0.01;
    int iter_max = 10000;

    // initialize and BC
    init(psi, psi_old, imax, jmax, dx);

    // run PointGS
    Point_GS(psi, psi_old, imax, jmax, error_max, iter_max);
    writeFile(psi, imax, jmax, dx, dy);


    return 0;
}
