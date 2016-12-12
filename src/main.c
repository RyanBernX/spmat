#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"

int eig_power(SpMat*, double*, double*, int, double, int);

int main(int argc, char **argv){
    SpMat A;
    double ev, *x;

    if (argc < 2){
        fprintf(stderr, "Usage: ./main.x <filename>\n");
        return 1;
    }

    if (read_csr(argv[1], &A)) return 1;

    x = (double*)malloc(A.n * sizeof(double));

    eig_power(&A, &ev, x, 10000, 1E-6, 1);

    free(x);
    SpMat_Destroy(&A);
    return 0;
}

