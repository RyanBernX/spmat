#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "mkl.h"

int eig_power(SpMat *A, double *ev, double *out, \
    int itr, double tol, int verbose){

    double *work, *x, *y, *t;
    double l1, res;
    int i, n;
    int iseed[] = {233, 233, 233, 233};

    char trans[] = {'N'};
    double alpha = 1.0, beta = 0.0;

    /* initialize */
    n = A->n;
    work = (double*)malloc(2 * n * sizeof(double));
    x = work; y = work + n;

    LAPACKE_dlarnv(3, iseed, n, x);

    for (i = 0; i < itr; i++){
        mkl_dcsrmv(trans, &n, &n, &alpha, A->matdesc, A->val, A->indx, \
            A->pntr, A->pntr + 1, x, &beta, y);
        l1 = y[cblas_idamax(n, y, 1)];
        cblas_dscal(n, 1 / l1, y, 1);
        /* check residual */
        cblas_daxpy(n, -1, y, 1, x, 1);
        res = cblas_dnrm2(n, x, 1);
        /* swap the pointer */
        t = x; x = y; y = t;
        if (verbose)
            printf("itr: %d\tl1: %e\tres: %e\n", i + 1, l1, res);
        if (res < tol) break;
    }
    cblas_dcopy(n, x, 1, out, 1);
    *ev = l1;

    free(work);
    if (verbose && res > tol)
        printf("Warning. Residual too large at %e.\n", res);

    return res > tol;
}
