#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include "mkl.h"

int eig_power(SpMat *A, double *ev, double *out, \
    int itr, double tol, int verbose){

    /* local variables */


    /* initialize */

    /* initial value of X */    

    for (i = 0; i < itr; i++){
        /* call SPMV */

        /* compute l-infty */

        /* scale the vector */

        /* check residual */

        /* swap the pointer */

        /* ??? */

        /* verbose mode */
        if (verbose)
            printf("itr: %d\tl1: %e\tres: %e\n", i + 1, );
        if (res < tol) break;
    }

    /* output */

    free(work);
    if (verbose && res > tol)
        printf("Warning. Residual too large at %e.\n", res);

    return res > tol;
}
