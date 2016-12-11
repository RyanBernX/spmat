#ifndef __SPMAT_H_
#define __SPMAT_H_

#include "mkl_types.h"

typedef struct{
    double *val;
    MKL_INT *indx;
    MKL_INT *pntr;
    int nnz;
    int m;
    int n;
    char label[9];
    char desc[73];
    char matdesc[4];
} SpMat;

int read_csr(const char *, SpMat *);
void SpMat_Destroy(SpMat *);

#endif

