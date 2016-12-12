#include <stdlib.h>
#include "mex.h"
#include "spmat.h"

int eig_power(SpMat*, double*, double*, int, double, int);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  double *ev, *x;
  SpMat A;
  size_t n;
  int i;
  mwIndex *Ir, *Jc;

  if (nrhs == 0){
    mexErrMsgTxt("One input is needed.");
    return;
  }

  if (!mxIsSparse(prhs[0])){
    mexErrMsgTxt("A must be sparse.");
    return;
  }

  n = mxGetM(prhs[0]);
  mexPrintf("%d\n", sizeof(mwIndex));

  if (nlhs != 2){
    mexErrMsgTxt("Usage: [ev, x] = eigpow(A).");
    return;
  }

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  ev = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
  x = mxGetPr(plhs[1]);

  A.m = A.n = n;
  A.matdesc[0] = 'G';
  A.matdesc[1] = 'L';
  A.matdesc[2] = 'N';
  A.matdesc[3] = 'C';
  A.val = mxGetPr(prhs[0]);
  Ir = mxGetIr(prhs[0]);
  Jc = mxGetJc(prhs[0]);

  A.nnz = Jc[n];

  A.indx = (int*)malloc(A.nnz * sizeof(int));
  A.pntr = (int*)malloc((n + 1) * sizeof(int));

  for (i = 0; i <= n; i++) A.pntr[i] = Jc[i];
  for (i = 0; i < A.nnz; i++) A.indx[i] = Ir[i];

  eig_power(&A, ev, x, 10000, 1E-6, 0);

  free(A.indx);
  free(A.pntr);

  return;
}
