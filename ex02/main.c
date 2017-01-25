#include <stdio.h>
#include <stdlib.h>
#include "sparse_matrix.h"

#define N 10

int main(int argc, char **argv){
  int i, j, k, ri, nnz;
  sparse_matrix A, A2;
  double *x, *b;

  initialize_matrix(&A, N*N*N, N*N*N, 7);

  for (i = 0; i < N; ++i){
    for (j = 0; j < N; ++j){
      for (k = 0; k < N; ++k){
        ri = (i * N + j) * N + k;
        matrix_add_entry(&A, ri, ri, 6.0);
        if (i > 0) matrix_add_entry(&A, ri, ri - N * N, -1.0);
        if (i < N - 1) matrix_add_entry(&A, ri, ri + N * N, -1.0);
        if (j > 0) matrix_add_entry(&A, ri, ri - N, -1.0);
        if (j < N - 1) matrix_add_entry(&A, ri, ri + N, -1.0);
        if (k > 0) matrix_add_entry(&A, ri, ri - 1, -1.0);
        if (k < N - 1) matrix_add_entry(&A, ri, ri + 1, -1.0);
      }
    }
  }

  matrix_compress(&A);

  nnz = mat_mult_max_nnz(&A, &A);
  printf("%d\n", nnz);
  initialize_matrix(&A2, A.m, A.n, nnz);

  mat_mat_mult(&A, &A, &A2);

  x = (double *)calloc(N*N*N, sizeof(double));
  b = (double *)calloc(N*N*N, sizeof(double));

  for (i = 0; i < N*N*N; i++) b[i] = 1.0;

  gs_iterate(&A2, x, b, 1E-6);

  for (i = 0; i < N; ++i){
    if (i != 5) continue;
    for (j = 0; j < N; ++j){
      for (k = 0; k < N; ++k){
        printf("%f ", x[(i*N+j)*N+k]);
      }
      printf("\n");
    }
    //printf("\n");
  }

  free(x); free(b);
  destroy_matrix(&A);
  destroy_matrix(&A2);

  return 0;
}
