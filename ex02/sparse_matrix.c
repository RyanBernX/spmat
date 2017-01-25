#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sparse_matrix.h"

void initialize_matrix(sparse_matrix * M, int m, int n, int n_max_row_nnz){
  int i, j, *p_rs;
  M->m = m;
  M->n = n;

  M->pntr = (int*)malloc((m + 1) * sizeof(int));
  p_rs = M->pntr;
  for (i = 0; i <= m; ++i){
    (*p_rs) = i * n_max_row_nnz;
    p_rs ++;
  }
  M->indx = (int*)malloc(m * n_max_row_nnz * sizeof(int));
  for (i = 0; i < m; ++i){
    j = M->pntr[i];
    M->indx[j] = i;
    for (; j < M->pntr[i+1]; ++j){
      M->indx[j] = -1;
    }
  }
  M->val = (double*)calloc(m * n_max_row_nnz, sizeof(double));
}

int matrix_add_entry(sparse_matrix *M, int ri, int ci, double ev){
  int i, *p_ci;
  for (i = M->pntr[ri]; i < M->pntr[ri + 1]; ++i){
    p_ci = M->indx + i;
    if ( *p_ci == ci ){
      M->val[i] += ev;
      return 0;
    } else if ( *p_ci == -1){
      (*p_ci) = ci;
      M->val[i] = ev;
      return 0;
    }
  }
  fprintf(stderr, "Warning. Row %d is full.\n", ri);
  return 1;
}

void matrix_compress(sparse_matrix *M){
  int i, j, k;
  int *rs, *ci;
  double *ev;

  rs = (int *)malloc((M->m + 1) * sizeof(int));
  rs[0] = 0;
  for (i = 0; i < M->m; ++i){
    for (j = M->pntr[i]; j < M->pntr[i + 1]; ++j){
      if (M->indx[j] == -1) break;
    }
    rs[i + 1] = rs[i] + j - M->pntr[i];
  }

  ci = (int *)malloc(rs[M->m] * sizeof(int));
  ev = (double *)malloc(rs[M->m] * sizeof(double));

  k = 0;
  for (i = 0; i < M->m; ++i){
    for (j = M->pntr[i]; j < M->pntr[i + 1]; ++j){
      if (M->indx[j] != -1){
        ci[k] = M->indx[j];
        ev[k] = M->val[j];
        ++ k;
      }
    }
  }

  printf("nnz: %d\n", rs[M->m]);

  free(M->pntr); M->pntr = rs;
  free(M->indx); M->indx = ci;
  free(M->val);  M->val = ev;
}

double lp_residual(sparse_matrix *M, double *x, double *b){
  int i, j;
  double res = 0, ri;

  for (i = 0; i < M->m; i++){
    ri = b[i];
    for (j = M->pntr[i]; j < M->pntr[i + 1]; j++){
      ri -= M->val[j] * x[M->indx[j]];
    }
    res += ri * ri;
  }
  return sqrt(res);
}

void gs_iterate(sparse_matrix *M, double *x, double *b, double tol){
  int i, j, k;
  double res, init_res, *r, ri;

  init_res = lp_residual(M, x, b);
  init_res = sqrt(init_res);
  do {
    for (k = 0; k < 10; ++k){
      for (i = 0; i < M->m; ++i){
        ri = b[i];
        for (j = M->pntr[i] + 1; j < M->pntr[i + 1]; ++j){
          ri -= M->val[j] * x[M->indx[j]];
        }
        x[i] = ri / M->val[M->pntr[i]];
      }
    }
    res = lp_residual(M, x, b) / init_res;
    fprintf(stderr, "residual = %e\n", res);
  } while (res > tol);
}

int mat_mult_max_nnz(sparse_matrix *A, sparse_matrix *B){
  int nnz = 0, row_nnz = 0;
  int i, j, j1, j2, k;
  int *buf, *idx;

  buf = (int *)calloc(B->n, sizeof(int));
  idx = (int *)malloc(B->n * sizeof(int));


  for (i = 0; i < A->m; ++i){
    row_nnz = 0;
    for (j1 = A->pntr[i]; j1 < A->pntr[i + 1]; ++j1){
      k = A->indx[j1];
      for (j2 = B->pntr[k]; j2 < B->pntr[k + 1]; ++ j2){
        j = B->indx[j2];
        if(buf[j] == 0){
          idx[row_nnz] = j;
          row_nnz += 1;
          buf[j] = 1;
        }
      }
    }
    nnz = row_nnz > nnz ? row_nnz : nnz;
    for (j1 = 0; j1 < row_nnz; ++j1){
      buf[idx[j1]] = 0;
    }
  }


  free(buf);
  free(idx);

  return nnz;
}

void mat_mat_mult(sparse_matrix *A, sparse_matrix *B, sparse_matrix *C){
  int row_nnz;
  int i, j, j1, j2, k;
  int *idx, *ibuf;
  double *buf;

  idx = (int *)malloc(C->n * sizeof(int));
  ibuf = (int *)calloc(C->n, sizeof(int));
  buf = (double *)calloc(C->n, sizeof(double));


  for (i = 0; i < A->m; ++i){
    row_nnz = 0;
    for (j1 = A->pntr[i]; j1 < A->pntr[i + 1]; ++j1){
      k = A->indx[j1];
      // A->val[j1] -> a(i, k)
      for (j2 = B->pntr[k]; j2 < B->pntr[k + 1]; ++j2){
        j = B->indx[j2];
        if (ibuf[j] == 0){
          idx[row_nnz] = j;
          row_nnz += 1;
          ibuf[j] = 1;
          
        }
        buf[j] += A->val[j1] * B->val[j2];
        // B->val[j2] -> b(k, j)
      }
    }

    for (j1 = 0; j1 < row_nnz; ++j1){
      matrix_add_entry(C, i, idx[j1], buf[idx[j1]]);
      ibuf[idx[j1]] = 0;
      buf[idx[j1]] = 0;
    }
  }

  free(idx);
  free(ibuf);
  free(buf);

  matrix_compress(C);
}
void destroy_matrix(sparse_matrix *M){
  free(M->pntr);
  free(M->indx);
  free(M->val);
}
