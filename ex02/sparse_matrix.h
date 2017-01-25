#ifndef __SPARSE_MATRIX_H_
#define __SPARSE_MATRIX_H_

typedef struct{
  int m, n;

  int *pntr;
  int *indx;

  double *val;
} sparse_matrix;


void initialize_matrix(sparse_matrix *, int, int, int);
int matrix_add_entry(sparse_matrix *, int, int, double);
void matrix_compress(sparse_matrix *);
void destroy_matrix(sparse_matrix *);

void gs_iterate(sparse_matrix*, double*, double*, double);
double lp_residual(sparse_matrix*, double*, double*);

int mat_mult_max_nnz(sparse_matrix *, sparse_matrix *);
void mat_mat_mult(sparse_matrix *, sparse_matrix *, sparse_matrix *);
#endif

