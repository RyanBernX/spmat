#include <stdio.h>
#include <string.h>
#include "spmat.h"

/** @fn int read_csr(const char *filename, SpMat *A)
 *  @brief Read the sparse matrix from file. The sparse matrix should have
 *  Rutherford Boeing format.
 *  @param filename The matrix filename.
 *  @param A Pointer to the sparse matrix object.
 *  @return (int) The function returns 0 if it ends with no errors. Otherwise,
 *  either it fails to open the file or the file is illegal.
 */

int read_csr(const char *filename, SpMat* A){
    FILE *fp = fopen(filename, "r");
    int i, nnz, row, col;
    char headerl1[82], headerl2[82], headerl3[82], headerl4[82], which[4];
    MKL_INT *pntr, *indx;
    double *val;

    if(fp == NULL){
        printf("Cannot open file %s. \n", filename);
        return -1;
    }

    fgets(headerl1, 82, fp);
    fgets(headerl2, 82, fp);
    fgets(headerl3, 82, fp);
    fgets(headerl4, 82, fp);
    sscanf(headerl3, "%3s %d %d %d", which, &row, &col, &nnz);

    if(strcmp(which, "rsa") == 0 || strcmp(which, "rua") == 0){
        strncpy(A->desc, headerl1, 72);
        strncpy(A->label, headerl1 + 72, 8);
        A->m = row; A->n = col; A->nnz = nnz;

        /* matrix descriptor */
        A->matdesc[0] = which[1] == 's' ? 'S' : 'G';
        A->matdesc[1] = 'U';
        A->matdesc[2] = 'N';
        A->matdesc[3] = 'F';

        /* allocating memory */
        A->val = (double*)malloc(nnz * sizeof(double));
        A->indx = (MKL_INT*)malloc(nnz * sizeof(MKL_INT));
        A->pntr = (MKL_INT*)malloc((col + 1) * sizeof(MKL_INT));
        val = A->val;
        pntr = A->pntr;
        indx = A->indx;

        /* input PNTR */
        for(i = 0; i <= col; i++) 
            fscanf(fp, "%d", pntr + i);
        /* input INDX */
        for(i = 0; i < nnz; i++)
            fscanf(fp, "%d", indx + i);
        /* input VAL */
        for(i = 0; i < nnz; i++)
            fscanf(fp, "%le", val+i);
    }
    else{
        printf("Format error. Rutherford Boeing CSR format matrix needed.\n");
        return -1;
    }
    fclose(fp);
    return 0;
}

void SpMat_Destroy(SpMat *A){
    free(A->val);
    free(A->pntr);
    free(A->indx);
}
