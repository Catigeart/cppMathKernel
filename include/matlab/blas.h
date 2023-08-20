//
// Created by Catigeart on 2023/6/21.
//

#ifndef CUBEGPA6SAB_BLAS_H
#define CUBEGPA6SAB_BLAS_H

#include "matrix.h"

void dscal(int n, double a, double *x, int incx);

void dscal2v(int n, double a, const double *x, int incx, double* y, int incy);

void dgecoo(int m, int n, int lda, const double* a, double* val, int* rowind, int* colind, int* nnz);

void dcooge(const double* vala, const int* rowinda, const int* colinda, int nnza, int m, int n, int ldb, double* b);

void cgecoo(int m, int n, int lda, const double complex* a, double complex* val, int* rowind, int* colind, int* nnz);

void ccooge(const double complex* vala, const int* rowinda, const int* colinda, int nnza, int m, int n, int ldb, double complex* b);

void dcoocoomm2ge(SparseTrans transa, const double* vala, const int* rowinda, const int* colinda, int nnza,
                  SparseTrans transb, const double* valb, const int* rowindb, const int* colindb, int nnzb,
                  int m, int n, int ldc, double* c);

void ccoocoomm2ge(SparseTrans transa, const double complex* vala, const int* rowinda, const int* colinda, int nnza,
                  SparseTrans transb, const double complex* valb, const int* rowindb, const int* colindb, int nnzb,
                  int m, int n, int ldc, double complex* c);

void cdcoocoomm2ge(SparseTrans transa, const double complex* vala, const int* rowinda, const int* colinda, int nnza,
                   SparseTrans transb, const double* valb, const int* rowindb, const int* colindb, int nnzb,
                   int m, int n, int ldc, double complex* c);

void dgecoomm2ge(int m, int n, int k, DenseTrans transa, int lda, const double* a,
                 SparseTrans transb, const double* valb, const int* rowindb, const int* colindb, int nnzb,
                 int ldc, double* c);

void cgecoomm2ge(int m, int n, DenseTrans transa, int lda, const double complex* a,
                 SparseTrans transb, const double complex* valb, const int* rowindb, const int* colindb, int nnzb,
                 int ldc, double complex* c);

void dgeadd2coo(int m, int n, double alpha, DenseTrans transa, int lda, const double* a,
                 double beta, DenseTrans transb, int ldb, const double* b,
                 double* valc, int* rowindc, int* colindc, int* nnzc);

void cgeadd2coo(int m, int n, double complex alpha, DenseTrans transa, int lda, const double complex* a,
                double complex beta, DenseTrans transb, int ldb, const double complex* b,
                double complex* valc, int* rowindc, int* colindc, int* nnzc);

void dgecopy(int m, int n, int lda, double* a, int ldb, double* b);

void dgeadd(int m, int n, double alpha, int lda, double* a, double beta, int ldb, double* b);

// double-type matrix adding
void dmadd(int m, int n, const double* a, const double* b, double* c);

void dgemm(int m, int n, int k, double alpha, const double* a, const double* b, double beta, double* c);


#endif //CUBEGPA6SAB_BLAS_H
