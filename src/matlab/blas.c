//
// Created by Catigeart on 2023/6/21.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <complex.h>
#include "matlab/matrix.h"
#include "matlab/util.h"



void dscal(const int n, const double a, double *x, const int incx) {
    for (int i = 0; i < n; i++) {
        x[i * incx] *= a;
    }
}

void dscal2v(const int n, const double a, const double *x, const int incx, double* y, const int incy) {
    for (int i = 0; i < n; i++) {
        y[i * incy] = x[i * incx] * a;
    }
}

/**
 * 普通矩阵转csc格式稀疏矩阵
 * @param m
 * @param n
 * @param lda
 * @param a
 * @param val size=m*n
 * @param ia size=m*n
 * @param ja size=n
 * @param nz
 */
void dgecsc(int m, int n, int lda, const double* a, double* val, int* ia, int* ja, int nz) {
    nz = 0;
    for (int j = 0; j < n; j++) {
        ja[j] = nz;
        for (int i = 0; i < m; i++) {
            val[nz] = a[i + j * lda];
            ia[nz] = i;
        }
        nz++;
    }
}

void dgecoo(int m, int n, int lda, const double* a, double* val, int* rowind, int* colind, int* nnz) {
    *nnz = 0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            if (fabs(a[i + j * lda]) < 1e-6) {
                continue;
            }
            val[*nnz] = a[i + j * lda];
            rowind[*nnz] = i;
            colind[*nnz] = j;
            (*nnz)++;
        }
    }
}

void dcooge(const double* vala, const int* rowinda, const int* colinda, int nnza, int m, int n, int ldb, double* b) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            b[i + j * ldb] = 0;
        }
    }

    for (int i = 0; i < nnza; i++) {
        b[rowinda[i] + colinda[i] * ldb] = vala[i];
    }
}

void cgecoo(int m, int n, int lda, const double complex* a, double complex* val, int* rowind, int* colind, int* nnz) {
    *nnz = 0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            if (cabs(a[i + j * lda]) < 1e-6) {
                continue;
            }
            val[*nnz] = a[i + j * lda];
            rowind[*nnz] = i;
            colind[*nnz] = j;
            (*nnz)++;
        }
    }
}

void ccooge(const double complex* vala, const int* rowinda, const int* colinda, int nnza, int m, int n, int ldb, double complex* b) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            b[i + j * ldb] = 0;
        }
    }

    for (int i = 0; i < nnza; i++) {
        b[rowinda[i] + colinda[i] * ldb] = vala[i];
    }
}

void dcoocoomm2ge(SparseTrans transa, const double* vala, const int* rowinda, const int* colinda, int nnza,
                  SparseTrans transb, const double* valb, const int* rowindb, const int* colindb, int nnzb,
                  int m, int n, int ldc, double* c) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            c[i + j * ldc] = 0;
        }
    }

    if (transa == SPARSE_NON_TRANS && transb == SPARSE_NON_TRANS) {
        for (int ap = 0; ap < nnza; ap++) {
            for (int bp = 0; bp < nnzb; bp++) {
                if (colinda[ap] == rowindb[bp]) {
                    c[rowinda[ap] + colindb[bp] * ldc] += vala[ap] * valb[bp];
                }
            }
        }
    }
    else if (transa == SPARSE_NON_TRANS && transb == SPARSE_TRANS) {
        for (int ap = 0; ap < nnza; ap++) {
            for (int bp = 0; bp < nnzb; bp++) {
                if (colinda[ap] == colindb[bp]) {
                    c[rowinda[ap] + rowindb[bp] * ldc] += vala[ap] * valb[bp];
                }
            }
        }
    }
    else {
        fprintf(stdout, "Unsupported operation of dcoocoomm2ge!");
    }
}

void ccoocoomm2ge(SparseTrans transa, const double complex* vala, const int* rowinda, const int* colinda, int nnza,
                  SparseTrans transb, const double complex* valb, const int* rowindb, const int* colindb, int nnzb,
                  int m, int n, int ldc, double complex* c) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            c[i + j * ldc] = 0;
        }
    }

    if (transa == SPARSE_NON_TRANS && transb == SPARSE_NON_TRANS) {
        for (int ap = 0; ap < nnza; ap++) {
            for (int bp = 0; bp < nnzb; bp++) {
                if (colinda[ap] == rowindb[bp]) {
                    c[rowinda[ap] + colindb[bp] * ldc] += vala[ap] * valb[bp];
                }
            }
        }
    }
    else if (transa == SPARSE_NON_TRANS && transb == SPARSE_TRANS) {
        for (int ap = 0; ap < nnza; ap++) {
            for (int bp = 0; bp < nnzb; bp++) {
                if (colinda[ap] == colindb[bp]) {
                    c[rowinda[ap] + rowindb[bp] * ldc] += vala[ap] * valb[bp];                }
            }
        }
    }
    else if (transa == SPARSE_NON_TRANS && transb == SPARSE_CONJ_TRANS) {
        for (int ap = 0; ap < nnza; ap++) {
            for (int bp = 0; bp < nnzb; bp++) {
                if (colinda[ap] == colindb[bp]) {
// #define DEBUG
#ifdef DEBUG
                    int idx = rowinda[ap] + rowindb[bp] * ldc;
                    printf("%+lf%+lfi", creal(c[idx]), cimag(c[idx]));
                    c[idx] += vala[ap] * (creal(valb[bp]) - cimag(valb[bp]));
                    printf("+(%+lf%+lfi)*(%+lf%+lfi)=%+lf%+lfi=c[%d]\n",
                           creal(vala[ap]), cimag(vala[ap]),
                           creal(valb[bp]), -cimag(valb[bp]),
                           creal(c[idx]), cimag(c[idx]),
                           idx);
                    printf("%+lf%+lfi\t", creal(valb[bp]), -cimag(valb[bp]));
                    double complex new_val = creal(valb[bp]) - cimag(valb[bp]);
                    printf("(%+lf%+lfi)*(%+lf%+lfi)=%+lf%+lfi\n",
                           creal(vala[ap]), cimag(vala[ap]),
                           creal(new_val), cimag(new_val),
                           creal(vala[ap] * new_val), cimag(vala[ap] * new_val)
                    );
#else
                    // cimag()返回的是double，要求共轭复数的话需要再手动乘一个I
                    c[rowinda[ap] + rowindb[bp] * ldc] += vala[ap] * (creal(valb[bp]) - cimag(valb[bp]) * I);
#endif
// #undef DEBUG
                }
            }
        }
    }
    else {
        fprintf(stdout, "Unsupported operation of ccoocoomm2ge!");
    }
}

void cdcoocoomm2ge(SparseTrans transa, const double complex* vala, const int* rowinda, const int* colinda, int nnza,
                  SparseTrans transb, const double* valb, const int* rowindb, const int* colindb, int nnzb,
                  int m, int n, int ldc, double complex* c) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            c[i + j * ldc] = 0;
        }
    }

    if (transa == SPARSE_NON_TRANS && transb == SPARSE_NON_TRANS) {
        for (int ap = 0; ap < nnza; ap++) {
            for (int bp = 0; bp < nnzb; bp++) {
                if (colinda[ap] == rowindb[bp]) {
                    c[rowinda[ap] + colindb[bp] * ldc] += vala[ap] * valb[bp];
                }
            }
        }
    }
    else {
        fprintf(stderr, "Unsupported operation of dccoocoomm2ge!");
    }
}

/// C := AB
void dgecoomm2ge(int m, int n, int k, DenseTrans transa, int lda, const double* a,
                 SparseTrans transb, const double* valb, const int* rowindb, const int* colindb, int nnzb,
                 int ldc, double* c) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            c[i + j * ldc] = 0;
        }
    }

    if (transa == DENSE_NON_TRANS && transb == SPARSE_NON_TRANS) {
        for (int bp = 0; bp < nnzb; bp++) {
            for (int i = 0; i < m; i++) {
                c[i + colindb[bp] * ldc] += a[i + rowindb[bp] * lda] * valb[bp];
            }
        }
    }
    else if (transa == DENSE_NON_TRANS && transb == SPARSE_TRANS) {
        for (int bp = 0; bp < nnzb; bp++) {
            for (int i = 0; i < m; i++) {
                c[i + rowindb[bp] * ldc] += a[i + colindb[bp] * lda] * valb[bp];
            }
        }
    }
    else {
        fprintf(stderr, "Unsupported operation of dgecoomm2ge!");
    }
}

void cgecoomm2ge(int m, int n, DenseTrans transa, int lda, const double complex* a,
                 SparseTrans transb, const double complex* valb, const int* rowindb, const int* colindb, int nnzb,
                 int ldc, double complex* c) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            c[i + j * ldc] = 0;
        }
    }

    if (transa == DENSE_NON_TRANS && transb == SPARSE_NON_TRANS) {
        for (int bp = 0; bp < nnzb; bp++) {
            for (int i = 0; i < m; i++) {
                c[i + colindb[bp] * ldc] += a[i + rowindb[bp] * lda] * valb[bp];
            }
        }
    }
    else if (transa == DENSE_NON_TRANS && transb == SPARSE_TRANS) {
        for (int bp = 0; bp < nnzb; bp++) {
            for (int i = 0; i < m; i++) {
                c[i + rowindb[bp] * ldc] += a[i + colindb[bp] * lda] * valb[bp];
            }
        }
    }
    else if (transa == DENSE_NON_TRANS && transb == SPARSE_CONJ_TRANS) {
        for (int bp = 0; bp < nnzb; bp++) {
            for (int i = 0; i < m; i++) {
                c[i + rowindb[bp] * ldc] += a[i + colindb[bp] * lda] * (creal(valb[bp]) - cimag(valb[bp]) * I);
            }
        }
    }
    else {
        fprintf(stderr, "Unsupported operation of dgecoomm2ge!");
    }
}

/// C = AB
void dcoogemm2ge(SparseTrans transa, const double* vala, const int* rowinda, const int* colinda, int nnza,
                 int m, int n, int k, DenseTrans transb, int ldb, const double* b, int ldc, double* c) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            c[i + j * ldc] = 0;
        }
    }

    fprintf(stderr, "Unsupported operation of dcoogemm2ge!");
}

/// C := alpha * A + beta * b
void dgeadd2coo(int m, int n, double alpha, DenseTrans transa, int lda, const double* a,
                 double beta, DenseTrans transb, int ldb, const double* b,
                 double* valc, int* rowindc, int* colindc, int* nnzc) {
    *nnzc = 0;
    if (transa == DENSE_NON_TRANS && transb == DENSE_NON_TRANS) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                double absum = alpha * a[i + j * lda] + beta * b[i + j * ldb];
                if (fabs(absum) > 1e-6) {
                    valc[*nnzc] = absum;
                    rowindc[*nnzc] = i;
                    colindc[*nnzc] = j;
                    (*nnzc)++;
                }
            }
        }
    }
}

/// C := alpha * A + beta * b
void cgeadd2coo(int m, int n, double complex alpha, DenseTrans transa, int lda, const double complex* a,
                double complex beta, DenseTrans transb, int ldb, const double complex* b,
                double complex* valc, int* rowindc, int* colindc, int* nnzc) {
    *nnzc = 0;
    if (transa == DENSE_NON_TRANS && transb == DENSE_NON_TRANS) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                double complex absum = alpha * a[i + j * lda] + beta * b[i + j * ldb];
                if (cabs(absum) > 1e-6) {
                    valc[*nnzc] = absum;
                    rowindc[*nnzc] = i;
                    colindc[*nnzc] = j;
                    (*nnzc)++;
                }
            }
        }
    }
}

// B := A
void dgecopy(int m, int n, int lda, double* a, int ldb, double* b) {
    for (int j = 0; j < n; j++) {
        memcpy(&b[0 + j * ldb], &a[0 + j * lda], m * sizeof(double));
    }
}

// B := alpha * A + beta * B
void dgeadd(int m, int n, double alpha, int lda, double* a, double beta, int ldb, double* b) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            b[i + j * ldb] = alpha * a[i + j * lda] + beta * b[i + j * ldb];
        }
    }
}

// double-type matrix adding
void dmadd(int m, int n, const double* a, const double* b, double* c) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            c[i + j * m] = a[i + j * m] + b[i + j * m];
        }
    }
}

// double-type C := alpha * A * B + beta * C
void dgemm(int m, int n, int k, double alpha, const double* a, const double* b, double beta, double* c) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            c[i + j * m] = beta * c[i + j * m];
        }
    }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            for (int kk = 0; kk < k; kk++) {
                c[i + j * m] += alpha * a[i + kk * m] * b[kk + j * k];
            }
        }
    }
}

void daxpy(int n, double a, const double* x, int incx, double* y, int incy) {
    for (int i = 0; i < n; i++) {
        y[i * incy] += a * x[i * incx];
    }
}
/*
void dspadd(int m, int n, SparseTrans transA, SparseTrans transB, double alpha, dSpMatrix* spA, double beta, dSpMatrix* spB, dSpMatrix* spC) {
    if (transA == SPARSE_NON_TRANS && transB == SPARSE_NON_TRANS) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {

            }
        }
    }
}

void dsp2m(SparseTrans transA, SparseTrans transB, double alpha, dSpMatrix* spA, dSpMatrix* spB, double beta, dSpMatrix* spC) {
    for (int i = 0; i < spC->len; i++) {
        spC->data[i].val *= beta;
    }
    realloc(spC->data, spA->len * spB->len * sizeof(dCOO)); // 预留空间

    if (transA == SPARSE_NON_TRANS && transB == SPARSE_NON_TRANS) {
        for (int i = 0; i < spA->len; i++) {
            for (int j = 0; j < spB->len; j++) {
                if (spA->data[i].col == spB->data[i].row) {
                    bool has_val = false;
                    for (int k = 0; k < spC->len; k++) {
                        if (spC->data[k].row == spA->data[i].row && spC->data[k].col == spB->data[i].col) {
                            spC->data[k].val += alpha * spA->data[i].val * spB->data[i].val;
                            has_val = true;
                            break;
                        }
                    }
                    if (!has_val) {
                        spC->data[spC->len].row = spA->data[i].row;
                        spC->data[spC->len].col = spB->data[i].col;
                        spC->data[spC->len].val = alpha * spA->data[i].val * spB->data[i].val;
                        spC->len++;
                    }
                }
            }
        }
    }
    else if (transA == SPARSE_NON_TRANS && transB == SPARSE_TRANS) {
        for (int i = 0; i < spA->len; i++) {
            for (int j = 0; j < spB->len; j++) {
                if (spA->data[i].col == spB->data[i].col) {
                    bool has_val = false;
                    for (int k = 0; k < spC->len; k++) {
                        if (spC->data[k].row == spA->data[i].row && spC->data[k].col == spB->data[i].row) {
                            spC->data[k].val += alpha * spA->data[i].val * spB->data[i].val;
                            has_val = true;
                            break;
                        }
                    }
                    if (!has_val) {
                        spC->data[spC->len].row = spA->data[i].row;
                        spC->data[spC->len].col = spB->data[i].row;
                        spC->data[spC->len].val = alpha * spA->data[i].val * spB->data[i].val;
                        spC->len++;
                    }
                }
            }
        }
    }
    else {
        // TODO
    }

}

// 稀疏矩阵相乘，结果存储为稠密矩阵
void dsp2md(SparseTrans transA, SparseTrans transB, int m, int n, double alpha, dSpMatrix* spA, dSpMatrix* spB,
            double beta, int ldc, double* c) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            c[i + j * ldc] *= beta;
        }
    }
    if (transA == SPARSE_NON_TRANS && transB == SPARSE_TRANS) {
        for (int i = 0; i < spA->len; i++) {
            for (int j = 0; j < spB->len; j++) {
                if (spA->data[i].col == spB->data[i].col) {
                    c[spA->data[i].row + spB->data[i].row * ldc] += alpha * spA->data[i].val * spB->data[i].val;
                }
            }
        }
    }
    else {
        // TODO
    }
}*/