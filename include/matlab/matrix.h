//
// Created by Catigeart on 2023/6/21.
//

#ifndef CUBEGPA6SAB_MATRIX_H
#define CUBEGPA6SAB_MATRIX_H

#include <complex.h>

#define complex _Complex

#define DCOO_MALLOC(coo, n) coo.val=(double*)malloc(n*sizeof(double));coo.rowind=(int*)malloc(n*sizeof(int));coo.colind=(int*)malloc(n*sizeof(int))
#define CCOO_MALLOC(coo, n) coo.val=(double complex*)malloc(n*sizeof(double complex));coo.rowind=(int*)malloc(n*sizeof(int));coo.colind=(int*)malloc(n*sizeof(int))

#define DCSR_MALLOC(csr, m, n) csr.val=(double*)malloc(m*n*sizeof(double));csr.ptrB=(int*)malloc(m*sizeof(int));csr.colind=(int*)malloc(m*n*sizeof(int))
#define CCSR_MALLOC(csr, m, n) csr.val=(double complex*)malloc(m*n*sizeof(double complex));csr.ptrB=(int*)malloc(m*sizeof(int));csr.colind=(int*)malloc(m*n*sizeof(int))

#define COO_FREE(coo) free(coo.val);free(coo.rowind);free(coo.colind)

typedef enum {
    SPARSE_NON_TRANS,
    SPARSE_TRANS,
    SPARSE_CONJ_TRANS
} SparseTrans;

typedef enum {
    DENSE_NON_TRANS,
    DENSE_TRANS
} DenseTrans;

typedef struct {
    int nnz;
    int* rowind;
    int* colind;
    double* val;
} dCOOMatrix;

typedef struct {
    int nnz;
    int* rowind;
    int* colind;
    double complex* val;
} cCOOMatrix;

typedef struct {
    int nnz; // 非零元素数
    int idx; // base-idx，指示从0开始计数还是从1开始计数
    int* ptrB; // 索引指向每行的Begin位置，size=行数, 指向的val[i]的下标=ptrB[i]-idx
    int* colind; // size=nnz，列坐标
    double* val; // size=nnz，变量
} dCSRMatrix;

typedef struct {
    int m; // 行数
    int nnz; // 非零元素数
    int idx; // base-idx，指示从0开始计数还是从1开始计数
    int* ptrB; // 索引指向每行的Begin位置，size=行数, 指向的val[i]的下标=ptrB[i]-idx
    int* colind; // size=nnz，列坐标
    double complex* val; // size=nnz，变量
} cCSRMatrix;

#endif //CUBEGPA6SAB_MATRIX_H
