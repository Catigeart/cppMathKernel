//
// Created by Catigeart on 2023/6/21.
//

#ifndef CUBEGPA6SAB_UTIL_H
#define CUBEGPA6SAB_UTIL_H

#include <stdlib.h>

void* safe_malloc(size_t size);

void safe_free(void *p);

double get_time();

void print_dmat(char* name, int m, int n, int lda, double* a);

void print_cmat(char* name, int m, int n, int lda, double complex* a);

void print_imat(char* name, int m, int n, int lda, int* a);

void print_dcoo(char* name, double* val, int* rowind, int* colind, int nnz, int base_idx);

void print_ccoo(char* name, double complex* val, int* rowind, int* colind, int nnz, int base_idx);

void print_dcsr(char* name, int m, int nnz, double* a, int* ia, int* ja);

void print_ccsr(char* name, int m, int nnz, double complex* a, int* ia, int* ja);

#endif //CUBEGPA6SAB_UTIL_H
