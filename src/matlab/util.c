//
// Created by Catigeart on 2023/6/21.
//
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include <time.h>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef WIN32
int gettimeofday(struct timeval *tp, void *tzp)
{
    time_t clock;
    struct tm tm;
    SYSTEMTIME wtm;

    GetLocalTime(&wtm);
    tm.tm_year   = wtm.wYear  - 1900;
    tm.tm_mon    = wtm.wMonth - 1;
    tm.tm_mday   = wtm.wDay;
    tm.tm_hour   = wtm.wHour;
    tm.tm_min    = wtm.wMinute;
    tm.tm_sec    = wtm.wSecond;
    tm.tm_isdst  = -1;

    clock = mktime(&tm);
    tp->tv_sec   = clock;
    tp->tv_usec  = wtm.wMilliseconds * 1000;
    return (0);
}
#endif

void* safe_malloc(size_t size) {
    void *p = malloc(size);
    if (!p) {
        fprintf(stderr, "Malloc failed!");
        exit(0);
    }
    return p;
}

void safe_free(void *p) {
    if (p)
        free(p);
}

double get_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec * 1.e-6;
}

void print_dmat(char* name, int m, int n, int lda, double* a) {
    printf("Matrix %s:\n", name);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("\t%.2lf", a[i + j * lda]);
        }
        printf("\n");
    }
}

void print_cmat(char* name, int m, int n, int lda, double complex* a) {
    printf("Matrix %s:\n", name);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("\t%+.2lf%+.2lfi", creal(a[i + j * lda]), cimag(a[i + j * lda]));
        }
        printf("\n");
    }
}

void print_imat(char* name, int m, int n, int lda, int* a) {
    printf("Matrix %s:\n", name);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("\t%d", a[i + j * lda]);
        }
        printf("\n");
    }
}

void print_dcoo(char* name, double* val, int* rowind, int* colind, int nnz, int base_idx) {
    printf("Sparse matrix %s:\n", name);
    for (int i = 0; i < nnz; i++) {
        printf("\t{(%d, %d) : %.2lf}\n", rowind[i] + base_idx, colind[i] + base_idx, val[i]);
    }
}

void print_ccoo(char* name, double complex* val, int* rowind, int* colind, int nnz, int base_idx) {
    printf("Sparse matrix %s:\n", name);
    for (int i = 0; i < nnz; i++) {
        printf("\t{(%d, %d) : %+.2lf%+.2lfi}\n", rowind[i] + base_idx, colind[i] + base_idx, creal(val[i]), cimag(val[i]));
    }
}

void print_dcsr(char* name, int m, int nnz, double* a, int* ia, int* ja) {
    printf("CSR matrix %s:\n", name);
    printf("\tCSR base-idx = %d\n", ia[0]);
    printf("\tValues = [%.2lf", a[0]);
    for (int i = 1; i < nnz; i++) {
        printf(", %.2lf", a[i]);
    }
    printf("]\n\tCols = [%d", ja[0]);
    for (int i = 1; i < nnz; i++) {
        printf(", %d", ja[i]);
    }
    printf("]\n\tRow ptr = [%d", ia[1]);
    for (int i = 2; i <= m; i++) {
        printf(", %d", ia[i]);
    }
    printf("]\n");
}

void print_ccsr(char* name, int m, int nnz, double complex* a, int* ia, int* ja) {
    printf("CSR matrix %s:\n", name);
    printf("\tCSR base-idx = %d\n", ia[0]);
    printf("\tValues = [%+.2lf%+.2lfi", creal(a[0]), cimag(a[0]));
    for (int i = 1; i < nnz; i++) {
        printf(", %+.2lf%+.2lfi", creal(a[i]), cimag(a[i]));
    }
    printf("]\n\tCols = [%d", ja[0]);
    for (int i = 1; i < nnz; i++) {
        printf(", %d", ja[i]);
    }
    printf("]\n\tRow ptr = [%d", ia[1]);
    for (int i = 2; i <= m; i++) {
        printf(", %d", ia[i]);
    }
    printf("]\n");
}