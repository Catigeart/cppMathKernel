//
// Created by Catigeart on 2022/11/23.
//
#include <cstdio>
#include <cstdlib>
#include "factor/utils.h"

void* safe_malloc(size_t size) {
    void *p = malloc(size);
    if (!p) {
        fprintf(stderr, "Malloc failed");
        exit(0);
    }
    return p;
}

void safe_free(void *p) {
    if (p)
        free(p);
}

void read_matrix(const int& m, const int& n, double *a, const int& lda, FILE *f) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            fscanf(f, "%lf", &A(i, j));
        }
    }
}

void input_matrix(const int& m, const int& n, double *a, const int& lda) {
    // a = DMALLOC(lda * n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            scanf("%lf", &A(i, j));
        }
    }
}

void print_matrix(const int& m, const int& n, const double *a, const int& lda) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%.2lf\t", A(i, j));
        }
        putchar('\n');
    }
}

void read_vector(const int& n, double* x, const int& incx, FILE *f) {
    for (int i = 0; i < n; ++i) {
        fscanf(f, "%lf", &x(i));
    }
}

void input_vector(const int& n, double* x, const int& incx) {
    // x = DMALLOC(n * incx);
    for (int i = 0; i < n; ++i) {
        scanf("%lf", &x(i));
    }
}

void print_vector(const int& n, const double* x, const int& incx) {
    for (int i = 0; i < n; ++i) {
        printf("%.2lf\t", x(i));
    }
    putchar('\n');
}