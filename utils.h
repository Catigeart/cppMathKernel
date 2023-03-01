//
// Created by Catigeart on 2022/11/23.
//

#ifndef MATRIX_V2_UTILS_H
#define MATRIX_V2_UTILS_H

#include <cstdio>
#include <cstdlib>

/// 宏语法糖，用于简化矩阵和向量的坐标写法
#define A(x, y) a[x + (y) * lda]
#define B(x, y) b[x + (y) * ldb]
#define C(x, y) c[x + (y) * ldc]
#define AA(x, y) aa[x + (y) * lda]
#define AT(x, y) at[x + (y) * ldat]
#define L(x, y) l[x + (y) * ldl]
#define U(x, y) u[x + (y) * ldu]
#define Q(x, y) q[x + (y) * ldq]
#define R(x, y) r[x + (y) * ldr]
#define P(x, y) p[x + (y) * ldp]
#define T(x, y) t[x + (y) * ldt]
#define Rn(x, y) rn[x + (y) * ldrn]
#define Pn(x, y) pn[x + (y) * ldpn]
#define N(x, y) n[x + (y) * ldn]
#define V(x, y) v[x + (y) * ldv]
#define RANGE(x, y) range[x + (y) * ldr]
#define NULLSP(x, y) nullsp[x + (y) * ldn]

#define x(n) x[(n) * incx]
#define y(n) y[(n) * incy]
#define u(n) u[(n) * incu]

/// 内存分配语法糖
#define DMALLOC(size) (double *) safe_malloc(size * sizeof(double))

/// malloc()封装
void* safe_malloc(size_t size);

/// free()封装
void safe_free(void *p);

/// 从文件中读取【已分配内存】的矩阵
void read_matrix(const int& m, const int& n, double *a, const int& lda, FILE *f);

/// 从控制台中读取【已分配内存】的矩阵
void input_matrix(const int& m, const int& n, double *a, const int& lda);

/// 输出矩阵
void print_matrix(const int& m, const int& n, const double *a, const int& lda);

/// 从文件中读取【已分配内存】的向量
void read_vector(const int& n, double* x, const int& incx, FILE *f);

/// 从控制台中读取【已分配内存】的向量
void input_vector(const int& n, double* x, const int& incx);

/// 输出向量
void print_vector(const int& n, const double* x, const int& incx);

#endif //MATRIX_V2_UTILS_H
