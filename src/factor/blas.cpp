//
// Created by Catigeart on 2022/11/23.
//
#include <cstring>
#include <cmath>
#include <cfloat>
#include "factor/blas.h"
#include "factor/utils.h"

int damax(const int& n, const double* x, const int& incx) {
    auto maxVal = -DBL_MAX;
    int idx = -INT_MAX;
    for (int i = 0; i < n; ++i) {
        if(x(i) > maxVal) {
            maxVal = x(i);
            idx = i;
        }
    }
    return idx;
}

int daabsmax(const int& n, const double* x, const int& incx) {
    auto maxVal = DBL_MIN;
    int idx = -INT_MAX;
    for (int i = 0; i < n; ++i) {
        if ((fabs(x(i)) < 1e-6)) continue;
        if (fabs(x(i)) > maxVal) {
            maxVal = fabs(x(i));
            idx = i;
        }
    }
    return idx;
}

void dswap(const int& n, double* x, const int& incx, double* y, const int& incy) {
    for (int i = 0; i < n; ++i) {
        double tmp = x(i);
        x(i) = y(i);
        y(i) = tmp;
    }
}

/**
 * y := x
 * @param n
 * @param x
 * @param incx
 * @param y
 * @param incy
 */
void dcopy(const int& n, const double* x, const int& incx, double* y, const int& incy) {
    for (int i = 0; i < n; ++i) {
        y(i) = x(i);
    }
}

double ddot(const int& n, const double* x, const int& incx, const double* y, const int& incy) {
    double ret = 0;
    for (int i = 0; i < n; ++i) {
        ret += x(i) * y(i);
    }
    return ret;
}

/**
 * norm for vector
 * @param n
 * @param x
 * @param inc
 * @return
 */
double dnormv(const int& n, const double* x, const int& incx) {
    return sqrt(ddot(n, x, incx, x, incx));
}

/**
 * m, n是源尺寸，不是转置后的目标尺寸。mn->nm
 * @param m
 * @param n
 * @param a_dst
 * @param lda
 * @param b_src
 * @param ldb
 */
void dtrans(const int& m, const int& n, double* a_dst, const int& lda, const double* b_src, const int& ldb) {
    double *&a = a_dst;
    const double *&b = b_src;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A(j, i) = B(i, j);
        }
    }
}

void dtransin(const int& n, double *a, const int& lda) {
    double tmp;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            tmp = A(i, j);
            A(i, j) = A(j, i);
            A(j, i) = tmp;
        }
    }
}

/**
 * A := B
 * @param m
 * @param n
 * @param a_dst
 * @param lda
 * @param b_src
 * @param ldb
 */
void dmatcpy(const int& m, const int& n, double* a, const int& lda, const double* b, const int& ldb) {
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            A(i, j) = B(i, j);
        }
    }
}

void dempty(const int& m, const int& n, double* a, const int& lda) {
    for (int j = 0; j < n; ++j) {
        memset(a + j * lda, 0, m * sizeof(double));
    }
}

/**
 * y=op(A)x
 */
void dgmv(BLAS_TRANS trans, const int& m, const int& n,
          const double* a, const int& lda, const double* x, const int& incx, double* y, const int& incy) {

    /// op(A)=A^{T}, A^{T}_{n,m},x_{m},y_{n}
    if (trans == BlasTrans) {
        memset(y, 0, n * incy * sizeof(double));

            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                y(j) += A(i, j) * x(i);
            }
        }
    }
    /// op(A)=A, A_{m, n} ,x_{n}, y_{m}
    else {
        memset(y, 0, m * incy * sizeof(double));
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                y(i) += A(i, j) * x(j);
            }
        }
    }
}

/**
 * C_{m,n} := A_{m,k}*B_{k,n}
 * @param m
 * @param n
 * @param k
 * @param a
 * @param lda
 * @param b
 * @param ldb
 * @param c
 * @param ldc
 */
void dgmm(BLAS_TRANS trans, const int& m, const int& n, const int& k,
         const double* a, const int& lda, const double* b, const int& ldb, double* c, const int& ldc) {
    if (trans == BlasNoTrans) {
        const int &ldat = k;
        auto *at = DMALLOC(ldat * m);
        dtrans(m, k, at, ldat, a, lda);
        for (int im = 0; im < m; ++im) {
            for (int in = 0; in < n; ++in) {
                C(im, in) = 0; // set empty
                for (int ik = 0; ik < k; ++ik) {
                    C(im, in) += AT(ik, im) * B(ik, in);
                }
                //printf("C[%d][%d]=%.2lf\n", im, in, C(im, in));
            }
        }

        safe_free(at);
    }
    else { // BlasTrans

        for (int im = 0; im < m; ++im) {
            for (int in = 0; in < n; ++in) {
                C(im, in) = 0; // set empty
                for (int ik = 0; ik < k; ++ik) {
                    C(im, in) += A(ik, im) * B(ik, in);
                }
            }
        }
    }

}

/**
 * Ax=b, while A is a triangle matrix.
 * @param uplo
 * @param n
 * @param a
 * @param lda
 * @param x A vector contains b's elements. It contains x's elements after being updated.
 * @param incx
 */
void dtrsv(const BLAS_UPLO uplo, const int& n, const double* a, const int& lda, double* x, const int& incx) {
    const int &m = n;
    if (uplo == BlasUpper) {
        x(m - 1) /= A(m - 1, n - 1);
        for (int i = m - 2; i >= 0; --i) {
            for (int j = n - 1; j > i; --j) {
                x(i) -= A(i, j) * x(j);
            }
            x(i) /= A(i, i);
        }
    } else { // uplo == blas_lower
        x(0) /= A(0, 0);
        for (int i = 1; i < m; ++i) {
            for (int j = 0; j < i; ++j) {
                x(i) -= A(i, j) * x(j);
            }
            x(i) /= A(i, i);
        }
    }
}

double ddiagdot(const int& n, const double* a, const int& lda) {
    double ret = A(0, 0);
    for (int j = 1; j < n; ++j)
        ret *= A(j, j);
    return ret;
}