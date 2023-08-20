//
// Created by Catigeart on 2023/6/21.
//
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matlab/blas.h"
#include "matlab/matrix.h"
#include "matlab/util.h"

// matlab函数，向靠近零的方向舍入
void fix(double* X, int* Y, int len) {
    for (int i = 0; i < len; i++) {
        if (X[i] > 0) {
            Y[i] = floor(X[i]);
        }
        else {
            Y[i] = ceil(X[i]);
        }
    }
}

void deye(int n, double* a) {
    memset(a, 0, n * sizeof(double));
    for (int i = 0; i < n; i++) {
        a[i + i * n] = 1;
    }
}

void dkron(int xm, int xn, const double* x, int ym, int yn, const double* y, double* z) {
    for (int j = 0; j < xn * yn; j++) {
        int x_col = j / yn;
        int y_col = j % yn;
        for (int i = 0; i < xm; i++) {
            dscal2v(ym, x[i + x_col * xm], &y[0 + y_col * ym], 1, &z[i*ym + j * xm*ym], 1);
        }
    }
}

void cnorm2(int n, const double complex* x, int incx, double complex* norm_) {
    *norm_ = 0;
    for (int i = 0; i < n; i++) {
        *norm_ += sqrt(x[i * incx] * x[i * incx]);
    }
}

