//
// Created by Catigeart on 2022/11/23.
//
#include <cstdio>
#include "factor/linpack.h"
#include "factor/blas.h"
#include "factor/utils.h"
#include <cmath>
/*
void dgmv_test() {
    FILE *f = fopen("dgmv_test.txt", "r");
    int t, m, n, &lda = m, incx = 1, &incy = incx;
    fscanf(f, "%d", &t);
    while (t--) {
        fscanf(f, "%d %d", &m, &n);
        auto *a = (double *) safe_malloc(lda * n * sizeof(double));
        auto *x = (double *) safe_malloc(incx * n * sizeof(double));
        auto *y = (double *) safe_malloc(incy * m * sizeof(double));
        read_matrix(m, n, a, lda, f);
        read_vector(n, x, incx, f);
        dgmv(BlasTrans, m, n, a, lda, x, incx, y, incy);
        print_vector(n, y, incy);
        dgmv(BlasNoTrans, m, n, a, lda, x, incx, y, incy);
        print_vector(n, y, incy);
        safe_free(a);
        safe_free(x);
        safe_free(y);
    }
    fclose(f);
}

void Householder_test() {
    FILE *f = fopen("hh_test.txt", "r");
    int t, n, &ldr = n, incx = 1;
    fscanf(f, "%d", &t);
    while (t--) {
        fscanf(f, "%d", &n);
        auto *r = DMALLOC(ldr * n);
        auto *x = DMALLOC(n * incx);
        read_vector(n, x, incx, f);
        Householder(n, r, ldr, x, incx);
        print_matrix(n, n, r, ldr);
        safe_free(r);
        safe_free(x);
    }
}

void GaussJordan_RNTest() {
    FILE *f = fopen("gjrn_test.txt", "r");
    int t, m, n, &ldr = m, &ldn = n, &lda = m, &ldb = m, rank = 0;
    fscanf(f, "%d", &t);
    while (t--) {
        fscanf(f, "%d %d", &m, &n);
        auto *a = (double*) safe_malloc(lda * n * sizeof(double));
        auto *b = DMALLOC(ldb * n);
        read_matrix(m, n, a, lda, f);
        GaussJordan(m, n, a, lda, b, ldb, rank);
        printf("A:\n");
        print_matrix(m, n, a, lda);
        printf("Modified A:\n");
        print_matrix(m, n, b, ldb);


        auto *range = DMALLOC(ldr * rank);
        auto *nullsp = DMALLOC(ldn * (n - rank));
        Gauss2RankRangeNull(m, n, a, lda, b, ldb, range, ldr, nullsp, ldn, rank);
        printf("range:\n");
        print_matrix(m, rank, range, ldr);
        printf("nullsp:\n");
        print_matrix(n, n - rank, nullsp, ldn);



        safe_free(a);
        safe_free(b);
        safe_free(range);
        safe_free(nullsp);
    }
    fclose(f);
}


void gj(const int& m, const int& n, const double* a, const int& lda, double *b, const int& ldb, int& rank) {
    if (m == 1) {
        dcopy(n, a, lda, b, ldb);
        return;
    }

    rank = 0;
    dmatcpy(m, n, b, ldb, a, lda);

    // print_matrix(m, n, a, lda);
    // print_matrix(m, n, b, ldb);

    for (int i = 0, j = 0; i < m && j < n; ++i) { // 行阶梯形的情况下可能一次移动多个列，因此j需要手动加
        int idx;

        while (j < n) {
            idx = daabsmax(m - i, b + i + j * ldb, 1); // 找下一行的列最大主元
            // 没有找到或找到的最小元素是零，说明为非基本列
            if (idx == -1)  {
                ++j;
                continue;
            }
            idx += i;
            if (fabs(B(idx, j)) < 1e-6) {
                ++j; // 找下一列
                continue;
            }
            break;
        }
        if (idx != i) {
            dswap(n, b + i, ldb, b + idx, ldb); // 交换行
        }
        ++rank;



        for (int ii = 0; ii < m; ++ii) { // Jordan,每一行都要把主元位置消掉
            [[unlikely]] if (i == ii)  { // 当前行跳过，最后再把主元处理成1
                continue;
            }

            double coef = B(ii, j) / B(i, j); // 求得比值
            B(ii, j) = 0; // 直接给当前位置赋0

            for (int jj = j + 1; jj < n; ++jj) { // 从下一个位置开始逐个开始减
                B(ii, jj) -= B(i, jj) * coef;
            }
        }

        // 使主元为1
        for (int jj = j + 1; jj < n; ++jj) {
            B(i, jj) /= B(i, j);
        }
        B(i, j) = 1;

        ++j; //  该列计算完成，前进到下一列

        //printf("row %d: \n", i);
        //print_matrix(m, n, b, ldb);

    }

}

void memTest() {
    FILE *f = fopen("gjrn_test.txt", "r");
    int t, m, n, &ldr = m, &ldn = n, &lda = m, &ldb = m, rank = 0;
    fscanf(f, "%d", &t);
    while (t--) {
        fscanf(f, "%d %d", &m, &n);
        auto *a = (double*) safe_malloc(lda * n * sizeof(double));
        auto *b = DMALLOC(ldb * n);

        read_matrix(m, n, a, lda, f);
        gj(m, n, a, lda, b, ldb, rank);
        printf("A:\n");
        print_matrix(m, n, a, lda);
        printf("Modified A:\n");
        print_matrix(m, n, b, ldb);

        safe_free(a);
        auto *range = DMALLOC(ldr * rank);
        auto *nullsp = DMALLOC(ldn * (n - rank));
        Gauss2RankRangeNull(m, n, b, ldb, range, ldr, nullsp, ldn, rank);
        printf("range:\n");
        print_matrix(m, rank, range, ldr);
        printf("nullsp:\n");
        print_matrix(n, n - rank, nullsp, ldn);

        safe_free(a);
        safe_free(b);
        // safe_free(range);
        // safe_free(nullsp);
    }
    fclose(f);
}


void dgmmTest() {
    FILE *f = fopen("dgmm_test.txt", "r");
    int t, m, n, k, &lda = m, &ldb = k, &ldc = m;
    fscanf(f, "%d", &t);
    while (t--) {
        fscanf(f,"%d %d %d", &m, &n, &k);
        auto *a = DMALLOC(lda * k);
        auto *b = DMALLOC(ldb * n);
        auto *c = DMALLOC(ldc * n);
        read_matrix(m, k, a, lda, f);
        read_matrix(k, n, b, ldb, f);
        dgmm(BlasTrans, m, n, k, a, lda, b, ldb, c, ldc);
        print_matrix(m, n, c, ldc);

        safe_free(a);
        safe_free(b);
        safe_free(c);
    }
}*/