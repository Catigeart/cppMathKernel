//
// Created by Catigeart on 2022/11/22.
//

/// ##ColMajor## ///

#include <cstring>
#include <cassert>
#include <cmath>
#include "factor/linpack.h"
#include "factor/blas.h"
#include "factor/utils.h"

void lup(const int& n, const double* a, const int& lda,
               double* l, const int& ldl, double* u, const int& ldu, double* p, const int& ldp, int& swap_cnt) {
    swap_cnt = 0;
    int pFlag[n];
    for (int j = 0; j < n; ++j) {
        pFlag[j] = j;
        for (int i = 0; i < j; ++i) {
            L(i, j) = 0;
        }
        L(j, j) = 1;
        /*
        for (int i = j + 1; i < n; ++i) {
            U(i, j) = 0;
        }
        */
    }

    // printf("A:\n");
    // print_matrix(n, n, a, lda);

    // 第一次求解，并初始化
    int idx = daabsmax(n, a, 1);
    std::swap(pFlag[idx], pFlag[0]);
    dcopy(n, a + idx, lda, u, ldu);
    for (int i = 1; i < n; ++i) {
        U(i, 0) = 0;
        if (i != idx) {
            swap_cnt++;
            L(i, 0) = A(i, 0) / U(0, 0);
            for (int j = 1; j < n; ++j) {
                U(i, j) = A(i, j) - U(0, j) * L(i, 0);
            }
        }
        else {
            L(i, 0) = A(0, 0) / U(0, 0);
            for (int j = 1; j < n; ++j) {
                U(i, j) = A(0, j) - U(0, j) * L(i, 0);
            }
        }
    }

     // printf("L0:\n");
     // print_matrix(n, n, l, ldl);
     // printf("U0:\n");
     // print_matrix(n, n, u, ldu);

    for (int i = 1; i < n; ++i) {
        idx = daabsmax(n - i, u + i + i * ldu, 1) + i;
        if (i != idx) {
            dswap(i, l + i, ldl, l + idx, ldl);
            dswap(n - i, u + i + i * ldu, ldu, u + idx + i * ldu, ldu);
            std::swap(pFlag[idx], pFlag[i]);
            swap_cnt++;
        }
        for (int ii = i + 1; ii < n; ++ii) {
            L(ii, i) = U(ii, i) / U(i, i);
            U(ii, i) = 0;
            for (int j = i + 1; j < n; ++j) {
                U(ii, j) = U(ii, j) - U(i, j) * L(ii, i);
            }
        }

         // printf("L%d:\n", i);
         // print_matrix(n, n, l, ldl);
         // printf("U%d:\n", i);
         // print_matrix(n, n, u, ldu);
    }

    for (int j = 0; j < n; ++j) {
        memset(p + j * ldp, 0, n * sizeof(double));
    }
    for (int i = 0; i < n; ++i) {
        P(i, pFlag[i]) = 1;
    }


    // printf("P:\n");
    // print_matrix(n, n, p, ldp);
}

bool lu(const int& n, const double* a, const int& lda, double* l, const int& ldl, double* u, const int& ldu) {
    // TODO: opt
    memset(l, 0, ldl * n * sizeof(double));
    memset(u, 0, ldu * n * sizeof(double));

    for (int i = 0; i < n; ++i) {
        // L
        L(i, i) = 1;
        for (int j = 0; j < i; ++j) {
            if (U(j, j) == 0) {
                return false; // 此时l,u是被修改过的，原有的数据已被破坏
            }
            L(i, j) = A(i, j);
            for (int k = 0; k < j; ++k) {
                L(i, j) -= L(i, k) * U(k, j);
            }
            L(i, j) /= U(j, j);
        }
        // U
        for (int j = i; j < n; ++j) {
            U(i, j) = A(i, j);
            for (int k = 0; k < i; k++) {
                U(i, j) -= L(i, k) * U(k, j);
            }
        }
    }

    return true;
}



/**
 * Q_{m,n}, R{n,n}
 * @param m
 * @param n
 * @param a
 * @param lda
 * @param q
 * @param ldq
 * @param r
 * @param ldr
 */
void gs_qr(const int& m, const int& n, const double *a, const int& lda, double* q, const int& ldq, double* r, const int& ldr) {
    auto *aa = (double *) safe_malloc(lda * n * sizeof(double));
    memcpy(aa, a, lda * n * sizeof(double));
    memset(r, 0, ldr * n * sizeof(double));
    for (int k = 0; k < n; ++k) {
        double norm2 = 0;
        for (int i = 0; i < n; ++i) {
            norm2 += AA(i, k) * AA(i, k);
        }
        R(k, k) = sqrt(norm2);
        for (int i = 0; i < m; ++i) {
            Q(i, k) = AA(i, k) / R(k, k);
        }
        for (int j = k + 1; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                R(k, j) += AA(i, j) * Q(i, k);
            }
            for (int i = 0; i < m; ++i) {
                AA(i, j) -= R(k, j) * Q(i, k);
            }
        }
    }
    safe_free(aa);
}

void householder(int n, double* r, const int& ldr, const double* x, const int& incx) {
   // printf("-----Householder-----\n");
   // printf("x:\n");
  // print_vector(n, x, incx);

    int incu = 1;
    auto *u = DMALLOC(n * incu);
    for (int i = 0; i < n; ++i) {
        u(i) = x(i);
    }
    u(0) -= dnormv(n, x, incx);
    dgmm(BlasNoTrans, n, n, 1, u, 1, u, 1, r, ldr);
    double unorm2 = ddot(n, u, incu, u, incu);
    //print_matrix(n, n, r, ldr);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            R(i, j) *= 2;
            R(i, j) /= unorm2;
            R(i, j) = -R(i, j);
        }
        R(j, j) += 1;
    }
    safe_free(u);

   // printf("R:\n");
   // print_matrix(n, n, r, ldr);
   // printf("----------\n");
}

/**
 * A_{m,n}=Q_{m,m}R_{m,n}
 * @param m
 * @param n
 * @param a
 * @param lda
 * @param q
 * @param ldq
 * @param r
 * @param ldr
 */
void hh_qr(const int& m, const int& n, const double *a, const int& lda, double* q, const int& ldq, double* r, const int& ldr) {
    // PA=T,所以T就是R。用别名提高可读性
    // Rn_{m,m}, P_{m,m}, T_{m,n}
    double *&t = r;
    const int& ldt = ldr;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            T(i, j) = A(i, j); // 初始化T
        }
    }
    const int &ldrn = m, &ldp = m;
    // R1
    auto *rn = (double *) safe_malloc(ldrn * m * sizeof(double));
    householder(m, rn, ldrn, t, 1); // 第一列

    //printf("col 0:\n");
    //print_matrix(m, m, rn, ldrn);

    // 初始化 P,T
    auto *p = (double *) safe_malloc(ldp * m * sizeof(double));
    dmatcpy(m, m, p, ldp, rn, ldrn);
    // dmatcpy(m, n, t, ldt, a, lda);

    // two tmp matrix
    // b临时保存P \gets RnP的结果，维度为{m,m}；//{{{c临时保存T \gets PA的结果，维度为{m,n}}}}
    const int &ldb = m, &ldc = m;
    auto *b = (double *) safe_malloc(ldb * m * sizeof(double));
    auto *c = (double *) safe_malloc(ldc * n * sizeof(double));

    // 更新t
    dgmm(BlasNoTrans, m, n, m, rn, ldrn, t, ldt, c, ldc);
    dmatcpy(m, n, t, ldt, c, ldc);
    //printf("R_0A:\n");
    //print_matrix(m, n, t, ldt);

    for (int j = 1; j < n - 1; ++j) { // 最后一列不用再求了
        // 利用指针偏移传入子矩阵
        householder(m - j, rn + j + j * ldrn, ldrn, t + ldt * j + j, 1); // 第j列
        // 补全rn的左上角部分
        for (int k = 0; k < j; ++k) {
            Rn(k, k) = 1;
            // memset(rn + (k + 1) + k * ldrn, 0, n - k - 1); // 置列为0
            for (int kk = k + 1; kk < n; ++kk) {
                Rn(k, kk) = 0;
                Rn(kk, k) = 0;
            }
        }

        //printf("col %d:\n", j);
        //print_matrix(m, m, rn, ldrn);

        // 更新p
        dgmm(BlasNoTrans, m, m, m, rn, ldrn, p, ldp, b, ldb);
        dmatcpy(m, n, p, ldp, b, ldb);
        // 更新t
        dgmm(BlasNoTrans, m, n, m, rn, ldrn, t, ldt, c, ldc);
        dmatcpy(m, n, t, ldt, c, ldc);

        //printf("R_%dA:\n", j);
        //print_matrix(m, n, t, ldt);
    }
    // Q是P的转置
    dtrans(m, m, q, ldq, p, ldp);

    safe_free(rn);
    safe_free(p);
    safe_free(b);
    safe_free(c);
}

void givens(const int& n, double* p, const int& ldp, const double* x, const int& incx, const int& xi, const int& xj) {
    //printf("x:\n");
    //print_vector(n, x, incx);
    dempty(n, n, p, ldp);
    double sqrtij = sqrt(x(xi) * x(xi) + x(xj) * x(xj));
    P(xi, xi) = x(xi) / sqrtij;
    P(xj, xj) = P(xi, xi);
    P(xi, xj) = x(xj) / sqrtij;
    P(xj, xi) = -P(xi, xj);
    for (int i = 0; i < n; ++i) {
        if (i == xi || i == xj) [[unlikely]] {
            continue;
        }
        P(i, i) = 1;
    }
}

/**
 *
 * @param m
 * @param n
 * @param a
 * @param lda
 * @param q
 * @param ldq
 * @param r
 * @param ldr
 */
void gv_qr(const int& m, const int& n, const double *a, const int& lda, double* q, const int& ldq, double* r, const int& ldr) {
    // PA=T,所以T就是R。用别名提高可读性
    double *&t = r;
    const int &ldt = ldr, &ldp = m, &ldpn = m;
    dmatcpy(m, n, t, ldt, a, lda); // 初始化t
    auto *p = DMALLOC(ldp * m);
    // Pn
    auto *pn = DMALLOC(ldpn * m);
    // gemm临时矩阵
    const int &ldb = m, &ldc = m;
    auto *b = (double *) safe_malloc(ldb * m * sizeof(double));
    auto *c = (double *) safe_malloc(ldc * n * sizeof(double));

    givens(m, p, ldp, t, 1, 0, 1);
    //printf("P10:\n");
    //print_matrix(m, m, p, ldp);
    // 更新t
    dgmm(BlasNoTrans, m, n, m, p, ldp, t, ldt, c, ldc);
    dmatcpy(m, n, t, ldt, c, ldc);
    //printf("T:\n");
    //print_matrix(m, n, t, ldt);

    for (int j = 0; j < m - 1; ++j) {
        for (int i = j + 1; i < m; ++i) {
            if (j == 0 && i == 1) [[unlikely]] {
                continue; // 第一个Givens已经计算存入P了
            }
            givens(m, pn, ldpn, t + j * ldt, 1, j, i);
            //printf("P%d%d:\n", i, j);
            //print_matrix(m, m, pn, ldpn);
            // 更新p
            dgmm(BlasNoTrans, m, m, m, pn, ldpn, p, ldp, b, ldb);
            dmatcpy(m, n, p, ldp, b, ldb);
            //printf("P:\n");
            //print_matrix(m, m, p, ldp);
            // 更新t
            dgmm(BlasNoTrans, m, n, m, pn, ldpn, t, ldt, c, ldc);
            dmatcpy(m, n, t, ldt, c, ldc);
            //printf("T:\n");
            //print_matrix(m, n, t, ldt);
        }
    }

    // Q是P的转置
    dtrans(m, m, q, ldq, p, ldp);

    safe_free(pn);
    safe_free(p);
    safe_free(b);
    safe_free(c);
}

 /*
void GaussJordan(const int& m, const int& n, const double* a, const int& lda, double *b, const int& ldb, int& rank) {
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
*/

/**
 * range_{m, rank},nullspace_{n, n-rank}
 */
void gauss2RankRangeNull(const int& m, const int& n, const double* a, const int& lda,
                         const double* b, const int& ldb, double* range, const int& ldr, double* nullsp, const int& ldn, const int& rank) {

    // rank = 0;
    bool isPivot[n];
    int ip = 0, jp = 0;
    while (jp < n) {
        if(fabs(B(ip, jp)) > 1e-6) {
            isPivot[jp] = true;
            // ++rank;
            ++ip;
        }
        else {
            isPivot[jp] = false;
        }
        ++jp;
    }

    // 先置零。
    // TODO: 可优化
    for (int jj = 0; jj < n - rank; ++jj) {
        for (int ii = 0; ii < n; ++ii) {
            NULLSP(ii, jj) = 0;
        }
    }

    int rj = 0, nj = 0; // rj -> range 当前等待输入数据的列，nj -> nullspace 当前等待输入的列
    for (int j = 0; j < n; ++j) {
        if (isPivot[j]) {
            dcopy(m, a + j * lda, 1, range + rj * ldr, 1);
            // assert(rj < m);
            ++rj;
        }

        else {
            for (int i = 0; i < rank; ++i) {
                // assert((i + nj * ldn) < ldn * n);
                NULLSP(i, nj) = -B(i, j);
            }
            // dcopy(rank, a + j * lda, 1, nullsp + nj * ldn, 1);
            ++nj;
            // assert(nj < n);
        }

    }

    for (int j = 0; j < n - rank; ++j) {
        NULLSP(j + rank, j) = 1;
    }
}

void gaussJordan(const int& m, const int& n, const double* a, const int& lda, double *b, const int& ldb, int& rank) {
    rank = 0;
    dmatcpy(m, n, b, ldb, a, lda);

    for (int i = 0, j = 0; i < m && j < n; ++j) {
        int idx = daabsmax(m - i, b + i + j * ldb, 1) + i;
        if (idx < 0 || fabs(B(idx, j)) < 1e-6) {
            continue;
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

        // printf("row %d: \n", i);
        // print_matrix(m, n, b, ldb)b;

        ++i;
    }
}

/**
 *  U_{m,m}, R_{m,n} ,V_{n,n}
 */
void urv(const int& m, const int& n, const double* a, const int& lda,
         double* u, const int& ldu, double* r, const int& ldr, double* v, const int& ldv) {
    memset(u, 1, ldu * m * sizeof(double));
    memset(v, 1, ldv * n * sizeof(double));

    int rank = 0, rankT = 0;
    const int &ldb = m, &ldat = n, &ldbt = n;

    auto *b = DMALLOC(ldb * n);
    gaussJordan(m, n, a, lda, b, ldb, rank);
    //printf("A:\n");
    //print_matrix(m, n, a, lda);
    // printf("Modified A:\n");
    // print_matrix(m, n, b, ldb);
    gauss2RankRangeNull(m, n, a, lda, b, ldb, u, ldu, v + rank * ldv, ldv, rank);
    //printf("U:\n");
    //print_matrix(m, m, u, ldu);
    //printf("V:\n");
    //print_matrix(n, n, v, ldv);

    auto *at = DMALLOC(ldat * m);
    auto *bt = DMALLOC(ldbt * m);
    dtrans(m, n, at,ldat, a, lda);
    gaussJordan(n, m, at, ldat, bt, ldbt, rankT);
    //printf("A:\n");
    //print_matrix(m, n, a, lda);
    //printf("Modified A:\n");
    //print_matrix(m, n, b, ldb);

    gauss2RankRangeNull(n, m, at, ldat, bt, ldbt, v, ldv, u + rank * ldu, ldu, rankT);

    //printf("U:\n");
    //print_matrix(m, m, u, ldu);
    //printf("V:\n");
    //print_matrix(n, n, v, ldv);
/*
    const int &ldmm = m, &ldnn = n;
    auto *mm = DMALLOC(ldmm * m);
    auto *nn = DMALLOC(ldnn * n);
    int rk1, rk2;
    GaussJordan(m, m, u, ldu, mm, ldmm, rk1);
    GaussJordan(n, n, v, ldv, nn, ldnn, rk2);
    printf("rk1=%d,rk2=%d\n", rk1, rk2);
    safe_free(mm);
    safe_free(nn);
*/
    // 将u和v的向量长度标准化
    for (int j = 0; j < m; ++j) {
        double norm = dnormv(m, u + j * ldu, 1);
        for (int i = 0; i < m; ++i) {
            U(i, j) /= norm;
        }
    }
    for (int j = 0; j < n; ++j) {
        double norm = dnormv(n, v + j * ldv, 1);
        for (int i = 0; i < n; ++i) {
            V(i, j) /= norm;
        }
    }

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            if (i < rank && j < rank) {
                int incy = 1;
                auto *y = DMALLOC(incy * m);
                dgmv(BlasNoTrans, m, n, a, lda, v + j * ldv, 1, y, incy);
                R(i, j) = ddot(m, u + i * ldu, 1, y, incy);
            }
            else [[likely]] {
                R(i, j) = 0;
            }
        }
    }
    /*
    printf("v1:\n");
    print_matrix(m, n, r, ldr);
    const int &ldc = m;
    auto *c = DMALLOC(ldc * n); // 存放计算过程的临时矩阵
    dgmm(BlasTrans, m, n, m, u, ldu, a, lda, c, ldc);
    dgmm(BlasNoTrans, m, n, n, c, ldc, v, ldv, r, ldr);
    safe_free(c);
    printf("v2:\n");
    print_matrix(m, n, r, ldr);
     */
}