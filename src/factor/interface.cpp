//
// Created by Catigeart on 2022/11/27.
//
#include "factor/interface.h"
#include "factor/linpack.h"
#include "factor/blas.h"
#include "factor/utils.h"

int Rank(const int &m, const int &n, const double *a, const int &lda) {
    const int &ldb = m;
    auto *b = DMALLOC(ldb * n);
    int rank;
    gaussJordan(m, n, a, lda, b, ldb, rank);
    safe_free(b);
    return rank;
}
double Det(const int& n, const double* a, const int& lda) {
    if (Rank(n, n, a, lda) < n)
        return 0;
    const int &ldl = n, ldu = n, ldp = n;
    auto *l = DMALLOC(ldl * n);
    auto *u = DMALLOC(ldu * n);
    auto *p = DMALLOC(ldp * n);
    int swap_cnt;
    lup(n, a, lda, l, ldl, u, ldu, p, ldp, swap_cnt);
    double det = ddiagdot(n, u, ldu);
    if (swap_cnt % 2 == 1)
        det = -det;
    safe_free(l);
    safe_free(u);
    safe_free(p);
    return det;
}

void GaussJordan(const int& m, const int& n, const double* a, const int&lda, double* b, const int& ldb) {
    int rank;
    gaussJordan(m, n, a, lda, b, ldb, rank);
}

EXEC_STATE LU(const int& n, const double* a, const int& lda,
                           double* l, const int& ldl, double* u, const int&ldu, double& det) {
    if (Rank(n, n, a, lda) < n)
        return ExecNoFullRank;
    if (!lu(n, a, lda, l, ldl, u, ldu))
        return ExecLUZeroPivot;
    det = ddiagdot(n, u, ldu);
    return ExecSucc;
}

EXEC_STATE LUP(const int& n, const double* a, const int& lda,
               double* l, const int& ldl, double* u, const int&ldu, double* p, const int& ldp, double& det) {
    if (Rank(n, n, a, lda) < n)
        return ExecNoFullRank;
    int swap_cnt;
    lup(n, a, lda, l, ldl, u, ldu, p, ldp, swap_cnt);
    det = ddiagdot(n, u, ldu);
    if (swap_cnt % 2 == 1)
        det = -det;
    return ExecSucc;
}

EXEC_STATE QR(const int& n, const double* a, const int& lda,
              double* q, const int& ldq, double* r, const int&ldr, QRType type) {
    if (Rank(n, n, a, lda) < n)
        return ExecNoFullRank;
    switch (type) {
        case SchimidtType:
            gs_qr(n, n, a, lda, q, ldq, r, ldr);
            break;
        case HouseholderType:
            hh_qr(n, n, a, lda, q, ldq, r, ldr);
            break;
        case GivensType:
            gv_qr(n, n, a, lda, q, ldq, r, ldr);
    }
    return ExecSucc;
}

void URV(const int& m, const int& n, const double* a, const int& lda,
               double* u, const int& ldu, double* r, const int& ldr, double* v, const int& ldv) {
    urv(m, n, a, lda, u, ldu, r, ldr, v, ldv);
}

void LUSolve(const int& n, const double* l, const int& ldl, const double* u, const int&ldu, double* x, const int& incx) {
    dtrsv(BlasLower, n, l, ldl, x, incx);
    dtrsv(BlasUpper, n, u, ldu, x, incx);
}

void LUPSolve(const int& n, const double* l, const int& ldl, const double* u, const int&ldu,
              const double* p, const int& ldp, double* b, const int& incb, double* x, const int& incx) {
    dgmv(BlasNoTrans, n, n, p, ldp, b, incb, x, incx);
    // printf("Pb=b': ");
    // print_vector(n, x, incx);
    dtrsv(BlasLower, n, l, ldl, x, incx);
    // printf("Ly=b': ");
    // print_vector(n, x, incx);
    dtrsv(BlasUpper, n, u, ldu, x, incx);
    // printf("Ux=y:  ");
}

void QRSolve(const int& n, const double* q, const int& ldq, const double* r, const int& ldr,
             const double* b, const int& incb, double* x, const int& incx) {
    dgmv(BlasTrans, n, n, q, ldq, b, incb, x, incx);
    dtrsv(BlasUpper, n, r, ldr, x, incx);
}