//
// Created by Catigeart on 2022/11/23.
//

#ifndef MATRIX_V2_LINPACK_H
#define MATRIX_V2_LINPACK_H

#include "utils.h"

/**
 * 行列式求解函数。
 * @param n n*n矩阵A的维度大小
 * @param a 矩阵A
 * @param lda 矩阵A的主维度
 * @return Det(A_{n,n})
 */
double Det(const int& n, const double* a, const int& lda);

/**
 * LUP分解函数。
 * @param n 矩阵A、L、U、P的行数/列数
 * @param a 矩阵A
 * @param lda A的主维度
 * @param l 矩阵L
 * @param ldl L的主维度
 * @param u 矩阵U
 * @param ldu U的主维度
 * @param p 矩阵P
 * @param ldp P的主维度
 * @param swap_cnt LUP分解过程中行交换的次数，为上层函数计算行列式提供依据。
 */
void lup(const int& n, const double* a, const int& lda,
         double* l, const int& ldl, double* u, const int& ldu, double* p, const int& ldp, int& swap_cnt);

/**
 * LU分解函数
 * @param n 矩阵A、L、U的行数/列数
 * @param a 矩阵A
 * @param lda A的主维度
 * @param l 矩阵L
 * @param ldl L的主维度
 * @param u 矩阵U
 * @param ldu U的主维度
 * @return 布尔值，若为真则成功分解，若为假则运算过程中出现0主元，分解失败。
 */
bool lu(const int& n, const double* a, const int& lda, double* l, const int& ldl, double* u, const int& ldu);

/**
 * 基于施密特正交化的QR分解函数，其中A_{m,n}，Q_{m,n}，R{n,n}
 * @param m A、Q的行数
 * @param n A、Q、R的列数、R的行数
 * @param a 矩阵A
 * @param lda A的主维度
 * @param q 矩阵Q的存放位置
 * @param ldq Q的主维度
 * @param r 矩阵R的存放位置
 * @param ldr R的主维度
 */
void gs_qr(const int& m, const int& n, const double *a, const int& lda, double* q, const int& ldq, double* r, const int& ldr);

/**
 * 对向量进行Householder分解。
 * @param n 向量x的长度、矩阵R的行数/列数
 * @param r 矩阵R的存放位置
 * @param ldr R的主维度
 * @param x 向量x
 * @param incx x的步长
 */
void householder(int n, double* r, const int& ldr, const double* x, const int& incx);

/**
 * 基于Householder的QR分解函数，其中A_{m,n}、Q_{n,n}、R_{m,n}
 * @param m A、R的行数
 * @param n A、Q、R的列数、Q的行数
 * @param a 矩阵A
 * @param lda A的主维度
 * @param q 矩阵Q
 * @param ldq Q的主维度
 * @param r 矩阵R
 * @param ldr R的主维度
 */
void hh_qr(const int& m, const int& n, const double *a, const int& lda, double* q, const int& ldq, double* r, const int& ldr);

/**
 * 基于Givens的QR分解函数，其中A_{m,n}、Q_{n,n}、R_{m,n}
 * @param m A、R的行数
 * @param n A、Q、R的列数、Q的行数
 * @param a 矩阵A
 * @param lda A的主维度
 * @param q 矩阵Q
 * @param ldq Q的主维度
 * @param r 矩阵R
 * @param ldr R的主维度
 */
void gv_qr(const int& m, const int& n, const double *a, const int& lda, double* q, const int& ldq, double* r, const int& ldr);

/**
 * 高斯-约旦分解函数，保存化简后的函数和秩的计算结果。
 * @param m A的行数
 * @param n A的列数
 * @param a 矩阵A
 * @param lda A的主维度
 * @param b 化简后的矩阵B存放位置
 * @param ldb B的主维度
 * @param rank 秩
 */
void gaussJordan(const int& m, const int& n, const double* a, const int& lda, double *b, const int& ldb, int& rank);

/**
 * 基于原矩阵以及高斯-约旦化简的行阶梯形式求解值空间和零空间的一组基。
 * @param m A、B的行数
 * @param n A、B的列数
 * @param a 矩阵A
 * @param lda A的主维度
 * @param b 高斯-约旦化简后的矩阵B
 * @param ldb B的主维度
 * @param range 由值空间的一组基组成的矩阵的存放位置
 * @param ldr range的主维度
 * @param nullsp 由零空间的一组基组成的矩阵的存放位置
 * @param ldn nullsp的主维度
 * @param rank 传入的秩，用于指示range和nullsp的求解
 */
void gauss2RankRangeNull(const int& m, const int& n, const double* a, const int& lda,
                         const double* b, const int& ldb, double* range, const int& ldr, double* nullsp, const int& ldn, const int& rank);

/**
 * URV分解函数。
 * @param m A、R的行数，U的行/列数
 * @param n A、R的列数，V的行/列数
 * @param a 矩阵A
 * @param lda A的主维度
 * @param u 矩阵U的存放位置
 * @param ldu U的主维度
 * @param r 矩阵R的存放位置
 * @param ldr R的主维度
 * @param v 矩阵V的存放位置
 * @param ldv V的主维度
 */
void urv(const int& m, const int& n, const double* a, const int& lda,
         double* u, const int& ldu, double* r, const int& ldr, double* v, const int& ldv);

#endif //MATRIX_V2_LINPACK_H
