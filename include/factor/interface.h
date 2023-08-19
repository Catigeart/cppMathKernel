//
// Created by Catigeart on 2022/11/27.
//

#ifndef MATRIX_V2_INTERFACE_H
#define MATRIX_V2_INTERFACE_H

/**
 * 接口函数的执行结果。
 */
enum EXEC_STATE {
    ExecSucc, /// 成功执行
    ExecNoFullRank, /// 不满秩
    ExecLUZeroPivot /// LU分解过程中出现0主元
};

/**
 * QR分解的类型。
 */
enum QRType {
    SchimidtType, /// 施密特分解
    HouseholderType, /// Householder分解
    GivensType /// 吉文斯分解
};

/**
 * 秩计算接口，计算A_{m,n}的秩。
 * @param m A的行数
 * @param n A的列数
 * @param a 矩阵A
 * @param lda A的主维度
 * @return 矩阵A_{m,n}的秩。
 */
int Rank(const int &m, const int &n, const double *a, const int &lda);

/**
 * 行列式计算接口，计算A_{n,n}的行列式
 * @param n A的维度大小
 * @param a 矩阵A
 * @param lda A的主维度
 * @return det(A_{n,n})
 */
double Det(const int& n, const double* a, const int& lda);

/**
 * 高斯-约旦化简接口，对A_{m,n}进行化简。
 * @param m A的行数
 * @param n A的列数
 * @param a 原矩阵A
 * @param lda 原矩阵A的主维度
 * @param b 修改后的矩阵存放位置
 * @param ldb 修改后的矩阵的主维度
 */
void GaussJordan(const int& m, const int& n, const double* a, const int&lda, double* b, const int& ldb);

/**
 * LU分解接口，对A_{n,n}进行分解。
 * @param n A的行数/列数
 * @param a 矩阵A
 * @param lda 矩阵A的主维度
 * @param l 矩阵L的存放位置
 * @param ldl 矩阵L的主维度
 * @param u 矩阵U的存放位置
 * @param ldu 矩阵U的主维度
 * @param det 利用U对角线求得的行列式值
 * @return 执行状态，包括成功执行、矩阵不满秩和运算过程中出现0主元三种。
 */
EXEC_STATE LU(const int& n, const double* a, const int& lda, double* l, const int& ldl, double* u, const int&ldu, double& det);

/**
 * LUP分解接口，对A_{n,n}进行分解。
 * @param n A的行数/列数
 * @param a 矩阵A
 * @param lda 矩阵A的主维度
 * @param l 矩阵L的存放位置
 * @param ldl 矩阵L的主维度
 * @param u 矩阵U的存放位置
 * @param ldu 矩阵U的主维度
 * @param p 矩阵P的存放位置
 * @param ldp 矩阵P的主维度
 * @param det 利用U对角线求得的行列式值
 * @return 执行状态，包括成功执行、矩阵不满秩两种。
 */
EXEC_STATE LUP(const int& n, const double* a, const int& lda,
               double* l, const int& ldl, double* u, const int&ldu, double* p, const int& ldp, double& det);

/**
 * QR分解接口，对Q_{n,n}进行分解。
 * @param n A的行数/列数
 * @param a 矩阵A
 * @param lda 矩阵A的主维度
 * @param q 矩阵Q的存放位置
 * @param ldq 矩阵Q的主维度
 * @param r 矩阵R的存放位置
 * @param ldr 矩阵R的主维度
 * @param type QR分解的类型，包括施密特正交化、Householder和Givens三种。
 * @return 执行状态，包括成功执行、矩阵不满秩两种。
 */
EXEC_STATE QR(const int& n, const double* a, const int& lda,
              double* q, const int& ldq, double* r, const int&ldr, QRType type);

/**
 * URV分解接口，对矩阵A_{m,n}进行分解
 * @param m A的行数
 * @param n A的列数
 * @param a 矩阵A
 * @param lda 矩阵A的主维度
 * @param u 矩阵U的存放位置
 * @param ldu 矩阵U的主维度
 * @param r 矩阵R的存放位置
 * @param ldr 矩阵R的主维度
 * @param v 矩阵V的存放位置
 * @param ldv 矩阵V的主维度
 */
void URV(const int& m, const int& n, const double* a, const int& lda,
               double* u, const int& ldu, double* r, const int& ldr, double* v, const int& ldv);

/**
 * 根据LU分解进行求解的接口。
 * @param n 矩阵L、U的行数/列数以及向量x/b的长度
 * @param l 矩阵L
 * @param ldl 矩阵L的主维度
 * @param u 矩阵U
 * @param ldu 矩阵U的主维度
 * @param x 向量x/b，运算前存储b的值，运算后存储计算后的x的结果
 * @param incx 向量x的步长
 */
void LUSolve(const int& n, const double* l, const int& ldl, const double* u, const int&ldu, double* x, const int& incx);

/**
 * 根据LUP分解进行求解的接口。
 * @param n 矩阵L、U、P的行数/列数以及向量x、b的长度
 * @param l 矩阵L
 * @param ldl 矩阵L的主维度
 * @param u 矩阵U
 * @param ldu 矩阵U的主维度
 * @param p 矩阵P
 * @param ldp 矩阵P的主维度
 * @param b 向量b
 * @param incb 向量b的步长
 * @param x 向量x的存放位置
 * @param incx 向量x的步长
 */
void LUPSolve(const int& n, const double* l, const int& ldl, const double* u, const int&ldu,
              const double* p, const int& ldp, double* b, const int& incb, double* x, const int& incx);

/**
 * 根据QR分解进行求解的接口。
 * @param n 矩阵L、U、P的行数/列数以及向量x、b的长度
 * @param q 矩阵Q的存放位置
 * @param ldq 矩阵Q的主维度
 * @param r 矩阵R的存放位置
 * @param ldr 矩阵R的主维度
 * @param b 向量b
 * @param incb 向量b的步长
 * @param x 向量x的存放位置
 * @param incx 向量x的步长
 */
void QRSolve(const int& n, const double* q, const int& ldq, const double* r, const int& ldr,
             const double* b, const int& incb, double* x, const int& incx);

#endif //MATRIX_V2_INTERFACE_H
