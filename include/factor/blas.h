//
// Created by Catigeart on 2022/11/23.
//

#ifndef MATRIX_V2_BLAS_H
#define MATRIX_V2_BLAS_H

/**
 * 指示三角矩阵是上三角矩阵还是下三角矩阵。
 */
enum BLAS_UPLO {
    BlasUpper,
    BlasLower
};

/**
 * 指示该矩阵按原矩阵计算还是按转置矩阵计算。
 */
enum BLAS_TRANS {
    BlasNoTrans,
    BlasTrans,
};

/**
 * 求一个向量中的最大值。
 * @param n 向量的长度
 * @param x 向量x
 * @param incx x的步长
 * @return 最大值的位置索引
 */
int damax(const int& n, const double* x, const int& incx);

/**
 * 找一个向量中的最大绝对值。
 * @param n 向量的长度
 * @param x 向量x
 * @param incx x的步长
 * @return 最大绝对值的位置索引
 */
int daabsmax(const int& n, const double* x, const int& incx);

/**
 * 交换两个向量的值。
 * @param n 向量的长度
 * @param x 向量x
 * @param incx x的步长
 * @param y 向量y
 * @param incy y的步长
 */
void dswap(const int& n, double* x, const int& incx, double* y, const int& incy);

/**
 *  y := x。
 * @param n 向量的长度
 * @param x 向量x
 * @param incx x的步长
 * @param y 向量y
 * @param incy y的步长
 */
void dcopy(const int& n, const double* x, const int& incx, double* y, const int& incy);

/**
 * 计算两向量的点乘（内积）。
 * @param n 向量的长度
 * @param x 向量x
 * @param incx x的步长
 * @param y 向量y
 * @param incy y的步长
 * @return x^{T}y
 */
double ddot(const int& n, const double* x, const int& incx, const double* y, const int& incy);

/**
 * 计算向量的模长。
 * @param n 向量的长度
 * @param x 向量x
 * @param incx x的步长
 * @return 向量的模长
 */
double dnormv(const int& n, const double* x, const int& incx);

/**
 * 计算矩阵的转置。
 * @param m B的行数、A的列数
 * @param n B的列数、A的行数
 * @param a_dst 转置后矩阵A的存放位置
 * @param lda A的主维度
 * @param b_src 原矩阵B
 * @param ldb B的主维度
 */
void dtrans(const int& m, const int& n, double* a_dst, const int& lda, const double* b_src, const int& ldb);

/**
 * n*n矩阵的原地转置。
 * @param n A的行数/列数
 * @param a 矩阵A，计算后为A^T
 * @param lda A的主维度
 */
void dtransin(const int& n, double *a, const int& lda);

/**
 * A := B。
 * @param m B的行数、A的列数
 * @param n B的列数、A的行数
 * @param a_dst 矩阵A的存放位置
 * @param lda A的主维度
 * @param b_src 原矩阵B
 * @param ldb B的主维度
 */
void dmatcpy(const int& m, const int& n, double* a_dst, const int& lda, const double* b_src, const int& ldb);

/**
 * 将矩阵置零。
 * @param m A的行数
 * @param n A的列数
 * @param a 矩阵A
 * @param lda A的主维度
 */
void dempty(const int& m, const int& n, double* a, const int& lda);

/**
 * y := op(A)x
 * @param trans 若trans=BlasNoTrans，则op(A)=A；若为BlasTrans，则op(A)=A^T
 * @param m A的行数
 * @param n A的列数
 * @param a 矩阵A
 * @param lda A的主维度
 * @param x 向量x
 * @param incx x的步长
 * @param y 向量y的存放位置
 * @param incy y的步长
 */
void dgmv(BLAS_TRANS trans, const int& m, const int& n,
          const double* a, const int& lda, const double* x, const int& incx, double* y, const int& incy);

/**
 * C := op(A)B
 * @param trans 若trans=BlasNoTrans，则op(A)=A；若为BlasTrans，则op(A)=A^T
 * @param m A、C的行数
 * @param n B、C的列数
 * @param k A的列数、B的行数
 * @param a 矩阵A
 * @param lda A的主维度
 * @param b 矩阵B
 * @param ldb B的主维度
 * @param c 矩阵C
 * @param ldc C的主维度
 */
void dgmm(BLAS_TRANS trans, const int& m, const int& n, const int& k,
          const double* a, const int& lda, const double* b, const int& ldb, double* c, const int& ldc);

/**
 * 对Ax=b进行求解，其中A为三角矩阵。
 * @param uplo 指示A的类型，若为BlasUpper，则A为上三角矩阵；若为BlasLower，则A为下三角矩阵
 * @param n A的行数/列数、x的长度
 * @param a 矩阵A
 * @param lda A的主维度
 * @param x 计算前为向量b，计算后存储求解得到的向量x
 * @param incx x的步长
 */
void dtrsv(BLAS_UPLO uplo, const int& n, const double* a, const int& lda, double* x, const int& incx);

/**
 * 求解矩阵的迹。
 * @param n A的行数/列数
 * @param a 矩阵A
 * @param lda 主维度A
 * @return 矩阵的迹
 */
double ddiagdot(const int& n, const double* a, const int& lda);

#endif //MATRIX_V2_BLAS_H
