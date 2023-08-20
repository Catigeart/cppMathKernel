# cppMathKernel
该项目为本人手搓过的各类数学内核函数的整合，包括矩阵分解、BLAS等。目前已整合矩阵分解部分函数和个别matlab函数。项目中的函数大部分都是平台通用的c/cpp未优化函数，主要为个人参考使用和积累学习，也供大家参考。

## 1 矩阵分解器

**警告：本项目为课堂作业/练习项目，不能保证正确性！**

该项目为参考BLAS（Basic Linear Algebra Subprograms）接口和内存设计风格的开发的C++项目，所有矩阵、向量运算全部自底层从零实现。现已实现的功能有（仅支持double类型）：

- 对n*n矩阵的LU分解及其行列式计算、线性方程组求解
- 对n*n矩阵的LUP分解及其行列式计算、线性方程组求解
- 对n*n矩阵的QR分解（施密特正交化、Householder、Givens）及其线性方程组求解
- 对m*n矩阵的URV分解
- 对m*n矩阵的高斯-约旦化简
- 对m*n矩阵的秩的计算
- 对n*n矩阵的行列式的计算

### 项目结构

该项目自底向上可划分为四个层次：

#### blas.cpp/blas.h（底层函数层）

即基础线性代数运算层。参考BLAS（Basic Linear Algebra Subprograms）的接口风格和内存设计，对基础的矩阵、向量运算予以实现，包括：

- 求向量的最大值
  - int damax(const int& n, const double* x, const int& incx)
- 求向量的最大绝对值
  - int daabsmax(const int& n, const double* x, const int& incx)
- 交换两向量
  - void dswap(const int& n, double* x, const int& incx, double* y, const int& incy)
- 复制向量
  - void dcopy(const int& n, const double* x, const int& incx, double* y, const int& incy)
- 向量内积
  - double ddot(const int& n, const double* x, const int& incx, const double* y, const int& incy)
- 求向量的模
  - double dnormv(const int& n, const double* x, const int& incx)
- 矩阵转置
  - void dtrans(const int& m, const int& n, double* a_dst, const int& lda, const double* b_src, const int& ldb)
- 矩阵原地转置
  - void dtransin(const int& n, double *a, const int& lda)
- 复制矩阵
  - void dmatcpy(const int& m, const int& n, double* a_dst, const int& lda, const double* b_src, const int& ldb)
- 0矩阵
  - void dempty(const int& m, const int& n, double* a, const int& lda)
- 矩阵-向量乘
  - void dgmv(BLAS_TRANS trans, const int& m, const int& n,
        const double* a, const int& lda, const double* x, const int& incx, double* y, const int& incy)
- 矩阵-矩阵乘
  - void dgmm(BLAS_TRANS trans, const int& m, const int& n, const int& k,
        const double* a, const int& lda, const double* b, const int& ldb, double* c, const int& ldc)
- 三角矩阵的线性方程组求解
  - void dtrsv(BLAS_UPLO uplo, const int& n, const double* a, const int& lda, double* x, const int& incx)
- 矩阵的迹
  - double ddiagdot(const int& n, const double* a, const int& lda)

#### linpack.cpp/linpack.h（高层函数层）

借用线性系统软件包(Linear system package) 的缩写，主要是上层函数的实现，包括：

- 求矩阵的秩
  - int Rank(const int &m, const int &n, const double *a, const int &lda)
- 求矩阵的行列式
  - double Det(const int& n, const double* a, const int& lda)
- 高斯-约旦化简
  - void GaussJordan(const int& m, const int& n, const double* a, const int&lda, double* b, const int& ldb)
- LU分解
  - EXEC_STATE LU(const int& n, const double* a, const int& lda, double* l, const int& ldl, double* u, const int&ldu, double& det)
- LUP分解
  - EXEC_STATE LUP(const int& n, const double* a, const int& lda,
        double* l, const int& ldl, double* u, const int&ldu, double* p, const int& ldp, double& det)
- 三种QR分解
  - EXEC_STATE QR(const int& n, const double* a, const int& lda,
        double* q, const int& ldq, double* r, const int&ldr, QRType type)
- URV分解
  - void URV(const int& m, const int& n, const double* a, const int& lda,
        double* u, const int& ldu, double* r, const int& ldr, double* v, const int& ldv)
  - **该函数计算可能存在错误。**
- 基于LU分解的线性方程组求解
  - void LUSolve(const int& n, const double* l, const int& ldl, const double* u, const int&ldu, double* x, const int& incx)
- 基于LUP分解的线性方程组求解
  - void LUPSolve(const int& n, const double* l, const int& ldl, const double* u, const int&ldu,
        const double* p, const int& ldp, double* b, const int& incb, double* x, const int& incx)
- 基于QR分解的线性方程组求解
  - void QRSolve(const int& n, const double* q, const int& ldq, const double* r, const int& ldr,
        const double* b, const int& incb, double* x, const int& incx)

#### interface.cpp/interface.h*（接口层）

对上层函数的封装，为应用层调用提供参数友好的接口。

#### main.cpp（应用层）

面向用户的函数调用和测试，拥有完善的提示文字，交互友好，支持文件输入和控制台输入两种方式。

### 项目参数与内存设计解释

参考高性能计算的开发习惯，本项目采用列主序的方式存储，即元素按列连续存储，在C/C++中以一维动态数组实现。为方便理解，以下对向量和矩阵的参数作简要解释，各函数的具体参数解释可参考源代码注释。

#### 矩阵与主维度（Leading dimension）

在本项目中，一个矩阵A包括以下参数：

- 维度大小，m\*n矩阵有m、n两个参数，n\*n矩阵有n一个参数
- 矩阵指针a
- 主维度lda。在列主序中，lda>=m。

列主序下的 Leading dimension ，即矩阵从某个元素移动到相邻列的另一个元素所需要移动的指针的距离。因此lda>=m。当lda>m时，说明当前的矩阵在内存中实际上是更大矩阵的子矩阵。

对矩阵A的第i行，第j列元素，其坐标为i+j*lda。在该项目中，绝大部分情况下均取ld?=m，部分涉及子矩阵计算情况时，存在ld?>m的情况（如Givens分解）。

#### 向量

向量x存储的逻辑相对简单，包括以下参数：

- 向量长度n
- 向量指针x
- 向量步长incx，即访问下一个向量元素间隔的距离。

该项目的实现中未涉及inc? != 1的情况。

## 2 Matlab函数重写
项目整合了一部分我在一个改写matlab代码为c语言代码工作中所发明的简陋轮子。
