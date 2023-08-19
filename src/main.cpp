#include <cstdio>
#include "factor/interface.h"
#include "factor/utils.h"

// 要求完成课堂上讲的关于矩阵分解的LU、QR（Gram-Schmidt）、Orthogonal Reduction (Householder reduction 和Givens reduction)和 URV程序实现，要求如下：
//
//        1、一个综合程序，根据选择参数的不同，实现不同的矩阵分解；在此基础上，实现Ax=b方程组的求解，以及计算A的行列式；
//
//         2、可以用matlab、Python等编写程序，需附上简单的程序说明，比如参数代表什么意思，输入什么，输出什么等等，附上相应的例子；
//
//         3、一定是可执行文件，例如 .m文件等,不能是word或者txt文档。附上源代码，不能为直接调用matlab等函数库;

/**
 * 用于指示函数调用时读入数据的类型。
 */
enum DataSourceType {
    FileType,
    InputType
};

/**
 * LU分解及其求解测试，支持对满秩的，运算过程中主元不出现零的n*n矩阵的LU分解，以及利用LU分解结果求解。
 * 输入数据应包括一个n*n的矩阵A和一个长度为n的向量b。
 * @param type 读取数据的方式。
 */
void LUTest(DataSourceType type) {
    int t, n, &lda = n, &ldl = n, &ldu = n, incx = 1;
    double *a, *l, *u, *x;
    FILE *f;
    if (type == FileType) {
        f = fopen("test/lu_test.txt", "r");
        fscanf(f, "%d", &t);
    }
    else {
        printf("请输入一个n*n的矩阵和一个长度为n的向量，先输入维度，再输入矩阵和向量。示例：\n");
        printf("3\n"
               "2 2 2\n"
               "4 7 7\n"
               "6 18 22\n"
               "12 24 12\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("测试用例%d:\n", i);
            fscanf(f, "%d", &n);
            a = (double *) safe_malloc(lda * n * sizeof(double));
            l = (double *) safe_malloc(ldl * n * sizeof(double));
            u = (double *) safe_malloc(ldu * n * sizeof(double));
            x = (double *) safe_malloc(incx * n * sizeof(double));
            read_matrix(n, n, a, lda, f);
            read_vector(n, x, incx, f);
        }
        else {
            scanf("%d", &n);
            a = DMALLOC(lda * n);
            x = DMALLOC(incx * n);
            input_matrix(n, n, a, lda);
            input_vector(n, x, incx);
            l = DMALLOC(ldl * n);
            u = DMALLOC(ldu * n);
        }
        printf("测试矩阵A：\n");
        print_matrix(n, n, a, lda);
        printf("测试向量b：\n");
        print_vector(n, x, incx);
        double det = 0;
        EXEC_STATE state = LU(n, a, lda, l, ldl, u, ldu, det);
        switch (state) {
            case ExecSucc:
                printf("L:\n");
                print_matrix(n, n, l, ldl);
                printf("U:\n");
                print_matrix(n, n, u, ldu);
                printf("利用LU分解，可求得行列式值：%.2lf\n", det);
                printf("利用LU分解，对Ax=b <=> L(Ux)=b进行求解，求得x：\n");
                LUSolve(n, l, ldl, u, ldu, x, incx);
                print_vector(n, x, incx);
                break;
            case ExecNoFullRank:
                printf("该矩阵不满秩，无法进行LU分解！\n");
                break;
            case ExecLUZeroPivot:
                printf("主元在运算过程中出现0，无法执行LU分解！\n");
        }
        /*        if () {


                        printf("利用LU分解，对Ax=b <=> L(Ux)=b进行求解。\n");
                printf("对Ly=b求解，求得y:\n");
                dtrsv(BlasLower, n, l, ldl, x, incx);
                print_vector(n, x, incx);
                printf("对Ux=y求解，求得x:\n");
                dtrsv(BlasUpper, n, u, ldu, x, incx);
                print_vector(n, x, incx);
                printf("A的行列式的值：%.2lf\n", ddiagdot(n, u, ldu));

            printf("求解完成。\n");
        } else {
            printf("主元在运算过程中出现0，无法执行LU分解！\n");
        }    */
        safe_free(a);
        safe_free(l);
        safe_free(u);
        safe_free(x);
    }

    if (type == FileType)
        fclose(f);
}

/**
 * LUP分解及其求解测试，支持对满秩的的n*n矩阵的LU分解，以及利用LUP分解结果求解。
 * 输入数据应包括一个n*n的矩阵A和一个长度为n的向量b。
 * @param type 读取数据的方式。
 */
void LUPTest(DataSourceType type) {
    FILE *f;
    int t, n, &lda = n, &ldl = n, &ldu = n, &ldp = n, incx = 1, &incb = incx;
    double *a, *l, *u, *p, *b, *x;
    if (type == FileType) {
        f = fopen("test/lup_test.txt", "r");
        fscanf(f, "%d", &t);
    } else {
        printf("请输入一个n*n的矩阵和一个长度为n的向量，先输入维度，再输入矩阵和向量。示例：\n");
        printf("3\n"
               "2 2 2\n"
               "4 7 7\n"
               "6 18 22\n"
               "12 24 12\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("测试用例%d:\n", i);
            fscanf(f, "%d", &n);
            a = (double *) safe_malloc(lda * n * sizeof(double));
            l = (double *) safe_malloc(ldl * n * sizeof(double));
            u = (double *) safe_malloc(ldu * n * sizeof(double));
            p = DMALLOC(ldp * n);
            b = DMALLOC(incb * n);
            x = (double *) safe_malloc(incx * n * sizeof(double));
            read_matrix(n, n, a, lda, f);
            read_vector(n, b, incb, f);
        }
        else {
            scanf("%d", &n);
            a = DMALLOC(lda * n);
            b = DMALLOC(incb * n);
            input_matrix(n, n, a, lda);
            input_vector(n, b, incb);
            l = DMALLOC(ldl * n);
            u = DMALLOC(ldu * n);
            p = DMALLOC(ldp * n);
            x = DMALLOC(incx * n);
        }
        printf("测试矩阵A：\n");
        print_matrix(n, n, a, lda);
        printf("测试向量b：\n");
        print_vector(n, b, incb);
        double det = 0;
        EXEC_STATE state = LUP(n, a, lda, l, ldl, u, ldu, p, ldp, det);
        switch (state) {
            case ExecSucc:
                printf("L:\n");
                print_matrix(n, n, l, ldl);
                printf("U:\n");
                print_matrix(n, n, u, ldu);
                printf("P:\n");
                print_matrix(n, n, p, ldp);
                printf("利用LUP分解，可求得行列式值：%.2lf\n", det);
                printf("利用LUP分解，对Ax=b <=> PL(Ux)=b进行求解，求得x：\n");
                LUPSolve(n, l, ldl, u, ldu, p, ldp, b, incb, x, incx);
                print_vector(n, x, incx);
                break;
            case ExecNoFullRank:
                printf("该矩阵不满秩，无法进行LUP分解！\n");
        }
        safe_free(a);
        safe_free(l);
        safe_free(u);
        safe_free(p);
        safe_free(b);
        safe_free(x);
    }
    if (type == FileType)
        fclose(f);
}

/**
 * QR分解测试，包括3种QR分解的方式。支持对满秩的的n*n矩阵的QR分解，以及利用QR分解结果求解。
 * @param type QR分解的类型，支持施密特正交化、Householder正交化和Givens正交化。
 * @param dataType 读取数据的方式。
 */
void QRTest(QRType type, DataSourceType dataType) {
    FILE *f = nullptr;
    int t, n, &m = n;
    if (dataType == FileType) {
        switch (type) {
            case SchimidtType:
                f = fopen("test/gs_qr_test.txt", "r");
                break;
            case HouseholderType:
                f = fopen("test/hh_qr_test.txt", "r");
                break;
            case GivensType:
                f = fopen("test/gv_qr_test.txt", "r");
        }
        fscanf(f, "%d", &t);
    }
    else {
        printf("请输入一个n*n的矩阵和一个长度为n的向量，先输入维度，再输入矩阵和向量。示例：\n");
        printf("3\n"
               "2 2 2\n"
               "4 7 7\n"
               "6 18 22\n"
               "12 24 12\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        int &lda = m, &ldq = m, ldr, incx = 1, &incb = incx;
        /// Sch-QR: A_{m,n}=Q_{m,n}R_{n,n}
        /// Hh/Gv-QR: A_{m,n}=Q_{m,m}R_{m,n}
        double *a, *q, *r, *x, *b;
        if (dataType == FileType) {
            printf("测试用例%d:\n", i);
            fscanf(f, "%d", &n);
            switch (type) {
                case SchimidtType:
                    q = DMALLOC(ldq * n);
                    ldr = n;
                    break;
                case HouseholderType:
                case GivensType:
                    q = DMALLOC(ldq * m);
                    ldr = m;
            }
            a = (double *) safe_malloc(lda * n * sizeof(double));
            r = (double *) safe_malloc(ldr * n * sizeof(double));
            x = (double *) safe_malloc(incx * n * sizeof(double));
            b = (double *) safe_malloc(incb * m * sizeof(double));
            read_matrix(m, n, a, lda, f);
            read_vector(m, b, incb, f);
        }
        else {
            scanf("%d", &n);
            switch (type) {
                case SchimidtType:
                    q = DMALLOC(ldq * n);
                    ldr = n;
                    break;
                case HouseholderType:
                case GivensType:
                    q = DMALLOC(ldq * m);
                    ldr = m;
            }
            a = DMALLOC(lda  * n);
            b = DMALLOC(incb * n);
            input_matrix(n, n, a, lda);
            input_vector(n, b, incb);
            r = DMALLOC(ldr * n);
            x = DMALLOC(incx * n);
        }

        printf("测试矩阵A：\n");
        print_matrix(m, n, a, lda);
        printf("测试向量b：\n");
        print_vector(m, b, incb);
        EXEC_STATE state = QR(n, a, lda, q, ldq, r, ldr, type);
        if (state == ExecNoFullRank) {
            printf("该矩阵不满秩，无法进行QR分解！\n");
        }
        else if (state == ExecSucc) {
            printf("Q:\n");
            print_matrix(n, n, q, ldq);
            printf("R:\n");
            print_matrix(n, n, r, ldr);
            printf("利用QR分解，对Ax=b <=> Rx=Q^{T}b进行求解，解得x：\n");
            QRSolve(n, q, ldq, r, ldr, b, incb, x, incx);
            print_vector(n, x, incx);
        }


/*
        if (m != n) {
            printf("暂不支持A不为方阵时的非施密特正交化下的线性方程组求解！\n");
            return;
        }

        printf("利用QR分解，对Ax=b <=> Rx=Q^{T}b进行求解。\n");
        /// 注意：x实际上是b
        printf("对Q^{T}b=y求解，求得y：\n");
        dgmv(BlasTrans, m, n, q, ldq, x, incx, y, incy);
        print_vector(n, y, incy);
        /// 注意：y实际上是b->x
        printf("对Rx=y求解，求得x：\n");
        dtrsv(BlasUpper, n, r, ldr, y, incy);
        print_vector(n, y, incy);

        printf("A的行列式的值：%.2lf\n", Det(n, a, lda));
*/
        safe_free(a);
        safe_free(q);
        safe_free(r);
        safe_free(x);
        safe_free(b);
    }
    if (dataType == FileType)
        fclose(f);
}

/**
 * URV分解测试，支持对m*n矩阵的URV分解。输入应包括一个m*n的矩阵A。
 * @param type 输入数据的方式。
 */
void URVTest(DataSourceType type) {
    FILE *f;
    int t, m, n, &lda = m, &ldu =  m, &ldr = m, &ldv = n;
    double *a, *u, *r, *v;
    if (type == FileType) {
        f = fopen("test/urv_test.txt", "r");
        fscanf(f, "%d", &t);
    }
    else {
        printf("请输入一个m*n的矩阵，先输入维度，再输入矩阵。示例：\n");
        printf("3 4\n"
               "-4 -2 -4 -2\n"
               "2 -2 2 1\n"
               "-4 1 -4 -2\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("测试用例%d:\n", i);
            fscanf(f, "%d %d", &m, &n);
            a = DMALLOC(lda * n);

            read_matrix(m, n, a, lda, f);
        }
        else {
            scanf("%d %d", &m, &n);
            a = DMALLOC(lda * n);
            input_matrix(m, n, a, lda);
        }
        u = DMALLOC(ldu * m);
        r = DMALLOC(ldr * n);
        v = DMALLOC(ldv * n);

        printf("A:\n");
        print_matrix(m, n, a, lda);
        URV(m, n, a, lda, u, ldu, r, ldr, v, ldv);
        printf("U:\n");
        print_matrix(m, m, u, ldu);
        printf("R:\n");
        print_matrix(m, n, r, ldr);
        printf("V:\n");
        print_matrix(n, n, v, ldv);
/*
        int &ldb = m, &ldc = n;
        auto *b = DMALLOC(ldb * n);
        auto *c = DMALLOC(ldc * n);
        printf("验算：\n");
        dgmm(BlasNoTrans, m, n, m, u, ldu, a, lda, b, ldb);
        dgmm(BlasTrans, m, n, n, b, ldb, v, ldv, c, ldc);
        print_matrix(m, n, c, ldc);
        printf("验算2：\n");
        print_matrix(m, n, c, ldc);
        dtransin(n, v, ldv);
        dgmm(BlasNoTrans, m, n, n, b, ldb, v, ldv, c, ldc);
*/
        safe_free(a);
        safe_free(u);
        safe_free(r);
        safe_free(v);
        // safe_free(b);
        // safe_free(c);
    }
    if (type == FileType)
        fclose(f);
}

/**
 * 高斯-约旦分解测试，输入数据应包括一个m*n的矩阵A。
 * @param type 输入数据的方式。
 */
void GJTest(DataSourceType type) {
    FILE *f;
    int t, m, n, &lda = m, &ldb = m;
    double* a, *b;
    if (type == FileType) {
        f = fopen("test/gj_test.txt", "r");
        fscanf(f, "%d", &t);
    }
    else {
        printf("请输入一个m*n的矩阵，先输入维度，再输入矩阵。示例：\n");
        printf("3 4\n"
               "-4 -2 -4 -2\n"
               "2 -2 2 1\n"
               "-4 1 -4 -2\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("测试用例%d:\n", i);
            fscanf(f, "%d %d", &m, &n);
            a = DMALLOC(lda * n);
            b = DMALLOC(ldb * n);
            read_matrix(m, n, a, lda, f);
        }
        else {
            scanf("%d %d", &m, &n);
            a = DMALLOC(lda * n);
            input_matrix(m, n, a, lda);
            b = DMALLOC(ldb * n);
        }
        GaussJordan(m, n, a, lda, b, ldb);
        printf("A:\n");
        print_matrix(m, n, a, lda);
        printf("Modified A:\n");
        print_matrix(m, n, b, ldb);

        safe_free(a);
        safe_free(b);
    }

    if (type == FileType)
        fclose(f);
}

/**
 * 行列式求解。输入数据包括一个n*n的矩阵A。
 * @param type 输入数据的方式。
 */
void DetTest(DataSourceType type) {
    FILE *f;
    int t, n, &lda = n;
    double *a;
    if (type == FileType) {
        f = fopen("test/det_test.txt", "r");
        fscanf(f, "%d", &t);
    }
    else {
        printf("请输入一个n*n的矩阵，先输入维度，再输入矩阵。示例：\n");
        printf("3\n"
               "2 2 2\n"
               "4 7 7\n"
               "6 18 22\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("测试用例%d:\n", i);
            fscanf(f, "%d", &n);
            a = DMALLOC(lda * n);
            read_matrix(n, n, a, lda, f);
            printf("A:\n");
            print_matrix(n, n, a, lda);
        }
        else {
            scanf("%d", &n);
            a = DMALLOC(lda * n);
            input_matrix(n, n, a, lda);
        }
        printf("A的行列式：%.2lf\n", Det(n, a, lda));

        safe_free(a);
    }
    if (type == FileType)
        fclose(f);
}

/**
 * 秩求解测试，输入数据包括一个m*n的矩阵A。
 * @param type 输入数据的方式。
 */
void rankTest(DataSourceType type) {
    FILE *f;
    int t, m, n, &lda = m;
    double *a;
    if (type == FileType) {
        f = fopen("test/rank_test.txt", "r");
        fscanf(f, "%d", &t);
    }
    else {
        printf("请输入一个m*n的矩阵，先输入维度，再输入矩阵。示例：\n");
        printf("3 4\n"
               "-4 -2 -4 -2\n"
               "2 -2 2 1\n"
               "-4 1 -4 -2\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("测试用例%d:\n", i);
            fscanf(f, "%d %d", &m, &n);
            a = DMALLOC(lda * n);
            read_matrix(m, n, a, lda, f);
            printf("A:\n");
            print_matrix(m, n, a, lda);
        }
        else {
            scanf("%d %d", &m, &n);
            a = DMALLOC(lda * n);
            input_matrix(m, n, a, lda);
        }
        printf("秩：%d\n", Rank(m, n, a, lda));

        safe_free(a);
    }

    if (type == FileType)
        fclose(f);
}

int main() {
    int choice;
    // bool out = false;
    DataSourceType type = FileType;
    printf("------矩阵分解器------\n");
    printf("现支持对n*n矩阵的LU分解、LUP分解和3种QR分解及其求解，以及行列式和秩的计算；支持对m*n矩阵的URV分解和Gauss-Jordan化简。\n");
    printf("默认输入方式为从文件中读取。\n");
    while (true) {
        printf("1. LU分解 2. LUP分解 3. QR分解（施密特正交化） 4. QR分解（Householder） 5. QR分解（Givens）\n");
        printf("6. URV分解 7. Gauss-Jordan化简 8. 计算行列式 9. 计算秩 10. 修改输入方式 0.关闭程序\n");
        printf("请输入您选择的序号：");
        scanf("%d", &choice);
        switch (choice) {
            case 1:
                LUTest(type);
                break;
            case 2:
                LUPTest(type);
                break;
            case 3:
                // GS_QRTest();
                QRTest(SchimidtType, type);
                break;
            case 4:
                // Hh_QRTest();
                QRTest(HouseholderType, type);
                break;
            case 5:
                QRTest(GivensType, type);
                break;
            case 6:
                URVTest(type);
                break;
            case 7:
                GJTest(type);
                break;
            case 8:
                DetTest(type);
                break;
            case 9:
                rankTest(type);
                break;
            case 10:
                if (type == FileType) {
                    type = InputType;
                    printf("输入模式已修改为从控制台输入！\n");
                }
                else {
                    type = FileType;
                    printf("输入模式已修改为从文件输入！\n");
                }
                break;
            case 0:
                return 0;
            //case 11:
           //    dgmv_test();
           //     break;
            default:
                printf("输入序号有误，请重新输入！\n");
        }
    }
}

/*
void GS_QRTest() {
    FILE *f = fopen("qr_test.txt", "r");
    int t, m, n, &lda = m, &ldq = m, &ldr = m, incx = 1, &incy = incx;
    fscanf(f, "%d", &t);
    while (t--) {
        fscanf(f, "%d %d", &m, &n);
        auto *a = (double *) safe_malloc(lda * n * sizeof(double));
        auto *q = (double *) safe_malloc(ldq * n * sizeof(double));
        auto *r = (double *) safe_malloc(ldr * n * sizeof(double));
        auto *x = (double *) safe_malloc(incx * n * sizeof(double));
        auto *y = (double *) safe_malloc(incy * m * sizeof(double));
        printf("测试矩阵A：\n");
        read_matrix(m, n, a, lda, f);
        print_matrix(m, n, a, lda);
        printf("测试向量b：\n");
        read_vector(m, x, incx, f);
        print_vector(m, x, incx);
        GS_QR(m, n, a, lda, q, ldq, r, ldr);
        printf("Q:\n");
        print_matrix(n, n, q, ldq);
        printf("R:\n");
        print_matrix(n, n, r, ldr);
        printf("利用QR分解，对Ax=b <=> Rx=Q^{T}b进行求解。\n");
        /// 注意：x实际上是b
        printf("对Q^{T}b=y求解，求得y：\n");
        dgmv(BlasTrans, m, n, q, ldq, x, incx, y, incy);
        print_vector(n, y, incy);
        /// 注意：y实际上是b->x
        printf("对Rx=y求解，求得x：\n");
        dtrsv(BlasUpper, n, r, ldr, y, incy);
        print_vector(n, y, incy);
        safe_free(a);
        safe_free(q);
        safe_free(r);
        safe_free(x);
        safe_free(y);
    }
    fclose(f);
}

void Hh_QRTest() {
    FILE *f = fopen("hh_qr_test.txt", "r");
    int t, m, n, &lda = m, &ldq = m, &ldr = m, incx = 1, &incy = incx;
    fscanf(f, "%d", &t);
    while (t--) {
        fscanf(f, "%d %d", &m, &n);
        auto *a = DMALLOC(lda * n);
        auto *q = DMALLOC(ldq * m);
        auto *r = DMALLOC(ldr * n);
        auto *x = DMALLOC(n * incx);
        auto *y = DMALLOC(m * incy);
        printf("测试矩阵A：\n");
        read_matrix(m, n, a, lda, f);
        print_matrix(m, n, a, lda);
        printf("测试向量b：\n");
        read_vector(m, x, incx, f);
        print_vector(m, x, incx);
        Hh_QR(m, n, a, lda, q, ldq, r, ldr);
        printf("Q:\n");
        print_matrix(n, n, q, ldq);
        printf("R:\n");
        print_matrix(n, n, r, ldr);
        printf("利用QR分解，对Ax=b <=> Rx=Q^{T}b进行求解。\n");

        if (m != n) {
            printf("暂不支持A不为方阵时的线性方程组求解！\n");
            return;
        }

        /// 注意：x实际上是b
        printf("对Q^{T}b=y求解，求得y：\n");
        dgmv(BlasTrans, m, n, q, ldq, x, incx, y, incy);
        print_vector(n, y, incy);
        /// 注意：y实际上是b->x
        printf("对Rx=y求解，求得x：\n");
        dtrsv(BlasUpper, n, r, ldr, y, incy);
        print_vector(n, y, incy);
        safe_free(a);
        safe_free(q);
        safe_free(r);
        safe_free(x);
        safe_free(y);
    }

    fclose(f);
}
*/

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