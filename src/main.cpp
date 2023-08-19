#include <cstdio>
#include "factor/interface.h"
#include "factor/utils.h"

// Ҫ����ɿ����Ͻ��Ĺ��ھ���ֽ��LU��QR��Gram-Schmidt����Orthogonal Reduction (Householder reduction ��Givens reduction)�� URV����ʵ�֣�Ҫ�����£�
//
//        1��һ���ۺϳ��򣬸���ѡ������Ĳ�ͬ��ʵ�ֲ�ͬ�ľ���ֽ⣻�ڴ˻����ϣ�ʵ��Ax=b���������⣬�Լ�����A������ʽ��
//
//         2��������matlab��Python�ȱ�д�����踽�ϼ򵥵ĳ���˵���������������ʲô��˼������ʲô�����ʲô�ȵȣ�������Ӧ�����ӣ�
//
//         3��һ���ǿ�ִ���ļ������� .m�ļ���,������word����txt�ĵ�������Դ���룬����Ϊֱ�ӵ���matlab�Ⱥ�����;

/**
 * ����ָʾ��������ʱ�������ݵ����͡�
 */
enum DataSourceType {
    FileType,
    InputType
};

/**
 * LU�ֽ⼰�������ԣ�֧�ֶ����ȵģ������������Ԫ���������n*n�����LU�ֽ⣬�Լ�����LU�ֽ�����⡣
 * ��������Ӧ����һ��n*n�ľ���A��һ������Ϊn������b��
 * @param type ��ȡ���ݵķ�ʽ��
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
        printf("������һ��n*n�ľ����һ������Ϊn��������������ά�ȣ�����������������ʾ����\n");
        printf("3\n"
               "2 2 2\n"
               "4 7 7\n"
               "6 18 22\n"
               "12 24 12\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("��������%d:\n", i);
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
        printf("���Ծ���A��\n");
        print_matrix(n, n, a, lda);
        printf("��������b��\n");
        print_vector(n, x, incx);
        double det = 0;
        EXEC_STATE state = LU(n, a, lda, l, ldl, u, ldu, det);
        switch (state) {
            case ExecSucc:
                printf("L:\n");
                print_matrix(n, n, l, ldl);
                printf("U:\n");
                print_matrix(n, n, u, ldu);
                printf("����LU�ֽ⣬���������ʽֵ��%.2lf\n", det);
                printf("����LU�ֽ⣬��Ax=b <=> L(Ux)=b������⣬���x��\n");
                LUSolve(n, l, ldl, u, ldu, x, incx);
                print_vector(n, x, incx);
                break;
            case ExecNoFullRank:
                printf("�þ������ȣ��޷�����LU�ֽ⣡\n");
                break;
            case ExecLUZeroPivot:
                printf("��Ԫ����������г���0���޷�ִ��LU�ֽ⣡\n");
        }
        /*        if () {


                        printf("����LU�ֽ⣬��Ax=b <=> L(Ux)=b������⡣\n");
                printf("��Ly=b��⣬���y:\n");
                dtrsv(BlasLower, n, l, ldl, x, incx);
                print_vector(n, x, incx);
                printf("��Ux=y��⣬���x:\n");
                dtrsv(BlasUpper, n, u, ldu, x, incx);
                print_vector(n, x, incx);
                printf("A������ʽ��ֵ��%.2lf\n", ddiagdot(n, u, ldu));

            printf("�����ɡ�\n");
        } else {
            printf("��Ԫ����������г���0���޷�ִ��LU�ֽ⣡\n");
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
 * LUP�ֽ⼰�������ԣ�֧�ֶ����ȵĵ�n*n�����LU�ֽ⣬�Լ�����LUP�ֽ�����⡣
 * ��������Ӧ����һ��n*n�ľ���A��һ������Ϊn������b��
 * @param type ��ȡ���ݵķ�ʽ��
 */
void LUPTest(DataSourceType type) {
    FILE *f;
    int t, n, &lda = n, &ldl = n, &ldu = n, &ldp = n, incx = 1, &incb = incx;
    double *a, *l, *u, *p, *b, *x;
    if (type == FileType) {
        f = fopen("test/lup_test.txt", "r");
        fscanf(f, "%d", &t);
    } else {
        printf("������һ��n*n�ľ����һ������Ϊn��������������ά�ȣ�����������������ʾ����\n");
        printf("3\n"
               "2 2 2\n"
               "4 7 7\n"
               "6 18 22\n"
               "12 24 12\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("��������%d:\n", i);
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
        printf("���Ծ���A��\n");
        print_matrix(n, n, a, lda);
        printf("��������b��\n");
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
                printf("����LUP�ֽ⣬���������ʽֵ��%.2lf\n", det);
                printf("����LUP�ֽ⣬��Ax=b <=> PL(Ux)=b������⣬���x��\n");
                LUPSolve(n, l, ldl, u, ldu, p, ldp, b, incb, x, incx);
                print_vector(n, x, incx);
                break;
            case ExecNoFullRank:
                printf("�þ������ȣ��޷�����LUP�ֽ⣡\n");
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
 * QR�ֽ���ԣ�����3��QR�ֽ�ķ�ʽ��֧�ֶ����ȵĵ�n*n�����QR�ֽ⣬�Լ�����QR�ֽ�����⡣
 * @param type QR�ֽ�����ͣ�֧��ʩ������������Householder��������Givens��������
 * @param dataType ��ȡ���ݵķ�ʽ��
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
        printf("������һ��n*n�ľ����һ������Ϊn��������������ά�ȣ�����������������ʾ����\n");
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
            printf("��������%d:\n", i);
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

        printf("���Ծ���A��\n");
        print_matrix(m, n, a, lda);
        printf("��������b��\n");
        print_vector(m, b, incb);
        EXEC_STATE state = QR(n, a, lda, q, ldq, r, ldr, type);
        if (state == ExecNoFullRank) {
            printf("�þ������ȣ��޷�����QR�ֽ⣡\n");
        }
        else if (state == ExecSucc) {
            printf("Q:\n");
            print_matrix(n, n, q, ldq);
            printf("R:\n");
            print_matrix(n, n, r, ldr);
            printf("����QR�ֽ⣬��Ax=b <=> Rx=Q^{T}b������⣬���x��\n");
            QRSolve(n, q, ldq, r, ldr, b, incb, x, incx);
            print_vector(n, x, incx);
        }


/*
        if (m != n) {
            printf("�ݲ�֧��A��Ϊ����ʱ�ķ�ʩ�����������µ����Է�������⣡\n");
            return;
        }

        printf("����QR�ֽ⣬��Ax=b <=> Rx=Q^{T}b������⡣\n");
        /// ע�⣺xʵ������b
        printf("��Q^{T}b=y��⣬���y��\n");
        dgmv(BlasTrans, m, n, q, ldq, x, incx, y, incy);
        print_vector(n, y, incy);
        /// ע�⣺yʵ������b->x
        printf("��Rx=y��⣬���x��\n");
        dtrsv(BlasUpper, n, r, ldr, y, incy);
        print_vector(n, y, incy);

        printf("A������ʽ��ֵ��%.2lf\n", Det(n, a, lda));
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
 * URV�ֽ���ԣ�֧�ֶ�m*n�����URV�ֽ⡣����Ӧ����һ��m*n�ľ���A��
 * @param type �������ݵķ�ʽ��
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
        printf("������һ��m*n�ľ���������ά�ȣ����������ʾ����\n");
        printf("3 4\n"
               "-4 -2 -4 -2\n"
               "2 -2 2 1\n"
               "-4 1 -4 -2\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("��������%d:\n", i);
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
        printf("���㣺\n");
        dgmm(BlasNoTrans, m, n, m, u, ldu, a, lda, b, ldb);
        dgmm(BlasTrans, m, n, n, b, ldb, v, ldv, c, ldc);
        print_matrix(m, n, c, ldc);
        printf("����2��\n");
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
 * ��˹-Լ���ֽ���ԣ���������Ӧ����һ��m*n�ľ���A��
 * @param type �������ݵķ�ʽ��
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
        printf("������һ��m*n�ľ���������ά�ȣ����������ʾ����\n");
        printf("3 4\n"
               "-4 -2 -4 -2\n"
               "2 -2 2 1\n"
               "-4 1 -4 -2\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("��������%d:\n", i);
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
 * ����ʽ��⡣�������ݰ���һ��n*n�ľ���A��
 * @param type �������ݵķ�ʽ��
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
        printf("������һ��n*n�ľ���������ά�ȣ����������ʾ����\n");
        printf("3\n"
               "2 2 2\n"
               "4 7 7\n"
               "6 18 22\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("��������%d:\n", i);
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
        printf("A������ʽ��%.2lf\n", Det(n, a, lda));

        safe_free(a);
    }
    if (type == FileType)
        fclose(f);
}

/**
 * �������ԣ��������ݰ���һ��m*n�ľ���A��
 * @param type �������ݵķ�ʽ��
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
        printf("������һ��m*n�ľ���������ά�ȣ����������ʾ����\n");
        printf("3 4\n"
               "-4 -2 -4 -2\n"
               "2 -2 2 1\n"
               "-4 1 -4 -2\n");
        t = 1;
    }
    for (int i = 1; i <= t; ++i) {
        if (type == FileType) {
            printf("��������%d:\n", i);
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
        printf("�ȣ�%d\n", Rank(m, n, a, lda));

        safe_free(a);
    }

    if (type == FileType)
        fclose(f);
}

int main() {
    int choice;
    // bool out = false;
    DataSourceType type = FileType;
    printf("------����ֽ���------\n");
    printf("��֧�ֶ�n*n�����LU�ֽ⡢LUP�ֽ��3��QR�ֽ⼰����⣬�Լ�����ʽ���ȵļ��㣻֧�ֶ�m*n�����URV�ֽ��Gauss-Jordan����\n");
    printf("Ĭ�����뷽ʽΪ���ļ��ж�ȡ��\n");
    while (true) {
        printf("1. LU�ֽ� 2. LUP�ֽ� 3. QR�ֽ⣨ʩ������������ 4. QR�ֽ⣨Householder�� 5. QR�ֽ⣨Givens��\n");
        printf("6. URV�ֽ� 7. Gauss-Jordan���� 8. ��������ʽ 9. ������ 10. �޸����뷽ʽ 0.�رճ���\n");
        printf("��������ѡ�����ţ�");
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
                    printf("����ģʽ���޸�Ϊ�ӿ���̨���룡\n");
                }
                else {
                    type = FileType;
                    printf("����ģʽ���޸�Ϊ���ļ����룡\n");
                }
                break;
            case 0:
                return 0;
            //case 11:
           //    dgmv_test();
           //     break;
            default:
                printf("��������������������룡\n");
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
        printf("���Ծ���A��\n");
        read_matrix(m, n, a, lda, f);
        print_matrix(m, n, a, lda);
        printf("��������b��\n");
        read_vector(m, x, incx, f);
        print_vector(m, x, incx);
        GS_QR(m, n, a, lda, q, ldq, r, ldr);
        printf("Q:\n");
        print_matrix(n, n, q, ldq);
        printf("R:\n");
        print_matrix(n, n, r, ldr);
        printf("����QR�ֽ⣬��Ax=b <=> Rx=Q^{T}b������⡣\n");
        /// ע�⣺xʵ������b
        printf("��Q^{T}b=y��⣬���y��\n");
        dgmv(BlasTrans, m, n, q, ldq, x, incx, y, incy);
        print_vector(n, y, incy);
        /// ע�⣺yʵ������b->x
        printf("��Rx=y��⣬���x��\n");
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
        printf("���Ծ���A��\n");
        read_matrix(m, n, a, lda, f);
        print_matrix(m, n, a, lda);
        printf("��������b��\n");
        read_vector(m, x, incx, f);
        print_vector(m, x, incx);
        Hh_QR(m, n, a, lda, q, ldq, r, ldr);
        printf("Q:\n");
        print_matrix(n, n, q, ldq);
        printf("R:\n");
        print_matrix(n, n, r, ldr);
        printf("����QR�ֽ⣬��Ax=b <=> Rx=Q^{T}b������⡣\n");

        if (m != n) {
            printf("�ݲ�֧��A��Ϊ����ʱ�����Է�������⣡\n");
            return;
        }

        /// ע�⣺xʵ������b
        printf("��Q^{T}b=y��⣬���y��\n");
        dgmv(BlasTrans, m, n, q, ldq, x, incx, y, incy);
        print_vector(n, y, incy);
        /// ע�⣺yʵ������b->x
        printf("��Rx=y��⣬���x��\n");
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

    for (int i = 0, j = 0; i < m && j < n; ++i) { // �н����ε�����¿���һ���ƶ�����У����j��Ҫ�ֶ���
        int idx;

        while (j < n) {
            idx = daabsmax(m - i, b + i + j * ldb, 1); // ����һ�е��������Ԫ
            // û���ҵ����ҵ�����СԪ�����㣬˵��Ϊ�ǻ�����
            if (idx == -1)  {
                ++j;
                continue;
            }
            idx += i;
            if (fabs(B(idx, j)) < 1e-6) {
                ++j; // ����һ��
                continue;
            }
            break;
        }
        if (idx != i) {
            dswap(n, b + i, ldb, b + idx, ldb); // ������
        }
        ++rank;



        for (int ii = 0; ii < m; ++ii) { // Jordan,ÿһ�ж�Ҫ����Ԫλ������
            [[unlikely]] if (i == ii)  { // ��ǰ������������ٰ���Ԫ�����1
                continue;
            }

            double coef = B(ii, j) / B(i, j); // ��ñ�ֵ
            B(ii, j) = 0; // ֱ�Ӹ���ǰλ�ø�0

            for (int jj = j + 1; jj < n; ++jj) { // ����һ��λ�ÿ�ʼ�����ʼ��
                B(ii, jj) -= B(i, jj) * coef;
            }
        }

        // ʹ��ԪΪ1
        for (int jj = j + 1; jj < n; ++jj) {
            B(i, jj) /= B(i, j);
        }
        B(i, j) = 1;

        ++j; //  ���м�����ɣ�ǰ������һ��

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