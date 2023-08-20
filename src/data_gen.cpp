//
// Created by Catigeart on 2023/8/20.
//
#include <algorithm>
#include <cstdlib>
#include <ctime>

void gen_d_full_rank_matrix(int n, double* A, int lda) {
    srand((unsigned int)time(nullptr));
    // 随机生成标准正交矩阵（相当于单位阵做1类初等变换）
    int j_arr[n];
    for (int i = 0; i < n; i++) {
        j_arr[i] = i;
    }
    std::random_shuffle(j_arr, j_arr + n);
    for (int i = 0; i < n; i++) {
        A[i + j_arr[i] * lda] = 1;
    }
    // 随机进行n次3类初等变换
    // 2类初等变换可以视为3类初等变换加回本行，因此忽略
    for (int p = 0; p < n; p++) {
        int src_row = 0 + rand() % n;
        int dst_row = 0 + rand() % n;
        double k = (double) (-50 + rand() % 100) / 10.0;
        for (int j = 0; j < n; j++) {
            A[dst_row + j * lda] += k * A[src_row + j * lda];
        }
    }
}

void gen_d_positive_definite_matrix(int n, double* A, int lda) {
    // 可逆矩阵的转置相乘必为正定矩阵
    double *M = (double*) malloc(n * n * sizeof(double));
    gen_d_full_rank_matrix(n, M, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i + j * lda] = M[i + j * n] * M[j + i * n];
        }
    }
}