//AVX优化 对齐
// #include <iostream>
// #include <chrono>
// #include <immintrin.h>
// #include <cstdlib>
// using namespace std;
// using namespace std::chrono;

// // 使用一维数组表示矩阵：A[i][j] = A[i * n + j]
// void gaussian_elimination(double *A, double *b, double *x, int n)
// {
//     // 消去过程优化（AVX）
//     for (int k = 0; k < n; k++) {
//         // 选择主元（列主元法）
//         int max_row = k;
//         double max_val = std::abs(A[k * n + k]);
//         for (int i = k + 1; i < n; i++) {
//             double val = std::abs(A[i * n + k]);
//             if (val > max_val) {
//                 max_val = val;
//                 max_row = i;
//             }
//         }
//         // 交换行
//         if (max_row != k) {
//             for (int j = 0; j < n; j++) {
//                 swap(A[k * n + j], A[max_row * n + j]);
//             }
//             swap(b[k], b[max_row]);
//         }

//         // 消去当前列
//         double pivot = A[k * n + k];
//         for (int i = k + 1; i < n; i++) {
//             double factor = A[i * n + k] / pivot;
//             __m256d factor_vec = _mm256_set1_pd(factor);

//             // 向量化处理j循环
//             int j = k + 1;
//             for (; j <= n - 4; j += 4) {
//                 double *a_i = &A[i * n + j];
//                 double *a_k = &A[k * n + j];
//                 __m256d a_i_vec = _mm256_loadu_pd(a_i);
//                 __m256d a_k_vec = _mm256_loadu_pd(a_k);
//                 __m256d res = _mm256_sub_pd(a_i_vec, _mm256_mul_pd(factor_vec, a_k_vec));
//                 _mm256_storeu_pd(a_i, res);
//             }

//             // 处理剩余元素
//             for (; j < n; j++) {
//                 A[i * n + j] -= factor * A[k * n + j];
//             }

//             b[i] -= factor * b[k];
//         }
//     }

//     // 回代过程优化（AVX）
//     x[n - 1] = b[n - 1] / A[(n - 1) * n + (n - 1)];
//     for (int i = n - 2; i >= 0; i--) {
//         double sum = b[i];
//         int j = i + 1;

//         // 向量化处理j循环
//         while (j <= n - 4) {
//             double *a_row = &A[i * n + j];
//             double *x_ptr = &x[j];
//             __m256d a_vec = _mm256_loadu_pd(a_row);
//             __m256d x_vec = _mm256_loadu_pd(x_ptr);
//             __m256d product = _mm256_mul_pd(a_vec, x_vec);

//             // 合并向量元素到标量sum
//             __m128d part1 = _mm256_castpd256_pd128(product);
//             __m128d part2 = _mm256_extractf128_pd(product, 1);
//             part1 = _mm_add_pd(part1, part2);
//             __m128d total = _mm_hadd_pd(part1, part1);
//             double temp[2];
//             _mm_store_pd(temp, total);
//             sum -= temp[0] + temp[1]; // 合并两个部分的总和

//             j += 4;
//         }

//         // 处理剩余元素
//         for (; j < n; j++) {
//             sum -= A[i * n + j] * x[j];
//         }

//         x[i] = sum / A[i * n + i];
//     }
// }

// int main()
// {
//     int n;
//     cout << "Enter matrix dimension: ";
//     cin >> n;

//     size_t align = 32; // AVX要求32字节对齐
//     size_t matrix_size = sizeof(double) * n * n;
//     size_t vector_size = sizeof(double) * n;

//     // 对齐内存分配
//     double *A = (double*)_aligned_malloc(matrix_size, align);
//     double *b = (double*)_aligned_malloc(vector_size, align);
//     double *x = (double*)_aligned_malloc(vector_size, align);

//     if (!A || !b || !x) {
//         cout << "Memory allocation failed" << endl;
//         return 1;
//     }

//     // 初始化矩阵（对角线偏移避免除零）
//     for (int i = 0; i < n; ++i) {
//         b[i] = 1.0;
//         for (int j = 0; j < n; ++j) {
//             A[i * n + j] = 1.0 + (i == j); 
//         }
//     }

//     auto start = high_resolution_clock::now();
//     gaussian_elimination(A, b, x, n);
//     auto end = high_resolution_clock::now();
//     auto duration = duration_cast<microseconds>(end - start);
//     cout << "AVX:" << duration.count() << "us" << endl;

//     // 释放内存
//     _aligned_free(A);
//     _aligned_free(b);
//     _aligned_free(x);

//     return 0;
// }


#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <cstdlib>
using namespace std;
using namespace std::chrono;

// 使用一维数组表示矩阵：A[i][j] = A[i * n + j]
void gaussian_elimination(double *A, double *b, double *x, int n)
{
    for (int k = 0; k < n; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = A[i * n + k] / A[k * n + k];
            __m256d factor_vec = _mm256_set1_pd(factor);

            int j = k + 1;
            for (; j + 3 < n; j += 4)
            {
                // 确保地址对齐
                double *addr_i = &A[i * n + j];
                double *addr_k = &A[k * n + j];

                if (((uintptr_t)addr_i % 32 == 0) && ((uintptr_t)addr_k % 32 == 0))
                {
                    __m256d row_i = _mm256_load_pd(addr_i);
                    __m256d row_k = _mm256_load_pd(addr_k);
                    __m256d res = _mm256_sub_pd(row_i, _mm256_mul_pd(factor_vec, row_k));
                    _mm256_store_pd(addr_i, res);
                }
                else
                {
                    // 如果未对齐，使用 _mm256_loadu_pd
                    __m256d row_i = _mm256_loadu_pd(addr_i);
                    __m256d row_k = _mm256_loadu_pd(addr_k);
                    __m256d res = _mm256_sub_pd(row_i, _mm256_mul_pd(factor_vec, row_k));
                    _mm256_storeu_pd(addr_i, res);
                }
            }

            for (; j < n; ++j)
            {
                A[i * n + j] -= factor * A[k * n + j];
            }

            b[i] -= factor * b[k];
        }
    }

    // 回代过程优化
    x[n - 1] = b[n - 1] / A[(n - 1) * n + (n - 1)];
    for (int i = n - 2; i >= 0; i--) {
        double sum = b[i];
        int j = i + 1;

        // 向量化处理j的循环
        while (j <= n - 4) { // 每次处理4个元素
            double *A_ptr = &A[i * n + j];
            double *x_ptr = &x[j];
            __m256d a_vec = _mm256_loadu_pd(A_ptr);
            __m256d x_vec = _mm256_loadu_pd(x_ptr);
            __m256d product = _mm256_mul_pd(a_vec, x_vec);

            // 计算四个元素的总和
            __m128d a_low = _mm256_castpd256_pd128(product);
            __m128d a_high = _mm256_extractf128_pd(product, 1);
            a_low = _mm_add_pd(a_low, a_high); // a0+a2, a1+a3
            __m128d sum1 = _mm_hadd_pd(a_low, a_low); // (a0+a1+a2+a3, ...)
            double temp[2];
            _mm_store_pd(temp, sum1);
            sum -= temp[0];

            j += 4;
        }

        // 处理剩下的j
        for (; j < n; j++) {
            sum -= A[i * n + j] * x[j];
        }

        x[i] = sum / A[i * n + i];
    }
}

int main()
{
    int n;
    cout << "Enter matrix dimension: ";
    cin >> n;

    size_t align = 32;
    size_t matrix_size = sizeof(double) * n * n;
    size_t vector_size = sizeof(double) * n;

    double *A = (double *)_aligned_malloc(matrix_size, align);
    double *b = (double *)_aligned_malloc(vector_size, align);
    double *x = (double *)_aligned_malloc(vector_size, align);

    if (!A || !b || !x)
    {
        cout << "Memory allocation failed" << endl;
        return 1;
    }

    for (int i = 0; i < n; ++i)
    {
        b[i] = 1.0;
        for (int j = 0; j < n; ++j)
        {
            A[i * n + j] = 1.0 + (i == j); // 加一个对角线偏移避免除0
        }
    }

    auto start = high_resolution_clock::now();
    gaussian_elimination(A, b, x, n);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "AVX:" << duration.count() << "us" << endl;
    _aligned_free(A);
    _aligned_free(b);
    _aligned_free(x);

    return 0;
}