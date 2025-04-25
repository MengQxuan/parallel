// #include <iostream>
// #include <chrono>
// #include <immintrin.h>
// #include <cstdlib>
// using namespace std;
// using namespace std::chrono;

// void gaussian_elimination_avx512(double *A, double *b, double *x, int n)
// {
//     // 消去过程（AVX512 8路向量化）
//     for (int k = 0; k < n; ++k)
//     {
//         for (int i = k + 1; i < n; ++i)
//         {
//             double factor = A[i * n + k] / A[k * n + k];
//             __m512d factor_vec = _mm512_set1_pd(factor);

//             int j = k + 1;
//             for (; j + 7 < n; j += 8) // 每次处理8个元素
//             {
//                 // 直接使用未对齐的加载/存储
//                 __m512d row_i = _mm512_loadu_pd(&A[i * n + j]);
//                 __m512d row_k = _mm512_loadu_pd(&A[k * n + j]);
//                 __m512d result = _mm512_sub_pd(row_i, _mm512_mul_pd(factor_vec, row_k));
//                 _mm512_storeu_pd(&A[i * n + j], result);
//             }

//             // 处理剩余元素（非8的倍数）
//             for (; j < n; ++j)
//             {
//                 A[i * n + j] -= factor * A[k * n + j];
//             }

//             b[i] -= factor * b[k];
//         }
//     }

//     // 回代过程（AVX512 8路向量化）
//     x[n - 1] = b[n - 1] / A[(n - 1) * n + (n - 1)];
//     for (int i = n - 2; i >= 0; --i)
//     {
//         int j = i + 1;
//         __m512d sum_vec = _mm512_setzero_pd();
//         double sum = b[i]; // 初始化为b[i]

//         // 向量化处理j的循环（每次处理8个元素）
//         for (; j + 7 < n; j += 8)
//         {
//             // 加载未对齐的8个元素
//             __m512d a_vec = _mm512_loadu_pd(&A[i * n + j]);
//             __m512d x_vec = _mm512_loadu_pd(&x[j]);
//             __m512d product = _mm512_mul_pd(a_vec, x_vec);

//             // 累加到sum_vec中
//             sum_vec = _mm512_add_pd(sum_vec, product);
//         }

//         // 将向量结果合并到标量sum中
//         alignas(64) double temp[8]; // 8个元素对齐到64字节
//         _mm512_storeu_pd(temp, sum_vec);
//         sum -= (temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7]);

//         // 处理剩余元素
//         for (; j < n; ++j)
//         {
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

//     size_t align = 64; // AVX512建议64字节对齐，但此处允许内存不对齐
//     size_t matrix_size = sizeof(double) * n * n;
//     size_t vector_size = sizeof(double) * n;

//     double *A = (double *)_aligned_malloc(matrix_size, align);
//     double *b = (double *)_aligned_malloc(vector_size, align);
//     double *x = (double *)_aligned_malloc(vector_size, align);

//     if (!A || !b || !x)
//     {
//         cout << "Memory allocation failed." << endl;
//         return 1;
//     }

//     // 初始化矩阵（对角线加1避免除零）
//     for (int i = 0; i < n; ++i)
//     {
//         b[i] = 1.0;
//         for (int j = 0; j < n; ++j)
//         {
//             A[i * n + j] = 1.0 + (i == j);
//         }
//     }

//     auto start = high_resolution_clock::now();
//     gaussian_elimination_avx512(A, b, x, n);
//     auto end = high_resolution_clock::now();
//     auto duration = duration_cast<microseconds>(end - start);

//     cout << "AVX512" << duration.count() << "us" << endl;

//     _aligned_free(A);
//     _aligned_free(b);
//     _aligned_free(x);

//     return 0;
// }

// #include <iostream>
// #include <chrono>
// #include <immintrin.h>
// #include <malloc.h>
// using namespace std;
// using namespace std::chrono;

// void gaussian_elimination_avx512(double *A, double *b, double *x, int n)
// {
//     for (int k = 0; k < n; ++k)
//     {
//         for (int i = k + 1; i < n; ++i)
//         {
//             double factor = A[i * n + k] / A[k * n + k];
//             __m512d factor_vec = _mm512_set1_pd(factor);

//             int j = k + 1;
//             for (; j + 7 < n; j += 8)
//             {
//                 __m512d row_i = _mm512_loadu_pd(&A[i * n + j]);
//                 __m512d row_k = _mm512_loadu_pd(&A[k * n + j]);
//                 __m512d result = _mm512_sub_pd(row_i, _mm512_mul_pd(factor_vec, row_k));
//                 _mm512_storeu_pd(&A[i * n + j], result);
//             }

//             for (; j < n; ++j)
//                 A[i * n + j] -= factor * A[k * n + j];
//             b[i] -= factor * b[k];
//         }
//     }

//     x[n - 1] = b[n - 1] / A[(n - 1) * n + (n - 1)];
//     for (int i = n - 2; i >= 0; --i)
//     {
//         int j = i + 1;
//         __m512d sum_vec = _mm512_setzero_pd();
//         double sum = b[i]; // 初始化为b[i]

//         for (; j + 7 < n; j += 8)
//         {
//             __m512d a_vec = _mm512_loadu_pd(&A[i * n + j]);
//             __m512d x_vec = _mm512_loadu_pd(&x[j]);
//             __m512d product = _mm512_mul_pd(a_vec, x_vec);
//             sum_vec = _mm512_add_pd(sum_vec, product);
//         }

//         alignas(64) double temp[8];
//         _mm512_storeu_pd(temp, sum_vec);
//         sum -= (temp[0] + temp[1] + temp[2] + temp[3] 
//               + temp[4] + temp[5] + temp[6] + temp[7]);

//         for (; j < n; ++j)
//             sum -= A[i * n + j] * x[j];

//         x[i] = sum / A[i * n + i]; // 确保分母非零
//     }
// }

// int main()
// {
//     int n;
//     cout << "Enter matrix dimension: ";
//     cin >> n;

//     size_t align = 64; // AVX512要求64字节对齐
//     size_t matrix_size = sizeof(double) * n * n;
//     size_t vector_size = sizeof(double) * n;

//     double *A = (double*)_aligned_malloc(matrix_size, align);
//     double *b = (double*)_aligned_malloc(vector_size, align);
//     double *x = (double*)_aligned_malloc(vector_size, align);

//     if (!A || !b || !x)
//     {
//         cout << "Memory allocation failed." << endl;
//         return 1;
//     }

//     for (int i = 0; i < n; ++i)
//     {
//         b[i] = 1.0;
//         for (int j = 0; j < n; ++j)
//             A[i * n + j] = 1.0 + (i == j ? 1.0 : 0.0); // 对角线为2.0，其他为1.0
//     }

//     auto start = high_resolution_clock::now();
//     gaussian_elimination_avx512(A, b, x, n);
//     auto end = high_resolution_clock::now();
//     auto duration = duration_cast<microseconds>(end - start);
//     cout << "AVX512 Time: " << duration.count() << "us" << endl;

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

void gaussian_elimination_avx512(double *A, double *b, double *x, int n)
{
    for (int k = 0; k < n; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = A[i * n + k] / A[k * n + k];
            __m512d factor_vec = _mm512_set1_pd(factor); // 512-bit 向量

            int j = k + 1;
            for (; j + 7 < n; j += 8)
            {
                double *addr_i = &A[i * n + j];
                double *addr_k = &A[k * n + j];

                __m512d row_i = _mm512_loadu_pd(addr_i);
                __m512d row_k = _mm512_loadu_pd(addr_k);
                __m512d res = _mm512_sub_pd(row_i, _mm512_mul_pd(factor_vec, row_k));
                _mm512_storeu_pd(addr_i, res);
            }

            for (; j < n; ++j)
            {
                A[i * n + j] -= factor * A[k * n + j];
            }

            b[i] -= factor * b[k];
        }
    }

    // 回代过程
    x[n - 1] = b[n - 1] / A[(n - 1) * n + (n - 1)];
    for (int i = n - 2; i >= 0; i--)
    {
        double sum = b[i];
        int j = i + 1;

        // 每次处理8个double元素
        for (; j + 7 < n; j += 8)
        {
            double *A_ptr = &A[i * n + j];
            double *x_ptr = &x[j];

            __m512d a_vec = _mm512_loadu_pd(A_ptr);
            __m512d x_vec = _mm512_loadu_pd(x_ptr);
            __m512d prod = _mm512_mul_pd(a_vec, x_vec);

            // 累加prod中8个double
            __m256d lo = _mm512_castpd512_pd256(prod);                     // 前4个
            __m256d hi = _mm512_extractf64x4_pd(prod, 1);                  // 后4个
            __m256d sum256 = _mm256_add_pd(lo, hi);                        // 加起来成4个

            __m128d low = _mm256_castpd256_pd128(sum256);                 // 前2个
            __m128d high = _mm256_extractf128_pd(sum256, 1);              // 后2个
            __m128d sum128 = _mm_add_pd(low, high);                       // 加成2个
            __m128d final_sum = _mm_hadd_pd(sum128, sum128);              // 横向相加

            double temp[2];
            _mm_store_pd(temp, final_sum);
            sum -= temp[0];
        }

        for (; j < n; j++)
        {
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

    size_t align = 64;
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
            A[i * n + j] = 1.0 + (i == j);
        }
    }

    auto start = high_resolution_clock::now();
    gaussian_elimination_avx512(A, b, x, n);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "AVX-512" << duration.count() << "us" << endl;

    _aligned_free(A);
    _aligned_free(b);
    _aligned_free(x);

    return 0;
}
