#include <iostream>
#include <chrono>
#include <cstdlib>
#include <malloc.h>
#include <immintrin.h>
using namespace std;
using namespace std::chrono;

void gaussian_elimination_avx_no_alignment(double *A, double *b, double *x, int n)
{
    for (int k = 0; k < n; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = A[i * n + k] / A[k * n + k];
            __m256d factor_vec = _mm256_set1_pd(factor); // factor 向量 4路向量化
            //_mm256_set1_pd(factor) 创建了一个包含 4 个 factor 值的 AVX 向量，类型为 __m256d。
            // 每个 __m256d 类型包含 4 个 double 类型的元素（每个元素占 8 字节，总共 32 字节，256 位

            int j = k + 1;
            for (; j + 3 < n; j += 4)
            {
                double *addr_i = &A[i * n + j];
                double *addr_k = &A[k * n + j];

                // 无论是否对齐，都使用 _mm256_loadu_pd 和 _mm256_storeu_pd
                __m256d vec_i = _mm256_loadu_pd(addr_i); // 非对齐加载
                __m256d vec_k = _mm256_loadu_pd(addr_k); // 非对齐加载
                __m256d product = _mm256_mul_pd(factor_vec, vec_k);
                __m256d result = _mm256_sub_pd(vec_i, product);
                _mm256_storeu_pd(addr_i, result); // 非对齐存储
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
        for (int j = i + 1; j < n; j++)
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

    size_t matrix_size = n * n * sizeof(double);
    size_t vector_size = n * sizeof(double);

    // 内存不对齐分配
    double *A = static_cast<double *>(malloc(matrix_size));
    double *b = static_cast<double *>(malloc(vector_size));
    double *x = static_cast<double *>(malloc(vector_size));

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
            A[i * n + j] = 1.0 + (i == j); // 避免除0
        }
    }

    auto start = high_resolution_clock::now();
    gaussian_elimination_avx_no_alignment(A, b, x, n);
    auto end = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(end - start);
    cout << "AVX Gaussian Elimination Done, Time: " << duration.count() << "us" << endl;

    free(A);
    free(b);
    free(x);

    return 0;
}
