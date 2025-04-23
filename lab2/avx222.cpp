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
    cout << "AVX Gaussian Elimination Done, Time: " << duration.count() << "us" << endl;

    _aligned_free(A);
    _aligned_free(b);
    _aligned_free(x);

    return 0;
}