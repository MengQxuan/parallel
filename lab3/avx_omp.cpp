//AVX + OpenMP
#include <iostream>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <immintrin.h>
#include <omp.h>
using namespace std;
using namespace std::chrono;

int n;
double *A, *b, *x;

void initialize_matrix()
{
    srand(time(0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            A[i * n + j] = 0.0;

        A[i * n + i] = rand() % 100 + 1000;
        for (int j = i + 1; j < n; ++j)
            A[i * n + j] = rand() % 100 + 1;

        b[i] = rand() % 100 + 1;
    }

    for (int k = 0; k < n / 2; ++k)
    {
        int row1 = rand() % n;
        int row2 = rand() % n;
        double factor = (rand() % 100) / 100.0;
        for (int j = 0; j < n; ++j)
        {
            A[row1 * n + j] += factor * A[row2 * n + j];
        }
        b[row1] += factor * b[row2];
    }
}

void gaussian_elimination_avx_openmp(int num_threads)
{
    omp_set_num_threads(num_threads);

#pragma omp parallel shared(A, b, n)
    {
        for (int k = 0; k < n; ++k)
        {
            // 串行标准化第 k 行
#pragma omp single
            {
                double pivot = A[k * n + k];
                for (int j = k + 1; j < n; ++j)
                    A[k * n + j] /= pivot;
                b[k] /= pivot;
                A[k * n + k] = 1.0;
            }

            // 并行消元（k+1 到 n 行）
#pragma omp for schedule(static)
            for (int i = k + 1; i < n; ++i)
            {
                double factor = A[i * n + k];
                __m256d factor_vec = _mm256_set1_pd(factor);

                int j = k + 1;
                for (; j + 3 < n; j += 4)
                {
                    __m256d row_i = _mm256_loadu_pd(&A[i * n + j]);
                    __m256d row_k = _mm256_loadu_pd(&A[k * n + j]);
                    __m256d result = _mm256_sub_pd(row_i, _mm256_mul_pd(factor_vec, row_k));
                    _mm256_storeu_pd(&A[i * n + j], result);
                }

                for (; j < n; ++j)
                    A[i * n + j] -= factor * A[k * n + j];

                b[i] -= factor * b[k];
                A[i * n + k] = 0.0;
            }
        }
    }

    // 回代阶段（串行）
    x[n - 1] = b[n - 1] / A[(n - 1) * n + (n - 1)];
    for (int i = n - 2; i >= 0; i--)
    {
        int j = i + 1;
        __m256d sum_vec = _mm256_setzero_pd();
        double sum = 0.0;

        for (; j + 3 < n; j += 4)
        {
            __m256d a_vec = _mm256_loadu_pd(&A[i * n + j]);
            __m256d x_vec = _mm256_loadu_pd(&x[j]);
            sum_vec = _mm256_add_pd(sum_vec, _mm256_mul_pd(a_vec, x_vec));
        }

        alignas(32) double temp[4];
        _mm256_storeu_pd(temp, sum_vec);
        sum = temp[0] + temp[1] + temp[2] + temp[3];

        for (; j < n; ++j)
            sum += A[i * n + j] * x[j];

        x[i] = (b[i] - sum) / A[i * n + i];
    }
}

int main()
{
    cout << "Enter matrix dimension: ";
    cin >> n;

    int num_threads = 2; // 可配置线程数

    A = (double *)malloc(sizeof(double) * n * n);
    b = (double *)malloc(sizeof(double) * n);
    x = (double *)malloc(sizeof(double) * n);

    if (!A || !b || !x)
    {
        cout << "Memory allocation failed." << endl;
        return 1;
    }

    initialize_matrix();

    auto start = high_resolution_clock::now();
    gaussian_elimination_avx_openmp(num_threads);
    auto end = high_resolution_clock::now();

    cout << "AVX+OpenMP用时：" << duration_cast<microseconds>(end - start).count() << " us" << endl;

    free(A);
    free(b);
    free(x);
    return 0;
}
