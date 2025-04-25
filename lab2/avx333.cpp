//4路 优化回代过程 不对齐
#include <iostream>
#include <chrono>
#include <cstdlib>
#include <immintrin.h>
using namespace std;
using namespace std::chrono;

void gaussian_elimination_back_substitution_avx(double *A, double *b, double *x, int n)
{
    // 消去过程
    for (int k = 0; k < n; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = A[i * n + k] / A[k * n + k];
            for (int j = k + 1; j < n; ++j)
            {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    // 回代过程AVX优化
    x[n - 1] = b[n - 1] / A[(n - 1) * n + (n - 1)];
    for (int i = n - 2; i >= 0; i--)
    {
        __m256d sum_vec = _mm256_setzero_pd(); // 用于累加的向量
        int j = i + 1;
        for (; j + 3 < n; j += 4)
        {
            __m256d a_vec = _mm256_loadu_pd(&A[i * n + j]);
            __m256d x_vec = _mm256_loadu_pd(&x[j]);
            __m256d prod = _mm256_mul_pd(a_vec, x_vec);
            sum_vec = _mm256_add_pd(sum_vec, prod);
        }

        // 将 sum_vec 的 4 个 double 加在一起
        double sum_array[4];
        _mm256_storeu_pd(sum_array, sum_vec);
        double sum = sum_array[0] + sum_array[1] + sum_array[2] + sum_array[3];

        // 剩余部分串行累加
        for (; j < n; j++)
        {
            sum += A[i * n + j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i * n + i];
    }
}

int main()
{
    int n;
    cout << "Enter matrix dimension: ";
    cin >> n;

    size_t matrix_size = n * n * sizeof(double);
    size_t vector_size = n * sizeof(double);

    double *A = static_cast<double *>(malloc(matrix_size));
    double *b = static_cast<double *>(malloc(vector_size));
    double *x = static_cast<double *>(malloc(vector_size));

    if (!A || !b || !x)
    {
        cout << "Memory allocation failed" << endl;
        return 1;
    }

    // 初始化矩阵
    for (int i = 0; i < n; ++i)
    {
        b[i] = 1.0;
        for (int j = 0; j < n; ++j)
        {
            A[i * n + j] = 1.0 + (i == j); // 避免除0
        }
    }

    auto start = high_resolution_clock::now();
    gaussian_elimination_back_substitution_avx(A, b, x, n);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "Back of AVX333:" << duration.count() << "us" << endl;
    free(A);
    free(b);
    free(x);
    return 0;
}
