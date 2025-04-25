//4路 优化回代过程 对齐
#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <cstdlib>
using namespace std;
using namespace std::chrono;

void gaussian_elimination(double *A, double *b, double *x, int n)
{
    // 消去过程
    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            double factor = A[i * n + k] / A[k * n + k];
            for (int j = k + 1; j < n; j++) {
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

    if (!A || !b || !x) {
        cout << "Memory allocation failed" << endl;
        return 1;
    }

    for (int i = 0; i < n; ++i) {
        b[i] = 1.0;
        for (int j = 0; j < n; ++j) {
            A[i * n + j] = 1.0 + (i == j); // 加一个对角线偏移避免除0
        }
    }

    auto start = high_resolution_clock::now();
    gaussian_elimination(A, b, x, n);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "Back of AVX:" << duration.count() << "us" << endl;

    _aligned_free(A);
    _aligned_free(b);
    _aligned_free(x);
    return 0;
}
