#include <pmmintrin.h> //sse3
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <chrono>
using namespace std;

// 消去优化
void sse1(float *A, float *b, float *x, int n)
{
    for (int k = 0; k < n; ++k)
    {
        for (int i = k + 1; i < n; ++i)
        {
            float factor = A[i * n + k] / A[k * n + k];
            __m128 factor_vec = _mm_set1_ps(factor);
            int j = k + 1;
            for (; j + 3 < n; j += 4)
            {
                // 预取后续数据
                if (j + 16 < n)
                    _mm_prefetch((const char *)&A[i * n + j + 16], _MM_HINT_T0);

                __m128 a_i = _mm_loadu_ps(&A[i * n + j]);
                __m128 a_k = _mm_loadu_ps(&A[k * n + j]);
                __m128 res = _mm_sub_ps(a_i, _mm_mul_ps(factor_vec, a_k));
                _mm_storeu_ps(&A[i * n + j], res);
            }

            for (; j < n; ++j)
                A[i * n + j] -= factor * A[k * n + j];
            b[i] -= factor * b[k];
        }
    }

    for (int i = n - 1; i >= 0; --i)
    {
        float sum = b[i];
        for (int j = i + 1; j < n; ++j)
            sum -= A[i * n + j] * x[j];
        x[i] = sum / A[i * n + i];
    }
}

// 回代优化
void sse2(float *A, float *b, float *x, int n)
{
    for (int k = 0; k < n; ++k)
    {
        for (int i = k + 1; i < n; ++i)
        {
            float factor = A[i * n + k] / A[k * n + k];
            for (int j = k + 1; j < n; ++j)
                A[i * n + j] -= factor * A[k * n + j];
            b[i] -= factor * b[k];
        }
    }

    for (int i = n - 1; i >= 0; --i)
    {
        float sum = b[i];
        int j = i + 1;
        __m128 sum_vec = _mm_setzero_ps();
        for (; j + 3 < n; j += 4)
        {
            __m128 a_vec = _mm_loadu_ps(&A[i * n + j]);
            __m128 x_vec = _mm_loadu_ps(&x[j]);
            sum_vec = _mm_add_ps(sum_vec, _mm_mul_ps(a_vec, x_vec));
        }

        // 使用SSE3的水平加法指令合并向量元素
        sum_vec = _mm_hadd_ps(sum_vec, sum_vec); // 合并为两个元素
        sum_vec = _mm_hadd_ps(sum_vec, sum_vec); // 合并为一个元素
        sum -= _mm_cvtss_f32(sum_vec);           // 提取单精度浮点数

        for (; j < n; ++j)
            sum -= A[i * n + j] * x[j];

        x[i] = sum / A[i * n + i];
    }
}

void sse3(float *A, float *b, float *x, int n)
{
    for (int k = 0; k < n; ++k)
    {
        for (int i = k + 1; i < n; ++i)
        {
            float factor = A[i * n + k] / A[k * n + k];
            __m128 factor_vec = _mm_set1_ps(factor);
            int j = k + 1;
            for (; j + 3 < n; j += 4)
            {
                // 预取后续数据
                if (j + 16 < n)
                    _mm_prefetch((const char *)&A[i * n + j + 16], _MM_HINT_T0);

                __m128 a_i = _mm_loadu_ps(&A[i * n + j]);
                __m128 a_k = _mm_loadu_ps(&A[k * n + j]);
                __m128 res = _mm_sub_ps(a_i, _mm_mul_ps(factor_vec, a_k));
                _mm_storeu_ps(&A[i * n + j], res);
            }

            for (; j < n; ++j)
                A[i * n + j] -= factor * A[k * n + j];
            b[i] -= factor * b[k];
        }
    }

    for (int i = n - 1; i >= 0; --i)
    {
        float sum = b[i];
        int j = i + 1;
        __m128 sum_vec = _mm_setzero_ps();
        for (; j + 3 < n; j += 4)
        {
            __m128 a_vec = _mm_loadu_ps(&A[i * n + j]);
            __m128 x_vec = _mm_loadu_ps(&x[j]);
            sum_vec = _mm_add_ps(sum_vec, _mm_mul_ps(a_vec, x_vec));
        }

        // 使用SSE3的水平加法指令合并向量元素
        sum_vec = _mm_hadd_ps(sum_vec, sum_vec); // 合并为两个元素
        sum_vec = _mm_hadd_ps(sum_vec, sum_vec); // 合并为一个元素
        sum -= _mm_cvtss_f32(sum_vec);           // 提取单精度浮点数

        for (; j < n; ++j)
            sum -= A[i * n + j] * x[j];

        x[i] = sum / A[i * n + i];
    }
}

int main()
{
    int n;
    cout << "Enter matrix size: ";
    cin >> n;

    float *A = (float *)malloc(sizeof(float) * n * n);
    float *b = (float *)malloc(sizeof(float) * n);
    float *x = (float *)malloc(sizeof(float) * n);

    float *A_copy = (float *)malloc(sizeof(float) * n * n);
    float *b_copy = (float *)malloc(sizeof(float) * n);
    float *x_copy = (float *)malloc(sizeof(float) * n);

    if (!A || !b || !x)
    {
        cout << "Memory allocation failed!" << endl;
        return 1;
    }

    for (int i = 0; i < n; ++i)
    {
        b[i] = 1.0f;
        for (int j = 0; j < n; ++j)
        {
            A[i * n + j] = (i == j) ? 2.0f : 1.0f;
        }
    }

    memcpy(A_copy, A, sizeof(float) * n * n);
    memcpy(b_copy, b, sizeof(float) * n);

    auto start = chrono::high_resolution_clock::now();
    sse1(A, b, x, n);
    auto end = chrono::high_resolution_clock::now();
    cout << "SSE1:" << chrono::duration_cast<chrono::microseconds>(end - start).count() << "us" << endl;

    memcpy(A, A_copy, sizeof(float) * n * n);
    memcpy(b, b_copy, sizeof(float) * n);
    start = chrono::high_resolution_clock::now();
    sse2(A, b, x, n);
    end = chrono::high_resolution_clock::now();
    cout << "SSE2:" << chrono::duration_cast<chrono::microseconds>(end - start).count() << "us" << endl;

    memcpy(A, A_copy, sizeof(float) * n * n);
    memcpy(b, b_copy, sizeof(float) * n);
    start = chrono::high_resolution_clock::now();
    sse3(A, b, x, n);
    end = chrono::high_resolution_clock::now();
    cout << "SSE3:" << chrono::duration_cast<chrono::microseconds>(end - start).count() << "us" << endl;

    free(A);
    free(b);
    free(x);
    free(A_copy);
    free(b_copy);
    free(x_copy);

    return 0;
}
