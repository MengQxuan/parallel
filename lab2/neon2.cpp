#include <arm_neon.h>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <chrono>
using namespace std;

// 串行消去过程
void neon0(float *A, float *b, float *x, int n)
{
    for (int k = 0; k < n; ++k)
    {
        for (int i = k + 1; i < n; ++i)
        {
            float factor = A[i * n + k] / A[k * n + k];
            for (int j = k + 1; j < n; ++j)
            {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    for (int i = n - 1; i >= 0; --i)
    {
        float sum = b[i];
        for (int j = i + 1; j < n; ++j)
        {
            sum -= A[i * n + j] * x[j];
        }
        x[i] = sum / A[i * n + i];
    }
}

// 优化消去过程
void neon1(float *A, float *b, float *x, int n)
{
    for (int k = 0; k < n; ++k)
    {
        for (int i = k + 1; i < n; ++i)
        {
            float factor = A[i * n + k] / A[k * n + k];
            float32x4_t factor_vec = vdupq_n_f32(factor);
            int j = k + 1;
            for (; j + 3 < n; j += 4)
            {
                float32x4_t a_i = vld1q_f32(&A[i * n + j]);
                float32x4_t a_k = vld1q_f32(&A[k * n + j]);
                float32x4_t res = vsubq_f32(a_i, vmulq_f32(factor_vec, a_k));
                vst1q_f32(&A[i * n + j], res);
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
        {
            sum -= A[i * n + j] * x[j];
        }
        x[i] = sum / A[i * n + i];
    }
}

// 优化回代过程
void neon2(float *A, float *b, float *x, int n)
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
        float32x4_t sum_vec = vmovq_n_f32(0.0f);
        for (; j + 3 < n; j += 4)
        {
            float32x4_t a_vec = vld1q_f32(&A[i * n + j]);
            float32x4_t x_vec = vld1q_f32(&x[j]);
            sum_vec = vmlaq_f32(sum_vec, a_vec, x_vec);
        }

        float32x2_t sum_low = vget_low_f32(sum_vec);
        float32x2_t sum_high = vget_high_f32(sum_vec);
        float32x2_t total = vpadd_f32(sum_low, sum_high);
        float final_sum = vget_lane_f32(vpadd_f32(total, total), 0);
        sum -= final_sum;

        for (; j < n; ++j)
            sum -= A[i * n + j] * x[j];

        x[i] = sum / A[i * n + i];
    }
}

// 消去和回代过程都优化
void neon3(float *A, float *b, float *x, int n)
{
    for (int k = 0; k < n; ++k)
    {
        for (int i = k + 1; i < n; ++i)
        {
            float factor = A[i * n + k] / A[k * n + k];
            float32x4_t factor_vec = vdupq_n_f32(factor);
            int j = k + 1;
            for (; j + 3 < n; j += 4)
            {
                float32x4_t a_i = vld1q_f32(&A[i * n + j]);
                float32x4_t a_k = vld1q_f32(&A[k * n + j]);
                float32x4_t res = vsubq_f32(a_i, vmulq_f32(factor_vec, a_k));
                vst1q_f32(&A[i * n + j], res);
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
        float32x4_t sum_vec = vmovq_n_f32(0.0f);
        for (; j + 3 < n; j += 4)
        {
            float32x4_t a_vec = vld1q_f32(&A[i * n + j]);
            float32x4_t x_vec = vld1q_f32(&x[j]);
            sum_vec = vmlaq_f32(sum_vec, a_vec, x_vec);
        }

        float32x2_t sum_low = vget_low_f32(sum_vec);
        float32x2_t sum_high = vget_high_f32(sum_vec);
        float32x2_t total = vpadd_f32(sum_low, sum_high);
        float final_sum = vget_lane_f32(vpadd_f32(total, total), 0);
        sum -= final_sum;

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

    // aligned_alloc
    float *A = (float *)aligned_alloc(16, sizeof(float) * n * n);
    float *b = (float *)aligned_alloc(16, sizeof(float) * n);
    float *x = (float *)aligned_alloc(16, sizeof(float) * n);

    if (!A || !b || !x)
    {
        cout << "Allocation failed\n";
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

    auto start = chrono::high_resolution_clock::now();
    neon0(A, b, x, n);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "neon0:" << duration.count() << "us" << endl;

    start = chrono::high_resolution_clock::now();
    neon1(A, b, x, n);
    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "neon1:" << duration.count() << "us" << endl;

    start = chrono::high_resolution_clock::now();
    neon2(A, b, x, n);
    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "neon2:" << duration.count() << "us" << endl;

    start = chrono::high_resolution_clock::now();
    neon3(A, b, x, n);
    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "neon3:" << duration.count() << "us" << endl;

    free(A);
    free(b);
    free(x);
    return 0;
}
