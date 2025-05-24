#include <arm_neon.h>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <chrono>
using namespace std;

float *A, *b, *x;
int n;

void initialize_matrix()
{
    srand(time(0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            A[i * n + j] = 0.0f;

        A[i * n + i] = rand() % 100 + 1000;
        for (int j = i + 1; j < n; ++j)
            A[i * n + j] = rand() % 100 + 1;

        b[i] = rand() % 100 + 1;
    }

    for (int k = 0; k < n / 2; ++k)
    {
        int row1 = rand() % n;
        int row2 = rand() % n;
        float factor = (rand() % 100) / 100.0f;
        for (int j = 0; j < n; ++j)
            A[row1 * n + j] += factor * A[row2 * n + j];
        b[row1] += factor * b[row2];
    }
}

void neon(float *A, float *b, float *x, int n)
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
    // cout << "Enter matrix size: ";
    // cin >> n;
    n = 2048;

    A = (float *)malloc(sizeof(float) * n * n);
    b = (float *)malloc(sizeof(float) * n);
    x = (float *)malloc(sizeof(float) * n);

    if (!A || !b || !x)
    {
        cout << "Allocation failed\n";
        return 1;
    }

    initialize_matrix();

    auto start = chrono::high_resolution_clock::now();
    neon(A, b, x, n);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "NEON，用时：" << duration.count() << "us" << endl;

    free(A);
    free(b);
    free(x);
    return 0;
}
