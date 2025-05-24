#include <xmmintrin.h> // SSE
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <cstring>
using namespace std;

int n;
float *A, *b, *x;

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
        {
            A[row1 * n + j] += factor * A[row2 * n + j];
        }
        b[row1] += factor * b[row2];
    }
}

// sse优化
void sse(float *A, float *b, float *x, int n)
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

        float temp[4];
        _mm_storeu_ps(temp, sum_vec);
        sum -= temp[0] + temp[1] + temp[2] + temp[3];

        for (; j < n; ++j)
            sum -= A[i * n + j] * x[j];

        x[i] = sum / A[i * n + i];
    }
}

int main()
{
    cout << "Enter matrix size: ";
    cin >> n;

    A = (float *)malloc(sizeof(float) * n * n);
    b = (float *)malloc(sizeof(float) * n);
    x = (float *)malloc(sizeof(float) * n);

    float *A_copy = (float *)malloc(sizeof(float) * n * n);
    float *b_copy = (float *)malloc(sizeof(float) * n);

    if (!A || !b || !x || !A_copy || !b_copy)
    {
        cout << "Memory allocation failed!" << endl;
        return 1;
    }

    initialize_matrix();
    memcpy(A_copy, A, sizeof(float) * n * n);
    memcpy(b_copy, b, sizeof(float) * n);

    auto start = chrono::high_resolution_clock::now();
    sse(A, b, x, n);
    auto end = chrono::high_resolution_clock::now();

    cout << "SSE用时: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << "us" << endl;

    free(A);
    free(b);
    free(x);
    free(A_copy);
    free(b_copy);

    return 0;
}
