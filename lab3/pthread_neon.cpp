// neon+静态线程+barrier同步
#include <arm_neon.h>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <chrono>
#include <pthread.h>
using namespace std;

float *A, *b, *x;
int n;
int num_threads;
pthread_barrier_t barrier;

struct ThreadParam
{
    int thread_id;
};

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

void *threadFunc(void *arg)
{
    ThreadParam *param = (ThreadParam *)arg;
    int t_id = param->thread_id;

    for (int k = 0; k < n; ++k)
    {
        pthread_barrier_wait(&barrier);

        for (int i = k + 1 + t_id; i < n; i += num_threads)
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

        pthread_barrier_wait(&barrier);
    }

    pthread_exit(NULL);
}

void back_substitution()
{
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
    // cout << "Enter number of threads: ";
    // cin >> num_threads;
    n = 2048;
    num_threads = 8;

    A = (float *)malloc(sizeof(float) * n * n);
    b = (float *)malloc(sizeof(float) * n);
    x = (float *)malloc(sizeof(float) * n);

    if (!A || !b || !x)
    {
        cout << "Allocation failed\n";
        return 1;
    }

    initialize_matrix();
    pthread_barrier_init(&barrier, NULL, num_threads);

    pthread_t *threads = new pthread_t[num_threads];
    ThreadParam *params = new ThreadParam[num_threads];
    for (int t_id = 0; t_id < num_threads; ++t_id)
    {
        params[t_id].thread_id = t_id;
        pthread_create(&threads[t_id], NULL, threadFunc, &params[t_id]);
    }

    auto start = chrono::high_resolution_clock::now();

    for (int t_id = 0; t_id < num_threads; ++t_id)
        pthread_join(threads[t_id], NULL);

    back_substitution();

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "静态线程+barrier+NEON优化，用时：" << duration.count() << "us" << endl;

    pthread_barrier_destroy(&barrier);
    delete[] threads;
    delete[] params;
    free(A);
    free(b);
    free(x);
    return 0;
}