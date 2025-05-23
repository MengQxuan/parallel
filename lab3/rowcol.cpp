#include <iostream>
#include <pthread.h>
#include <chrono>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace std;
using namespace std::chrono;

#define MAX_THREADS 8

double **A, *b, *x;
int n;
int num_threads;

pthread_barrier_t barrier_Division;
pthread_barrier_t barrier_Elimination;

struct ThreadData
{
    int thread_id;
};

void initialize_matrix()
{
    srand(time(0));
    for (int i = 0; i < n; ++i)
    {
        fill(A[i], A[i] + n, 0.0);
        A[i][i] = rand() % 100 + 1000;
        for (int j = i + 1; j < n; ++j)
        {
            A[i][j] = rand() % 100 + 1;
        }
        b[i] = rand() % 100 + 1;
    }
    for (int k = 0; k < n / 2; ++k)
    {
        int row1 = rand() % n;
        int row2 = rand() % n;
        double factor = (rand() % 100) / 100.0;
        for (int j = 0; j < n; ++j)
        {
            A[row1][j] += factor * A[row2][j];
        }
        b[row1] += factor * b[row2];
    }
}

void *threadFunc(void *param)
{
    ThreadData *p = (ThreadData *)param;
    int t_id = p->thread_id;

    for (int k = 0; k < n; ++k)
    {
        // 线程 0 执行除法操作
        if (t_id == 0)
        {
            for (int j = k + 1; j < n; ++j)
            {
                A[k][j] = A[k][j] / A[k][k];
            }
            b[k] = b[k] / A[k][k];
            A[k][k] = 1.0;
        }

        pthread_barrier_wait(&barrier_Division);

        // // 支持水平划分（按行）和垂直划分（按列）
        // bool row_partition = false; // false为列划分

        // if (row_partition)
        // {
        //     for (int i = k + 1 + t_id; i < n; i += num_threads)
        //     {
        //         double factor = A[i][k];
        //         for (int j = k + 1; j < n; ++j)
        //         {
        //             A[i][j] -= factor * A[k][j];
        //         }
        //         A[i][k] = 0.0;
        //         b[i] -= factor * b[k];
        //     }
        // }
        // else
        // {
        //     for (int i = k + 1; i < n; ++i)
        //     {
        //         double factor = A[i][k];
        //         for (int j = k + 1 + t_id; j < n; j += num_threads)
        //         {
        //             A[i][j] -= factor * A[k][j];
        //         }
        //         if (t_id == 0)
        //         {
        //             A[i][k] = 0.0;
        //             b[i] -= factor * b[k];
        //         }
        //     }
        // }

        for (int i = k + 1; i < n; ++i) // 列划分
        {
            double factor = A[i][k];
            for (int j = k + 1 + t_id; j < n; j += num_threads)
            {
                A[i][j] -= factor * A[k][j];
            }
            if (t_id == 0)
            {
                A[i][k] = 0.0;
                b[i] -= factor * b[k];
            }
        }

        pthread_barrier_wait(&barrier_Elimination);
    }

    pthread_exit(NULL);
    return NULL;
}

void back_substitution()
{
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        double sum = b[i];
        for (int j = i + 1; j < n; j++)
        {
            sum -= A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }
}

int main()
{
    cout << "输入矩阵维度: ";
    cin >> n;
    cout << "输入线程数 (建议 <= " << MAX_THREADS << "): ";
    cin >> num_threads;
    if (num_threads > MAX_THREADS)
        num_threads = MAX_THREADS;

    A = new double *[n];
    for (int i = 0; i < n; ++i)
    {
        A[i] = new double[n];
    }
    b = new double[n];
    x = new double[n];

    initialize_matrix();

    pthread_barrier_init(&barrier_Division, NULL, num_threads);
    pthread_barrier_init(&barrier_Elimination, NULL, num_threads);

    pthread_t handles[MAX_THREADS];
    ThreadData param[MAX_THREADS];
    for (int t_id = 0; t_id < num_threads; ++t_id)
    {
        param[t_id].thread_id = t_id;
        pthread_create(&handles[t_id], NULL, threadFunc, (void *)&param[t_id]);
    }

    auto start = high_resolution_clock::now();

    for (int t_id = 0; t_id < num_threads; ++t_id)
    {
        pthread_join(handles[t_id], NULL);
    }

    back_substitution();

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "（行/列划分）静态线程+barrier同步优化，用时：" << duration.count() << "us" << endl;

    pthread_barrier_destroy(&barrier_Division);
    pthread_barrier_destroy(&barrier_Elimination);

    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}
