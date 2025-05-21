// 静态线程+信号量同步优化
#include <iostream>
#include <pthread.h>
#include <semaphore.h>
#include <chrono>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace std;
using namespace std::chrono;

#define MAX_THREADS 8

// 全局变量
double **A, *b, *x;
int n;
int num_threads;

sem_t sem_main;
sem_t sem_workerstart[MAX_THREADS];
sem_t sem_workerend[MAX_THREADS];

struct ThreadData
{
    int thread_id;
};

// 初始化上三角矩阵并进行线性组合
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
        sem_wait(&sem_workerstart[t_id]);

        for (int i = k + 1 + t_id; i < n; i += num_threads)
        {
            double factor = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; ++j)
            {
                A[i][j] -= factor * A[k][j];
            }
            A[i][k] = 0.0;
            b[i] -= factor * b[k];
        }

        sem_post(&sem_main);
        sem_wait(&sem_workerend[t_id]);
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

    sem_init(&sem_main, 0, 0);
    for (int i = 0; i < num_threads; ++i)
    {
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
    }

    pthread_t handles[MAX_THREADS];
    ThreadData param[MAX_THREADS];
    for (int t_id = 0; t_id < num_threads; ++t_id)
    {
        param[t_id].thread_id = t_id;
        pthread_create(&handles[t_id], NULL, threadFunc, (void *)&param[t_id]);
    }

    auto start = high_resolution_clock::now();

    for (int k = 0; k < n; ++k)
    {
        for (int j = k + 1; j < n; ++j)
        {
            A[k][j] /= A[k][k];
        }
        b[k] /= A[k][k];
        A[k][k] = 1.0;

        for (int t_id = 0; t_id < num_threads; ++t_id)
        {
            sem_post(&sem_workerstart[t_id]);
        }

        for (int t_id = 0; t_id < num_threads; ++t_id)
        {
            sem_wait(&sem_main);
        }

        for (int t_id = 0; t_id < num_threads; ++t_id)
        {
            sem_post(&sem_workerend[t_id]);
        }
    }

    for (int t_id = 0; t_id < num_threads; ++t_id)
    {
        pthread_join(handles[t_id], NULL);
    }

    back_substitution();

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "静态线程+信号量同步优化，用时：" << duration.count() << "us" << endl;

    sem_destroy(&sem_main);
    for (int i = 0; i < num_threads; ++i)
    {
        sem_destroy(&sem_workerstart[i]);
        sem_destroy(&sem_workerend[i]);
    }

    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}
