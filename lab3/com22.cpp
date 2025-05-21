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
pthread_barrier_t barrier;
sem_t sem_main, sem_worker;
double **A, *b, *x;
int n;
int num_threads;
int current_k = 0;

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

    // 线性组合操作
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

void *worker_thread(void *arg)
{
    ThreadData *data = (ThreadData *)arg;
    int tid = data->thread_id;

    for (int k = 0; k < n; k++)
    {
        sem_wait(&sem_worker);

        for (int i = k + 1 + tid; i < n; i += num_threads)
        {
            double factor = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; j++)
            {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }

        sem_post(&sem_main);
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

    pthread_t threads[MAX_THREADS];
    ThreadData thread_data[MAX_THREADS];

    sem_init(&sem_main, 0, 0);
    sem_init(&sem_worker, 0, 0);

    auto start = high_resolution_clock::now();

    for (int i = 0; i < num_threads; i++)
    {
        thread_data[i].thread_id = i;
        pthread_create(&threads[i], NULL, worker_thread, (void *)&thread_data[i]);
    }

    for (int k = 0; k < n; k++)
    {
        current_k = k;

        // 主线程执行除法操作
        double divisor = A[k][k];
        for (int j = k + 1; j < n; j++)
        {
            A[k][j] /= divisor;
        }
        b[k] /= divisor;

        // 唤醒工作线程进行消去操作
        for (int i = 0; i < num_threads; i++)
        {
            sem_post(&sem_worker);
        }

        // 等待工作线程完成消去操作
        for (int i = 0; i < num_threads; i++)
        {
            sem_wait(&sem_main);
        }
    }

    for (int i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    back_substitution();

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "信号量同步优化后的静态线程高斯消元完成，用时：" << duration.count() << "us" << endl;

    sem_destroy(&sem_main);
    sem_destroy(&sem_worker);

    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}
