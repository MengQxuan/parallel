// #include <iostream>
// #include <pthread.h>
// #include <chrono>
// #include <vector>
// using namespace std;
// using namespace std::chrono;

// #define MAX_THREADS 8

// // 全局变量
// pthread_barrier_t barrier;
// double **A, *b, *x;
// int n;
// int num_threads;

// struct ThreadData
// {
//     int thread_id;
// };

// void *gaussian_elimination_row(void *arg)
// {
//     ThreadData *data = (ThreadData *)arg;
//     int tid = data->thread_id;

//     for (int k = 0; k < n; k++)
//     {
//         pthread_barrier_wait(&barrier);

//         for (int i = k + 1 + tid; i < n; i += num_threads)
//         {
//             double factor = A[i][k] / A[k][k];
//             for (int j = k + 1; j < n; ++j)
//             {
//                 A[i][j] -= factor * A[k][j];
//             }
//             b[i] -= factor * b[k];
//         }

//         pthread_barrier_wait(&barrier);
//     }
//     pthread_exit(NULL);
//     return NULL;
// }

// void *gaussian_elimination_col(void *arg)
// {
//     ThreadData *data = (ThreadData *)arg;
//     int tid = data->thread_id;

//     for (int k = 0; k < n; k++)
//     {
//         pthread_barrier_wait(&barrier);

//         for (int j = k + 1 + tid; j < n; j += num_threads)
//         {
//             for (int i = k + 1; i < n; i++)
//             {
//                 double factor = A[i][k] / A[k][k];
//                 A[i][j] -= factor * A[k][j];
//             }
//         }

//         pthread_barrier_wait(&barrier);
//     }
//     pthread_exit(NULL);
//     return NULL;
// }

// void back_substitution()
// {
//     x[n - 1] = b[n - 1] / A[n - 1][n - 1];
//     for (int i = n - 2; i >= 0; i--)
//     {
//         double sum = b[i];
//         for (int j = i + 1; j < n; j++)
//         {
//             sum -= A[i][j] * x[j];
//         }
//         x[i] = sum / A[i][i];
//     }
// }

// int main()
// {
//     cout << "输入矩阵维度: ";
//     cin >> n;
//     cout << "输入线程数 (建议 <= " << MAX_THREADS << "): ";
//     cin >> num_threads;

//     if (num_threads > MAX_THREADS)
//         num_threads = MAX_THREADS;

//     A = new double *[n];
//     for (int i = 0; i < n; ++i)
//     {
//         A[i] = new double[n];
//     }
//     b = new double[n];
//     x = new double[n];

//     srand(time(0));
//     for (int i = 0; i < n; ++i)
//     {
//         fill(A[i], A[i] + n, 0.0); // 行初始化为0
//         A[i][i] = 1.0;             // 对角线置1
//         for (int j = i + 1; j < n; ++j)
//         {
//             A[i][j] = rand() % 100 + 1; // 生成1~100的随机数
//         }
//     }

//     pthread_t threads[MAX_THREADS];
//     ThreadData thread_data[MAX_THREADS];
//     pthread_barrier_init(&barrier, NULL, num_threads);

//     auto start = high_resolution_clock::now();

//     cout << "选择任务划分方式 (1: 按行划分, 2: 按列划分): ";
//     int partition_choice;
//     cin >> partition_choice;

//     for (int i = 0; i < num_threads; i++)
//     {
//         thread_data[i].thread_id = i;
//         if (partition_choice == 1)
//         {
//             pthread_create(&threads[i], NULL, gaussian_elimination_row, (void *)&thread_data[i]);
//         }
//         else
//         {
//             pthread_create(&threads[i], NULL, gaussian_elimination_col, (void *)&thread_data[i]);
//         }
//     }

//     for (int i = 0; i < num_threads; i++)
//     {
//         pthread_join(threads[i], NULL);
//     }

//     back_substitution();

//     auto end = high_resolution_clock::now();
//     auto duration = duration_cast<microseconds>(end - start);
//     cout << "静态线程高斯消元完成，用时：" << duration.count() << "us" << endl;

//     pthread_barrier_destroy(&barrier);

//     for (int i = 0; i < n; ++i)
//     {
//         delete[] A[i];
//     }
//     delete[] A;
//     delete[] b;
//     delete[] x;

//     return 0;
// }

#include <iostream>
#include <pthread.h>
#include <chrono>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace std;
using namespace std::chrono;

#define MAX_THREADS 8

// 全局变量
pthread_barrier_t barrier;
double **A, *b, *x;
int n;
int num_threads;

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
        A[i][i] = rand() % 100 + 1;
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

void *gaussian_elimination_row(void *arg)
{
    ThreadData *data = (ThreadData *)arg;
    int tid = data->thread_id;

    for (int k = 0; k < n; k++)
    {
        pthread_barrier_wait(&barrier);

        for (int i = k + 1 + tid; i < n; i += num_threads)
        {
            double factor = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; ++j)
            {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }

        pthread_barrier_wait(&barrier);
    }
    pthread_exit(NULL);
    return NULL;
}

void *gaussian_elimination_col(void *arg)
{
    ThreadData *data = (ThreadData *)arg;
    int tid = data->thread_id;

    for (int k = 0; k < n; k++)
    {
        pthread_barrier_wait(&barrier);

        for (int j = k + 1 + tid; j < n; j += num_threads)
        {
            for (int i = k + 1; i < n; i++)
            {
                double factor = A[i][k] / A[k][k];
                A[i][j] -= factor * A[k][j];
            }
        }

        pthread_barrier_wait(&barrier);
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
    pthread_barrier_init(&barrier, NULL, num_threads);

    auto start = high_resolution_clock::now();

    cout << "选择任务划分方式 (1: 按行划分, 2: 按列划分): ";
    int partition_choice;
    cin >> partition_choice;

    for (int i = 0; i < num_threads; i++)
    {
        thread_data[i].thread_id = i;
        if (partition_choice == 1)
        {
            pthread_create(&threads[i], NULL, gaussian_elimination_row, (void *)&thread_data[i]);
        }
        else
        {
            pthread_create(&threads[i], NULL, gaussian_elimination_col, (void *)&thread_data[i]);
        }
    }

    for (int i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    back_substitution();

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "静态线程高斯消元完成，用时：" << duration.count() << "us" << endl;

    pthread_barrier_destroy(&barrier);

    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}
