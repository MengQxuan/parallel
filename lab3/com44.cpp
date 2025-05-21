// 动态线程
#include <iostream>
#include <pthread.h>
#include <chrono>
#include <cstdlib>
#include <ctime>
using namespace std;
using namespace std::chrono;

struct ThreadParam
{
    int k;
    int t_id;
};

double **A, *b, *x;
int n;

void initialize_matrix(double **A, double *b, int n)
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
    ThreadParam *p = (ThreadParam *)param;
    int k = p->k;
    int t_id = p->t_id;
    int i = k + t_id + 1;

    if (i >= n)
        pthread_exit(NULL);

    double factor = A[i][k];
    for (int j = k + 1; j < n; ++j)
    {
        A[i][j] -= factor * A[k][j];
    }
    A[i][k] = 0.0;
    b[i] -= factor * b[k];

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

    A = new double *[n];
    for (int i = 0; i < n; ++i)
    {
        A[i] = new double[n];
    }
    b = new double[n];
    x = new double[n];

    initialize_matrix(A, b, n);

    auto start = high_resolution_clock::now();

    for (int k = 0; k < n; ++k)
    {
        for (int j = k + 1; j < n; ++j)
        {
            A[k][j] /= A[k][k];
        }
        b[k] /= A[k][k];
        A[k][k] = 1.0;

        int worker_count = n - 1 - k;
        pthread_t *handles = new pthread_t[worker_count];
        ThreadParam *params = new ThreadParam[worker_count];

        for (int t_id = 0; t_id < worker_count; ++t_id)
        {
            params[t_id].k = k;
            params[t_id].t_id = t_id;
            pthread_create(&handles[t_id], NULL, threadFunc, (void *)&params[t_id]);
        }

        for (int t_id = 0; t_id < worker_count; ++t_id)
        {
            pthread_join(handles[t_id], NULL);
        }

        delete[] handles;
        delete[] params;
    }

    back_substitution();

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "动态线程，用时：" << duration.count() << "us" << endl;

    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}
