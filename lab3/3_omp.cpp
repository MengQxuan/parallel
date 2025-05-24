// static dynamic guided
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <vector>
using namespace std;
using namespace std::chrono;

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

long long gaussian_elimination_schedule(double **A, double *b, double *x, int n, int num_threads, string schedule_type, int chunk_size)
{
    int i, j, k;
    omp_set_num_threads(num_threads);
    auto t1 = high_resolution_clock::now();

#pragma omp parallel shared(A, b, n) private(i, j, k)
    {
        for (k = 0; k < n; ++k)
        {
#pragma omp single
            {
                double pivot = A[k][k];
                for (j = k + 1; j < n; ++j)
                    A[k][j] /= pivot;
                b[k] /= pivot;
                A[k][k] = 1.0;
            }

            if (schedule_type == "static")
            {
#pragma omp for schedule(static)
                for (i = k + 1; i < n; ++i)
                {
                    double factor = A[i][k];
                    for (j = k + 1; j < n; ++j)
                        A[i][j] -= factor * A[k][j];
                    b[i] -= factor * b[k];
                    A[i][k] = 0.0;
                }
            }
            else if (schedule_type == "dynamic")
            {
#pragma omp for schedule(dynamic, chunk_size)
                for (i = k + 1; i < n; ++i)
                {
                    double factor = A[i][k];
                    for (j = k + 1; j < n; ++j)
                        A[i][j] -= factor * A[k][j];
                    b[i] -= factor * b[k];
                    A[i][k] = 0.0;
                }
            }
            else if (schedule_type == "guided")
            {
#pragma omp for schedule(guided, chunk_size)
                for (i = k + 1; i < n; ++i)
                {
                    double factor = A[i][k];
                    for (j = k + 1; j < n; ++j)
                        A[i][j] -= factor * A[k][j];
                    b[i] -= factor * b[k];
                    A[i][k] = 0.0;
                }
            }
        }
    }

    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        double sum = b[i];
        for (int j = i + 1; j < n; j++)
            sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
    }

    auto t2 = high_resolution_clock::now();
    return duration_cast<microseconds>(t2 - t1).count();
}

int main()
{
    int n = 512;

    int num_threads = 8;

    vector<int> chunk_sizes = {1, 4, 8, 16, 32, 64, 128, 256};
    vector<string> schedules = {"static", "dynamic", "guided"};

    double **A = new double *[n];
    double **backup = new double *[n];
    for (int i = 0; i < n; ++i)
    {
        A[i] = new double[n];
        backup[i] = new double[n];
    }

    double *b = new double[n];
    double *b_backup = new double[n];
    double *x = new double[n];

    // 初始化一次数据
    initialize_matrix(A, b, n);
    for (int i = 0; i < n; ++i)
        copy(A[i], A[i] + n, backup[i]);
    copy(b, b + n, b_backup);

    for (auto schedule : schedules)
    {
        // static只执行一次，chunk_size = 0
        if (schedule == "static")
        {
            for (int i = 0; i < n; ++i)
                copy(backup[i], backup[i] + n, A[i]);
            copy(b_backup, b_backup + n, b);

            long long t = gaussian_elimination_schedule(A, b, x, n, num_threads, schedule, 0);
            cout << "已完成: " << schedule << ", chunk=N/A, 用时=" << t << "us" << endl;
        }
        else
        {
            for (int chunk : chunk_sizes)
            {
                for (int i = 0; i < n; ++i)
                    copy(backup[i], backup[i] + n, A[i]);
                copy(b_backup, b_backup + n, b);

                long long t = gaussian_elimination_schedule(A, b, x, n, num_threads, schedule, chunk);
                cout << "已完成: " << schedule << ", chunk=" << chunk << ", 用时=" << t << "us" << endl;
            }
        }
    }

    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
        delete[] backup[i];
    }
    delete[] A;
    delete[] backup;
    delete[] b;
    delete[] b_backup;
    delete[] x;

    return 0;
}
