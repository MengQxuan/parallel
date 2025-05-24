// 串行算法的openmp优化
#include <iostream>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <omp.h>
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

void gaussian_elimination_openmp(double **A, double *b, double *x, int n, int num_threads)
{
    omp_set_num_threads(num_threads);
    int i, j, k;

#pragma omp parallel shared(A, b, n) default(none) private(i, j, k) firstprivate(num_threads)
    {
        for (int k = 0; k < n; ++k)
        {
#pragma omp single
            {
                double pivot = A[k][k];
                for (int j = k + 1; j < n; ++j)
                    A[k][j] /= pivot;
                b[k] /= pivot;
                A[k][k] = 1.0;
            }

#pragma omp for
            for (int i = k + 1; i < n; ++i)
            {
                double factor = A[i][k];
                for (int j = k + 1; j < n; ++j)
                    A[i][j] -= factor * A[k][j];
                b[i] -= factor * b[k];
                A[i][k] = 0.0;
            }
        }
    }

    // 回代阶段（串行）
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        double sum = b[i];
        for (int j = i + 1; j < n; j++)
            sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
    }
}

int main()
{
    int n;
    cout << "输入矩阵维度: ";
    cin >> n;

    int num_threads = 8; // 改线程数

    double **A = new double *[n];
    for (int i = 0; i < n; ++i)
        A[i] = new double[n];

    double *b = new double[n];
    double *x = new double[n];

    initialize_matrix(A, b, n);

    auto start = high_resolution_clock::now();
    gaussian_elimination_openmp(A, b, x, n, num_threads);
    auto end = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(end - start);
    cout << "OpenMP并行优化用时: " << duration.count() << " us" << endl;

    for (int i = 0; i < n; ++i)
        delete[] A[i];
    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}
