#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <malloc.h>
using namespace std;
using namespace std::chrono;

// 使用一维数组模拟二维矩阵，访问 A[i][j] == A[i * n + j]
void gaussian_elimination(double *A, double *b, double *x, int n)
{
    // 消去过程
    for (int k = 0; k < n; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = A[i * n + k] / A[k * n + k];
            for (int j = k + 1; j < n; ++j)
            {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    // 回代过程
    x[n - 1] = b[n - 1] / A[(n - 1) * n + (n - 1)];
    for (int i = n - 2; i >= 0; i--)
    {
        double sum = b[i];
        for (int j = i + 1; j < n; j++)
        {
            sum -= A[i * n + j] * x[j];
        }
        x[i] = sum / A[i * n + i];
    }
}

int main()
{
    int n;
    cout << "输入矩阵维度: ";
    cin >> n;

    // 每个double是8字节，要求32字节对齐，所以内存大小必须是32的倍数
    size_t align = 32;
    size_t matrix_size = n * n * sizeof(double);
    size_t vector_size = n * sizeof(double);

    // 对齐分配矩阵和向量
    double *A = static_cast<double *>(_aligned_malloc(matrix_size, align));
    double *b = static_cast<double *>(_aligned_malloc(vector_size, align));
    double *x = static_cast<double *>(_aligned_malloc(vector_size, align));

    if (!A || !b || !x)
    {
        cerr << "内存分配失败" << endl;
        return 1;
    }

    for (int i = 0; i < n; ++i)
    {
        b[i] = 1.0;
        for (int j = 0; j < n; ++j)
        {
            A[i * n + j] = 1.0;
        }
    }

    auto start = high_resolution_clock::now();
    gaussian_elimination(A, b, x, n);
    auto end = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(end - start);
    cout << "内存对齐的普通高斯消元完成，用时：" << duration.count() << " 微秒" << endl;

    _aligned_free(A);
    _aligned_free(b);
    _aligned_free(x);

    return 0;
}
