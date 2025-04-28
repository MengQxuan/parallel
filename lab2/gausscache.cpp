#include <iostream>
#include <chrono>
#include <cmath>
using namespace std;
using namespace std::chrono;

void gaussian_elimination(double *A, double *b, double *x, int n)
{
    // 消去过程（调整循环顺序优化缓存）
    for (int k = 0; k < n; ++k)
    {
        const int pivot = k * n;
        const double pivot_val = A[pivot + k];

        for (int i = k + 1; i < n; ++i)
        {
            const int row = i * n;
            const double factor = A[row + k] / pivot_val;

            // 对j循环展开+连续内存访问
            int j = k + 1;
            for (; j <= n - 4; j += 4)
            { // 4路循环展开
                A[row + j] -= factor * A[pivot + j];
                A[row + j + 1] -= factor * A[pivot + j + 1];
                A[row + j + 2] -= factor * A[pivot + j + 2];
                A[row + j + 3] -= factor * A[pivot + j + 3];
            }
            // 处理剩余元素
            for (; j < n; ++j)
            {
                A[row + j] -= factor * A[pivot + j];
            }
            b[i] -= factor * b[k];
        }
    }

    // 回代过程（优化连续访问）
    x[n - 1] = b[n - 1] / A[(n - 1) * n + (n - 1)];
    for (int i = n - 2; i >= 0; --i)
    {
        const int row = i * n;
        double sum = b[i];

        int j = i + 1;
        for (; j <= n - 4; j += 4)
        { // 4路循环展开
            sum -= A[row + j] * x[j] + A[row + j + 1] * x[j + 1] + A[row + j + 2] * x[j + 2] + A[row + j + 3] * x[j + 3];
        }
        // 处理剩余元素
        for (; j < n; ++j)
        {
            sum -= A[row + j] * x[j];
        }
        x[i] = sum / A[row + i];
    }
}

int main()
{
    int n;
    cout << "输入矩阵维度: ";
    cin >> n;

    // 使用单个连续内存块存储矩阵（行优先）
    double *A = new double[n * n];
    double *b = new double[n];
    double *x = new double[n];

    // 初始化矩阵（对角线占优）
    for (int i = 0; i < n; ++i)
    {
        b[i] = 1.0;
        double row_sum = 0.0;
        for (int j = 0; j < n; ++j)
        {
            A[i * n + j] = (i == j) ? (n * 100.0) : (rand() % 10 + 1); // 确保对角线占优
            row_sum += abs(A[i * n + j]);
        }
        A[i * n + i] = row_sum * 1.1; // 强化对角线优势
    }

    auto start = high_resolution_clock::now();
    gaussian_elimination(A, b, x, n);
    auto end = high_resolution_clock::now();

    cout << "cache优化用时：" << duration_cast<microseconds>(end - start).count() << "us" << endl;

    delete[] A;
    delete[] b;
    delete[] x;
    return 0;
}