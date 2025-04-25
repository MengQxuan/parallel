//串行 不对齐
#include <iostream>
#include <chrono>
using namespace std;
using namespace std::chrono;

void gaussian_elimination(double **A, double *b, double *x, int n)
{
    // 消去过程
    for (int k = 0; k < n; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; ++j)
            {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // 回代过程
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
    int n;
    cout << "输入矩阵维度: ";
    cin >> n;

    // 堆上动态分配二维数组 A 和一维数组 b、x
    double **A = new double *[n];
    for (int i = 0; i < n; ++i)
    {
        A[i] = new double[n];
    }
    double *b = new double[n];
    double *x = new double[n];

    for (int i = 0; i < n; ++i)
    {
        b[i] = 1.0;
        for (int j = 0; j < n; ++j)
        {
            A[i][j] = 1.0;
        }
    }

    // // 初始化矩阵A（参考m_reset逻辑）
    // srand(time(0));
    // for (int i = 0; i < n; ++i)
    // {
    //     fill(A[i], A[i] + n, 0.0); // 行初始化为0
    //     A[i][i] = 1.0;             // 对角线置1
    //     for (int j = i + 1; j < n; ++j)
    //     {
    //         A[i][j] = rand() % 100 + 1; // 生成1~100的随机数
    //     }
    // }

    auto start = high_resolution_clock::now();
    gaussian_elimination(A, b, x, n);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "内存不对齐的普通高斯消元完成，用时：" << duration.count() << " 微秒" << endl;

    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;
    return 0;
}