#include <iostream>
#include <vector>
#include <chrono>
using namespace std;
using namespace std::chrono;

void gaussian_elimination(vector<vector<double>> &A, vector<double> &b, vector<double> &x)
{
    int n = A.size();

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
    cout << "请输入矩阵的维度：";
    cin >> n;
    // 初始化矩阵 A 和向量 b，所有元素设为 1
    vector<vector<double>> A(n, vector<double>(n, 1.0));
    vector<double> b(n, 1.0);
    vector<double> x(n, 0.0); // 解向量

    auto start = high_resolution_clock::now();
    gaussian_elimination(A, b, x);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "普通高斯消元完成，用时：" << duration.count() << " 微秒" << endl;
    return 0;
}
