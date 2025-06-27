#include <iostream>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cuda_runtime.h>
using namespace std;
using namespace std::chrono;

// 除法核函数
__global__ void division_kernel(double *A, double *b, int k, int n)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid > k && tid < n)
    {
        double factor = A[tid * n + k] / A[k * n + k];
        for (int j = k + 1; j < n; ++j)
        {
            A[tid * n + j] -= factor * A[k * n + j];
        }
        b[tid] -= factor * b[k];
        A[tid * n + k] = 0.0;
    }
}

// 初始化矩阵
void initialize_matrix(double *A, double *b, int n)
{
    srand(time(0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            A[i * n + j] = rand() % 100 + 1;
        }
        b[i] = rand() % 100 + 1;
    }
}

// 反向替代（在 CPU 上执行）
void back_substitution(double *A, double *b, double *x, int n)
{
    x[n - 1] = b[n - 1] / A[n * n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j)
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

    // 使用一维数组表示 A[n][n]
    double *A = new double[n * n];
    double *b = new double[n];
    double *x = new double[n];

    initialize_matrix(A, b, n);

    double *d_A, *d_b;
    cudaMalloc(&d_A, sizeof(double) * n * n);
    cudaMalloc(&d_b, sizeof(double) * n);

    cudaMemcpy(d_A, A, sizeof(double) * n * n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, sizeof(double) * n, cudaMemcpyHostToDevice);

    auto start = high_resolution_clock::now();

    for (int k = 0; k < n - 1; ++k)
    {
        int threadsPerBlock = 256;
        int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
        division_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_b, k, n);
        cudaDeviceSynchronize();
    }

    // 把 A 和 b 拷贝回 CPU 执行回代
    cudaMemcpy(A, d_A, sizeof(double) * n * n, cudaMemcpyDeviceToHost);
    cudaMemcpy(b, d_b, sizeof(double) * n, cudaMemcpyDeviceToHost);

    back_substitution(A, b, x, n);

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "GPU高斯消元完成，用时：" << duration.count() << "us" << endl;

    cudaFree(d_A);
    cudaFree(d_b);
    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}
