#include <iostream>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cuda_runtime.h>
using namespace std;

// 策略B：每个Block处理一行，每个线程处理该行的一列
__global__ void division_kernel_row_col(double *A, double *b, int k, int n)
{
    int row = k + 1 + blockIdx.x;  // 每个Block对应一行
    int col = k + 1 + threadIdx.x; // 每个线程处理该行的一个元素

    if (row < n && col < n)
    {
        double Aik = A[row * n + k];
        double Akk = A[k * n + k];
        A[row * n + col] -= Aik * A[k * n + col] / Akk;
    }

    __syncthreads(); // 块内同步

    // 线程0更新b和A[row][k]
    if (row < n && threadIdx.x == 0)
    {
        b[row] -= A[row * n + k] * b[k] / A[k * n + k];
        A[row * n + k] = 0.0;
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

// 反向替代
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

    double *A = new double[n * n];
    double *b = new double[n];
    double *x = new double[n];

    initialize_matrix(A, b, n);

    double *d_A, *d_b;
    cudaMalloc(&d_A, sizeof(double) * n * n);
    cudaMalloc(&d_b, sizeof(double) * n);

    cudaMemcpy(d_A, A, sizeof(double) * n * n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, sizeof(double) * n, cudaMemcpyHostToDevice);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float elapsedTime = 0.0;

    cudaEventRecord(start, 0);

    for (int k = 0; k < n - 1; ++k)
    {
        dim3 threadsPerBlock(256);
        dim3 blocksPerGrid(n - k - 1); // 每个block对应一行
        division_kernel_row_col<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_b, k, n);
        cudaDeviceSynchronize();
    }

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cout << "[策略B] GPU高斯消元完成，用时：" << elapsedTime << " ms" << endl;

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    cudaMemcpy(A, d_A, sizeof(double) * n * n, cudaMemcpyDeviceToHost);
    cudaMemcpy(b, d_b, sizeof(double) * n, cudaMemcpyDeviceToHost);
    back_substitution(A, b, x, n);

    cudaFree(d_A);
    cudaFree(d_b);
    delete[] A;
    delete[] b;
    delete[] x;
    return 0;
}
