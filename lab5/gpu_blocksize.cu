// #include <iostream>
// #include <chrono>
// #include <cstdlib>
// #include <ctime>
// #include <cuda_runtime.h>
// using namespace std;

// // 除法核函数  每个线程处理一整行
// __global__ void division_kernel(double* A, double* b, int k, int n) {
//     int tid = blockIdx.x * blockDim.x + threadIdx.x;
//     if (tid > k && tid < n) {
//         double factor = A[tid * n + k] / A[k * n + k];
//         for (int j = k + 1; j < n; ++j) {
//             A[tid * n + j] -= factor * A[k * n + j];
//         }
//         b[tid] -= factor * b[k];
//         A[tid * n + k] = 0.0;
//     }
// }

// // 初始化矩阵
// void initialize_matrix(double* A, double* b, int n) {
//     srand(time(0));
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             A[i * n + j] = rand() % 100 + 1;
//         }
//         b[i] = rand() % 100 + 1;
//     }
// }

// int main() {
//     int n;
//     cout << "输入矩阵维度: ";
//     cin >> n;

//     // 尝试不同的 block size 配置
//     int block_sizes[] = {64, 128, 256, 512, 1024};
//     int num_tests = sizeof(block_sizes) / sizeof(int);

//     for (int t = 0; t < num_tests; ++t) {
//         int threadsPerBlock = block_sizes[t];

//         // 主机内存分配
//         double* A = new double[n * n];
//         double* b = new double[n];
//         initialize_matrix(A, b, n);

//         // 设备内存分配
//         double *d_A, *d_b;
//         cudaMalloc(&d_A, sizeof(double) * n * n);
//         cudaMalloc(&d_b, sizeof(double) * n);

//         cudaMemcpy(d_A, A, sizeof(double) * n * n, cudaMemcpyHostToDevice);
//         cudaMemcpy(d_b, b, sizeof(double) * n, cudaMemcpyHostToDevice);

//         // 初始化计时器
//         cudaEvent_t start, stop;
//         cudaEventCreate(&start);
//         cudaEventCreate(&stop);
//         float elapsedTime = 0.0;

//         cudaEventRecord(start, 0);

//         // GPU 高斯消元主循环
//         for (int k = 0; k < n - 1; ++k) {
//             int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
//             division_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_b, k, n);
//             cudaDeviceSynchronize();
//         }

//         cudaEventRecord(stop, 0);
//         cudaEventSynchronize(stop);
//         cudaEventElapsedTime(&elapsedTime, start, stop);

//         cout << "[block size = " << threadsPerBlock << "] GPU高斯消元用时：" << elapsedTime << " ms" << endl;

//         // 释放资源
//         cudaEventDestroy(start);
//         cudaEventDestroy(stop);
//         cudaFree(d_A);
//         cudaFree(d_b);
//         delete[] A;
//         delete[] b;
//     }

//     return 0;
// }

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cuda_runtime.h>
using namespace std;

// 策略B：每个Block处理一行，每个线程处理该行的一列
__global__ void division_kernel_row_col(double *A, double *b, int k, int n)
{
    int row = k + 1 + blockIdx.x;
    int col = k + 1 + threadIdx.x;

    if (row < n && col < n)
    {
        double Aik = A[row * n + k];
        double Akk = A[k * n + k];
        A[row * n + col] -= Aik * A[k * n + col] / Akk;
    }

    __syncthreads(); // 块内同步

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

    int block_sizes[] = {64, 128, 256, 512, 1024};
    int num_tests = sizeof(block_sizes) / sizeof(int);

    for (int t = 0; t < num_tests; ++t)
    {
        int threadsPerBlock = block_sizes[t];

        // 主机内存
        double *A = new double[n * n];
        double *b = new double[n];
        double *x = new double[n];
        initialize_matrix(A, b, n);

        // 设备内存
        double *d_A, *d_b;
        cudaMalloc(&d_A, sizeof(double) * n * n);
        cudaMalloc(&d_b, sizeof(double) * n);
        cudaMemcpy(d_A, A, sizeof(double) * n * n, cudaMemcpyHostToDevice);
        cudaMemcpy(d_b, b, sizeof(double) * n, cudaMemcpyHostToDevice);

        // 计时器
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        float elapsedTime = 0.0;
        cudaEventRecord(start, 0);

        // 消元主循环
        for (int k = 0; k < n - 1; ++k)
        {
            int needed_cols = n - (k + 1); // 每行剩余列数
            int actual_threads = min(threadsPerBlock, needed_cols);
            dim3 blockDim(actual_threads);
            dim3 gridDim(n - k - 1); // 每行一个 block
            division_kernel_row_col<<<gridDim, blockDim>>>(d_A, d_b, k, n);
            cudaDeviceSynchronize();
        }

        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&elapsedTime, start, stop);
        cout << "[block size = " << threadsPerBlock << "] GPU高斯消元用时：" << elapsedTime << " ms" << endl;

        // 反代
        cudaMemcpy(A, d_A, sizeof(double) * n * n, cudaMemcpyDeviceToHost);
        cudaMemcpy(b, d_b, sizeof(double) * n, cudaMemcpyDeviceToHost);
        back_substitution(A, b, x, n);

        // 清理资源
        cudaFree(d_A);
        cudaFree(d_b);
        delete[] A;
        delete[] b;
        delete[] x;
        cudaEventDestroy(start);
        cudaEventDestroy(stop);
    }
    return 0;
}
