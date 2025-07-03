#include <iostream>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cuda_runtime.h>
using namespace std;

// 自定义 atomicAdd 支持 double 类型
__device__ double atomicAddDouble(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

// 策略B：一个Block负责一行，线程负责列元素
__global__ void division_kernel_row_col(double* A, double* b, int k, int n) {
    int row = k + 1 + blockIdx.x;
    int col = k + 1 + threadIdx.x;

    if (row < n && col < n) {
        double Aik = A[row * n + k];
        double Akk = A[k * n + k];
        A[row * n + col] -= Aik * A[k * n + col] / Akk;
    }

    __syncthreads();

    if (row < n && threadIdx.x == 0) {
        b[row] -= A[row * n + k] * b[k] / A[k * n + k];
        A[row * n + k] = 0.0;
    }
}

// 回代阶段：每次求解 x[k]
__global__ void back_substitution_kernel(double* A, double* b, double* x, int k, int n) {
    int j = threadIdx.x;

    if (j > k && j < n) {
        atomicAddDouble(&b[k], -A[k * n + j] * x[j]);
    }

    __syncthreads();

    if (j == 0) {
        x[k] = b[k] / A[k * n + k];
    }
}

// 初始化矩阵
void initialize_matrix(double* A, double* b, int n) {
    srand(time(0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i * n + j] = rand() % 100 + 1;
        }
        b[i] = rand() % 100 + 1;
    }
}

int main() {
    int n;
    cout << "输入矩阵维度: ";
    cin >> n;

    double* A = new double[n * n];
    double* b = new double[n];
    double* x = new double[n];

    initialize_matrix(A, b, n);

    double *d_A, *d_b, *d_x;
    cudaMalloc(&d_A, sizeof(double) * n * n);
    cudaMalloc(&d_b, sizeof(double) * n);
    cudaMalloc(&d_x, sizeof(double) * n);

    cudaMemcpy(d_A, A, sizeof(double) * n * n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, sizeof(double) * n, cudaMemcpyHostToDevice);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float elapsedTime = 0.0;

    cudaEventRecord(start, 0);

    // 前向消元
    for (int k = 0; k < n - 1; ++k) {
        dim3 threadsPerBlock(256);
        dim3 blocksPerGrid(n - k - 1);
        division_kernel_row_col<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_b, k, n);
        cudaDeviceSynchronize();
    }

    // 回代求解（倒序）
    for (int k = n - 1; k >= 0; --k) {
        int threads = (n > 1024) ? 1024 : n; // 最大线程数限制
        back_substitution_kernel<<<1, threads>>>(d_A, d_b, d_x, k, n);
        cudaDeviceSynchronize();
    }

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cout << "[策略B + GPU回代] 完整GPU高斯消元完成，用时：" << elapsedTime << " ms" << endl;

    cudaMemcpy(x, d_x, sizeof(double) * n, cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_b);
    cudaFree(d_x);
    delete[] A;
    delete[] b;
    delete[] x;
    return 0;
}
