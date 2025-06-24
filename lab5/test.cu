// test.cu
#include <iostream>
#include <vector>

// CUDA核函数：向量加法
__global__ void addVectors(float *A, float *B, float *C, int N)
{
    // 一个GPU线程负责一个元素的加法
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N)
    {
        C[idx] = A[idx] + B[idx];
    }
}

int main()
{
    int N = 1 << 20; // 向量大小，2^20 = 1048576
    size_t size = N * sizeof(float);

    // 主机（CPU）端内存
    std::vector<float> h_A(N);
    std::vector<float> h_B(N);
    std::vector<float> h_C(N);

    // 初始化主机端数据
    for (int i = 0; i < N; ++i)
    {
        h_A[i] = static_cast<float>(i);
        h_B[i] = static_cast<float>(i * 2);
    }

    // 设备（GPU）端内存指针
    float *d_A, *d_B, *d_C;

    // 分配设备端内存
    cudaMalloc(&d_A, size);
    cudaMalloc(&d_B, size);
    cudaMalloc(&d_C, size);

    // 将主机端数据拷贝到设备端
    cudaMemcpy(d_A, h_A.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B.data(), size, cudaMemcpyHostToDevice);

    // 定义CUDA网格和块的维度：
    // 因为核函数中只计算一个元素，所以所有维度的线程总数需要等于N
    int threadsPerBlock = 256; // // 一个线程块最大线程数由GPU硬件决定，一般不超过1024
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

    // 调用CUDA核函数
    addVectors<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, N);

    // 将设备端结果拷贝回主机端
    // cudaMemcpy默认是隐式同步的，所以这里可以不需要cudaDeviceSynchronize()
    cudaMemcpy(h_C.data(), d_C, size, cudaMemcpyDeviceToHost);

    // 验证结果
    bool success = true;
    for (int i = 0; i < N; ++i)
    {
        if (h_C[i] != h_A[i] + h_B[i])
        {
            success = false;
            break;
        }
    }

    if (success)
    {
        std::cout << "向量加法成功！" << std::endl;
    }
    else
    {
        std::cout << "向量加法失败！" << std::endl;
    }

    // 释放设备端内存
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    return 0;
}