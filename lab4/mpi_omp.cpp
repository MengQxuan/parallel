#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <mpi.h>
#include <omp.h>
using namespace std;

// 查找行k所属的进程
int find_owner(int k, int n, int num_procs, int base_rows, int remainder)
{
    if (k < remainder * (base_rows + 1))
    {
        return k / (base_rows + 1);
    }
    return remainder + (k - remainder * (base_rows + 1)) / base_rows;
}

// 初始化矩阵（仅在0号进程执行）
void initialize_matrix(double **A, double *b, int n)
{
    srand(time(0));
    for (int i = 0; i < n; ++i)
    {
        fill(A[i], A[i] + n, 0.0);
        A[i][i] = rand() % 100 + 1000; // 对角线元素大值
        for (int j = i + 1; j < n; ++j)
        {
            A[i][j] = rand() % 100 + 1;
        }
        b[i] = rand() % 100 + 1;
    }

    // 随机线性组合增强数值稳定性
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

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int n = 0;
    double **A = nullptr, *b = nullptr;
    double *A_1D = nullptr; // 矩阵的连续存储

    // 设置 OpenMP 线程数
    omp_set_num_threads(8);

    // 0号进程读取矩阵维度
    if (rank == 0)
    {
        cout << "输入矩阵维度: ";
        cin >> n;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 计算负载均衡的行分配
    int base_rows = n / num_procs;
    int remainder = n % num_procs;
    int local_rows = base_rows + (rank < remainder ? 1 : 0);
    int start_row = rank * base_rows + min(rank, remainder);
    int end_row = start_row + local_rows - 1;

    // 0号进程初始化全局矩阵
    if (rank == 0)
    {
        A = new double *[n];
        b = new double[n];
        for (int i = 0; i < n; ++i)
        {
            A[i] = new double[n];
        }
        initialize_matrix(A, b, n);

        // 将矩阵转换为连续存储
        A_1D = new double[n * n];
        for (int i = 0; i < n; ++i)
        {
            copy(A[i], A[i] + n, A_1D + i * n);
        }
    }

    // 分配本地存储
    double *local_A = new double[local_rows * n];
    double *local_b = new double[local_rows];

    // 准备Scatterv参数
    int *sendcounts = new int[num_procs];
    int *displs = new int[num_procs];
    for (int i = 0; i < num_procs; ++i)
    {
        int i_rows = base_rows + (i < remainder ? 1 : 0);
        sendcounts[i] = i_rows * n;
        displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];
    }

    // 分发矩阵数据
    MPI_Scatterv(A_1D, sendcounts, displs, MPI_DOUBLE,
                 local_A, local_rows * n, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    // 准备b向量的分发参数
    int *b_counts = new int[num_procs];
    int *b_displs = new int[num_procs];
    for (int i = 0; i < num_procs; ++i)
    {
        b_counts[i] = base_rows + (i < remainder ? 1 : 0);
        b_displs[i] = (i == 0) ? 0 : b_displs[i - 1] + b_counts[i - 1];
    }

    // 分发b向量
    MPI_Scatterv(b, b_counts, b_displs, MPI_DOUBLE,
                 local_b, local_rows, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    // 分配行缓冲区
    double *row_buf = new double[n + 1]; // n个矩阵元素 + 1个向量元素

    const int num_runs = 10; // 运行10次取平均
    double total_time = 0.0;

    for (int run = 0; run < num_runs; ++run)
    {
        // 复制初始数据（每次运行前重置）
        double *run_local_A = new double[local_rows * n];
        double *run_local_b = new double[local_rows];
        copy(local_A, local_A + local_rows * n, run_local_A);
        copy(local_b, local_b + local_rows, run_local_b);

        MPI_Barrier(MPI_COMM_WORLD);
        double start_time = MPI_Wtime();

        // 主消元循环
        for (int k = 0; k < n; k++)
        {
            int owner = find_owner(k, n, num_procs, base_rows, remainder);

            // 行拥有者处理除法
            if (rank == owner)
            {
                int local_k = k - start_row;
                double pivot = run_local_A[local_k * n + k];

                // 除法操作
                for (int j = k + 1; j < n; j++)
                {
                    run_local_A[local_k * n + j] /= pivot;
                }
                run_local_A[local_k * n + k] = 1.0;

                // 打包行数据（矩阵行 + b元素）
                copy(run_local_A + local_k * n, run_local_A + local_k * n + n, row_buf);
                row_buf[n] = run_local_b[local_k];
            }

            // 广播行数据给需要它的进程
            if (rank >= owner)
            {
                if (rank == owner)
                {
                    // 发送给所有更高排名的进程
                    for (int dest = owner + 1; dest < num_procs; dest++)
                    {
                        MPI_Send(row_buf, n + 1, MPI_DOUBLE, dest, k, MPI_COMM_WORLD);
                    }
                }
                else
                {
                    MPI_Recv(row_buf, n + 1, MPI_DOUBLE, owner, k,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
#pragma omp parallel for
            // 本地消元操作
            for (int i = 0; i < local_rows; i++)
            {
                int global_i = start_row + i;
                if (global_i > k)
                {
                    // 获取消元因子
                    double factor = run_local_A[i * n + k];

                    // 消元操作
                    for (int j = k + 1; j < n; j++)
                    {
                        run_local_A[i * n + j] -= factor * row_buf[j];
                    }
                    run_local_b[i] -= factor * row_buf[n];

                    // 设置已处理元素为0
                    run_local_A[i * n + k] = 0.0;
                }
            }
        }

        // 收集结果到0号进程
        double *gathered_A = nullptr;
        double *gathered_b = nullptr;

        if (rank == 0)
        {
            gathered_A = new double[n * n];
            gathered_b = new double[n];
        }

        // 收集矩阵A
        MPI_Gatherv(run_local_A, local_rows * n, MPI_DOUBLE,
                    gathered_A, sendcounts, displs, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        // 收集向量b
        MPI_Gatherv(run_local_b, local_rows, MPI_DOUBLE,
                    gathered_b, b_counts, b_displs, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        // 回代求解（仅在0号进程执行）
        double *x = nullptr;
        if (rank == 0)
        {
            x = new double[n];
            // 回代过程
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = gathered_b[i];
                for (int j = i + 1; j < n; j++)
                {
                    sum -= gathered_A[i * n + j] * x[j];
                }
                x[i] = sum / gathered_A[i * n + i];
            }

            // 释放内存
            delete[] x;
            delete[] gathered_A;
            delete[] gathered_b;
        }

        double end_time = MPI_Wtime();
        double duration = end_time - start_time;

        if (rank == 0)
        {
            total_time += duration;
            cout << run + 1 << "用时: " << duration * 1000 << " ms" << endl;
        }

        delete[] run_local_A;
        delete[] run_local_b;
    }

    if (rank == 0)
    {
        cout << "平均时间: " << (total_time / num_runs) * 1000 << " ms" << endl;
    }

    delete[] local_A;
    delete[] local_b;
    delete[] row_buf;
    delete[] sendcounts;
    delete[] displs;
    delete[] b_counts;
    delete[] b_displs;

    if (rank == 0)
    {
        for (int i = 0; i < n; ++i)
        {
            delete[] A[i];
        }
        delete[] A;
        delete[] b;
        delete[] A_1D;
    }

    MPI_Finalize();
    return 0;
}