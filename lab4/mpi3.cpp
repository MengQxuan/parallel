// 本文件包含三种 MPI 优化通信方式的高斯消元程序：
// 1. 阻塞通信 blocking version (blocking_gauss())
// 2. 非阻塞双边通信 non-blocking version (nonblocking_gauss())
// 3. 单边通信 one-sided version (onesided_gauss())
#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
using namespace std;
void initialize_matrix(double *A, double *b, int n)
{
    srand(time(0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            A[i * n + j] = 0.0;
        A[i * n + i] = rand() % 100 + 1000;
        for (int j = i + 1; j < n; ++j)
            A[i * n + j] = rand() % 100 + 1;
        b[i] = rand() % 100 + 1;
    }
    for (int k = 0; k < n / 2; ++k)
    {
        int r1 = rand() % n, r2 = rand() % n;
        double factor = (rand() % 100) / 100.0;
        for (int j = 0; j < n; ++j)
            A[r1 * n + j] += factor * A[r2 * n + j];
        b[r1] += factor * b[r2];
    }
}

int find_owner(int k, int n, int num_procs, int base_rows, int remainder)
{
    return (k < remainder * (base_rows + 1)) ? k / (base_rows + 1)
                                             : remainder + (k - remainder * (base_rows + 1)) / base_rows;
}

void blocking_gauss(int n, int rank, int size, double *local_A, double *local_b, int base_rows, int remainder, int start_row, int local_rows)
{
    double *row_buf = new double[n + 1];
    for (int k = 0; k < n; ++k)
    {
        int owner = find_owner(k, n, size, base_rows, remainder);
        int local_k = k - start_row;
        if (rank == owner)
        {
            double pivot = local_A[local_k * n + k];
            for (int j = k + 1; j < n; ++j)
                local_A[local_k * n + j] /= pivot;
            local_A[local_k * n + k] = 1.0;
            for (int j = 0; j < n; ++j)
                row_buf[j] = local_A[local_k * n + j];
            row_buf[n] = local_b[local_k];
            for (int dest = owner + 1; dest < size; dest++)
                MPI_Send(row_buf, n + 1, MPI_DOUBLE, dest, k, MPI_COMM_WORLD);
        }
        else if (rank > owner)
        {
            MPI_Recv(row_buf, n + 1, MPI_DOUBLE, owner, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (int i = 0; i < local_rows; ++i)
        {
            int global_i = start_row + i;
            if (global_i > k)
            {
                double factor = local_A[i * n + k];
                for (int j = k + 1; j < n; ++j)
                    local_A[i * n + j] -= factor * row_buf[j];
                local_b[i] -= factor * row_buf[n];
                local_A[i * n + k] = 0.0;
            }
        }
    }
    delete[] row_buf;
}

void nonblocking_gauss(int n, int rank, int size, double *local_A, double *local_b, int base_rows, int remainder, int start_row, int local_rows)
{
    double *row_buf = new double[n + 1];
    MPI_Request request;
    for (int k = 0; k < n; ++k)
    {
        int owner = find_owner(k, n, size, base_rows, remainder);
        int local_k = k - start_row;
        if (rank == owner)
        {
            double pivot = local_A[local_k * n + k];
            for (int j = k + 1; j < n; ++j)
                local_A[local_k * n + j] /= pivot;
            local_A[local_k * n + k] = 1.0;
            for (int j = 0; j < n; ++j)
                row_buf[j] = local_A[local_k * n + j];
            row_buf[n] = local_b[local_k];
            if (rank + 1 < size)
                MPI_Isend(row_buf, n + 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD, &request);
        }
        else if (rank > owner)
        {
            MPI_Irecv(row_buf, n + 1, MPI_DOUBLE, rank - 1, k, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
            if (rank + 1 < size)
                MPI_Isend(row_buf, n + 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD, &request);
        }
        for (int i = 0; i < local_rows; ++i)
        {
            int global_i = start_row + i;
            if (global_i > k)
            {
                double factor = local_A[i * n + k];
                for (int j = k + 1; j < n; ++j)
                    local_A[i * n + j] -= factor * row_buf[j];
                local_b[i] -= factor * row_buf[n];
                local_A[i * n + k] = 0.0;
            }
        }
    }
    delete[] row_buf;
}

void onesided_gauss(int n, int rank, int size, double *local_A, double *local_b,
                    int base_rows, int remainder, int start_row, int local_rows)
{
    double *row_buf = new double[n + 1];
    MPI_Win win;
    MPI_Win_create(row_buf, (n + 1) * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    for (int k = 0; k < n; ++k)
    {
        int owner = find_owner(k, n, size, base_rows, remainder);
        int local_k = k - start_row;
        MPI_Win_fence(0, win);
        if (rank == owner)
        {
            double pivot = local_A[local_k * n + k];
            for (int j = k + 1; j < n; ++j)
                local_A[local_k * n + j] /= pivot;
            local_A[local_k * n + k] = 1.0;
            for (int j = 0; j < n; ++j)
                row_buf[j] = local_A[local_k * n + j];
            row_buf[n] = local_b[local_k];
        }
        MPI_Win_fence(0, win);
        if (rank != owner)
        {
            MPI_Get(row_buf, n + 1, MPI_DOUBLE, owner, 0, n + 1, MPI_DOUBLE, win);
        }
        MPI_Win_fence(0, win);

        for (int i = 0; i < local_rows; ++i)
        {
            int global_i = start_row + i;
            if (global_i > k)
            {
                double factor = local_A[i * n + k];
                for (int j = k + 1; j < n; ++j)
                    local_A[i * n + j] -= factor * row_buf[j];
                local_b[i] -= factor * row_buf[n];
                local_A[i * n + k] = 0.0;
            }
        }
    }

    MPI_Win_free(&win);
    delete[] row_buf;
}

extern void initialize_matrix(double *A, double *b, int n);
extern int find_owner(int k, int n, int num_procs, int base_rows, int remainder);
extern void blocking_gauss(int n, int rank, int size, double *local_A, double *local_b, int base_rows, int remainder, int start_row, int local_rows);
extern void nonblocking_gauss(int n, int rank, int size, double *local_A, double *local_b, int base_rows, int remainder, int start_row, int local_rows);
extern void onesided_gauss(int n, int rank, int size, double *local_A, double *local_b, int base_rows, int remainder, int start_row, int local_rows);

void run_and_time(void (*gauss_func)(int, int, int, double *, double *, int, int, int, int),
                  const char *name, int n, int rank, int size,
                  double *A_full, double *b_full, int base_rows, int remainder,
                  int start_row, int local_rows)
{
    double total_time = 0.0;
    for (int run = 0; run < 10; ++run)
    {
        double *local_A = new double[local_rows * n];
        double *local_b = new double[local_rows];

        // 分发 A 和 b 到每个进程
        for (int i = 0; i < local_rows; ++i)
        {
            int global_row = start_row + i;
            memcpy(&local_A[i * n], &A_full[global_row * n], n * sizeof(double));
            local_b[i] = b_full[global_row];
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();

        gauss_func(n, rank, size, local_A, local_b, base_rows, remainder, start_row, local_rows);

        MPI_Barrier(MPI_COMM_WORLD);
        double end = MPI_Wtime();

        if (rank == 0)
        {
            total_time += end - start;
        }

        delete[] local_A;
        delete[] local_b;
    }
    if (rank == 0)
    {
        cout << name << " 平均用时: " << (total_time / 10) * 1000 << " ms" << endl;
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n;
    if (rank == 0)
    {
        cout << "输入矩阵维度: ";
        cin >> n;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int base_rows = n / size;
    int remainder = n % size;
    int local_rows = base_rows + (rank < remainder ? 1 : 0);
    int start_row = rank * base_rows + min(rank, remainder);

    double *A_full = new double[n * n];
    double *b_full = new double[n];

    if (rank == 0)
    {
        initialize_matrix(A_full, b_full, n);
    }
    MPI_Bcast(A_full, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b_full, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    run_and_time(blocking_gauss, "阻塞通信", n, rank, size, A_full, b_full, base_rows, remainder, start_row, local_rows);
    run_and_time(nonblocking_gauss, "非阻塞通信", n, rank, size, A_full, b_full, base_rows, remainder, start_row, local_rows);
    run_and_time(onesided_gauss, "单边通信", n, rank, size, A_full, b_full, base_rows, remainder, start_row, local_rows);

    delete[] A_full;
    delete[] b_full;
    MPI_Finalize();
    return 0;
}
