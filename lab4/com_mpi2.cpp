#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <mpi.h>
using namespace std;

int find_owner(int k, int n, int num_procs, int base_rows, int remainder)
{
    if (k < remainder * (base_rows + 1))
        return k / (base_rows + 1);
    return remainder + (k - remainder * (base_rows + 1)) / base_rows;
}

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

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n;
    double *A_all = nullptr, *b_all = nullptr;
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

    // Scatter 分发矩阵
    double *A_flat = nullptr;
    double *b_vec = nullptr;
    if (rank == 0)
    {
        A_flat = new double[n * n];
        b_vec = new double[n];
        initialize_matrix(A_flat, b_vec, n);
    }

    int *sendcounts = new int[size], *displs = new int[size];
    int *b_counts = new int[size], *b_displs = new int[size];
    for (int i = 0; i < size; ++i)
    {
        int rows = base_rows + (i < remainder ? 1 : 0);
        sendcounts[i] = rows * n;
        b_counts[i] = rows;
        displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];
        b_displs[i] = (i == 0) ? 0 : b_displs[i - 1] + b_counts[i - 1];
    }

    double *local_A = new double[local_rows * n];
    double *local_b = new double[local_rows];
    MPI_Scatterv(A_flat, sendcounts, displs, MPI_DOUBLE, local_A, local_rows * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(b_vec, b_counts, b_displs, MPI_DOUBLE, local_b, local_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double *row_buf = new double[n + 1];
    const int num_runs = 10;
    double total_time = 0;

    for (int run = 0; run < num_runs; ++run)
    {
        double *A_run = new double[local_rows * n];
        double *b_run = new double[local_rows];
        copy(local_A, local_A + local_rows * n, A_run);
        copy(local_b, local_b + local_rows, b_run);

        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();

        for (int k = 0; k < n; ++k)
        {
            int owner = find_owner(k, n, size, base_rows, remainder);
            int local_k = k - start_row;

            if (rank == owner)
            {
                double pivot = A_run[local_k * n + k];
                for (int j = k + 1; j < n; ++j)
                    A_run[local_k * n + j] /= pivot;
                A_run[local_k * n + k] = 1.0;
                for (int j = 0; j < n; ++j)
                    row_buf[j] = A_run[local_k * n + j];
                row_buf[n] = b_run[local_k];
                if (rank + 1 < size)
                    MPI_Send(row_buf, n + 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Recv(row_buf, n + 1, MPI_DOUBLE, rank - 1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (rank + 1 < size)
                    MPI_Send(row_buf, n + 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD);
            }

            for (int i = 0; i < local_rows; ++i)
            {
                int global_i = start_row + i;
                if (global_i > k)
                {
                    double factor = A_run[i * n + k];
                    for (int j = k + 1; j < n; ++j)
                        A_run[i * n + j] -= factor * row_buf[j];
                    b_run[i] -= factor * row_buf[n];
                    A_run[i * n + k] = 0;
                }
            }
        }

        double *x = nullptr;
        if (rank == 0)
        {
            double *final_A = new double[n * n];
            double *final_b = new double[n];
            x = new double[n];

            MPI_Gatherv(A_run, local_rows * n, MPI_DOUBLE, final_A, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gatherv(b_run, local_rows, MPI_DOUBLE, final_b, b_counts, b_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            for (int i = n - 1; i >= 0; --i)
            {
                double sum = final_b[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= final_A[i * n + j] * x[j];
                x[i] = sum / final_A[i * n + i];
            }

            delete[] final_A;
            delete[] final_b;
            delete[] x;
        }
        else
        {
            MPI_Gatherv(A_run, local_rows * n, MPI_DOUBLE, nullptr, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gatherv(b_run, local_rows, MPI_DOUBLE, nullptr, b_counts, b_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        double end = MPI_Wtime();
        if (rank == 0)
        {
            total_time += end - start;
            cout << "第 " << run + 1 << " 次用时: " << (end - start) * 1000 << " ms" << endl;
        }

        delete[] A_run;
        delete[] b_run;
    }

    if (rank == 0)
    {
        cout << "平均用时: " << (total_time / num_runs) * 1000 << " ms" << endl;
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
        delete[] A_flat;
        delete[] b_vec;
    }

    MPI_Finalize();
    return 0;
}
