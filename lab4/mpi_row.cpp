#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace std;

void initialize_matrix(vector<vector<double>> &A, vector<double> &b, int n)
{
    srand(time(0));
    for (int i = 0; i < n; ++i)
    {
        fill(A[i].begin(), A[i].end(), 0.0);
        A[i][i] = rand() % 100 + 1000; // 避免除零
        for (int j = i + 1; j < n; ++j)
            A[i][j] = rand() % 100 + 1;
        b[i] = rand() % 100 + 1;
    }

    for (int k = 0; k < n / 2; ++k)
    {
        int row1 = rand() % n;
        int row2 = rand() % n;
        double factor = (rand() % 100) / 100.0;
        for (int j = 0; j < n; ++j)
            A[row1][j] += factor * A[row2][j];
        b[row1] += factor * b[row2];
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int n;
    if (rank == 0)
    {
        cout << "输入矩阵维度: ";
        cin >> n;
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    const int runs = 10;
    vector<vector<double>> A0(n, vector<double>(n));
    vector<double> b0(n);
    if (rank == 0)
        initialize_matrix(A0, b0, n);

    double total_time = 0.0;

    for (int r = 0; r < runs; ++r)
    {
        vector<vector<double>> A(n, vector<double>(n));
        vector<double> b(n);
        vector<double> x(n);

        if (rank == 0)
        {
            A = A0;
            b = b0;
        }

        // 广播矩阵和向量
        for (int i = 0; i < n; ++i)
            MPI_Bcast(A[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(b.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        int base_rows = n / numprocs;
        int extra = n % numprocs;
        int r1 = rank * base_rows + min(rank, extra);
        int rows = base_rows + (rank < extra ? 1 : 0);
        int r2 = r1 + rows - 1;

        MPI_Barrier(MPI_COMM_WORLD); // 保证同步
        double start_time = MPI_Wtime();

        // 消元过程
        for (int k = 0; k < n; ++k)
        {
            if (k >= r1 && k <= r2)
            {
                for (int j = k + 1; j < n; ++j)
                    A[k][j] /= A[k][k];
                b[k] /= A[k][k];
                A[k][k] = 1.0;

                for (int i = 0; i < numprocs; ++i)
                    if (i != rank)
                        MPI_Send(A[k].data(), n, MPI_DOUBLE, i, k, MPI_COMM_WORLD);
                for (int i = 0; i < numprocs; ++i)
                    if (i != rank)
                        MPI_Send(&b[k], 1, MPI_DOUBLE, i, n + k, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Recv(A[k].data(), n, MPI_DOUBLE, MPI_ANY_SOURCE, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&b[k], 1, MPI_DOUBLE, MPI_ANY_SOURCE, n + k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            for (int i = max(k + 1, r1); i <= r2; ++i)
            {
                double factor = A[i][k];
                for (int j = k + 1; j < n; ++j)
                    A[i][j] -= factor * A[k][j];
                b[i] -= factor * b[k];
                A[i][k] = 0.0;
            }
        }

        // 回代求解（仅 0 号进程）
        if (rank == 0)
        {
            for (int i = n - 1; i >= 0; --i)
            {
                x[i] = b[i];
                for (int j = i + 1; j < n; ++j)
                    x[i] -= A[i][j] * x[j];
            }

            double end_time = MPI_Wtime();
            double duration = (end_time - start_time) * 1e3; // 单位 ms
            total_time += duration;
            cout << r + 1 << "用时: " << duration << " ms" << endl;
        }
    }

    if (rank == 0)
    {
        cout << "平均时间: " << total_time / runs << " ms" << endl;
    }

    MPI_Finalize();
    return 0;
}
