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
        A[i][i] = rand() % 100 + 1000;
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

        for (int i = 0; i < n; ++i)
            MPI_Bcast(A[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(b.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        double start_time = MPI_Wtime();

        for (int k = 0; k < n; ++k)
        {
            double pivot;
            double bk;

            // 拥有 A[k][k] 的进程广播主元
            if (k % numprocs == rank)
            {
                pivot = A[k][k];
                bk = b[k];
            }
            MPI_Bcast(&pivot, 1, MPI_DOUBLE, k % numprocs, MPI_COMM_WORLD);
            MPI_Bcast(&bk, 1, MPI_DOUBLE, k % numprocs, MPI_COMM_WORLD);

            // 所有进程更新各自负责的列
            for (int j = k + 1; j < n; ++j)
            {
                if (j % numprocs == rank)
                {
                    A[k][j] /= pivot;
                }
            }
            b[k] = bk / pivot;
            A[k][k] = 1.0;

            // 更新后续行（本地处理）
            for (int i = k + 1; i < n; ++i)
            {
                double factor = A[i][k];
                for (int j = k + 1; j < n; ++j)
                {
                    if (j % numprocs == rank)
                    {
                        A[i][j] -= factor * A[k][j];
                    }
                }
                b[i] -= factor * b[k];
                A[i][k] = 0.0;
            }
        }

        // 回代（0号进程统一进行）
        if (rank == 0)
        {
            for (int i = n - 1; i >= 0; --i)
            {
                x[i] = b[i];
                for (int j = i + 1; j < n; ++j)
                    x[i] -= A[i][j] * x[j];
            }

            double end_time = MPI_Wtime();
            double duration = (end_time - start_time) * 1e3;
            total_time += duration;
            cout << r + 1 << "用时: " << duration << " ms" << endl;
        }
    }

    if (rank == 0)
    {
        cout << "平均运行时间: " << total_time / runs << " ms" << endl;
    }

    MPI_Finalize();
    return 0;
}
