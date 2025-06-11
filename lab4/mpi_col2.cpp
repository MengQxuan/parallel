// #include <mpi.h>
// #include <iostream>
// #include <vector>
// #include <cstdlib>
// #include <ctime>
// using namespace std;

// void initialize_matrix(vector<vector<double>> &A, vector<double> &b, int n)
// {
//     srand(time(0));
//     for (int i = 0; i < n; ++i)
//     {
//         fill(A[i].begin(), A[i].end(), 0.0);
//         A[i][i] = rand() % 100 + 1000;
//         for (int j = i + 1; j < n; ++j)
//             A[i][j] = rand() % 100 + 1;
//         b[i] = rand() % 100 + 1;
//     }

//     for (int k = 0; k < n / 2; ++k)
//     {
//         int row1 = rand() % n;
//         int row2 = rand() % n;
//         double factor = (rand() % 100) / 100.0;
//         for (int j = 0; j < n; ++j)
//             A[row1][j] += factor * A[row2][j];
//         b[row1] += factor * b[row2];
//     }
// }

// int main(int argc, char *argv[])
// {
//     MPI_Init(&argc, &argv);
//     int rank, numprocs;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

//     int n;
//     if (rank == 0)
//     {
//         cout << "输入矩阵维度: ";
//         cin >> n;
//     }

//     MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

//     const int runs = 10;
//     vector<vector<double>> A0(n, vector<double>(n));
//     vector<double> b0(n);
//     if (rank == 0)
//         initialize_matrix(A0, b0, n);

//     double total_time = 0.0;

//     for (int r = 0; r < runs; ++r)
//     {
//         vector<vector<double>> A(n, vector<double>(n));
//         vector<double> b(n);
//         vector<double> x(n);

//         if (rank == 0)
//         {
//             A = A0;
//             b = b0;
//         }

//         for (int i = 0; i < n; ++i)
//             MPI_Bcast(A[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//         MPI_Bcast(b.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//         MPI_Barrier(MPI_COMM_WORLD);
//         double start_time = MPI_Wtime();

//         for (int k = 0; k < n; ++k)
//         {
//             double pivot = 0.0;
//             double bk = 0.0;

//             if (k % numprocs == rank)
//             {
//                 pivot = A[k][k];
//                 bk = b[k];
//             }

//             // 使用非阻塞广播
//             MPI_Request reqs[2];
//             MPI_Ibcast(&pivot, 1, MPI_DOUBLE, k % numprocs, MPI_COMM_WORLD, &reqs[0]);
//             MPI_Ibcast(&bk, 1, MPI_DOUBLE, k % numprocs, MPI_COMM_WORLD, &reqs[1]);

//             // 这里可以插入 overlap 操作（当前没有），等待广播完成
//             MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

//             // 所有进程更新自己负责的列
//             for (int j = k + 1; j < n; ++j)
//             {
//                 if (j % numprocs == rank)
//                 {
//                     A[k][j] /= pivot;
//                 }
//             }
//             b[k] = bk / pivot;
//             A[k][k] = 1.0;

//             // 本地更新后续行
//             for (int i = k + 1; i < n; ++i)
//             {
//                 double factor = A[i][k];
//                 for (int j = k + 1; j < n; ++j)
//                 {
//                     if (j % numprocs == rank)
//                     {
//                         A[i][j] -= factor * A[k][j];
//                     }
//                 }
//                 b[i] -= factor * b[k];
//                 A[i][k] = 0.0;
//             }
//         }

//         // 回代（只在0号进程）
//         if (rank == 0)
//         {
//             for (int i = n - 1; i >= 0; --i)
//             {
//                 x[i] = b[i];
//                 for (int j = i + 1; j < n; ++j)
//                     x[i] -= A[i][j] * x[j];
//             }

//             double end_time = MPI_Wtime();
//             double duration = (end_time - start_time) * 1e3;
//             total_time += duration;
//             cout << r + 1 << "用时: " << duration << " ms" << endl;
//         }
//     }

//     if (rank == 0)
//     {
//         cout << "平均运行时间: " << total_time / runs << " ms" << endl;
//     }

//     MPI_Finalize();
//     return 0;
// }

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
            double pivot = 0.0;
            double bk = 0.0;

            if (k % numprocs == rank)
            {
                pivot = A[k][k];
                bk = b[k];
            }

            MPI_Request reqs[2];
            MPI_Ibcast(&pivot, 1, MPI_DOUBLE, k % numprocs, MPI_COMM_WORLD, &reqs[0]);
            MPI_Ibcast(&bk, 1, MPI_DOUBLE, k % numprocs, MPI_COMM_WORLD, &reqs[1]);

            // 通信重叠计算（先处理非主元列）
            for (int i = k + 1; i < n; ++i)
            {
                double factor = A[i][k];
                for (int j = k + 1; j < n; ++j)
                {
                    if (j % numprocs == rank && j != k)
                    {
                        A[i][j] -= factor * A[k][j];
                    }
                }
                // b[i] 暂不处理，等主元同步后再处理
            }

            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

            // 主元列做除法
            for (int j = k + 1; j < n; ++j)
            {
                if (j % numprocs == rank)
                {
                    A[k][j] /= pivot;
                }
            }
            b[k] = bk / pivot;
            A[k][k] = 1.0;

            // 本地更新 A[i][k] 与 b[i]
            for (int i = k + 1; i < n; ++i)
            {
                double factor = A[i][k];
                if (k % numprocs == rank)
                {
                    A[i][k] = 0.0;
                }
                b[i] -= factor * b[k];
            }
        }

        // 回代（只在0号进程）
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
