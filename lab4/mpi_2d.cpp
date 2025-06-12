// #include <mpi.h>
// #include <iostream>
// #include <vector>
// #include <cmath>
// #include <cstdlib>
// #include <ctime>
// #include <cassert>

// using namespace std;

// void initialize_matrix(vector<vector<double>> &A, vector<double> &b, int n)
// {
//     srand(time(0));
//     for (int i = 0; i < n; ++i)
//     {
//         fill(A[i].begin(), A[i].end(), 0.0);
//         A[i][i] = rand() % 100 + 1000; // 避免除零
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

//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     int n;
//     if (rank == 0)
//     {
//         cout << "输入矩阵维度（建议可被 2 整除）: ";
//         cin >> n;
//     }

//     MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

//     int q = sqrt(size);
//     if (q * q != size || n % q != 0)
//     {
//         if (rank == 0)
//             cerr << "进程数需为平方数，且 n 可被 sqrt(P) 整除。" << endl;
//         MPI_Finalize();
//         return -1;
//     }

//     const int block_size = n / q;
//     int row_rank = rank / q, col_rank = rank % q;

//     vector<vector<double>> A_block(block_size, vector<double>(block_size));
//     vector<double> b_block(block_size);

//     if (rank == 0)
//     {
//         vector<vector<double>> A(n, vector<double>(n));
//         vector<double> b(n);
//         initialize_matrix(A, b, n);

//         // 分发矩阵块
//         for (int i = 0; i < q; ++i)
//         {
//             for (int j = 0; j < q; ++j)
//             {
//                 int dest = i * q + j;
//                 for (int ii = 0; ii < block_size; ++ii)
//                 {
//                     if (dest == 0)
//                     {
//                         for (int jj = 0; jj < block_size; ++jj)
//                             A_block[ii][jj] = A[i * block_size + ii][j * block_size + jj];
//                         b_block[ii] = b[i * block_size + ii];
//                     }
//                     else
//                     {
//                         MPI_Send(&A[i * block_size + ii][j * block_size], block_size, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
//                     }
//                 }
//                 if (dest != 0)
//                 {
//                     MPI_Send(&b[i * block_size], block_size, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
//                 }
//             }
//         }
//     }
//     else
//     {
//         for (int i = 0; i < block_size; ++i)
//         {
//             MPI_Recv(A_block[i].data(), block_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         }
//         MPI_Recv(b_block.data(), block_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     }

//     vector<double> Arow(n);
//     vector<double> brow(1);
//     MPI_Barrier(MPI_COMM_WORLD);
//     double t0 = MPI_Wtime();

//     for (int k = 0; k < n; ++k)
//     {
//         int owner_row = k / block_size;
//         int owner_col = k / block_size;
//         int owner = owner_row * q + owner_col;

//         if (rank == owner)
//         {
//             int local_i = k % block_size;
//             fill(Arow.begin(), Arow.end(), 0.0);
//             for (int j = 0; j < block_size; ++j)
//                 Arow[owner_col * block_size + j] = A_block[local_i][j];
//             brow[0] = b_block[local_i];
//         }

//         MPI_Bcast(Arow.data(), n, MPI_DOUBLE, owner, MPI_COMM_WORLD);
//         MPI_Bcast(brow.data(), 1, MPI_DOUBLE, owner, MPI_COMM_WORLD);

//         for (int i = 0; i < block_size; ++i)
//         {
//             int global_i = row_rank * block_size + i;
//             if (global_i > k)
//             {
//                 double factor = A_block[i][k - col_rank * block_size] / Arow[k];
//                 for (int j = 0; j < block_size; ++j)
//                 {
//                     int global_j = col_rank * block_size + j;
//                     A_block[i][j] -= factor * Arow[global_j];
//                 }
//                 b_block[i] -= factor * brow[0];
//             }
//         }
//     }

//     MPI_Barrier(MPI_COMM_WORLD);
//     double t1 = MPI_Wtime();
//     if (rank == 0)
//         cout << "总时间: " << (t1 - t0) * 1e3 << " ms" << endl;

//     MPI_Finalize();
//     return 0;
// }

#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cassert>
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

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int sqrt_p = sqrt(size);
    assert(sqrt_p * sqrt_p == size); // must be perfect square

    int n;
    if (rank == 0)
    {
        cout << "输入矩阵维度（建议可被 " << sqrt_p << " 整除）: ";
        cin >> n;
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    assert(n % sqrt_p == 0);
    int block_size = n / sqrt_p;

    // 获取当前进程的2D位置
    int row_id = rank / sqrt_p;
    int col_id = rank % sqrt_p;

    // 为每个进程分配对应的子块
    vector<vector<double>> A_block(block_size, vector<double>(n));
    vector<double> b_block(block_size);

    vector<vector<double>> A_full;
    vector<double> b_full;
    if (rank == 0)
    {
        A_full.resize(n, vector<double>(n));
        b_full.resize(n);
        initialize_matrix(A_full, b_full, n);
    }

    // 广播矩阵和b向量
    for (int i = 0; i < n; ++i)
    {
        if (rank == 0)
        {
            MPI_Bcast(A_full[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Bcast(A_block[0].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD); // dummy recv
        }
    }

    if (rank == 0)
    {
        MPI_Bcast(b_full.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else
    {
        vector<double> temp_b(n);
        MPI_Bcast(temp_b.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (int i = 0; i < block_size; ++i)
            b_block[i] = temp_b[row_id * block_size + i];
    }

    // 提取本地块数据（行方向）
    for (int i = 0; i < block_size; ++i)
    {
        int global_row = row_id * block_size + i;
        for (int j = 0; j < n; ++j)
            A_block[i][j] = rank == 0 ? A_full[global_row][j] : A_block[i][j];
        if (rank == 0)
            b_block[i] = b_full[global_row];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // 高斯消元
    for (int k = 0; k < n; ++k)
    {
        int owner_row = k / block_size;
        int owner_rank = owner_row * sqrt_p + col_id;
        vector<double> A_k_row(n);
        double b_k;

        if (row_id == owner_row)
        {
            int local_k = k % block_size;
            A_k_row = A_block[local_k];
            b_k = b_block[local_k];
        }

        MPI_Bcast(A_k_row.data(), n, MPI_DOUBLE, owner_row, MPI_COMM_WORLD);
        MPI_Bcast(&b_k, 1, MPI_DOUBLE, owner_row, MPI_COMM_WORLD);

        for (int i = 0; i < block_size; ++i)
        {
            int global_i = row_id * block_size + i;
            if (global_i > k)
            {
                double factor = A_block[i][k];
                for (int j = k + 1; j < n; ++j)
                    A_block[i][j] -= factor * A_k_row[j];
                b_block[i] -= factor * b_k;
                A_block[i][k] = 0.0;
            }
        }
    }

    double end_time = MPI_Wtime();
    if (rank == 0)
    {
        cout << "总时间: " << (end_time - start_time) * 1000 << " ms" << endl;
    }

    MPI_Finalize();
    return 0;
}
