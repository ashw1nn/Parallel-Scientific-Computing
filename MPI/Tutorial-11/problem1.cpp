#include <iostream>
#include <math.h>
#include <string.h>
#include <mpi.h>
using namespace std;
#define MAT_SIZE 700

void fill_matrix(double A[][MAT_SIZE], int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            A[i][j] = (i + j);
        }   
    }
}


void fill_array(double B[], int size)
{
    for (int i = 0; i < size; i++)
    {
        B[i] = i;
    }
}

void print_array(double* array, int n)
{
    for (int i = 0; i < n; i++)
    {
        cout << array[i] << "\n";
    }
    cout << endl;
}


void print_matrix(double mat[][MAT_SIZE], int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << mat[i][j] << " ";   
        }
        cout << endl;
    }
    cout << endl;
}


void Matrix_Vector_Mult(int size, double* matrix, double* vector, double* result)
{
    int nprocs, myid;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    double* local_result = new double[size/nprocs];
    
    double local_matrix [size][size];   //local matrix

    MPI_Scatter(matrix, (size*size)/nprocs, MPI_DOUBLE, local_matrix, (size*size)/nprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD); //Scatter the Matrix
    MPI_Bcast(vector, size, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Broadcast the Vector

    for (int i = 0; i < size/nprocs; i++)
    {
        for (int j = 0; j < size; j++)
        {
            local_result[i] = vector[j]*local_matrix[i][j];
        }
    }

    MPI_Gather(local_result, size/nprocs, MPI_DOUBLE, result, size/nprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}


int main(int argc, char** argv)
{
    int nprocs, myid;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (MAT_SIZE % nprocs){   //Valid Communicator Size?
        MPI_Finalize();
        return 0;
    }

    double A[MAT_SIZE][MAT_SIZE]; double B[MAT_SIZE];
    if (myid == 0)
    {
        fill_matrix(A, MAT_SIZE);
        // print_matrix(A, MAT_SIZE);
        fill_array(B, MAT_SIZE);
    }

    double* result = new double[MAT_SIZE];
    Matrix_Vector_Mult(MAT_SIZE, (double*) A, B, result);
    if (myid == 0)
        print_array(result, MAT_SIZE);

    delete[] result;
    MPI_Finalize();

    return 0;
}
