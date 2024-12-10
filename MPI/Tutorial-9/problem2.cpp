#include <iostream>
#include <math.h>
#include <string.h>
#include <mpi.h>

void printMatrix(double* matrix, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%.2lf\t", matrix[i*size + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printArray(double* array, int size)
{
    for (int i = 0; i < size; i++)
    {
        printf("%.2lf ", array[i]);
    }
    printf("\n");
}


void fill_matrix(double* matrix, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrix[i*size + j] = i + j;
        }
    }
}

void fill_vector(double* vec, int size)
{
    for (int i = 0; i < size; i++)
    {
        vec[i] = i;
    }
}


int main(int argc, char** argv)
{
    int comm_size, myid;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int n = 256;
    
    double* mat = new double [n*n];
    
    if(myid == 0) fill_matrix(mat, n);
    double* vec = new double[n];
    fill_vector(vec, n);
    // if(myid == 0)
    // {
    //     printMatrix(mat, n);
    //     printArray(vec, n);
    // }

    int local_m = (int) (n)/(comm_size);

    MPI_Scatter(mat, n*local_m, MPI_DOUBLE, mat, n*local_m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // printf("Process %d: \n", myid);
    // printMatrix(mat, n);

    double* y = new double[n];
    double* local_y = new double[local_m];

    for (int i = 0; i < local_m; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
        {
           sum += mat[i*local_m + j] * vec[j];
        }
        local_y[i] = sum;
    }

    MPI_Allgather(local_y, local_m, MPI_DOUBLE, y, local_m, MPI_DOUBLE, MPI_COMM_WORLD);

    if(myid == 0) printArray(y, n);
    MPI_Finalize();
    return 0;
}




