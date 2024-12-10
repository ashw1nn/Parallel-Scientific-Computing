#include<stdio.h>
#include<string.h>
#include<mpi.h> 		/* new and contains the function prototypes, macros, variables.. */
#include <math.h>

void print_array(double* array, int size)
{
    for (int i = 0; i < size; i++)
    {
        printf("%.2lf ", array[i]);
    }
    printf("\n");
}


void fill_array(double* array, int size)
{
    for (int i = 0; i < size; i++)
    {
        array[i] = pow(i,1);
    }
}


int main(int argc, char** argv)
{
    int myid, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    double arr[50];
    int recv_size = (int)50/size;
    
    if (myid == 0)
        fill_array(arr, 50);

    if (myid == 0)
        print_array(arr, 50);

    MPI_Scatter(arr, recv_size, MPI_DOUBLE, 
                arr, recv_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Doubles the of values for the elements not on rank 0 
    if (myid != 0)
    {
        recv_size = (int)50/size;
        for (int i = 0; i < recv_size; i++)
        {
            arr[i] = arr[i] * 2;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(arr, recv_size, MPI_DOUBLE,
                arr, recv_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myid == 0)
        print_array(arr, 50);
    
    MPI_Finalize();

    return 0;
}