#include <stdio.h>
#include <math.h>
#include<string.h>
#include <mpi.h>

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

    int size, myid;
    double receive_arr[5];
    double arr[5];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if(myid == 0)
    {
        int a = 10, b = 20;
        fill_array(arr, 5);
        print_array(arr, 5);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatter(arr, 1, MPI_DOUBLE, receive_arr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    print_array(receive_arr, 5);
    
    // MPI_Gather()
    // MPI_Allgather()

    MPI_Finalize();



    return 0;
}