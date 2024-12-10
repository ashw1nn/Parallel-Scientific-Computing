#include<stdio.h>
#include<string.h>
#include<mpi.h> 		/* new and contains the function prototypes, macros, variables.. */

int main(int argc, char** argv)
{
    int myid, size;
    int send, sum = 0;
    int receive_buf[50];
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (myid == 0)
    {
        for (int i = 0; i < size - 1; i++)
        {
            MPI_Recv(&receive_buf[i], 1, MPI_INT, MPI_ANY_SOURCE, 101, MPI_COMM_WORLD, &status);
            sum += receive_buf[i];
        }
        printf("The total sum is %d\n", sum);
        
    }
    else
    {
        send = myid + 1;
        MPI_Send(&send, 1, MPI_INT, 0, 101, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();

    return 0;
}