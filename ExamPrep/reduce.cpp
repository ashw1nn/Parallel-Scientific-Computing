#include <stdio.h>
#include <mpi.h>
#include <string.h>

int main(int argc, char *argv[])
{
    int i, myid, size, tag=100;
    // char message_send[100], message_recv[100];
    MPI_Status status;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (size != 3)
    {
        return 0;
    }

    int a, b, c, d;
    a = 1; b = 0; c = 2; d = 0;
    
    if (myid == 0)
    {
        MPI_Reduce(&a, &b, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&c, &d, 1, MPI_INT, MPI_SUM, 1, MPI_COMM_WORLD);
    }
    else if (myid == 1)
    {
        MPI_Reduce(&c, &d, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&a, &b, 1, MPI_INT, MPI_SUM, 1, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Reduce(&a, &b, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&c, &d, 1, MPI_INT, MPI_SUM, 1, MPI_COMM_WORLD);
    }
    
    printf("b = %d, d = %d\n", b, d);

    MPI_Finalize();
    return 0;
}
