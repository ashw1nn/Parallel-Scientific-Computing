#include <stdio.h>
#include <mpi.h>
#include <string.h>

int main(int argc, char *argv[])
{
    int i, myid, size, tag=100;
    char message_send[100], message_recv[100];
    MPI_Status status;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (myid != 0)
    {
        sprintf(message_send, " Hello from process %d\n", myid);
        MPI_Send(message_send, 50, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
    }
    else
    {
        for (int i = 1; i < size; i++)
        {
            MPI_Recv(message_recv, 50, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status);
            printf("\n %s", message_recv);
        }
        
        sprintf(message_send, " Hello from process %d\n", myid);
        printf("\n %s", message_send);
    }
    
    MPI_Finalize();

    return 0;
}
