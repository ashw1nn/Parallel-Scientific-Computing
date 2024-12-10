#include<stdio.h>
#include<string.h>
#include<mpi.h> 		/* new and contains the function prototypes, macros, variables.. */

int main(int argc, char** argv)
{
  int myid, size, tag=100;
  char message_send[50], message_recv[50];
  MPI_Status status; 		/* data type that is defined in mpi.h... */

  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* tells about the number of processes */
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* will return the rank or id of the process that called it.  */

  /* printf("%d", size); */

  if (myid == 0) 		/* not the designated process to receive..send the data.. */
  {
    for (int i = 1; i < size; i++)
    {
      sprintf(message_send, "Hello");
      MPI_Send(message_send, 50, MPI_CHAR, i, tag, MPI_COMM_WORLD);
      printf("Sending msg from %d to %d\n", myid, i);
    }
  }
  else
  {
    MPI_Recv(message_recv, 50, MPI_CHAR, 0, tag, MPI_COMM_WORLD, &status);
    printf("%s from %d \n", message_recv, myid);
  }

  MPI_Finalize();

  return 0;
}
