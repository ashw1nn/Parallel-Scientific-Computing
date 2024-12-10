#include <stdio.h>
#include <math.h>
#include<string.h>
#include <mpi.h>
#define PI 3.141592653589793


double f(double x)
{
    return (sin(x)/(2*pow(x, 3))) ;
}




int main(int argc, char** argv)
{
    int myid, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int n = 256;
    double x_start = 1, x_end = PI;
    double h = (x_end - x_start)/(n);
    
    int recv_size = (int) n/size;
    double x, local_sum = 0, ans = 0;

    for (int i = 0; i < recv_size; i++)
    {
        x = x_start + myid*recv_size*h + i*h;
        local_sum += f(x)*h;
    }
    printf("%d: Local sum %lf\n", myid, local_sum);
    MPI_Allreduce(&local_sum, &ans, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    ans += f(x_end);

    if (myid == 0)
    {
        printf("The value of Intergral is %lf\n", ans);
        printf("The error bw exact and numerical is %lf\n", ans - 0.198557298811);
    }
    MPI_Finalize();

    return 0;
}


