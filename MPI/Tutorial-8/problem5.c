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

    int n;
    double x_start, x_end;
 
    if (myid == 0)
    {
        printf("Enter the value of n (in multiples of 8): \n");
        scanf("%d", &n);
        printf("Enter the value of a: \n");
        scanf("%lf", &x_start);
        printf("Enter the value of b: \n");
        scanf("%lf", &x_end);
        printf("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double h = (x_end - x_start)/(n);

    int recv_size = (int) n/size;
    double x, even_sum = 0, odd_sum = 0, ans = 0;
    double total_odd_sum = 0, total_even_sum = 0;

    for (int i = 1; i < recv_size; i++)
    {
        x = x_start + myid*recv_size*h + i*h;
        if (i % 2 == 1)
            odd_sum += f(x);
        else
            even_sum += f(x);
    }

    MPI_Allreduce(&even_sum, &total_even_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&odd_sum, &total_odd_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    ans = (h * (f(x_start) + f(x_end) + 4*total_odd_sum + 2*total_even_sum)) / 3;

    if (myid == 0)
    {
        printf("The value of Intergral is %lf\n", ans);
        printf("The error bw exact and numerical is %lf\n", ans - 0.198557298811);
    }
    MPI_Finalize();

    return 0;
}


