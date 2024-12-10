#include <iostream>
#include <math.h>
#include <string.h>
#include <mpi.h>
using namespace std;


void print_array(double* array, int n)
{
    for (int i = 0; i < n; i++)
    {
        cout << array[i] << "\n";
    }
    cout << endl;
}


double f(double x)
{
    return x*tan(x);
}

double df(double x, double grid_size)
{
    double e = 0.01;
    double deri;
    if (x >= -1 - e && x <= -1 + e)
    {
        deri = (f(x + grid_size) - f(x))/(grid_size);
    }
    else if (x <= 1 + e && x >= 1 - e)
    {
        deri = (f(x) - f(x - grid_size))/(grid_size);
    }
    else
    {
        deri = (f(x + grid_size) - f(x - grid_size))/(2 * grid_size);
    }

    return deri;
}


int main(int argc, char** argv)
{

    int myid, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    double ll = -1, ul = 1; 
    double grid_size = 0.000001;
    int n = (int) ((ul - ll)/grid_size + 1) ;
    
    double* derivatives = nullptr;
    if (myid == 0)
    {
        derivatives = new double[n];
    }

    int local_n = (int) (n/size);
    int remainder = n % size;
    double* local_derivatives = new double[local_n];
    
    double x;
    for (int i = 0; i < local_n; i++)
    {
        x = ll + local_n*grid_size*myid + i*grid_size;
        local_derivatives[i] = df(x, grid_size);
    }

    MPI_Gather(local_derivatives, local_n, MPI_DOUBLE,
                derivatives, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myid == 0)
    {
        int flag = remainder;
        while (flag)
        {
            x = ul - (flag - 1)*grid_size;
            derivatives[n - flag] = df(x, grid_size);
            flag--;
        }
    }


    if (myid == 0)
    {
        print_array(derivatives, n);
    }

    MPI_Finalize();

    return 0;
}