#include <iostream>
#include <math.h>
#include <string.h>
#include <mpi.h>
using namespace std;
#define PI 3.141592653589793


void print_array(double* array, int x_size)
{
    for (int j = 0; j < x_size; j++)
    {
        cout << array[j] << "\n";
    }
    cout << endl;
}


void update_ghost_points_forward(double* &local_u, int local_grid_points)
{
    int myid, nprocs;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int leftProc = myid - 1; int rightProc = myid + 1;
    if (myid == 0)
    {
        leftProc = MPI_PROC_NULL;
        local_u[0] = 0;
        local_u[1] = 0;
    }
    if (myid == nprocs - 1)
    {
        rightProc = MPI_PROC_NULL;
    }
    
    MPI_Sendrecv(&local_u[local_grid_points], 2, MPI_DOUBLE, rightProc, 100, 
                 local_u, 2, MPI_DOUBLE, leftProc, 100, MPI_COMM_WORLD, &status);

}


void update_ghost_points_backward(double* &local_u, int local_grid_points)
{
    int myid, nprocs;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int leftProc = myid - 1; int rightProc = myid + 1;
    if (myid == 0)
    {
        leftProc = MPI_PROC_NULL;
    }
    if (myid == nprocs - 1) 
    {
        rightProc = MPI_PROC_NULL;
        local_u[local_grid_points - 1] = 0;
    }
    
    MPI_Sendrecv(&local_u[2], 1, MPI_DOUBLE, leftProc, 100, 
                 &local_u[local_grid_points + 2], 1, MPI_DOUBLE, rightProc, 100, MPI_COMM_WORLD, &status);

}


void Euler_QUICK(double* &local_uxn_prev, double* &local_uxn_new, double delta_x, double delta_t, int time_steps, int local_grid_points, int mod_grid_points, int grid_points)
{
    int myid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    for (int k = 0; k < time_steps; k++)
    {
        for (int i = 2; i < local_grid_points + 2; i++)
        {
            if ((myid == 0 && i == 2))
            {
                local_uxn_new[i] = ((-delta_t/delta_x) * (local_uxn_prev[i] - local_uxn_prev[i-1])) + local_uxn_prev[i];
            }
            else
            {
                local_uxn_new[i] = ((-delta_t/delta_x) * ((3.0/8.0)*local_uxn_prev[i] - (7.0/8.0)*local_uxn_prev[i-1] + (1.0/8.0)*local_uxn_prev[i-2] + (3.0/8.0)*local_uxn_prev[i+1])) + local_uxn_prev[i];
            }
        }
        update_ghost_points_forward(local_uxn_new, local_grid_points);
        update_ghost_points_backward(local_uxn_new, local_grid_points);
        local_uxn_prev = local_uxn_new;
    }
    
}


int main(int argc, char** argv)
{
    int myid, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    double delta_x = 0.002; double delta_t = 0.0001; double L = 2.0;
    int grid_points = (int) (L/delta_x) + 1; 
    int mod_grid_points = grid_points;
    while (mod_grid_points % nprocs) mod_grid_points++;

    // u(x, t)
    double* uxn;
    int local_grid_points = mod_grid_points/nprocs;
    double* local_uxn = new double[local_grid_points + 3];
    double* local_uxn_new = new double[local_grid_points + 3];
    
    // Constructing u(x, 0) and adding ghost points for local_u(x, t)
    if (myid == 0)
    {
        uxn = new double[mod_grid_points];
        double x = 0;
        // Base Case
        for (int i = 0; i < mod_grid_points; i++)
        {
            x = 0 + i*delta_x;
            if (x <= 0.5)
                uxn[i] = sin(4*PI*x);
            else
                uxn[i] = 0;    
        }
        uxn[0] = 0; uxn[grid_points-1] = 0; // Boundary Conditions
    }

    MPI_Scatter(uxn, local_grid_points, MPI_DOUBLE,
                    &local_uxn[2], local_grid_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    update_ghost_points_forward(local_uxn, local_grid_points);
    update_ghost_points_backward(local_uxn, local_grid_points);

    double sol_time = atof(argv[1]);

    Euler_QUICK(local_uxn, local_uxn_new, delta_x, delta_t, (int)(sol_time/delta_t), local_grid_points, mod_grid_points, grid_points);
    MPI_Gather(&local_uxn_new[2], local_grid_points, MPI_DOUBLE,
                uxn, local_grid_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myid == 0)
    {
        print_array(uxn, grid_points);
    }

    MPI_Finalize();

    return 0;
}
