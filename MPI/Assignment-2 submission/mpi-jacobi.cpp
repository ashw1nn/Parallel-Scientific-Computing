#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>
#include <mpi.h>
using namespace std;
#define PI 3.141592653589793


void swapDoublePointerArrays(double*** a, double*** b)
{
    double** temp = *a;
    *a = *b;
    *b = temp;
}


void print_matrix(double** mat, int rows, int cols)
{
    for (int i = 0; i < rows - 2; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            cout << mat[i][j] << ", ";
        }
        cout << endl;
    }
    cout << endl;
}


void update_ghost_points(double** &mesh, int grid_points, int local_grid_points)
{   
    int myid, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Status status;

    int leftProc = myid - 1;
    int rightProc = myid + 1;
    if(myid == 0) leftProc = MPI_PROC_NULL;
    if(myid == nprocs - 1) rightProc = MPI_PROC_NULL;
    // Filling start ghost points
    MPI_Sendrecv(mesh[local_grid_points], grid_points, MPI_DOUBLE, rightProc, 100,
                 mesh[0], grid_points, MPI_DOUBLE, leftProc, 100, MPI_COMM_WORLD, &status);
    
    // Filling end ghost points
    MPI_Sendrecv(mesh[1], grid_points, MPI_DOUBLE, leftProc, 101,
                 mesh[local_grid_points + 1], grid_points, MPI_DOUBLE, rightProc, 101, MPI_COMM_WORLD, &status);
    

}


void init_local_mesh(double** &matrix, int x_size, int y_size, int ll, int ul, double delta)
{
    int myid, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    double x = ll; double y = ll;
    double e = 0.0001;
    for (int i = 1; i < y_size; i++)
    {
        y = ll + (i-1)*delta + myid*(y_size - 2)*delta;
        for (int j = 0; j < x_size; j++)
        {
            x = ll + (j)*delta;

            if (y > ul)
            {
                matrix[i][j] = 0;
            }
            else if (x >= ll - e && x <= ll + e)
            {
                matrix[i][j] = sin(2*PI*y);
            }
            else
            {
                matrix[i][j] = 0;
            }
        }
    }
}


double q(int i, int j, int ll, double grid_size)
{
    double summ = 0;
    double x = ll + j * grid_size;
    double y = ll + i * grid_size;
    summ = pow(x, 2) + pow(y, 2);
    return summ;
}


double l_norm(double** &mesh, int rows, int cols)
{
    double summ = 0;
    for (int i = 1; i < rows + 1; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            summ += pow(mesh[i][j], 2);
        }
    }
    return summ;
}


double l_diff(double** &mesh, double** &new_mesh, int rows, int cols)
{
    double sum_diff = 0;
    for (int i = 1; i < rows + 1; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            sum_diff += pow((mesh[i][j] - new_mesh[i][j]), 2);
        }
    }
    return sum_diff;
}


bool convergence(double** &mesh, double** &new_mesh, int rows, int cols)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    double val = l_diff(mesh, new_mesh, rows, cols);
    double norm = l_norm(mesh, rows, cols);
    double aggregate_val = 0, aggregate_norm = 0;
    MPI_Allreduce(&val, &aggregate_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&norm, &aggregate_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double aggregate = pow((aggregate_val/aggregate_norm), 1.0/2);
    
    if (myid == 0) cout << "Aggregate error: " <<  aggregate << endl;
    if ( aggregate < 0.0001)
    {
        if (myid == 0) cout << "Error less than threshold, MESH CONVERGED!!" << endl;
        return true;
    }
    return false;
}


void Jacobi(double** &mesh, double** &new_mesh, int grid_points, int mod_grid_points, int local_grid_points, int ll, int ul, double delta)
{
    int myid, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    // Copying mesh_value into new_mesh
    for (int i = 0; i < local_grid_points + 2; i++)
    {
        for (int j = 0; j < grid_points; j++)
        {
            new_mesh[i][j] = mesh[i][j];
        }
    }

    // Values for Mesh points
    for (int i = 1; i < local_grid_points+1; i++)
    {
        for (int j = 1; j < grid_points-1; j++)
        {
            new_mesh[i][j] = (1.0/4) * ((mesh[i+1][j] + mesh[i-1][j] + mesh[i][j+1] + mesh[i][j-1]) + (pow(delta, 2)*q(i, j, ll, delta)));
        }
        
    }

    // Right Boundary Condition
    for (int i = 1; i < local_grid_points + 2; i++)
    {
        new_mesh[i][grid_points-1] = (1.0/3)*(4*new_mesh[i][grid_points-2] - new_mesh[i][grid_points-3]);
    }

    // Updating the ghost points again
    update_ghost_points(new_mesh, grid_points, local_grid_points);

}


int main(int argc, char** argv)
{
    int myid, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (argc != 2)
    {
        cout << "Check Arguments";
        return 0;
    }

    double delta = atof(argv[1]);
    int ul = 1; int ll = -1;
    int grid_points = (int) ((ul - ll)/delta) + 1;

    // Extra Padding to compensate for division
    int mod_grid_points = grid_points;
    while (mod_grid_points % nprocs) mod_grid_points++;


    /*****************************
     Plan is to slice the mesh along rows, 
     so create ghost points along rows only
    *****************************/

    // Per Process grid points 
    int local_grid_points = mod_grid_points/nprocs;
    double** local_mesh = new double* [local_grid_points + 2];
    for (int i = 0; i < local_grid_points + 2; i++)
    {
        local_mesh[i] = new double[grid_points];
    }
    
    // Local Mesh for each process is filled
    init_local_mesh(local_mesh, grid_points, local_grid_points+2, ll, ul, delta);
    
    // Transferring data to ghost points
    update_ghost_points(local_mesh, grid_points, local_grid_points);

    // Printing the data in each process
    // cout << "Process: " << myid << endl;
    // print_matrix(local_mesh, local_grid_points + 2, grid_points);


    int flag = 0; int iter = 0;
    double** new_mesh = new double* [local_grid_points + 2];
    for (int i = 0; i < local_grid_points + 2; i++)
    {
        new_mesh[i] = new double[grid_points];
    }

    double start_time, end_time, total_time_per_thread, max_time, min_time;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    
    do
    {
        if(flag)
        {
            // cout << myid << ":" << *local_mesh << ": Local Mesh" << " Old Iteration: " << iter  << endl;
            // cout << myid << ":" << *new_mesh << ": New Mesh" << " Old Iteration: " << iter <<  endl;
            // Free up old mesh and point it to the new one
            // for (int i = 0; i < local_grid_points + 2; i++)
            // {
            //     delete[] local_mesh[i];
            // }
            // delete[] local_mesh;

            // Pointing local_mesh to new one
            swapDoublePointerArrays(&local_mesh, &new_mesh);
            
            // cout << myid << ":" << *local_mesh << ": Local Mesh" << " New Iteration: " << iter << endl;
            // cout << myid << ":" << *new_mesh << ": New Mesh" << " New Iteration: " << iter << endl;
            
            // cout << myid << ":" << "Succesfully Deallocated Memory ; " << "Iter: " << iter << endl;
        }

        // cout << myid << ":" << "Entering Jacobi; " << "iter: " << iter << endl;
        Jacobi(local_mesh, new_mesh, grid_points, mod_grid_points, local_grid_points, ll, ul, delta);

        // ************************************************************
        // PRINTING NEW AND OLD MESH
        // for (int p = 0; p < nprocs; p++)
        // {
        //     if (myid == 0)
        //     {
        //         cout << myid << ":" << "Done with Jacobi " << "; iter: " << iter << endl;
        //         cout << "---------------------------------------------------------------------------------------------" << endl;
        //         cout << myid << ":" << "Mesh: " << endl;
        //         print_matrix(local_mesh, local_grid_points + 2, grid_points);
        //         cout << endl;
        //         cout << "---------------------------------------------------------------------------------------------" << endl;
        //         cout << myid << ":" << "New mesh: " << endl;
        //         print_matrix(new_mesh, local_grid_points + 2, grid_points);
        //         cout << "---------------------------------------------------------------------------------------------" << endl;
        //     }
        //     MPI_Barrier(MPI_COMM_WORLD);
        // }
        // ************************************************************
        flag = 1;
        iter++;

    } while (!convergence(local_mesh, new_mesh, local_grid_points, grid_points));

    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    total_time_per_thread = end_time - start_time;
    MPI_Reduce(&total_time_per_thread, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_time_per_thread, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    if (myid == 0) 
    {
        cout << "It took " << iter << " iterations." << endl;
        cout << "Max Time: " << max_time << endl;
        cout << "Min Time " << min_time << endl << endl;
    }

    // cout << "Final Mesh: " << endl;
    for (int p = 0; p < nprocs; p++)
    {
        if (myid == p)
        {
            print_matrix(new_mesh, local_grid_points + 2, grid_points);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    // Printing new mesh to file
    char filename[100];
    sprintf(filename, "parallel_jacobi_%.3f_%dt.csv", delta, nprocs);
    ofstream MyFile;
    MyFile.open(filename, std::ios_base::app);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int p = 0; p < nprocs; p++)
    {
        if (myid == p)
        {
            for (int i = 0; i < local_grid_points; i++)
            {
                double y = ll + i*delta + myid*local_grid_points*delta;
                for (int j = 0; j < grid_points; j++)
                {
                    if (y != ll) MyFile << new_mesh[i][j] << ", ";
                }
                if (y != ll) MyFile << endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    for (int i = 0; i < local_grid_points + 2; i++)
    {
        delete[] local_mesh[i];
        delete[] new_mesh[i];
    }
    delete[] local_mesh;
    delete[] new_mesh;

    MPI_Finalize();
    return 0;
}
