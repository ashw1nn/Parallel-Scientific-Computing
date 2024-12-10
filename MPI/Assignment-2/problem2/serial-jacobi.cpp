#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>
#include <time.h>
using namespace std;
#define PI 3.141592653589793


void print_matrix(double** mat, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}


void init_mesh(double** matrix, int size, int ll, int ul, double grid_size)
{
    double x = ll; double y = ll;
    double e = 0.0001;
    for (int i = 0; i < size; i++)
    {
        y = ll + i * grid_size;
        for (int j = 0; j < size; j++)
        {
            x = ll + j * grid_size;
            if (x >= ll - e && x <= ll + e)
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


double norm2(double** mesh, int mesh_size)
{
    double summ = 0;
    for (int i = 0; i < mesh_size; i++)
    {
        for (int j = 0; j < mesh_size; j++)
        {
            summ += pow(mesh[i][j], 2);
        }
    }
    return summ;
}


double sq_diff(double** mesh, double** new_mesh, int mesh_size)
{
    double sum_diff = 0;
    for (int i = 0; i < mesh_size; i++)
    {
        for (int j = 0; j < mesh_size; j++)
        {
            sum_diff += pow((mesh[i][j] - new_mesh[i][j]), 2);
        }
    }
    return sum_diff;
}


bool convergence(double** mesh, double** new_mesh, int mesh_size)
{
    double val = sqrt(sq_diff(mesh, new_mesh, mesh_size)/norm2(mesh, mesh_size));
    cout << "Error: " <<  val << endl;
    if ( val < 0.0001)
    {
        cout << "Error less than threshold, MESH CONVERGED!!" << endl;
        return true;
    }
    return false;
}


double** Jacobi(double** mesh, int grid_points, double grid_size, int ll, int ul)
{
    double** new_mesh = new double* [grid_points];
    for (int i = 0; i < grid_points; i++)
    {
        new_mesh[i] = new double[grid_points];
    }
    // 3 Boundary conditions
    init_mesh(new_mesh, grid_points, ll, ul, grid_size);

    // Values for Mesh points
    for (int i = 1; i < grid_points-1; i++)
    {
        for (int j = 1; j < grid_points-1; j++)
        {
            new_mesh[i][j] = (1.0/4) * ((mesh[i+1][j] + mesh[i-1][j] + mesh[i][j+1] + mesh[i][j-1]) + (pow(grid_size, 2)*q(i, j, ll, grid_size)));
        }
        
    }
    // Right Boundary Condition
    for (int i = 0; i < grid_points; i++)
    {
        new_mesh[i][grid_points-1] = (1.0/3)*(4*new_mesh[i][grid_points-2] - new_mesh[i][grid_points-3]);
    }

    return new_mesh;
}


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        cout << "Check Arguments!!" << endl;
        return 0;
    }
    
    double grid_size = atof(argv[1]);
    int ul = 1; int ll = -1;
    int grid_points = (int) ((ul - ll)/grid_size) + 1;

    double** mesh = new double* [grid_points];
    for (int i = 0; i < grid_points; i++)
    {
        mesh[i] = new double[grid_points];
    }
    
    init_mesh(mesh, grid_points, ll, ul, grid_size);

    cout << "Intial mesh: " << endl;
    print_matrix(mesh, grid_points);
    cout << endl;
    
    int flag = 0; int iter = 0;
    double** new_mesh;

    
    clock_t start_time, end_time;
    start_time = clock();
    do
    {
        if(flag)
        {
            // Free up old mesh and point it to the new one
            for (int i = 0; i < grid_points; i++)
            {
                delete[] mesh[i];
            }
            mesh = new_mesh;
        }

        new_mesh = Jacobi(mesh, grid_points, grid_size, ll, ul);
        flag = 1;
        iter++;

        /************************************************************
        // PRINTING NEW AND OLD MESH
        // cout << "Mesh: " << endl;
        // print_matrix(mesh, grid_points);
        // cout << endl;
        // cout << endl << endl;
        // cout << "New mesh: " << endl;
        // print_matrix(new_mesh, grid_points);
        ************************************************************/

    } while (!convergence(mesh, new_mesh, grid_points));
    end_time = clock();

    double time_taken = double(end_time - start_time)/double(CLOCKS_PER_SEC);

    cout << "It took " << iter << " iterations." << endl;
    cout << "It took " << time_taken << " secs." << endl << endl;
    cout << "Final Mesh: " << endl;
    print_matrix(new_mesh, grid_points);

    // Print Data Packing
    
    ofstream MyFile("serial_jacobi_phi_vs_y_005.csv");
    for (int i = 0; i < grid_points; i++)
        MyFile << ll+i*grid_size << ", " << new_mesh[i][grid_points/2] << endl;

    ofstream MyFile2("serial_jacobi_phi_vs_x_005.csv");
    for (int j = 0; j < grid_points; j++)
        MyFile2 << ll+j*grid_size << ", " << new_mesh[grid_points/2][j] << endl;

    ofstream MyFile3("serial_jacobi_surface_005.csv");
        for (int i = 0; i < grid_points; i++)
        {
            for (int j = 0; j < grid_points; j++)
            {
                MyFile3 << new_mesh[i][j] << ", ";
            }
            MyFile3 << endl;
        }
        MyFile3 << endl;

    // Free up mesh and new_mesh
    for (int i = 0; i < grid_points; i++)
    {
        delete[] mesh[i];
        delete[] new_mesh[i];
    }
    mesh = nullptr; new_mesh = nullptr;


    return 0;
}