#include <iostream>
#include <math.h>
#include <string.h>
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


void Euler_Upwind(double** result, double delta_x, double delta_t, int time_steps, int grid_points)
{
    for (int k = 1; k < time_steps; k++)
    {
        for (int i = 1; i < grid_points; i++)
        {
            result[k][i] = ((result[k-1][i-1] - result[k-1][i])/delta_x)*delta_t + result[k-1][i];
        }
        result[k][0] = 0; result[k][grid_points-1] = 0;
    }
    
}


void Euler_QUICK(double* u_prev, double* u_new, double delta_x, double delta_t, int time_steps, int grid_points)
{
    for (int k = 1; k <= time_steps; k++)
    {
        for (int i = 1; i < grid_points; i++)
        {
            if (i == 1 || i == grid_points - 1)
            {
                u_new[i] = ((-delta_t/delta_x) * (u_prev[i] - u_prev[i-1])) + u_prev[i];
            }
            else
            {
                u_new[i] = ((-delta_t/delta_x) * ((3.0/8.0)*u_prev[i] - (7.0/8.0)*u_prev[i-1] + (1.0/8.0)*u_prev[i-2] + (3.0/8.0)*u_prev[i+1])) + u_prev[i];
            }
        }
        u_new[0] = 0; u_new[grid_points-1] = 0;
        
        // Copying values for next iteration
        u_prev = u_new;
    }
}


void Analytical_Solution(double t, double** result, double delta_x, double delta_t, int time_steps, int grid_points)
{
    for (int k = 0; k < time_steps; k++)
    {
        for (int i = 0; i < grid_points; i++)
        {
            double x = 0 + i*delta_x;
            if (((x - t) <= 0.5) && ((x - t) >= 0))
                result[k][i] = sin(4*PI*x);
            else
                result[k][i] = 0;    
        }
        result[k][0] = 0; result[k][grid_points-1] = 0;
    }
}


int main(int argc, char** argv)
{
    double delta_x = 0.002; double delta_t = 0.0001; double L = 2.0;
    int grid_points = (int) (L/delta_x) + 1; int time_steps = (int) (1/delta_t) + 1;

    // u(x, t)
    double** uxt = new double*[time_steps];
    for (int i = 0; i < time_steps; i++)
    {
        uxt[i] = new double[grid_points];
    }
    
    // Base Case
    for (int i = 0; i < grid_points; i++)
    {
        double x = 0 + i*delta_x;
        if (x <= 0.5)
            uxt[0][i] = sin(4*PI*x);
        else
            uxt[0][i] = 0;    
    }
    uxt[0][0] = 0; uxt[0][grid_points-1] = 0;
    
    double sol_time = 1.0;
    /************************************** Upwind **************************************/
    // Euler_Upwind(uxt, delta_x, delta_t, time_steps, grid_points);
    // print_array(uxt[ (int)(sol_time/delta_t)], grid_points);
    /************************************** Upwind **************************************/

    /************************************** Quick **************************************/
    double* u_quick = new double[grid_points];
    // Base Case
    for (int i = 0; i < grid_points; i++)
    {
        double x = 0 + i*delta_x;
        if (x <= 0.5)
            u_quick[i] = sin(4*PI*x);
        else
            u_quick[i] = 0;    
    }
    u_quick[0] = 0; u_quick[grid_points-1] = 0;
    double* u_quick_new = new double[grid_points];
    Euler_QUICK(u_quick, u_quick_new, delta_x, delta_t, (int)((sol_time - 0.037)/delta_t), grid_points);
    print_array(u_quick_new, grid_points);
    /************************************** Quick **************************************/

    /************************************** Analytical **************************************/
    // Analytical_Solution(sol_time, uxt, delta_x, delta_t, time_steps, grid_points);
    // print_array(uxt[ (int)(sol_time/delta_t)], grid_points);
    /************************************** Analytical **************************************/

    


    return 0;
}



