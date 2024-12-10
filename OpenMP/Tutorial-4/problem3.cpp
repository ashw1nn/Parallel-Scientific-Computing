#include <iostream>
#include <math.h> /* pow */
#include <time.h> /* clock_t */
using namespace std;
#ifdef _OPENMP
#include <omp.h>
#endif

void print_array(int n, double* array){
    for (int i = 0; i < n; i++)
    {
        cout << array[i] << endl;
    }
    cout << endl;
}


double u(double x)
{
    return 7 - x * tan(x);
}


double du(double x, double grid_len, double ll, double ul)
{   
    double dudx = 0;
    if (x == ll)
        dudx = (u(x + grid_len) - u(x))/(grid_len);
    else if (x == ul)
        dudx = (u(x) - u(x - grid_len))/(grid_len);
    else
        dudx = (u(x - 2*grid_len) - 8*u(x - grid_len) + 8*u(x + grid_len) - u(x + 2*grid_len))/(12*grid_len);

    return dudx;
}

int main()
{
    double ll = -1;
    double ul = 1;
    double grid_length = 0.001;

    // Creating array for derivatives
    int arr_size = ((ul - ll)/grid_length) + 1;
    double *dudx = new double[arr_size];

    clock_t t;
    t = clock();
    #pragma omp parallel for num_threads(2)
    for (int i = 0; i < arr_size; i++)
    {
        double x = ll + i*grid_length;
        dudx[i] = du(x, grid_length, ll, ul);
    }
    t = clock() - t;
    cout << "2 Threads time taken: " << t << endl;

    t = clock();
    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < arr_size; i++)
    {
        double x = ll + i*grid_length;
        dudx[i] = du(x, grid_length, ll, ul);
    }
    t = clock() - t;
    cout << "4 Threads time taken: " << t << endl;

    t = clock();
    #pragma omp parallel for num_threads(8)
    for (int i = 0; i < arr_size; i++)
    {
        double x = ll + i*grid_length;
        dudx[i] = du(x, grid_length, ll, ul);
    }
    t = clock() - t;
    cout << "8 Threads time taken: " << t << endl;

    print_array(arr_size, dudx);
    
    delete[] dudx;
    return 0;
}