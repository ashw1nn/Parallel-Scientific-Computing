#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif 
using namespace std;


void printMatrix(double **mat, int m, int n)
{
	for (int i = 0; i < m; i ++)
    {
		for (int j = 0; j < n; j++)
			cout << mat[i][j] << " ";
		cout << endl;
	}
}

double phi(double x, double y)
{
    return ((pow(x, 2) - 1) * (pow(y, 2) - 1));
}

double q(double x, double y)
{
    return 2*(2 - pow(x, 2) - pow(y, 2));
}

double squared_error(double** matrix, double** true_matrix, int n)
{
    double diff = 0.0, den = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            diff += pow((matrix[i][j] - true_matrix[i][j]), 2);
            den += pow((true_matrix[i][j]), 2);
        }
    }

    return pow(diff/den, 0.5);
}


void Poisson_Gauss_Seidel(double** spatial_matrix, int n, double delta, double ll)
{
    double x = -1.0; double y = -1.0;
    
    #pragma omp parallel for collapse(2)
    // Red points update
    for (int i = 1; i < n - 1; i ++)
    {   
        for (int j = 1; j < n - 1; j++)
        {
            if (((i + j) % 2 == 1))
            {
                x = ll + i * delta;
                y = ll + j * delta;
                spatial_matrix[i][j] = 0.25 * (spatial_matrix[i + 1][j] + spatial_matrix[i - 1][j] + spatial_matrix[i][j + 1] + spatial_matrix[i][j - 1] + pow(delta, 2) * q(x, y));    
            }
        }
    }

    #pragma omp parallel for collapse(2)
    // Black points update
    for (int i = 1; i < n - 1; i ++)
    {
        for (int j = 1; j < n - 1; j ++)
        {
            if (((i + j) % 2 == 0))
            {
                x = ll + i * delta;
                y = ll + j * delta;
                spatial_matrix[i][j] = 0.25 * (spatial_matrix[i + 1][j] + spatial_matrix[i - 1][j] + spatial_matrix[i][j + 1] + spatial_matrix[i][j - 1] + pow(delta, 2) * q(x, y));    
            }
        }
    }
}


int main(int argc, char** argv)
{
    if (argc != 3)
    {
        cout << "Enter Grid Size, No of threads !!" << endl;
        return 0;
    }

    double ll = -1.0; double ul = 1.0; double delta = atof(argv[1]); int threads = atoi((argv[2]));
    int n = (int)((ul - ll)/delta) + 1;
    
    // Creating grid points
    double** spatial_matrix = new double*[n];
    for (int i = 0; i < n; i++)
        spatial_matrix[i] = new double[n];
    // Setting Initial Grid Value to 0 
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            spatial_matrix[i][j] = 0; 

    double** true_matrix = new double*[n];
    for (int i = 0; i < n; i++)
        true_matrix[i] = new double[n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            true_matrix[i][j] = phi(ll+i*delta, ll+j*delta);

    char filename1[100];
    sprintf(filename1, "out/p3_red_black_%f_%d.txt", delta, threads);
    ofstream MyFile1(filename1);
    int iter = 0;
    double t1, t2;
    t1 = omp_get_wtime();
    #pragma omp parallel num_threads(threads) shared(spatial_matrix, n, delta, ll)
    {
    while (squared_error(spatial_matrix, true_matrix, n) > 0.01)
    {
        Poisson_Gauss_Seidel(spatial_matrix, n, delta, ll);
        #pragma omp master
        {
        if(iter++ % 100 == 0)
        {
            cout << "Iter: " << iter << " Error : " << squared_error(spatial_matrix, true_matrix, n) << endl;
            MyFile1 << "Iter: " << iter << " Error : " << squared_error(spatial_matrix, true_matrix, n) << endl;
        }
        }
            
    }
    }
    t2 = omp_get_wtime();

    cout << "It took " << iter << " iterations and " << t2 - t1 << " secs" << endl;
    MyFile1 << "It took " << iter << " iterations and " << t2 - t1 << " secs" << endl;
    printMatrix(spatial_matrix, n, n);

    // Packing Data to plot
    char filename[100];
    sprintf(filename, "plot_p3_red_black_%f_%d.txt", delta, threads);
    ofstream MyFile(filename);
    for (int i = 0; i < n; i ++)
    {
		for (int j = 0; j < n; j++)
			MyFile << (ll + i*delta) << ", " << (ll + j*delta )<< ", "  << spatial_matrix[i][j] << endl;
	}

    return 0;
}

