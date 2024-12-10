#include <iostream>
#include <math.h>
#include <fstream>
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
	return;
}

double phi(double x, double y)
{
    return ((pow(x, 2) - 1) * (pow(y, 2) - 1));
}

double q(double x, double y)
{
    return 2*(2 - pow(x, 2) - pow(y, 2));
}

void Poisson_Gauss_Seidel(double** spatial_matrix, int n, double delta, double ll, double ul)
{
    double x = -1;
    double y = -1;
    for (int i = 1; i < n - 1; i++)
    {
        x = ll + i*delta;
        for (int j = 1; j < n - 1; j++)
        {
            y = ll + j*delta;
            spatial_matrix[i][j] = 0.25 * (spatial_matrix[i+1][j] + spatial_matrix[i-1][j] + spatial_matrix[i][j+1] + spatial_matrix[i][j-1] + pow(delta,2) * q(x, y));
        }
    }
}

double squared_error(double** matrix, double** true_matrix, int n)
{
    double diff = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            diff += pow((matrix[i][j] - true_matrix[i][j]), 2);
        }
    }

    return pow(diff, 0.5);
}

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        cout << "Check Arguments" << endl;
        return 0;
    }

    double ll = -1; double ul = 1; double delta = atof(argv[1]);
    int n = (ul - ll)/delta + 1;
    
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
    {
        for (int j = 0; j < n; j++)
        {
            true_matrix[i][j] = phi(ll+i*delta, ll+j*delta);
        }
    }

    double t1, t2;
    t1 = omp_get_wtime();
    int iter = 0;
    while (squared_error(spatial_matrix, true_matrix, n) > 0.01)
    {
        if(iter++ % 100 == 0)
            cout << "Iter: " << iter << " Error : " << squared_error(spatial_matrix, true_matrix, n) << endl;
        Poisson_Gauss_Seidel(spatial_matrix, n, delta, ll, ul);
    }
    t2 = omp_get_wtime();
    cout << "It took " << iter << " iterations and " << t2 - t1 << " secs" << endl;
    printMatrix(spatial_matrix, n, n);

    // // Packing Data to plot
    // ofstream MyFile("p3a_phi_VS_x.txt");
    // for (int i = 0; i < n; i++)
    //     MyFile << ll + (i * delta) << ", " << spatial_matrix[i][15] << endl;

    return 0;
}

