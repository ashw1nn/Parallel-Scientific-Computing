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
			printf("%0.3lf ", mat[i][j]);
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

void Poisson_Gauss_Seidel_Diagonal(double** spatial_matrix, int n, double delta, double ll, double** true_matrix)
{
    double x = -1; double y = -1;
    
    int iter = 0;
    while (squared_error(spatial_matrix, true_matrix, n) > 0.01)
    { 
        if(iter++ % 100 == 0)
        {
            cout << "Iter: " << iter << " Error : " << squared_error(spatial_matrix, true_matrix, n) << endl;
        }

        int nd = 2*n -1; int istart = 0; int iend = 0; int j = 0;
    
        for (int l = 1; l < nd; l++)
        {
            if (l <= n)
            {
                istart = 0;
                iend = l - 1;
            }
            else
            {
                istart = l - n;
                iend = n - 1;
            }

            // #pragma omp parallel for num_threads(2)
            for (int i = istart; i <= iend; i++)
            {
                j = l - i - 1;
                if ((i == 0) || ((j) == 0) || (i == n-1) || (j == n-1))
                {
                    continue;
                }
                x = ll + i*delta;
                y = ll + j*delta;
                spatial_matrix[i][j] = 0.25 * (spatial_matrix[i + 1][j] + spatial_matrix[i - 1][j] + spatial_matrix[i][j + 1] + spatial_matrix[i][j - 1] + pow(delta, 2) * q(x, y));
            }   
            
        }
    }
    cout << "It took " << iter << " iterations" << endl;
}

void Poisson_Gauss_Seidel_Diagonal_parallel(double** spatial_matrix, int n, double delta, double ll, double** true_matrix)
{
    double x = -1; double y = -1;
    int nd = 2*n -1; int istart = 0; int iend = 0; int j = 0;
    
    for (int l = 1; l < nd; l++)
    {
        if (l <= n)
        {
            istart = 0;
            iend = l - 1;
        }
        else
        {
            istart = l - n;
            iend = n - 1;
        }

        #pragma omp for
        for (int i = istart; i <= iend; i++)
        {
            j = l - i - 1;
            if ((i == 0) || ((j) == 0) || (i == n-1) || (j == n-1))
            {
                continue;
            }
            x = ll + i*delta;
            y = ll + j*delta;
            spatial_matrix[i][j] = 0.25 * (spatial_matrix[i + 1][j] + spatial_matrix[i - 1][j] + spatial_matrix[i][j + 1] + spatial_matrix[i][j - 1] + pow(delta, 2) * q(x, y));
        }   
    }
}



int main(int argc, char** argv)
{
    if (argc != 3)
    {
        cout << "Check Arguments" << endl;
        return 0;
    }

    double ll = -1; double ul = 1; double delta = atof(argv[1]); int threads = atoi(argv[2]);
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
        for (int j = 0; j < n; j++)
            true_matrix[i][j] = phi(ll+i*delta, ll+j*delta);


    double t_start; double t_end;

    char filename1[100];
    sprintf(filename1, "out/p3_diagonal_%f_%d.txt", delta, threads);
    ofstream MyFile1(filename1);

    
    int iter = 0;
    t_start = omp_get_wtime();
    #pragma omp parallel num_threads(threads) shared(spatial_matrix, delta, ll , n)
    {
    while (squared_error(spatial_matrix, true_matrix, n) > 0.01)
    { 
        #pragma omp master
        {
        if(iter++ % 100 == 0)
        {
            cout << "Iter: " << iter << " Error : " << squared_error(spatial_matrix, true_matrix, n) << endl;
            MyFile1 << "Iter: " << iter << " Error : " << squared_error(spatial_matrix, true_matrix, n) << endl;
        }
        }
        Poisson_Gauss_Seidel_Diagonal_parallel(spatial_matrix, n, delta, ll, true_matrix);
    }
    }
    t_end = omp_get_wtime();


    cout << "Took " << t_end - t_start << " secs" << endl;
    MyFile1 << "Took " << t_end - t_start << " secs" << endl;
    printMatrix(spatial_matrix, n, n);

    // Packing Data to plot
    char filename[100];
    sprintf(filename, "plot_p3_diagonal_%f_%d.txt", delta, threads);
    ofstream MyFile(filename);
    for (int i = 0; i < n; i ++)
    {
		for (int j = 0; j < n; j++)
			MyFile << (ll + i*delta) << ", " << (ll + j*delta )<< ", "  << spatial_matrix[i][j] << endl;
	}

    return 0;
}

