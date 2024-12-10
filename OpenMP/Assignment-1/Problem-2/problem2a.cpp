#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

void printMatrix(double **mat, int m, int n)
{
	for (int i = 0; i < m; i ++)
    {
		for (int j = 0; j < n; j++)
			printf("%0.2lf ", mat[i][j]);
		cout << endl;
	}
	return;
}

void printArray(double* arr, int n)
{
    for (int i = 0; i < n; i++)
        cout << arr[i] << endl;
    cout << endl;
}

double f(double x)
{
    if (x >= 0 && x <= 3.01)
        return sin(5*x);
    return MAXFLOAT;
}

double** create_pade_matrix(int n)
{
    double** matrix = new double*[n];
    for (int i = 0; i < n; i++)
        matrix[i] = new double[n];
    
    for (int i = 0; i < n; i++)
    {   
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                matrix[i][j] = 4;
            else if (j == i - 1 || j == i + 1)
                matrix[i][j] = 1;
            else
                matrix[i][j] = 0.00;
        }
    }
    matrix[0][0] = 1;
    matrix[n - 1][n - 1] = 1;
    matrix[0][1] = 2;
    matrix[n - 1][n - 2] = 2;

    return matrix;
}

double* create_zero_array(int n)
{
    double* arr = new double[n];
    for (int i = 0; i < n; i++)
        arr[i] = 0;
    return arr;
}

double* create_pade_array(int n, double ll, double ul)
{
    double h = (ul - ll)/(n - 1);
    double* arr = new double[n];
    double x = 0;
    for (int i = 0; i < n; i++){
        arr[i] = (3.0/h)*(f(x + h) - f(x - h));
        x += h;
    }
    arr[0] = ((-5.0/2) * f(0) + 2 * f(h) + 0.5 * f(2*h))/h; 
    arr[n - 1] = (1.0/h) * ((2.5*f((n - 1)*h)) - (2*f((n - 1)*h - h)) - (0.5*f((n - 1)*h - 2*h)));
    return arr;
}

void LU_decompose(double** matrix, int size)
{
    for (int i = 1; i < size; i++)
    {
        matrix[i][i - 1] = matrix[i][i - 1]/(matrix[i-1][i-1]); // l
        matrix[i][i] = (matrix[i][i] - matrix[i][i-1] * matrix[i-1][i]); // u
    }
}

void forward_substitute(double* z, double** L, double* y, int n)
{
    z[0] = y[0];
    for (int i = 1; i < n; i++)
        z[i] = y[i] - L[i][i - 1] * z[i - 1];
}

void back_substitute(double* x, double** U, double* z, int n)
{
    x[n - 1] = z[n - 1]/U[n-1][n-1];
    for (int i = n - 2; i >= 0; i--)
        x[i] = (z[i] - U[i][i + 1] * x[i + 1]) / (U[i][i]);
}

double* Derivative_Pade_Scheme(int size, double ll, double ul)
{
    double **A = create_pade_matrix(size);
    LU_decompose(A, size);
    cout << endl;

    double* x = create_zero_array(size);
    double* y = create_pade_array(size, ll, ul);
    double* z = create_zero_array(size);

    forward_substitute(z, A, y, size);
    back_substitute(x, A, z, size);

    delete[] y;
    delete[] z;
    for (int i = 0; i < size; i++)
        delete[] A[i];
    delete[] A;
    A = nullptr; y = nullptr; z = nullptr;

    return x;
}

int main(int argc, char** argv)
{
    if (argc != 2){
        cout << "Enter the Total no of Grid Points as parameter!!" << endl;
        return 0;
    }
    
    int n = atoi(argv[1]);
    int size = n;

    double ll = 0; double ul = 3;
    double* x = Derivative_Pade_Scheme(size, ll, ul);
    cout << "The Derivative of sin(5x) within the range " << ll << " to " << ul << " is :" << endl;
    printArray(x, size);
    cout << "Step size is " << n << endl;

    // Packing Data to plot
    ofstream MyFile("p2a_plot_data.txt");
    double h = (ul - ll)/(n - 1);
    for (int i = 0; i < n; i++)
        MyFile << (i * h) << ", " << x[i] << endl;

    delete[] x;
    x = nullptr;

    return 0;
}