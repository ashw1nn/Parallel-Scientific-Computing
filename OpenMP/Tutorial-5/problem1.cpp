#include <iostream>
#include <math.h>       /* pow */
#include <time.h>       /* clock_t */
using namespace std;


void create_dense_matrix(double (*matrix)[100], int n)
{
    for (int i = 0; i < n; i++)
    {   
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                *(*(matrix + i) + j) = (i + 1) + (j + 1);
            else if (i == 0 && j == n - 1)
                *(*(matrix + i) + j) = 1.0;
            else if (i == n - 1 && j == 0)
                *(*(matrix + i) + j) = 2.0*n - 1.0;
            else
                *(*(matrix + i) + j) = 1.0/n;
        }
    }
}


void populate_array(int n, double* array)
{
    srand(1);
    for (int i = 0; i < n; i++)
    {
        *(array + i) = (double) (rand() % 100);
    }
}

void print_array(double* array, int n)
{
    for (int i = 0; i < n; i++)
    {
        cout << array[i] << " ";
    }
    cout << endl;
}


double* sub_vec(double *vec1, double *vec2, int n, double* result)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}

void printMatrix (double mat[][100], int m, int n)
{
	
	for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j++) {
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}
	return;
}


double norm(double* vec, int n)
{
    double normm = 0;
    for (int i = 0; i < n; i++)
    {
        normm += pow(vec[i], 2);
    }
    normm = pow(normm, 0.5);
    return normm;
}

int main()
{
    int n = 100;
    double epsilon = 0.000001;

    // A
    double A[100][100];
    create_dense_matrix(A, n);
    printMatrix(A, 100, 100);
    
    // b
    double b[100] = {0};
    b[0] = 1;
    print_array(b, 100);

    // x
    double x[100] = {0};
    populate_array(100, x);
    print_array(x, 100);
    

    // Jacobi Iterative Method
    double diff[100] = {0};
    double tempX[100] = {0};
    tempX[0] = 1;
    while ((norm(sub_vec(x, tempX, 100, diff), 100)/norm(tempX, 100)) > epsilon)
    // for (int k = 0; k < 1000; k++)
    {
        for (int i = 0; i < n; i++)
        {
            double summ = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    summ = summ + A[i][j]*x[j]; 
            }
            tempX[i] = (b[i] - summ)/A[i][i];  
        }

        // Update the original x array after completing the iteration
        double temp = 0;
        for (int i = 0; i < n; i++)
        {
            temp = x[i];
            x[i] = tempX[i];
            tempX[i] = temp;
        }
    }

    print_array(x, 100);
    

    return 0;
}