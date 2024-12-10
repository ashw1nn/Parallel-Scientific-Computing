#include <iostream>
#include <math.h>       /* pow */
#include <time.h>       /* clock_t */
using namespace std;
#ifdef _OPENMP
#include <omp.h>
#endif

double** create_dense_matrix(int n)
{
    double** matrix = new double*[n];
    for (int i = 0; i < n; i++)
        matrix[i] = new double[n];
    
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
    return matrix;
}

void create_zero_array(double* arr, int n)
{
    for (int i = 0; i < n; i++)
        arr[i] = 0;
}

void populate_array(int n, double* array)
{
    srand(1);
    for (int i = 0; i < n; i++)
        *(array + i) = (double) (rand() % 100);
}

void print_array(double* array, int n)
{
    for (int i = 0; i < n; i++)
        cout << array[i] << " ";
    cout << endl;
}


double* sub_vec(double *vec1, double *vec2, int n, double* result)
{
    for (int i = 0; i < n; i++)
        result[i] = vec1[i] - vec2[i];
    return result;
}

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


double norm(double* vec, int n)
{
    double normm = 0;
    for (int i = 0; i < n; i++)
        normm += pow(vec[i], 2);
    normm = pow(normm, 0.5);
    return normm;
}


double* JacobiIter(double** A, double* b, double* x, double* prevX, double* diff, int n, int epsilon, int threads)
{
    #pragma omp parallel num_threads(threads)
    while ((norm(sub_vec(x, prevX, n, diff), n)/norm(prevX, n)) > epsilon)
    {
        #pragma parallel for collapse(2)
        for (int i = 0; i < n; i++)
        {
            double summ = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    summ = summ + A[i][j]*x[j]; 
            }
            prevX[i] = (b[i] - summ)/A[i][i];  
        }

        // Update the original x array after completing the iteration
        double temp = 0;
        #pragma parallel for
        for (int i = 0; i < n; i++)
        {
            temp = x[i];
            x[i] = prevX[i];
            prevX[i] = temp;
        }
    }

    return x;
}


int main(int argc, char* argv[])
{
    clock_t t;
    int n = atoi(argv[1]);
    int threads[3] = {2, 4, 8}; // Threads
    double epsilon = 0.000001;

    t = clock();
    // A
    double** A = create_dense_matrix(n);
    // printMatrix(A, n, n);
    t = clock() - t;
    cout << "Time taken to create Matrix A is " << t << endl;

    // b
    double* b = new double[n];
    create_zero_array(b, n);
    b[0] = 1;
    cout << "b : ";
    print_array(b, n);

    // x
    double *x = new double[n];
    populate_array(n, x);
    cout << "x_i : ";
    print_array(x, n);
    

    // Jacobi Iterative Method
    double* diff = new double[n]; 
    double* prevX = new double[n];
    create_zero_array(diff, n);
    create_zero_array(prevX, n);
    prevX[0] = 1;

    t = clock();
    x = JacobiIter(A, b, x, prevX, diff, n, epsilon, threads[0]);
    t = clock() - t;
    cout << "x_f : ";
    print_array(x, n);
    cout << "Time taken for " <<  threads[0] << " threads is " << t << endl;

    t = clock();
    x = JacobiIter(A, b, x, prevX, diff, n, epsilon, threads[1]);
    t = clock() - t;
    // print_array(x, n);
    cout << "Time taken for " << threads[1] << " threads is " << t << endl;

    t = clock();
    x = JacobiIter(A, b, x, prevX, diff, n, epsilon, threads[2]);
    t = clock() - t;
    // print_array(x, n);
    cout << "Time taken for " <<  threads[2] << " threads is " << t << endl;

    delete[] b;
    delete[] x;
    delete[] prevX;
    delete[] diff;
    for (int i = 0; i < n; i++)
        delete[] A[i];
    delete[] A;
    A = NULL;
    
    return 0;
}