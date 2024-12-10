#include <iostream>
#include <math.h>       /* pow */
#include <time.h>       /* clock_t */
using namespace std;


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


double* sub_vec(double *vec1, double *vec2, int n)
{
    double* result = new double[n];
    #pragma acc parallel loop create(result[0:n]) copyout(result[0:n])
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

    #pragma acc parallel loop reduction(+:normm)
    for (int i = 0; i < n; i++)
        normm += std::pow(vec[i], 2);
    normm = std::pow(normm, 0.5);
    return normm;
}

double sumOfVectorMultiplication(double* vec1, double* vec2, int n)
{
    double sum = 0;
    // #pragma acc parallel loop copy(sum) reduction(+:sum) present(vec1[0:n], vec2[0:n], n)
    for (int i = 0; i < n; i++)
        sum += vec1[i]*vec2[i];
    return sum;
}


double* JacobiIter(double** A, double* b, double* x, double* prevX, int n, double epsilon)
{
    double* diff = sub_vec(x, prevX, n);
    double error = norm(diff, n)/norm(prevX, n);
    
    #pragma acc enter data copyin(A[0:n][0:n], b[0:n], x[0:n], prevX[0:n], diff[0:n], n) create(error)

    while (error > epsilon)
    {
        double summ;
        cout << "Error: " << error << endl;
        #pragma acc parallel loop present(A[0:n][0:n], b[0:n], x[0:n], n)
        for (int i = 0; i < n; i++)
        {
            summ = sumOfVectorMultiplication(A[i], x, n) - A[i][i]*x[i];
            prevX[i] = (b[i] - summ)/A[i][i];  
        }

        // Update the original x array after completing the iteration

        #pragma acc parallel loop present(x[0:n], prevX[0:n], n) // ccopyout(x[0:n], prevX[0:n])
        for (int i = 0; i < n; i++)
        {
            double temp = x[i];
            x[i] = prevX[i];
            prevX[i] = temp;
        }

        #pragma acc serial present(x[0:n], prevX[0:n], n) copyout(error)
        {
        diff = sub_vec(x, prevX, n), n;
        error = norm(diff, n)/norm(prevX, n);
        }
    }
    

    #pragma acc exit data delete(A[0:n][0:n], b[0:n], prevX[0:n], n) copyout(x[0:n])

    return x;
}


int main(int argc, char* argv[])
{
    clock_t t;
    int n = atoi(argv[1]);
    double epsilon = 0.00001;

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
    double* prevX = new double[n];
    create_zero_array(prevX, n);
    prevX[0] = 1;

    t = clock();
    x = JacobiIter(A, b, x, prevX, n, epsilon);
    t = clock() - t;
    cout << "x_f : ";
    print_array(x, n);
    cout << "Time taken is " << t << endl;

    delete[] b;
    delete[] x;
    delete[] prevX;
    for (int i = 0; i < n; i++)
        delete[] A[i];
    delete[] A;
    A = NULL;
    
    return 0;
}