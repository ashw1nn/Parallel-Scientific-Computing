#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define N 600

int main()
{
    // Declaring Matrix, spacing, and limits
    clock_t start = clock();
    double Matrix_A[N][N], Function[N], Matrix_b[N], Matrix_LU[N][N], Matrix_g[N], Matrix_x[N], Matrix_Xplot[N];
    double x_start = 0.0, x_end = 3.0;
    int n = N;
    double h = (x_end - x_start) / (n - 1);
    int i, j;

#pragma acc parallel loop num_gangs(1000)
    for (i = 0; i < n; i++)
    {
        Function[i] = sin(5 * (i * h));
        Matrix_Xplot[i] = 5 * cos(5 * (i * h));
    }

    // Filling in Boundary Conditions for Matrix b
    Matrix_b[0] = (1.0 / h) * ((-5.0 / 2.0) * Function[0] + 2 * Function[1] + 0.5 * Function[2]);
    Matrix_b[n - 1] = (1.0 / h) * ((5.0 / 2.0) * Function[n - 1] - 2 * Function[n - 2] - 0.5 * Function[n - 3]);

    // Filling in Boundary Conditions for Matrix A
    Matrix_A[0][0] = 1.0;
    Matrix_A[0][1] = 2.0;
    Matrix_A[n - 1][n - 2] = 2.0;
    Matrix_A[n - 1][n - 1] = 1.0;

#pragma acc parallel loop num_gangs(1000)
    for (i = 0; i < n - 2; i++)
    {
        Matrix_A[0][i + 2] = 0.0;
        Matrix_A[n - 1][i] = 0.0;
    }

    // Populating Matrix A and b
    for (i = 1; i < n - 1; i++)
    {
        // #pragma acc parallel loop
        for (j = 0; j < n; j++)
        {
            if (j == i)
                Matrix_A[i][j] = 4.0;
            else if (j == i - 1 || j == i + 1)
                Matrix_A[i][j] = 1.0;
            else
                Matrix_A[i][j] = 0.0;
        }
        Matrix_b[i] = (3.0 / h) * (Function[i + 1] - Function[i - 1]);
    }

    // LU Decomposition

    // Filling in the starting values , row 1
    Matrix_LU[0][0] = Matrix_A[0][0];
    Matrix_LU[0][1] = Matrix_A[0][1];
    for (i = 2; i < n; i++)
        Matrix_LU[0][i] = 0.0;

    // Populating the rest of the matrix
    for (i = 1; i < n; i++)
    {
        // #pragma acc parallel loop
        for (j = 0; j < n; j++)
        {
            Matrix_LU[i][j] = 0.0;
        }
        if (i != (n - 1))
            Matrix_LU[i][i + 1] = Matrix_A[i][i + 1]; // c

        Matrix_LU[i][i - 1] = (Matrix_A[i][i - 1] / Matrix_LU[i - 1][i - 1]);          // l
        Matrix_LU[i][i] = Matrix_A[i][i] - (Matrix_LU[i][i - 1] * Matrix_A[i - 1][i]); // u
    }

    // Applying Lg = b equations to solve for yhat matrix
    Matrix_g[0] = Matrix_b[0];
    // #pragma acc parallel loop copyin(Matrix_LU[0:N][0:N], Matrix_b[0:N]) copy(Matrix_g[0:N])
    for (i = 1; i < n; i++)
    {
        Matrix_g[i] = Matrix_b[i] - (Matrix_LU[i][i - 1] * Matrix_g[i - 1]);
    }

    // Backtracing to find matrix x using Ux = g
    Matrix_x[n - 1] = Matrix_g[n - 1] / Matrix_LU[n - 1][n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        Matrix_x[i] = (Matrix_g[i] / Matrix_LU[i][i]) - ((Matrix_LU[i][i + 1] / Matrix_LU[i][i]) * Matrix_x[i + 1]);
    }

    // Printing
    for (int i = 0; i < n; i++)
    {
        printf("%f \n", Matrix_x[i]);
    }

    // Output
    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    printf("%f \n", seconds);

    return 0;
}

/*1. Plot the numerical vs analytical for n=100 and num_gangs = 10
2. Plot the time taken for n=1000 and num_gangs = 10,100,1000
*/