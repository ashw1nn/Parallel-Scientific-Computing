#include <iostream>
#include <math.h>       /* pow */
#include <time.h>       /* clock_t */
#include <openacc.h>
using namespace std;

void create_zero_matrix(double** matrix, int m, int n)
{
    #pragma acc parallel loop collapse(2) create(matrix[0:m][0:n])
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            matrix[i][j] = 0;
}

void populate_matrix(double** matrix, int m, int n)
{
    srand(1);
    // #pragma acc parallel loop create(matrix[0:m][0:n])
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            matrix[i][j] = (double) (rand() % 100);
}

void copy_matrix(double** mat1, double** mat2, int m, int n)
{
    #pragma acc parallel loop collapse(2)
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            mat2[i][j] = mat1[i][j];
}

void copy_array(double* mat1, double* mat2, int m)
{
    #pragma acc parallel loop present(mat1[0:m], mat2[0:m])
    for (int i = 0; i < m; i++)
        mat2[i] = mat1[i];
}

double norm(double* vec, int n)
{
    double normm = 0;
    #pragma acc parallel loop present(vec[0:n]) reduction(+:normm) firstprivate(n)
    for (int i = 0; i < n; i++)
        normm += std::pow(vec[i], 2);
    normm = std::pow(normm, 0.5);
    return normm;
}


double vec_dot(double* vec1, double* vec2, int n)
{
    double sum = 0;
    #pragma acc parallel loop present(vec1[0:n], vec2[0:n]) reduction(+:sum)
    for (int i = 0; i < n; i++)
        sum += vec1[i]*vec2[i];
    return sum;
}

double* scale_vec(double* vec, int scalar, int n)
{
    double* result = new double[n];
    #pragma acc parallel loop create(result[0:n]) copyout(result[0:n]) present(vec[0:n]) firstprivate(scalar)
    for (int i = 0; i < n; i++)
        result[i] = vec[i]*scalar;
    return result;
}

void sub_vec(double* vec1, double* vec2, int n)
{
    #pragma acc parallel loop present(vec1[0:n], vec2[0:n])
    for (int i = 0; i < n; i++)
        vec1[i] -= vec2[i];
}

double** multiply_matrix(double** mat1, double** mat2, int m, int n)
{
    double** result = new double*[m];
    for (int i = 0; i < m; i++)
    {
        result[i] = new double[n];
    }

    #pragma acc data present(mat1[0:m][0:n], mat2[0:m][0:n]) create(result[0:m][0:n])
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                result[i][j] = 0;
                #pragma acc parallel loop reduction(+:result[i][j])
                for (int k = 0; k < n; k++)
                {
                    result[i][j] += mat1[i][k] * mat2[k][j];
                }
            }
        }
    }
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

int check_same_matrix(double** mat1, double** mat2, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (mat1[i][j] != mat2[i][j])
                return 0;
        }
    }
    return 1;
}


int main(int argc, char *argv[])
{
    int M = 150;
    int N = 100;
    double **A = new double*[M];
    for (int i = 0; i < M; i++)
    {
        A[i] = new double[N];
    }
    // create_zero_matrix(A, M, N);
    populate_matrix(A, M, N);

    #pragma acc enter data copyin(A[0:M][0:N])

    double **Q = new double*[M];
    for (int i = 0; i < M; i++)
    {
        Q[i] = new double[N];
    }
    copy_matrix(A, Q, M, N);

    double **R = new double*[N];
    for (int i = 0; i < N; i++)
    {
        R[i] = new double[N];
    }
    create_zero_matrix(R, N, N);

    #pragma acc enter data copyin(Q[0:M][0:N], R[0:N][0:N])

    for (int i = 0; i < N; i++)
    {
        // cout << "Hi" << endl;
        R[i][i] = norm(Q[i], N);
        for (int j = 0; j < N; j++)
        {
            Q[i][j] = Q[i][j] / R[i][i];
        }

        for (int j = i+1; j < N; j++)
        {
            R[j][i] = vec_dot(Q[j], Q[i], N);
            double* scaled_Q = scale_vec(Q[i], R[i][j], N);
            sub_vec(Q[j], scaled_Q, N);
            delete[] scaled_Q;
        }
    }
    
    double** output = multiply_matrix(Q, R, M, N);
    int is_same_matrix = check_same_matrix(output, A, M, N);

    #pragma acc exit data copyout(A[0:M][0:N], output[0:M][0:N])

    printf("The matrices are %d\n", is_same_matrix);
    cout << "A: " << endl; 
    printMatrix(A, M, N);
    cout << "output: " << endl;
    printMatrix(output, M, N);

    return 0;
}