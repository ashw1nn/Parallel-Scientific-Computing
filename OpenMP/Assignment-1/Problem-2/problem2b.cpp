#include <iostream>
#include <fstream>
#include <math.h>
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

void printArray(double* arr, int n)
{
    for (int i = 0; i < n; i++)
        cout << arr[i] << " ";
    cout << endl;
}

double f(double x)
{
    return sin(5*x);
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

double* pade_scheme_rd_array(int n, int threads)
{
    // Doesn't create a matrix at all. Only non zero elements are stored to save memory
    double* a = new double[n + 1]; // appending with zeroes to solve index issue
    double* b = new double[n + 1];
    double* c = new double[n + 1];
    double* alpha = new double[n + 1];
    double* beta = new double[n + 1];
    double* y2 = create_pade_array(n, 0, 3);
    double* y = new double[n + 1];
    y[0] = 0;
    for (int i = 1; i < n+1; i++)
        y[i] = y2[i - 1];

    for (int i = 0; i < n + 1; i++)
    {
        b[i] = 4;
        a[i] = 1;
        c[i] = 1;
        alpha[i] = 0;
        beta[i] = 0;
    }
    
    a[0] = 0; b[0] = 0; c[0] = 0;

    a[1] = 0; a[n] = 2;
    b[1] = 1; b[n] = 1;
    c[1] = 2; c[n] = 0;

    double N = ceil(log2(n));
  #pragma omp parallel num_threads(threads) shared(alpha, beta, a, b, c, y)
  {
    for (int k = 0; k <= N; k++)
    {
      #pragma omp for 
        for (int i = 1; i < n + 1; i++)
        {
            if (i >= (int(pow(2, k-1)) + 1) && i <= n)
            {
                alpha[i] = -a[i]/b[i - int(pow(2, k-1))];
            }
            if (i >= 1 && i <= (n - int(pow(2, k-1))))
            {
                beta[i] = -c[i]/b[i + int(pow(2, k-1))];   
            }
            if (i >= (int(pow(2, k)) + 1) && i <= n)
            {
                a[i] = alpha[i] * a[i - int(pow(2, k-1))];
            }
            if (i >= 1 && i <= (n - int(pow(2, k))))
            {
                c[i] = beta[i] * c[i + int(pow(2, k-1))];
            }
            
            b[i] = alpha[i] * c[i - int(pow(2, k-1))] + b[i] + beta[i] * a[i + int(pow(2, k-1))];
            y[i] = alpha[i] * y[i - int(pow(2, k-1))] + y[i] + beta[i] * y[i + int(pow(2, k-1))];
        }
    }
    // cout << "a: ";printArray(a, n);
    // cout << "b: ";printArray(b, n);
    // cout << "c: ";printArray(c, n);
  }
    
    double* x = new double[n];
    for (int i = 1; i < n + 1; i++)
        x[i - 1] = y[i]/b[i];

    delete[] a; a = nullptr;
    delete[] b; b = nullptr;
    delete[] c; c = nullptr;
    delete[] alpha; alpha = nullptr;
    delete[] beta; beta = nullptr;
    delete[] y; y = nullptr;

    return x;
}


int main(int argc, char** argv)
{
    if (argc != 2){
        cout << "Enter the No. of Grid Points as parameters" << endl;
        return 0;
    }

    int size = atoi(argv[1]); double ul = 3; double ll = 0;

    int threads[3] = {2, 4, 8};

    clock_t t, t1, t2, t3;
    t = clock();
    double* derivatives = pade_scheme_rd_array(size, threads[0]);
    t1 = clock() - t;
    cout << "Time taken for " << threads[0] << " threads is " << t1 << " microsecs." << endl;
    
    t = clock();
    pade_scheme_rd_array(size, threads[1]);
    t2 = clock() - t;
    cout << "Time taken for " << threads[1] << " threads is " << t2 << " microsecs." << endl;

    t = clock();
    pade_scheme_rd_array(size, threads[2]);
    t3 = clock() - t;
    cout << "Time taken for " << threads[2] << " threads is " << t3 << " microsecs." << endl;

    clock_t times[3] = {t1, t2, t3};

    // Packing Data to plot
    ofstream MyFile("p2b_plot_data.txt");
    double h = (ul - ll)/(size - 1);
    for (int i = 0; i < size; i++)
        MyFile << (i * h) << ", " << derivatives[i] << endl;

    // Packing Time Data to plot
    ofstream TimeFile("p2b_time_data.txt");
    for (int i = 0; i < 3; i++)
        TimeFile << threads[i] << ", " << times[i] << endl;

    printArray(derivatives, size);

    delete[] derivatives; derivatives = nullptr;

    return 0;
}