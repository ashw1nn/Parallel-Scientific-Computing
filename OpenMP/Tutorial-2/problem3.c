/* OpenMP parallel version of trapezoidal rule */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.14159265358
#define ACTUALRESULT 0.198557

double func(double x)
{
  return (sin(x)) / (2 * (pow(x, 3)));
}

void simpson_rule(int n, double a, double b, double *result, int thread_count)
{
  double h, x, total; /* private  */
  int i;
  h = (b - a) / n;

  total = (func(a) + func(b));
  
  printf("Num of threads is %d, Thread no. is %d \n", omp_get_num_threads(), omp_get_thread_num());
  
  #pragma omp parallel for num_threads(thread_count) reduction(+:total)
  for (i = 1; i <= n - 1; i++)
  {
    printf("Num of threads is %d, Thread no. is %d \n", omp_get_num_threads(), omp_get_thread_num());
    x = a + i * h;
    if (i % 2 == 0)
        total += 2 * func(x);
    else
        total += 4 * func(x);
  }
//   printf("Num of threads is %d, Thread no. is %d \n", omp_get_num_threads(), omp_get_thread_num());
  total = total * (h/3.0);
  *result = total;
}

int main(int argc, char *argv[])
{
  double a, b, final_result;
  int n;
  int thread_count = 1;

  if (argc == 2)
  {
    thread_count = strtol(argv[1], NULL, 10);
  }
  else
  {
    printf("\n A command line argument other than name of the executable is required...exiting the program..");
    return 1;
  }

  n = 32; /* number of trapezoids.. */
  a = 1;  /* shared  */
  b = PI;
  final_result = 0.0; /* shared  */


  simpson_rule(n, a, b, &final_result, thread_count);

  printf("\n The area under the curve between 0 to PI is equal to %lf \n\n", final_result);

  printf("\n The error = %lf \n\n", final_result - ACTUALRESULT);

  return 0;
}



