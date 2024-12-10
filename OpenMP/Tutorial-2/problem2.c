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

void trapezoidal_rule(int n, double a, double b, double *result)
{
  double h, x, total; /* private  */
  int i;
  h = (b - a) / n;

  total = (func(a) + func(b)) / 2.0;
  
  printf("Num of threads is %d, Thread no. is %d \n", omp_get_num_threads(), omp_get_thread_num());
  
  #pragma omp parallel for reduction(+:total)
  for (i = 1; i <= n - 1; i++)
  {
    printf("Num of threads is %d, Thread no. is %d \n", omp_get_num_threads(), omp_get_thread_num());
    x = a + i * h;
    total += func(x);
  }
  printf("Num of threads is %d, Thread no. is %d \n", omp_get_num_threads(), omp_get_thread_num());
  total = total * h;
  *result = total;
}

int main(int argc, char *argv[])
{
  double a, b, final_result;
  int n;

  n = 32; /* number of trapezoids.. */
  a = 1;  /* shared  */
  b = PI;
  final_result = 0.0; /* shared  */

  trapezoidal_rule(n, a, b, &final_result);
  printf("\n The area under the curve between 0 to PI is equal to %lf \n\n", final_result);
  printf("\n The error = %lf \n\n", final_result - ACTUALRESULT);

  return 0;
}