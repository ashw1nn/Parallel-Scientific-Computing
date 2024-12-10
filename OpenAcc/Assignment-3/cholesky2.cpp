/* serial code for Cholesky decomposition */
/* make sure that the init function setups a  */
/* symmetric and positive definite matrix  */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TYPE float
#define N 10
#define SMALLVALUE 0.001

void initmult(TYPE mat[][N])
{
  for (int ii = 0; ii < N; ++ii)
    for (int jj = 0; jj < N && jj < ii; ++jj)
    {
      mat[ii][jj] = (ii + jj) / (float)N / N;
      mat[jj][ii] = (ii + jj) / (float)N / N;
    }

  for (int ii = 0; ii < N; ++ii)
    mat[ii][ii] = 1.0;
}

void printMat(TYPE a[][N])
{
  for (int ii = 0; ii < N; ++ii)
  {
    for (int jj = 0; jj < N; ++jj)
      printf("%.2f ", a[ii][jj]);
    printf("\n");
  }
}

void cholesky(TYPE a[][N])
{
  #pragma acc enter data copyin(a[0:N][0:N])

  #pragma acc parallel loop num_gangs(1000)
    for (int ii = 0; ii < N; ++ii)
    {
        #pragma acc loop
        for (int jj = 0; jj < ii; ++jj)
        {
            #pragma acc loop
            for (int kk = 0; kk < jj; ++kk)
                a[ii][jj] += -a[ii][kk] * a[jj][kk];
            a[ii][jj] /= (a[jj][jj] > SMALLVALUE ? a[jj][jj] : 1);
        }
        #pragma acc loop
        for (int kk = 0; kk < ii; ++kk)
            a[ii][ii] += -a[ii][kk] * a[ii][kk];
        a[ii][ii] = sqrt(a[ii][ii]);
    }

  #pragma acc exit data copyout(a[0:N][0:N])

  // for (int ii = 0; ii < N; ii++)
  // {
  //   a[ii][ii] = sqrt(a[ii][ii]);
  // }
  
  
}

int main()
{
  clock_t start = clock();
  TYPE a[N][N];

  initmult(a);
  // printMat(a);
  cholesky(a);

  printMat(a);

  clock_t end = clock();
  double total_time = (double)(end - start) / CLOCKS_PER_SEC;
  printf("%f us\n", total_time*1000000);

  
  return 0;
}
