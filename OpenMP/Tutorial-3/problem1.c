#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

void populate_array(int n, double* array){
    for (int i = 0; i < n; i++)
    {
        *(array + i) = (i + 1);
    }
}

void print_array(int n, double* array){
    for (int i = 0; i < n; i++)
    {
        printf("%f ", array[i]);
    }
    printf("\n");
}

int main()
{
    int n = 100;
    double arr[100];
    populate_array(n, arr);
    // print_array(n, arr);


    double total  = 1;
    for (int i = 1; i <= n; i++)
    {
        total = total+i;
    }
    printf("%f\n", total);

    double res;
    #pragma omp parallel for reduction(-:total)
    for (int i = 0; i < n; i++)
    {
        if (i == 0) printf("Num of threads %d\n", omp_get_num_threads());
        total = total - (i + 1);
    }

    printf("%f\n", total);
}


