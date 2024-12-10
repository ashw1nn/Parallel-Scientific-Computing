#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

void populate_array(int n, double* array){
    for (int i = 0; i < n; i++)
    {
        *(array + i) = (double) (rand() % 100);
    }
}

void print_array(int n, double* array){
    printf("\n");
    printf("[");
    for (int i = 0; i < n; i++)
    {
        printf("%.1f ", array[i]);
    }
    printf("]");
    printf("\n\n\n");
}

void swap(double* p1, double* p2){
    double temp;
    temp = *p1;
    *p1 = *p2;
    *p2 = temp;
}

void odd_even_sort(int n, double* arr, int thread_count){
    #pragma omp parallel num_threads(thread_count)
    for (int pass = 0; pass < n; pass++)
    {
        // printf("Thead no: %d, Num of threads %d\n", omp_get_thread_num(), omp_get_num_threads());
        if (pass%2 == 0)
        {
            #pragma omp for
            for (int i = 0; i < n - 1; i+=2)
            {
                if (pass == 0 && i == 0) printf("In Even Pass Num of threads %d\n", omp_get_num_threads());
                if (arr[i] > arr[i + 1])
                    swap(&arr[i], &arr[i+1]);
            }
        }
        else
        {
            #pragma omp for
            for (int i = 1; i < n - 1; i+=2)
            {
                if (pass == 1 && i == 1) printf("In Odd Pass Num of threads %d\n", omp_get_num_threads());
                if (arr[i] > arr[i + 1])
                    swap(&arr[i], &arr[i + 1]);
            }
        }
    }
}


int main(int argc, char *argv[]){

    int thread_count;
    int n = 100;
    double arr[100];
    populate_array(n, arr);
    print_array(n, arr);

    if (argc == 2)
    {
        thread_count = strtol(argv[1], NULL, 10);
    }
    else
    {
        printf("\n A command line argument other than name of the executable is required...exiting the program..");
        return 1;
    }

// Sorting
    odd_even_sort(n, arr, thread_count);
    print_array(n, arr);

    return 0;
}
