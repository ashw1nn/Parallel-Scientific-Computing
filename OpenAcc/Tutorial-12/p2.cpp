#include <stdio.h>
#include <stdlib.h>
#include <iostream>

int main(int argc, char **argv)
{
    int numbers[100];
    int N = 100;
    for (int i = 0; i < N; i++)
    {
        numbers[i] = i;
    }
    
    int mynumber = atoi(argv[1]);
    int flag = 0;
    int loc = -1;
#pragma acc parallel
    {
#pragma acc loop worker
        for (int i = 0; i < N; i++)
        {
            printf("Hello world, i = %d\n", i);
            if (numbers[i] == mynumber)
            {
                loc = i;
                flag = 1;
            }
        }
        if (!flag)
            printf("Number not found \n");
        else
            printf("Number found at location %d \n", loc);
    }
}