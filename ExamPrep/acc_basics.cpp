#include <iostream>

int main()
{

/***********************************************************************************************************************/
//     int a = 0;

//     #pragma acc parallel num_gangs(2) num_workers(20) vector_length(3) copy(a)
//     {
//         #pragma acc loop gang
//         for (int i = 0; i < 10; i++)
//         {
//             #pragma acc loop worker
//             for (int j = 0; j < 200; j++)
//             {
//                 #pragma acc loop vector
//                 for (int k = 0; k < 300; k++)
//                 {
//                     printf("Hi from C++ %d, %d, %d\n", i, j, k);   
//                 }
                
                
//             }
//         }
//     }

//     printf("CPUa = %d\n", a);

//     #pragma acc wait
/***********************************************************************************************************************/

    #pragma acc parallel num_gangs(10) async(1)
    {
    printf(" Hello \n");
    printf(" World \n");
    printf(" Bye \n");
    }

    printf(" Host \n");

    #pragma acc parallel num_gangs(3)
    {
    printf(" One \n");
    printf(" Two \n");
    }

    #pragma acc wait

    return 0;
/***********************************************************************************************************************/
/*
    int a[10] = {0};
    #pragma acc kernels
    {
    printf("One \n");
    for (int i = 0 ; i < 1200; ++i)
    {
    a[i] = i;
    printf("first loop %d \n", i);
    }
    printf("Two: a[9] = %d \n", a[9]);
    for (int i = 0; i < 1200; ++i)
    {
    a[i] *= 2;
    printf("second loop %d \n", i);
    }
    printf("Three: a[9] = %d \n", a[9]);
    }*/
/***********************************************************************************************************************/
/*
    int a[10] = {0};
    #pragma acc kernels
    {
    printf("One \n");
    for (int i = 0 ; i < 10; ++i)
    {
    // a[i] = i;
    printf("first loop %d \n", i);
    }
    printf("Two \n");
    for (int i = 0; i < 10; ++i)
    {
    // a[i] *= 2;
    printf("second loop %d \n", i);
    }
    printf("Three \n");
    }
    return 0;*/

}

