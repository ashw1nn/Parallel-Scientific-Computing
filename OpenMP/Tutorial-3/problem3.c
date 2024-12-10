#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int main()
{
    int n = 4;

    int a = 10;
    int b = 20;
    static int c = 40;
    int d = 60;

    # pragma omp threadprivate(c) 
    #pragma omp parallel for default(none) shared(n) firstprivate(a) lastprivate(b) private(d) num_threads(4)
    for (int i = 0; i < n; i++)
    {
        printf("Thead no: %d, Num of threads %d\n", omp_get_thread_num(), omp_get_num_threads());

        printf("a before %d\n", a);
        printf("b before %d\n", b);
        printf("c before %d\n", c);
        printf("d before %d\n", d);
        a = 100;
        b = 200;
        c = 400;
        d = 600;
        printf("a after %d\n", a);
        printf("b after %d\n", b);
        printf("c after %d\n", c);
        printf("d after %d\n", d);
        printf("Hello World\n");
    }

    printf("a outer %d\n", a);
    printf("b outer %d\n", b);
    printf("d outer %d\n", d);
    
    return 0;
}
