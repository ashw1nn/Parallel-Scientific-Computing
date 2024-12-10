#include <iostream>
#include <omp.h>
using namespace std;

void foo(int* ptr)
{
    (*ptr)++;
    printf("foo: %d\n", *ptr);
}


int main(int argc, char** argv)
{
    static int a = 0;
    int th = atoi(argv[1]);

    #pragma omp threadprivate(a)

    #pragma omp parallel num_threads(th)
    {
    
    // cout << "a(before): " << a << endl;
    a = 20;
    if (omp_get_thread_num() == 0)
    {
        a++;
    }
    cout << "a(after): " << a << " from thread " << omp_get_thread_num() << endl;
    
    foo(&a);
    cout << "a(after foo): " << a << " from thread " << omp_get_thread_num() << endl;

    // for (int i = 0; i < 100; i++)
    {
        // cout << "Hello from thread:" << omp_get_thread_num() << endl;
    }
    }
    

    return 0;
}
