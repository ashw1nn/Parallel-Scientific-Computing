#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define N 		100
#define PROBZERO	33

int main() {
	int a[N][N];	// cities X localities
	int x[N], 	// sum of infected people for all the localities of a city
	    y[N];	// number of localities in a city

	srand(1);

	// #pragma acc parallel loop
	for (int ii = 0; ii < N; ++ii) {
		x[ii] = 0;
		for (int jj = 0; jj < N; ++jj) {
			a[ii][jj] = rand() % N;	//ii * N + jj;
			if (rand() % 100 < PROBZERO) a[ii][jj] = 0;	// add some zeros to make life interesting.
		}
	}

	
	// #pragma acc serial
	// {
	#pragma acc parallel loop collapse(2) reduction(+:x, y) copy(x[0:N], y[0:N]) copyin(a[0:N][0:N])
	for (int ii = 0; ii < N; ++ii) {
		for (int jj = 0; jj < N; ++jj) {
			x[ii] += a[ii][jj];
			if (a[ii][jj] > 0) (y[ii])++;
		}
	}
	// }

	// for (int i = 0; i < N; i++)
	// {
	// 	printf("%d ", x[i]);
		
	// }
	// printf("\n");

	// for (int i = 0; i < N; i++)
	// {
	// 	printf("%d ", y[i]);
	// }
	// printf("\n");
	
	

	for (int ii = 0; ii < N; ++ii)
		printf("%.0f%% ", x[ii] * 100.0 / (y[ii] * N));
	printf("\n");

	return 0;
}
