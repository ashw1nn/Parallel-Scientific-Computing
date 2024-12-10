#include <iostream>
#include <math.h>       /* pow */
#include <time.h>       /* clock_t */
// #include <cstdio>
using namespace std;
#ifdef _OPENMP
#include <omp.h>
#endif

void MulTwoMatrices (double (&mat1)[100][100], double (&mat2)[100][100], double (&c)[100][100], int m, int n, int t) {
	
    for (int i = 1; i < m; i ++) {
		for (int j = 0; j < n; j ++) {
			c[i][j] = 0;
			for (int k = 0; k < m; k ++) {
				c[i][j] += mat1[i][k]*mat2[k][j];
			}
		}
	}
	return;
}

void printMatrix (double mat[][100], int m, int n) {
	
	for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j++) {
			printf("%f\t",mat[i][j]);
		}
		cout << "\n";
	}
	return;
}


int main () {

	const double PI = 3.14;
	double mat1[100][100], mat2[100][100];
	double sum[100][100], prod[100][100];
    
	// Get matrix size
	int m, n;
	cout << "Enter the number of rows (<100) and cols in matrix ";
    cin >> m >> n;

	for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j++) {
			mat1[i][j] = (double) (rand() % 100);
			mat2[j][i] = (double) (rand() % 100);
		}
	}
	

//2	
	clock_t t; 
	t = clock();
	MulTwoMatrices(mat1,mat2,sum,m,n, 2);
	t = clock()-t;
    printf ("2 Thread Addition took %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);


//4
	t = clock();
	MulTwoMatrices(mat1,mat2,sum,m,n, 4);
	t = clock()-t;
    printf ("4 Thread Addition took %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

//8
	t = clock();
	MulTwoMatrices(mat1,mat2,sum,m,n, 8);
	t = clock()-t;
    printf ("8 Thread Addition took %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

//20
	t = clock();
	MulTwoMatrices(mat1,mat2,sum,m,n, 20);
	t = clock()-t;
   	printf ("20 Thread Addition took %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

	// printMatrix(mat1,m,n);
	// cout << "\n\n";
	// printMatrix(mat2,m,n);
	// cout << "Sum of the matrices: \n\n";
	// printMatrix(sum,m,n);

	return 0;
}