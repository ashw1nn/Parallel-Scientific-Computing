#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

vector<double> multiplyTwoMatrices(vector<vector<double>> a, vector<double> b) {
    if (a[0].size() != b.size())
    {
        cout << "Cant multiply!!" << endl;
        return b;
    }

    vector<double> c(a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        for (int k = 0; k < a.size(); k++) {
            c[i] += a[i][k] * b[k];
        }
    }
    return c;
}

vector<vector<double>> createMatrix(int s1) {
    vector<vector<double>> a1;
    int num = 1;
    for (int i=0; i<s1; i++) {
        vector<double> vec;
        for (int j = 0; j < s1; j++)
        {
            vec.push_back(num);
            num++;
        }
        a1.push_back(vec);
	}
    return a1;
}

vector<double> create1DArr(int s1) {
    vector<double> a1;
    int num = 1;
    for (int i=0; i<s1; i++) {
        a1.push_back(num);
        num++;
	}
    return a1;
}


void printArr(vector<double> mat) {
    for (int i = 0; i < mat.size(); i++) { 
        cout << mat[i] << " ";  
    }
    cout << endl;
    return;
}


int main()
{
    int arr[4] = {256, 512, 1024, 2048};
    vector<double> times;
    for (int i = 0; i < 4; i++)
    {
        vector<vector<double>> a1 = createMatrix(256);
        vector<double> c1 = create1DArr(256);
        clock_t t;
        t = clock();
        vector<double> res = multiplyTwoMatrices(a1, c1);
        t = clock()-t;
        times.push_back(t);
    }

    printArr(times);

    return 0;
}