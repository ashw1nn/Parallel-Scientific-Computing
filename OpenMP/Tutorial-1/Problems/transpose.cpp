#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

vector<vector<double>> transposeMatrix(vector<vector<double>> matrix) {
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix.size(); j++)
        {
            if (i>j)
            {
                double temp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] =  temp;
            }
        }
    }
    return matrix;
}


void printMatrix(vector<vector<double>> mat) {
    for (int i = 0; i < mat.size(); i++) { 
        for (int j = 0; j < mat[i].size(); j++) 
            cout << mat[i][j] << " "; 
        cout << endl; 
    }
    return;
}

vector<vector<double>> multiplyTwoMatrices(vector<vector<double>> a, vector<vector<double>> b) {
    if (a.size() != b.size())
    {
        cout << "Cant mul matrices" << endl;
        return a;
    }

    vector<vector<double>> c(a.size(), vector<double>(a.size(), 0));
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < a.size(); j++) {
            for (int k = 0; k < a.size(); k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return c;
}


int main() {

    int s1;
    cout << "Enter the array 1 size ";
    cin >> s1;
    vector<vector<double>> a1;
    int num = 1;
    for (int i=0; i<s1; i++) {
        vector<double> vec;
        for (int j = 0; j < s1; j++)
        {
            vec.push_back(pow(num,2));
            num++;
        }
        a1.push_back(vec);
	}

    int s2;
    cout << "Enter the array 2 size ";
    cin >> s2;
    vector<vector<double>> a2;
    num = 2;
    for (int i=0; i<s2; i++) {
        vector<double> vec;
        for (int j = 0; j < s2; j++)
        {
            vec.push_back(pow(num,2));
            num++;
        }
        a2.push_back(vec);
	}
    
    vector<vector<double>> a3 = multiplyTwoMatrices(a1, a2);
    
    // Transposing Matrices
    vector<vector<double>> a1T = transposeMatrix(a1);
    vector<vector<double>> a2T = transposeMatrix(a2);
    vector<vector<double>> a3T = transposeMatrix(a3);
    vector<vector<double>> RHS = multiplyTwoMatrices(a2T, a1T);

    if (a3T == RHS)
    {
        cout << "LHS:" << endl;
        printMatrix(a3T);
        cout << "Mat 2 trans:" << endl;
        printMatrix(a2T);
        cout << "Mat 1 trans:" << endl;
        printMatrix(a1T);
        cout << "RHS:" << endl;
        printMatrix(RHS);
        cout << endl;
        cout << "Success!!!" << endl;
    }
    

    return 0;
}
