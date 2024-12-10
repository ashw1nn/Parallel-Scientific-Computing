#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

vector<vector<double>> addTwoMatrices(vector<vector<double>> matrix1, vector<vector<double>> matrix2) {
    vector<vector<double>> res;
    for (int i = 0; i < matrix1.size(); i++)
    {
        vector<double> vec;
        for (int j = 0; j < matrix1.size(); j++)
        {
            vec.push_back((matrix1[i][j] + matrix2[i][j]));
        }
        res.push_back(vec);
    }
    return res;
}


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
}

int main(){
    int s1;
    cout << "Enter the array 1 size ";
    cin >> s1;
    vector<vector<double>> A;
    int num = 1;
    for (int i=0; i<s1; i++) {
        vector<double> vec;
        for (int j = 0; j < s1; j++)
        {
            vec.push_back(pow(num,2));
            num++;
        }
        A.push_back(vec);
	}

    vector<vector<double>> AT = transposeMatrix(A);
    vector<vector<double>> summ = addTwoMatrices(A, AT);
    if (summ == transposeMatrix(summ))
    {
        cout << "Success" << endl;
    }
    
    return 0;
}
