#pragma once
#include<iostream>
#include"LinearEquation.h"
#include <Eigen/Dense>
using namespace std;

/*
    This program is used to measure the efficiency of my own 
    LEQ class and Eigen
*/

double coeff_for_matrix[3][6] = {
    {0, 1, 2, 3, 4, 0},
    {0, 2, 4, 0, 2, 1},
    {0, 3, 6, 1, 4, 1}
};

double coeff_for_result[3] = {10, 9, 15};

int main() {
    LinearEquation<double> lq(3, 6);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            lq.setMatrix(coeff_for_matrix[i][j], i, j);
        }
    }
    for (int i = 0; i < 3; i++) {
        lq.setResult(coeff_for_result[i], i);
    }
    cout << "Before elimination: " << endl;
    lq.show();
    lq.gaussian_elimination();
    cout << "After elimination: " << endl;
    lq.show();
    return 0;
}