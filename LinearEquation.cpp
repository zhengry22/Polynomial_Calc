#pragma once
#include<iostream>
#include"LinearEquation.h"
//#include <Eigen/Dense>
using namespace std;

/*
    This program is used to measure the efficiency of my own 
    LEQ class and Eigen
*/

double coeff_for_matrix[4][4] = {
    {2, -1, -1, 1},
    {1, 1, -2, 1},
    {4, 6, 2, -2},
    {3, 6, -9, 7}
};

double coeff_for_result[4] = {2, 4, 4, 9};

int main() {
    LinearEquation<double> lq(4, 4);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            lq.setMatrix(coeff_for_matrix[i][j], i, j);
        }
    }
    for (int i = 0; i < 4; i++) {
        lq.setResult(coeff_for_result[i], i);
    }
    cout << "Before elimination: " << endl;
    lq.show();
    lq.gaussian_elimination();
    cout << "After elimination: " << endl;
    lq.show();
    return 0;
}