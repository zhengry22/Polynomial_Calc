#pragma once
#include<cmath>
using namespace std;

inline double my_exp(double x) { 
    // This method is more faster than the one in cmath
    x = 1.0 + x / 1024;   
    x *= x; x *= x; x *= x; x *= x;   
    x *= x; x *= x; x *= x; x *= x;   
    x *= x; x *= x;   
    return x; 
}


double sigmoid(double x) {
    return 1.0 / (1.0 + my_exp(-x));
}


double silu(double x) {
    return x * sigmoid(x);
}


double silu_derivative(double x) {
    return sigmoid(x) + x * sigmoid(x) * (1 - sigmoid(x));
}


double silu_second_derivative(double x) {
    return 2 * sigmoid(x) * (1 - sigmoid(x)) + x * sigmoid(x) * (1 - sigmoid(x)) * (1 - 2 * sigmoid(x));
}