#pragma once
#include<cmath>
#define B 1
#define RATIO 0.1
#define RAN 2
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


double relu(double x) {
    return (x > 0) ? x : 0;
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


double abs_test(double x) {
    return (x > 0) ? x : -x;
}


double circle_relu(double x) {
    if (x < -RAN) {
        return 0;
    }
    if (x > RAN / (sqrt(2))) {
        return x;
    }
    else {
        double deg = 67.5 * M_PI / 180.0;
        double radius = RAN * tan(deg);
        double ret = radius - sqrt((RAN * RAN * deg * deg) - (x + RAN) * (x + RAN));
        return ret;
    }
}


double square_plus(double x) {
    return (x + sqrt(x * x + B)) / 2;
}

double mix_relu(double x) {
    return (square_plus(x) * RATIO) + (silu(x) * (1 - RATIO)); 
}


double gelu(double x) {
    return 0.5 * x * (1.0 + std::erf(x / std::sqrt(2.0)));
}

double gelu_and_sqplus(double x) {
    return (square_plus(x) * RATIO) + (gelu(x) * (1 - RATIO)); 
}