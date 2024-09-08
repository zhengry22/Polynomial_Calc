#pragma once
#include<cmath>
#define BETA 1
#define K 1.1
#define RATIO 0.5
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


double silu_with_k(double x) {
    return K * silu(x);
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

double silu_with_bias(double x) {
    if (x < 3 && x > 0.25) {return silu(x) + 0.25;}
    return silu(x);
}

double square_plus(double x) {
    return (x + sqrt(x * x + BETA)) / 2;
}

double mix_relu(double x) {
    return (square_plus(x) * RATIO) + (silu(x) * (1 - RATIO)); 
}


double gelu(double x) {
    return 0.5 * x * (1.0 + std::erf(x / std::sqrt(2.0)));
}


double better_gelu(double x) {
    if (abs(x) >= 1) {return gelu(x);}
    return relu(x);
}

double gelu_and_sqplus(double x) {
    return (square_plus(x) * RATIO) + (gelu(x) * (1 - RATIO)); 
}

double gelu_with_k(double x) {
    if (x > 0 && x < 1) {
        return (2 - x) * gelu(x);
    }
    else if (x > 1 && x < 1.5) {
        return (2.5 - x) * gelu(x);
    }
    else {
        return gelu(x);
    }
}


double tanh(double x) {
    return (exp(x) - exp(-x)) / (exp(x) + exp(-x));
}


double mix_relu_tanh(double x) {
    return mix_relu(x) * (tanh(4 * x - 1) * 0.5 + 0.5);
}