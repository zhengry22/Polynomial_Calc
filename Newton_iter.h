#pragma once
#include <iostream>
#include "SiLU.h"
#include "polynomial.h"
#include <cmath>
#include <functional>
#include <unordered_map>
#define NEWTON_MAX_ITERATION 50
using namespace std;

std::unordered_map<std::function<double(double)>, std::function<double(double)>> derivative_ = {
    {silu, silu_derivative}
};

std::unordered_map<std::function<double(double)>, std::function<double(double)>> derivative_second = {
    {silu, silu_second_derivative}
};


template<typename T, typename U> 
T Newton_iteration(T (*calc_func)(U), const Polynomial<T> &poly, T low, T high) {
    /* Look for 1st and 2nd derivative */
    if (derivative_.find(calc_func) == derivative_.end()) {
        throw std::invalid_argument("Derivative function not found for the given function.");
    }

    if (derivative_second.find(calc_func) == derivative_second.end()) {
        throw std::invalid_argument("2nd Derivative function not found for the given function.");
    }

    T max_point = (low + high) / 2; /* This is the point we choose at first */
    

    for (int i = 0; i < NEWTON_MAX_ITERATION; i++) {
        T derivative_1 = derivative_[calc_func](max_point);
        T derivative_2 = derivative_second[calc_func](max_point);
        Polynomial<T> first_der = poly.derivative();
        Polynomial<T> second_der = poly.second_derivative();
        T poly_result_1 = first_der.get_poly_value(max_point);
        T poly_result_2 = second_der.get_poly_value(max_point);
        max_point = max_point - (derivative_1 - poly_result_1) / (derivative_2 - poly_result_2);
    }
    
    return max_point;
}


template<typename T, typename U>
T 