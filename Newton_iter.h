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
T Newton_iteration(T (*calc_func)(U), Polynomial<T> &poly, T low, T high) {
    /* Look for 1st and 2nd derivative */
    if (derivative_.find(calc_func) != derivative_.end()) {
        return derivative_[calc_func];
    } else {
        throw std::invalid_argument("Derivative function not found for the given function.");
    }

    if (derivative_second.find(calc_func) != derivative_second.end()) {
        return derivative_second[calc_func];
    } else {
        throw std::invalid_argument("2nd Derivative function not found for the given function.");
    }

    T max_point = (low + high) / 2; /* This is the point we choose at first */

    poly.derivative()

}