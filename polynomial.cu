#pragma once
#include"polynomial.h"
#include"SiLU.h"
using namespace std;


int main() {
    Taylor<double, double> taylor(15, silu);
    Polynomial<double> poly = taylor.generate_approx(15, 0);
    poly.check();
    return 0;
}