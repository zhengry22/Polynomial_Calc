#pragma once
#include"polynomial.h"
#include"SiLU.h"
using namespace std;


int main() {
    Taylor<double, double> taylor(5, silu);
    Polynomial<double> poly = taylor.generate_approx(5, 0);
    poly.prune();
    poly.check();
    EncryptPolynomial my_enc = round_polynomial(poly);
    my_enc.show();
    // We can now get the polynomial that we need
    return 0;
}