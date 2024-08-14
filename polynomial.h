#pragma once
#include "../src/troy.h"
#include<vector>
#include"SiLU.h"
#include<iostream>
#include<cmath>
#include <thread>
#include <cassert>
using namespace std;
#define EPSILON 0.05
//#define MYDEBUG

template<typename T>
inline T abs(T x) {
    return (x > 0) ? x : -x;
}

template<typename T>
class Polynomial {
private:
    size_t degree;
    vector<T> coeffs;
public:
    Polynomial() {}
    Polynomial(vector<T> &coeff):degree(coeff.size() - 1), coeffs(coeff){}
    inline size_t get_degree() {
        return this->degree;
    }
    inline void check() {
        cout << "Coeffs: ";
        for (auto e: coeffs) {
            cout << e << " ";
        }
        cout << endl;
    }
    inline T get_coeff_by_rank(size_t rank) {
        assert(rank < this->coeffs.size() && "Error: visiting rank that is too big!");
        return this->coeffs[rank];
    }
    inline T get_last_nonzero_coeff() {
        int this_size = this->coeffs.size();
        for (int i = this_size - 1; i >= 0; i--) {
            if (this->coeffs[i] != 0) 
                return this->coeffs[i];
        }
        return 0;
    }
    void prune() {
        for (int i = 0; i < this->coeffs.size(); i++) {
            if (i != 0 && (abs(this->coeffs[i]) * 10 < abs(this->coeffs[i - 1]))) {
                this->coeffs[i] = 0;
            }
        }
    }
};


struct EncryptPolynomial {
    /*
        In homomorphic encryption, we need the polynomial to have integer coeffs.
        But the polynomial we use to approximate ReLU has coefficients in double.
        So, we need to convert it into a integer coeff polynomial with a factor of k.
    */
    Polynomial<long long> poly;
    long long k;
    void show() {
        poly.check();
        cout << "expansion coeff is: " << k << endl;
    }
};


template<typename T>
EncryptPolynomial round_polynomial(Polynomial<T> &poly) {
    EncryptPolynomial ret;
    
    // Calculate k
    T last_coeff = poly.get_last_nonzero_coeff();
    cout << "last coeff is: " << last_coeff << endl;
    long long k = (long long)((T)(1) / last_coeff);
    cout << "k is: " << k << endl;
    ret.k = k;
    
    vector<long long> new_coeff;
    size_t degree = poly.get_degree();
    for (int i = 0; i <= degree; i++) {
        T expand = poly.get_coeff_by_rank(i) * (T)(k);
        long long tmp = (long long)(expand);
        tmp = (((T)(tmp + 1) - expand) < (expand - (T)tmp)) ? tmp + 1 : tmp; 
        new_coeff.push_back(tmp);
    }

    Polynomial p(new_coeff);
    ret.poly = p;

    return ret;
}


template<typename T, typename U>
T get_coeff_taylor(const int deg, T (*func)(U), const U input) {
    /*
        Calculating the coeff for degree `deg`, using the approximation formula:
        $f^{(n)}(x) \approx \frac {(-1)^0C_{n}^{0}f(x + n\Delta{x}) + 
        (-1)^1C_{n}^{1}f(x + (n-2)\Delta{x}) + ... 
        + (-1)^{n-1}C_{n}^{n - 1}f(x - (n-2)\Delta{x}) + (-1)^(n)C_{n}^{n}f(x - n\Delta{x})} {(2\Delta{x})^n}$

        input stands for the x in f(x)

        \Delta{x} = EPSILON
    */ 
    vector<T> func_values;
    vector<long long> comb_coeffs;  
    int cnt = deg; // initially set to be deg
    long long comb = 1;
    for (int i = deg; i >= 0; i--) {
        // Push back f(x + cnt * \delta{x})
        T retval = func((U)(input + (U)(cnt) * EPSILON));
        func_values.push_back(retval);

        // Push back C_{n} ^ {i}
        comb_coeffs.push_back(comb);

        // Change the coeffs
        cnt -= 2;
        comb *= i;
        comb /= (deg - i + 1);
    }
#ifdef MYDEBUG
    cout << "Func values: ";
    for (auto e: func_values) {
        cout << e << " ";
    }
    cout << endl;
    cout << "comb_coeffs: ";
    for (auto e: comb_coeffs) {
        cout << e << " ";
    }
    cout << endl;
#endif
    assert((func_values.size() == comb_coeffs.size()) && "Error: vec not of the same size!");
    assert((func_values.size() == deg + 1) && "Error: vec size not equal to deg + 1!");
    long long denom = 1 << deg;
    T numerator = 0; 
    for (int i = 0; i <= deg; i++) {
        int sign = (i % 2 == 0) ? 1 : -1;
        numerator += (double)sign * (double)comb_coeffs[i] * func_values[i];
    }
    // cout << "numerator: " << numerator << endl;
    numerator /= denom;
    //cout << "numerator: " << numerator << endl;
    for (int i = 1; i <= deg; i++) {
        numerator /= EPSILON;
    }
    //cout << "numerator: " << numerator << endl;
    return numerator;
}


template<typename T, typename U>
class PolyApprox {
protected:
    size_t poly_degree = 0;
    T (*func)(U) = nullptr; // function pointer
public: 
    PolyApprox():poly_degree(0){}
    PolyApprox(size_t deg):poly_degree(deg){}
    PolyApprox(size_t deg, T (*func)(U)): poly_degree(deg), func(func){}
    virtual Polynomial<T> generate_approx(int deg, U input) = 0; 
    // Pure virtual function, since each method has its own way of generating an approximate polynomial
};


/*
    Class Taylor refers to approximating the function via the corresponding taylor series.
*/
template<typename T, typename U>
class Taylor: public PolyApprox<T, U>  {
public:
    Taylor(size_t deg):PolyApprox<T, U>(deg){}
    Taylor(size_t deg, T (*func)(U)): PolyApprox<T, U>(deg, func){}
    //friend get_coeff(int deg, T (*func)(U), U input);
    Polynomial<T> generate_approx(int deg, U input) {
        /*
            To calculate all the coefficients, we could use multi-threads.
            But first we may use single thread instead for conveniece.
        */
        vector<T> coefficients;
        for (int i = 0; i <= deg; i++) {
            T get_val = get_coeff_taylor(i, this->func, input);
            for (int j = 1; j <= i; j++) {
                get_val /= (T)j;
            }
            coefficients.push_back(get_val);
        }
#ifdef MYDEBUG
        cout << "coefficients: ";
        for (auto e : coefficients) {
            cout << e << " ";
        }
        cout << endl;
#endif
        Polynomial<T> ret(coefficients);
        return ret;
    }
};


/*
    Taylor approximation seems not that good, we can try Remez Algorithm instead
*/