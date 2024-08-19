# Polyomial_Calc
    * This is the lib we build for approximating the aim function denoted by the function pointer `func`.

## Taylor
* Using taylor's equation: 
$f^{(n)}(x) \approx \frac {(-1)^0C_{n}^{0}f(x + n\Delta{x}) + (-1)^1C_{n}^{1}f(x + (n-2)\Delta{x}) + ... + (-1)^{n-1}C_{n}^{n - 1}f(x - (n-2)\Delta{x}) + (-1)^(n)C_{n}^{n}f(x - n\Delta{x})} {(2\Delta{x})^n}$
this may not be very accurate, as the result suggests.

## Remez
* Using an exchange theorem which I don't really know how to implement yet;
