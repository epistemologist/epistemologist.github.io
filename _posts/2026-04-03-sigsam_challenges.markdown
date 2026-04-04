---
layout: post
title:  "SIGSAM Challenges"
categories: 
---

# Problem 1
What is a 4 significant digit approximation to the condition number of the 256x256 Hilbert matrix?

 - the $n \times n$ *Hilbert matrix* is the matrix $H = (h_{ij})_{\substack {1 \le i \le n \\ 1 \le j \le n}}$ where $h_{ij} = 1/(i+j-1)$
 - the *condition number* of a $n \times n$ matrix $A$ is defined to be $\|A\|_2 \cdot \|A^{-1}\|_2$ with $\| \cdot \|_2$ is the [spectral norm](https://mathworld.wolfram.com/SpectralNorm.html)

We first attempt to solve this problem with standard libraries:
```python
import numpy as np

N = 256
H = np.array([[1./(i+j-1) for i in range(1, N+1)] for j in range(1, N+1)])
np.linalg.cond(H, p=2)
```

With this code, we get that the condition number $c \; \dot= \; 1.9488855527466815 \cdot 10^{20}$ - however, we do not know how many digits are correct.

## Some theory

We have for a real $n \times n$ matrix $A$, we have that $$ \| A \|_2 = \sup_{|x|_2 \ne 0} \frac{|Ax|_2}{|x|_2} = (\lambda_{\text{max}}(A^T A))^{1/2}$$

(Here, for a vector $v = (v_1 \cdots v_n) \in \mathbb{R}^n$, $|v|_2 = (\sum_{i=1}^n v_i^2)^{1/2}$ is the standard Euclidean norm and $\lambda_{\text{max}}(M)$ represents the largest eigenvalue of $M$)

Our plan of attack is as follows:
 - we can calculate $\lambda_{\text{max}}$ via [power iteration](https://en.wikipedia.org/wiki/Power_iteration); this is useful as it is less expensive than other methods like calculating the characteristic polynomial
 - we have an [explicit closed form](https://mathoverflow.net/questions/47561/deriving-inverse-of-hilbert-matrix) of $H^{-1}$; we have that $(H^{-1})_{ij} = (-1)^{i+j}(i+j-1) {{n+i-1}\choose{n-j}}{{n+j-1}\choose{n-i}}{{i+j-2}\choose{i-1}}^2$
 - we calculate $c = \|H \|_2 \cdot \| H^{-1}\|_2$ with various working precisions and take agreeing digits as confirmation of precision

We first generate $H$ and $H^{-1}$:
```python
import sympy
from math import comb as C

N = 256
H = sympy.Matrix([[sympy.Rational(1, i+j-1) for i in range(1, N+1)] for j in range(1, N+1)])
H_inv = sympy.Matrix( [[pow(-1, i+j)*(i+j-1)*C(N+i-1,N-j)*C(N+j-1, N-i)*C(i+j-2,i-1)**2 for i in range(1, N+1)] for j in range(1, N+1)] )
I = H_inv @ H; assert I == sympy.eye(N)
```

Immediately, we can see that the condition number given by `numpy` is not accurate - note that the condition number of a matrix and it's inverse should be the same; however, we have:

```python
In [58]: np.linalg.cond( np.array(H, dtype=float) )
Out[58]: np.float64(1.9488855527466815e+20)

In [59]: np.linalg.cond( np.array(H_inv, dtype=float) )
 ** On entry to DLASCL parameter number  4 had an illegal value
 ** On entry to DLASCL parameter number  5 had an illegal value
Out[59]: np.float64(inf)

In [60]: np.linalg.cond( np.linalg.inv( np.array(H, dtype=float) ) )
Out[60]: np.float64(7.398387590891174e+18)
```

Similarly, `mpmath.cond` fails to calculate the condition number of this matrix:
```python
In [65]: mpmath.cond( mpmath.matrix(H) )
---------------------------------------------------------------------------

/usr/local/lib/python3.12/dist-packages/mpmath/matrices/linalg.py in cond(ctx, A, norm)
    572         if norm is None:
    573             norm = lambda x: ctx.mnorm(x,1)
--> 574         return norm(A) * norm(ctx.inverse(A))
    575
    576     def lu_solve_mat(ctx, a, b):

/usr/local/lib/python3.12/dist-packages/mpmath/matrices/linalg.py in inverse(ctx, A, **kwargs)
    300             n = A.rows
    301             # get LU factorisation
--> 302             A, p = ctx.LU_decomp(A)
    303             cols = []
    304             # calculate unit vectors and solve corresponding system to get columns

/usr/local/lib/python3.12/dist-packages/mpmath/matrices/linalg.py in LU_decomp(ctx, A, overwrite, use_cache)
    140             ctx.swap_row(A, j, p[j])
    141             if ctx.absmin(A[j,j]) <= tol:
--> 142                 raise ZeroDivisionError('matrix is numerically singular')
    143             # calculate elimination factors and add rows
    144             for i in xrange(j + 1, n):

ZeroDivisionError: matrix is numerically singular
```

The following code calculates $c$ using `mpmath` for extended precision.

```python
def calc_condition_number(prec: int):
    import mpmath
    mpmath.mp.dps = prec
    # Convert H, H^-1 to mpmath matrices with given precision
    H_mp = mpmath.matrix(H)
    H_inv_mp = mpmath.matrix(H_inv)

    # Calculate largest eigenvalue using power iteration
    import random
    random.seed(0)
    def power_iteration(M: mpmath.matrix, max_iter=1000):
        norm = lambda v: mpmath.sqrt(sum(map(lambda x: x**2, v)))
        v = mpmath.matrix([random.random() for _ in range(M.rows)])
        curr_eigenval = 0
        for _ in tqdm( range(max_iter) ):
            v_new = M @ v
            v = v_new / norm(v_new)
            # Use Rayleigh quotient to calculate max eigenvalue
            new_eigenval = ( v.T @ M @ v )[0] /(v.T @ v)[0]
            if abs(new_eigenval - curr_eigenval) < mpmath.mp.eps:
                break
            curr_eigenval = new_eigenval
        return curr_eigenval
    c = power_iteration(H_mp.T @ H_mp)
    c *= power_iteration(H_inv_mp.T @ H_inv_mp)
    return mpmath.sqrt(c)
```

Using this code, we have the following estimates of $c$:

|prec|c|
|-|-|
|5|1.7659529e+389|
|10|1.76595432877e+389|
|20|1.7659543287839574798389e+389|
|50|1.7659543287839574798421271894754308224698704011142166e+389|
|100|1.76595432878395747984212718947543082246987040111421779530875785293223497320682247670477517227e+389|

We therefore have with some confidence that $\boxed{c \approx 1.7659\cdot10^{389}}$

# Problem 2
What is $\int_1^6 x^{x^x} dx$ to 7 significant digits?

Note that $f(x) = x^{x^x}$ is superexponential and increases rapidly on this interval - for instance, $f(6) = 6^{6^6} \approx 2.659 \cdot 10^{36305} \implies \int_1^6 f(x) dx \le (6-1) f(6) \approx 1.33 \cdot 10^{36306} $

![Plot of f](./2_plot.png)

## First Attempt
As the integral is dominated by its behavior near $x=6$, we break the integral up into two parts, writing $$\int_1^6 f(x) = \int_1^t f(x) + \int_t^6 f(x)$$ with $\left| \int_1^t f(t) \right| \le (t-1) f(t) \le \varepsilon$. We evaluate the second integral with the trapezoid rule; we have that

$$ \int_a^b f(x) \approx \frac{b-a}{N}\cdot \left(\frac{f(a)+f(b)}{2} + \sum_{k=1}^{N-1} f\left(a + \frac{k(b-a)}{N} \right) \right) $$

with error $|E|\le \frac{(b-a)^3}{12N^2}\max_{t \in (a,b)}|f''(t)|$.
```python
from decimal import Decimal, getcontext
from typing import Callable
from tqdm import tqdm

def f(x: Decimal) -> Decimal:
        return x**(x**x)

# Generated with sympy.pycode
LN = lambda x: Decimal(x).ln()
def f_double_prime(x: Decimal) -> Decimal:
        return x**x*x**(x**x)*(x**x*((LN(x) + 1)*LN(x) + 1/x)**2 + (LN(x) + 1)**2*LN(x) + 2*(LN(x) + 1)/x + LN(x)/x - 1/x**2)

# Use trapezoid rule to integrate f on a,b
def integrate(f: Callable, a: Decimal, b: Decimal, N: int) -> Decimal:
        out = (b-a) / N * (f(a) + f(b) / 2)
        for k in tqdm( range(1, N) ):
                out += f(a + k*(b-a)/N) * (b-a)/N
        # Note f'' increasing on (1,6), so error is <= (b-a)^3/12N^2 f''(b)
        return (out, (b-a)**3 / (12*N**2) * f_double_prime(b) )

print( integrate(f, Decimal("5.97"), Decimal("6"), 10**8) )
```

```
100%|█████████| 999999/999999 [06:22<00:00, 2611.22it/s]
(Decimal('1.102669808844761440471708212E+36300'), Decimal('3.479564484930620929550513623E+36298'))
```

With this approach, we get ~2 accurate decimal places - this is due to the fact that $| f''(6) | \approx 10^{36316} \approx 10^{-11}f(6)$. With more subdivisions, we may be able to get 1-2 more digits, but we need another approach.

## Integration by Parts
We use a standard approach in analyzing the asymptotic behavior of integrals - namely, [integration by parts](https://www.youtube.com/watch?v=fz9lDRvgHt8). Let $f(x) = e^{g(x)} \implies g(x) = x^x \ln x$. We have, repeatedly integrating by parts that:

$$\begin{aligned} \int_A^B e^{f(x)} dx &= \int_A^B (e^{f(x)})' \cdot \frac{1}{f'(x)} dx &\text{(integrate by parts)} \\ &= \left.\left( \frac{e^{f(x)}}{f(x)}\right)\right|_{x=A}^{x=B} - \int_A^B e^{f(x)} \left( \frac{1}{f'(x)} \right) ' dx   \\ &= \left.\left( \frac{e^{f(x)}}{f(x)}\right)\right|_{x=A}^{x=B} + \int_A^B e^{f(x)} \left( \frac{f''(x)}{(f'(x))^2} \right)  dx \\ &= \left.\left( \frac{e^{f(x)}}{f(x)}\right)\right|_{x=A}^{x=B} + \int_A^B (e^{f(x)})' \left( \frac{f''(x)}{(f'(x))^3} \right) dx \\ &=  \left.\left( \frac{e^{f(x)}}{f(x)}\right)\right|_{x=A}^{x=B} +  \left. \left(\frac{e^{f(x)} f''(x)}{f'(x)^3} \right)\right|_{x=A}^{x=B} - \int_A^B e^{f(x)} \left( \frac{f''(x)}{(f'(x))^3} \right)' dx  \end{aligned}$$

Upon each application of integration by parts, each integral decreases in absolute value (see [here](https://dl.acm.org/doi/pdf/10.1145/274888.274890) for more details). We automate this process using `sympy` to calculate the derivatives:

```python
import sympy
from sympy import Symbol, Function

x = Symbol('x')
f = Function('f')

Integral = sympy.integrals.integrals.Integral

def integrate_by_parts(integral: Integral, v, depth):
    if depth == 0:
        return integral
    # We have \int u dv = uv - \int v du
    sym, a, b = integral.limits[0]
    expr = integral.args[0]
    dv = sympy.diff(v, sym)
    u = expr / dv
    du = sympy.diff(u, sym)
    return (u*v).subs({sym: b}).expand() - (u*v).subs({sym: a}).expand() - integrate_by_parts( Integral(v*du, (sym, a, b)), v, depth=depth-1)


def calc_integral(d):
    x = Symbol('x')
    f = Function('f')
    I = sympy.Integral(sympy.exp(f(x)), (x,1,6))
    integral_res = integrate_by_parts(I, sympy.exp(f(x)), depth=d)
    terms, error = [], None
    for term in integral_res.args:
        if type(term) == Integral or any([type(i) == Integral for i in term.args]):
            error = term
        else:
            terms.append(term)
    print(f"terms: {terms}, error:{error}")
    integral_estimate = sum([term.replace(f, lambda x: x**x*sympy.log(x)).doit().evalf() for term in terms])
    return integral_estimate.n(), error
```

We have the following estimates of the integral:
| number of applications of integration by parts | estimate of integral|
|-|-|
|1|1.10265158296803e+36300|
|2|1.10266499903801e+36300|
|3|1.10266499936193e+36300|
|4|1.10266499936194e+36300|
|5| 1.10266499936194e+36300|

We can therefore say with some confidence that the integral $\int_1^6 x^{x^x} \approx 1.10266449993619 \cdot 10^{36300}$

# Problem 3
What is $$\sum_{n=1}^\infty (n^\pi + n^2 + n^{\sqrt2} + 1)^{-1/3}$$ to 14 significant digits?

Note as the $n^{\text{th}}$ term is $O(n^{-\pi/3}) \approx O(n^{-1.047})$, this sum converges very slowly. Additionally, we have that

$$ S = \sum_{n=1}^{\infty} (n^\pi + n^2 + n^{\sqrt2} + 1)^{-1/3} \le \sum_{n=1}^\infty n^{-\pi/3} = \zeta(\pi/3) \approx 21.768$$

Using standard libraries, we can see that this sum is very slowly convergent - for example, [`mpmath.nsum`](https://mpmath.org/doc/current/calculus/sums_limits.html?highlight=nsum#mpmath.nsum) gives $S \approx 7.614$ with an error that convergence has not been achieved:

```python
import mpmath
from mpmath import mpf, pi, sqrt; mpmath.mp.dps = 20
mpmath.nsum(
    lambda n: (n**pi + n**2 + n**(sqrt(2)) + 1)**mpf("-1/3"),
    [1, mpmath.inf],
    method="r+s+e", # Richardson + Shanks + Euler-Maclaurin
    verbose=True
) 
"""
Ran out of precision for Richardson
Shanks error: 0.000837997
Euler-Maclaurin error: 0.1
Warning: failed to converge to target accuracy

mpf('7.6144829452690711377934')
"""
```

Similarly, [`scipy.nsum`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.nsum.html) gives $S \approx 21.193$ and a warning that convergence has not been reached.

```python
import numpy as np
from scipy.integrate import nsum

nsum(lambda n: (n**np.pi + n**2 + n**np.sqrt(2) + 1)**(-1./3), 1, np.inf)

"""
success: False
  status: -2 # Numerical integration reached its iteration limit; the sum may be divergent.
     sum: 21.19274566368473
   error: 2.201728029049385e-06
    nfev: 1081393
"""
```


We instead opt for PARI-GP's [`sumnum`](https://pari.math.u-bordeaux.fr/dochtml/html/Sums__products__integrals_and_similar_functions.html#se:sumnum) function which uses Euler-Maclaurin summation under the hood:
```sh
? calculate_sum(prec) = {
    default(realprecision, prec);
    sumnum(n=1, (n^Pi + n^2 + n^(sqrt(2)) + 1)^(-1/3)) * precision(1., 14)
};
? foreach(vector(12, i, 50*i), x, print(x, " ", calculate_sum(x)))
50 21.189305632800655988
100 21.193233488570487398
150 21.193240269998838078
200 21.193240377521559000
250 21.193240377708693550
300 21.193240377711535662
350 21.193240377711540934
400 21.193240377711540944
450 21.193240377711540944
500 21.193240377711540944
550 21.193240377711540944
600 21.193240377711540944
```

With the above code, we estimate the sum to be $\approx 21.193240377711540944$

## Digression: Complex Analysis
We state the following theorem without proof:

**Theorem**: Let $f(z)$ be a sufficiently nice analytic function with $f(z) = O(z^{-\alpha})$ for some $\alpha > 1$ as $|z| \to \infty$. We have that $$ \sum_{k=1}^\infty f(k) = \frac{1}{2 \pi i} \int_C f(z) \pi \cot{(\pi z)}$$ where the contour $C$ is such that $C$ runs from $\infty$ to $\infty$, starting in the upper half plane and crossing the real line between 0 and 1, having to its left all positive integers

A proof of this utilizes the residue theorem and can be found [here, Theorem 3.6](www-m3.ma.tum.de/bornemann/challengebook/)

It turns out that our summand is sufficiently nice to apply this theorem: it remains to choose a contour $C$ that satisfies the above theorem and yields good convergence.

## An Example
Note that we can model our sum via $S' = \sum_{n =1}^\infty n^{-\pi/3} = \zeta(\pi/3)$ - we illustrate the above theorem by choosing the contour $C$ parametrized by $Z(t) = \frac{1}{2} - it$; we therefore have that
$$ \begin{aligned} S' &= \frac{1}{2i} \int_C z^{-\pi/3} \cot(\pi z) dz \\ &= \frac{1}{2i} \int_{z = 1/2+i \infty}^{z = 1/2- i\infty } z^{-\pi/3} \cot(\pi z) dz & z=\frac{1}{2} - it \\ &= \frac{1}{2i}  \int_{-\infty}^\infty \left(\frac{1}{2} - it\right)^{-\pi/3} \cot\left( \pi (1/2-it)\right) \cdot (-i dt) \\  &= -\frac{1}{2} \underbrace{ \int_{-\infty}^\infty \left(\frac{1}{2} - it\right)^{-\pi/3} \cot\left( \pi (1/2-it)\right) dt }_{I(t)} \end{aligned}$$

We can evaluate this integral in PARI/GP:
```gp
? calc_integral(prec) = {
    default(realprecision, prec);
    intnum(x = -oo, oo, (1/2-I*x)^(-Pi/3)*cotan(Pi*(1/2-I*x))) / (2*I) * precision(1., 14);
}
? foreach( vector(10, i, 20*i), i, print(i, " ", abs(calc_integral(i))))    
20 21.736563013800747076
40 21.764696595221773520
60 21.767723975165699842
80 21.768128736195472930
100 21.768174858565164674
120 21.768180646396210492
140 21.768181381459356978
160 21.768181471008903874
180 21.768181482096866264
200 21.768181483428099030
```
This suggests that **increasing the working precision by 20 digits yields 1 more correct digit of the answer** - not very good convergence!

Note that we have $\cot(\pi(1/2 \pm it)) \to \mp i$ as $t \to \infty$ and therefore:
$$\begin{aligned} |I(t)| &= |1/2 - it|^{-\pi/3} \cdot |\cot(\pi(1/2 - it))| \\ &\sim |t|^{-\pi/3} \cdot |\cot(\pi(1/2 - it))| \\ &= |t|^{-\pi/3} \approx |t|^{-1.047}  \end{aligned}$$

which is the same order of convergence as the sum.

We instead choose the contour $Z(t) = 1 - 1/2 \cosh(t) - it$ - for $t \in (-\infty, \infty)$, this parametrizes a parabola-like curve with $Z(\infty) = Z(-\infty) = \infty$ and $Z(0) = 1/2$. The following code calculates both the original sum $S$ and $S'$ with this method to 10 digits accuracy with this contour:

```python
import mpmath
from mpmath import cosh, sinh, cot, pi, sqrt
mpmath.mp.dps = 20

Z = lambda t: 1 - 0.5*cosh(t) - 1j*t
dZ = lambda t: -0.5 * sinh(t) - 1j

# Test summation for zeta(pi/3)
f = lambda z: z**(-pi/3)
integrand = lambda t: 1/(2j) * f(Z(t)) * cot(pi*Z(t)) * dZ(t)
M = 10**5
print("Integral calculation of zeta(pi/3):", mpmath.quad(integrand, [-M, -sqrt(M), 0, sqrt(M), M] ))
print("Actual value of zeta(pi/3):", mpmath.zeta(pi/3))

# Calculate actual sum
g = lambda z: z**(-pi/3) * (1 + z**(2 - pi) + z**(sqrt(2)-pi) + z**(-pi))**(-1/3)
integrand = lambda t: 1/(2j) * g(Z(t)) * cot(pi*Z(t)) * dZ(t)

print("Value of sum S: ", mpmath.quad(integrand, [-M, -sqrt(M), 0, sqrt(M), M]))
```

```
$ python3 script.py 
Integral calculation of zeta(pi/3): (21.768181483791584034 - 3.059147123664685429e-29j)
Actual value of zeta(pi/3): 21.768181483610468891
Value of sum S:  (21.193240377892622022 - 3.059147123664685429e-29j)
```

# Problem 4
What is the coefficient of $x^{3000}$ in the polynomial $$(x+1)^{2000}(x^2+x+1)^{1000}(x^4+x^3+x^2+x+1)^{500}$$ to 13 significant digits?

Note that most computer algebra systems nowadays can calculate the given polynomial exactly quite quickly - the following PARI/GP command calculates the above coefficient in less than a second:

```sh
$ time printf "default(parisizemax, 500M); \n 1. * polcoeff((x+1)^2000 * (x^2+x+1)^1000 * (x^4+x^3+x^2+x+1)^500, 3000)" | gp -q 2>/dev/null
3.9739422655800430396966762632992860446 E1426

real    0m0.247s
user    0m0.224s
sys     0m0.024s

```

We can ask a generalization of this question - define $F(k)$ to be the coefficient of $x^{6k}$ in the polynomial $p(x)= (x+1)^{4k} (x^2+x+1)^{2k} (x^4+x^3+x^2+x+1)^k$ - how high can we calculate $F(k)$?

[Standard libraries](https://gist.github.com/epistemologist/51a98fa1ae61ff2372e6316d573d686d) are able to calculate up to $k = 5000$ in less than 10 seconds:
![](4_plot1.png)

The general shape of the above plot confirms the $O(n^2)$ time complexity of naive polynomial multiplication algorithm and the $O(n \log n)$ time complexity of FFT-based multiplication algorithms used by modern libraries. 

However, we also have an example of the more general phenomenon that a single time complexity is often not enough to yield important information - for example, comparing [Sympy's implementation of polynomial multiplication](https://github.com/sympy/sympy/blob/16fa855354eb7bcabd3fe10993841e03b1382692/sympy/polys/densearith.py#L735) and [flint's implementation of polynomial multiplication](https://github.com/flintlib/flint/blob/fe1d6c6966a163d658ed1afe6a69ac7142e91940/src/fmpz_poly/mul.c), we can see different algorithms are used depending on the polynomial degree. These design choices ultimately explain the fact that Sympy is slower than Flint for small $k$ but comparable for larger $k$.

We showcase two techniques for calculating $F(k)$ for higher $k$

## Reconstruction via Chinese Remainder Theorem
Calculating the above polynomial in $\mathbb{Z}[x]$ is slow mainly due to the size of the polynomial's coefficients. Note that since all of the coefficients of $x$ in $p(x)$ are positive, we have that 

$$F(k) \le \text{sum of all coefficients} = p(1) = 2^{4k} \cdot 3^{2k} \cdot 5^k = 720^k$$

We can remedy this by instead calculating the polynomial in $\mathbb{Z}_p[x]$ for a small prime $p$:

![alt text](4_plot2.png)

We then can calculate $F(k) \bmod p$ for several small primes in parallel and then use the above bound with the Chinese Remainder Theorem to reconstruct the actual value of $F(k)$.

![alt text](4_diagram.png)

We illustrate this by calculating $F(50000)$. Since we know that $F(50000) \le 720^{50000}$, we calculate $F(50000) \bmod p$ for $\frac{\log_{2}(720^{50000})}{64} \approx 7416$ 64-bit primes in parallel and then reconstruct the value of $F(50000)$ with the Chinese Remainder Theorem:

```python
from concurrent.futures import ProcessPoolExecutor
from math import log2
from sympy import nextprime, poly, GF
from sympy.abc import x
from sympy.ntheory.modular import crt
from tqdm import tqdm

import sys; sys.set_int_max_str_digits(10**6)

k = 50000
NUM_PRIMES = 1 + int( log2(720 ** k ) / 64 ) 

# Calculate first `NUM_PRIMES` 64-bit primes to use in CRT calculation
print(f"[+] Calculating {NUM_PRIMES} 64-bit primes")
primes = [nextprime(2**64)]
while len(primes) < NUM_PRIMES:
        primes.append(nextprime(primes[-1]))

# Note: We get massive speedup with python-flint being installed
# Otherwise, this is rather slow
def calc_f(k, MOD):
    p1 = x+1
    p2 = x**2+x+1
    p3 = x**4+x**3+x**2+x+1
    p = poly(p1, domain=GF(MOD))**(4*k) * poly(p2, domain=GF(MOD))**(2*k) * poly(p3, domain=GF(MOD))**k
    return ( p.coeff_monomial(x**(6*k)) % MOD, MOD)

# Worker function to pass to ProcessPoolExecutor
# Because lambdas can't be pickled?
def worker(p):
        return calc_f(k, p)

# Calculate F(k) mod p in parallel
print(f"[+] Calculating F({k}) mod primes")
results = []
with ProcessPoolExecutor(max_workers=32) as executor:
        futures = executor.map( worker, primes )
        for future in tqdm(futures, total=len(primes)):
                results.append(future)

# Reconstruct with CRT
print(f"[+] Reconstructing F({k})...")
res = str( crt([i[1] for i in results], [i[0] for i in results], check=False)[0] )
print(f"F({k}) = {res[:10]}...{res[-10:]} ({len(res)} digits)")
```

```sh
$ time python3 solve_crt.py 
[+] Calculating 7416 64-bit primes
[+] Calculating F(50000) mod primes
100%|███| 7416/7416 [09:31<00:00, 12.97it/s]
[+] Reconstructing F(50000)...
F(50000) = 3612727222...0386752768 (142864 digits)

real    10m7.741s
user    303m3.718s
sys     2m2.420s
```

For more on algorithms of this kind, see [*Modern Computer Algebra* by von zur Gathen and Gerhard](https://www.cambridge.org/core/books/modern-computer-algebra/DB3563D4013401734851CF683D2F03F0) Chapter 5.

## Guessing A Closed Form
Instead of calculating all of the coefficients of $p(x)$ for large $k$ - we instead calculate only the coefficient of $x^{6k}$. We use the standard notation for the *coefficient of* operator:
$$ f(x) = \sum_{k \ge 0} a_k x^k \iff [x^k]f(x) = a_k $$
We have:
$$\begin{aligned}
F(k) &= [x^{6k}] \left( (x+1)^{4k} (x^2+x+1)^{2k} (x^4+x^3+x^2+x+1)^k \right) \\
&= [x^{6k}] \left( \left( \frac{1-x^2}{1-x} \right)^{4k}  \left( \frac{1-x^3}{1-x} \right)^{2k}  \left( \frac{1-x^5}{1-x} \right)^{k} \right)
\end{aligned}$$

Using the identities $(1-x^k)^\alpha = \sum_{i} \binom\alpha i (-1)^i x^{ik}$ and $(1-x)^{-k} = \sum_i \binom{i+k-1}{i} x^i$, we have:

$$\begin{aligned}
F(k) &= [x^{6k}] \left( (1-x^2)^{4k} \cdot (1-x^3)^{2k} \cdot (1-x^5)^k \cdot (1-x)^{-7k} \right) \\
&= [x^{6k}] \left( \left( \sum_a \binom{4k}{a} (-1)^a x^{2a}  \right) \left( \sum_b \binom{2k}{b} (-1)^b x^{3b} \right) \left( \sum_{c} \binom{k}{c} (-1)^c x^{5c} \right) \left( \sum_d \binom{d+7k-1}{d} x^d\right) \right) \\
&= \sum_{2a+3b+5c+d = 6k} \binom{4k}{a} \binom{2k}{b} \binom{k}{c} \binom{d+7k-1}{d} (-1)^{a+b+c} \\
&= \sum_{a,b,c} \binom{4k}{a} \binom{2k}{b} \binom{k}{c} \binom{13k-5c-3b-2a-1}{6k-5c-3b-2a} (-1)^{a+b+c}
\end{aligned}$$

This representation as a summation suggests that $F(k)$ may have a closed form or at least a simple recurrence. We attempt to [guess such a closed form with SageMath](https://ask.sagemath.org/question/68486/code-for-guessing-formula-for-integer-sequence/).

```python
sage: R.<x> = PolynomialRing(QQ)
sage: def F(k):
....:     p = (x+1)**(4*k) * (x**2+x+1)**(2*k) * (x**4+x**3+x**2+x+1)**k
....:     return p.coefficients()[6*k]
....: 
sage: seq = [F(k) for k in range(1, 200)]
sage: C = CFiniteSequences(QQ)
sage: C.guess(seq)
0 # sequence is likely not a linear recurrence
sage: from ore_algebra import OreAlgebra, guess
sage: guess(seq, OreAlgebra(QQ['n'], 'Sn'))
(-14799649383664094776320*n^18 - 893418722765459777802240*n^17 - 25303276015022778224540544*n^16 - 446643747105724431860812416*n^15 - 5505131487708177396732824256*n^14 - 50302707434088834373891146144*n^13 - 353048933103657334648924383432*n^12 - 1945559669575179990564767745648*n^11 - 8531616194900294571565011262728*n^10 - 29985136718410847763005218836912*n^9 - 84638318375841503388654077669304*n^8 - 191356940734680945876224201246736*n^7 - 343838242573284787274385667270296*n^6 - 484185531204263313823356848446704*n^5 - 522389764729600243646463785481120*n^4 - 416506138766744488482698887723200*n^3 - 231001942653247713118369912512000*n^2 - 79494141176489236394912978880000*n - 12768833561365149181778592000000)*Sn^4 + (10574637192146663913281680*n^18 + 633077727020755349722412920*n^17 + 17786105932552033654844595546*n^16 + 311517391757721333283482229125*n^15 + 3810841686472640247836090704503*n^14 + 34569441387458156224700609118555*n^13 + 240934823267935284828144508236139*n^12 + 1318830360958250836249615637939425*n^11 + 5746077817330178802134937200227689*n^10 + 20070496098027461734162772290446825*n^9 + 56317961501869956151578076053081999*n^8 + 126609904186122016808629263922171830*n^7 + 226275031561584005999221216743326428*n^6 + 317006598963145241143115137017188200*n^5 + 340360969443972621490842549437752416*n^4 + 270127659596378974023771699999261120*n^3 + 149169165908252073561231390531993600*n^2 + 51124214259713168758195729996032000*n + 8180541032410806173197596172800000)*Sn^3 + (58398173013172057869656320*n^18 + 3466957042689058721344305920*n^17 + 96590400172111838003633846784*n^16 + 1677671710737815554878205251504*n^15 + 20353229710997674251347967705516*n^14 + 183111329700203435724280949090436*n^13 + 1265792085650926368892623993044512*n^12 + 6872714833724188095447466899730772*n^11 + 29705133421005492755275457787287208*n^10 + 102941476004912548625692594623233388*n^9 + 286622975121427160534479463209155024*n^8 + 639486050687516269351594466298810924*n^7 + 1134423037161820089706622844462617356*n^6 + 1577849673724025352948134909509093856*n^5 + 1682240676040629035477200206743471280*n^4 + 1326077082836724563798305324166899200*n^3 + 727508944054623423147210249013152000*n^2 + 247777189760112493690451068197120000*n + 39410782358696625425788184601600000)*Sn^2 + (928200918220647453075200*n^18 + 54640923284178326304393600*n^17 + 1509146355395251895735401440*n^16 + 25979950995010086449458906320*n^15 + 312328444839098654136129792000*n^14 + 2783952831002425364128712541720*n^13 + 19063821850845675269549999705600*n^12 + 102523347973162090447534120721880*n^11 + 438869191432713748818720144882960*n^10 + 1506208364195816419837795938321240*n^9 + 4153384483714401059532296366152320*n^8 + 9178079038310490152558594637207000*n^7 + 16128387142219032127535679608857040*n^6 + 22226963184194261584173266950649040*n^5 + 23488124738338018640414868140341440*n^4 + 18359966834750828694285798294859200*n^3 + 9993982070093155276517161539552000*n^2 + 3379692571500045640851525893760000*n + 534240922865020313596415193600000)*Sn - 18271172078597647872000*n^18 - 1066443733331026652160000*n^17 - 29187192926698356341318400*n^16 - 497582520806954905994016000*n^15 - 5919788784242169915806412000*n^14 - 52179712883448861784377600000*n^13 - 353057180452312841604866935800*n^12 - 1874448733476776108690127186000*n^11 - 7913881495389601307512780875000*n^10 - 26760580203200407866383087610000*n^9 - 72624302482928912035338405328200*n^8 - 157750788490836646176485159310000*n^7 - 272128552299042841758581966061000*n^6 - 367619484864671334095707608030000*n^5 - 380204049615786105676151842317600*n^4 - 290363527673780474678425000488000*n^3 - 154131861228572325908223341280000*n^2 - 50724620693765637030651830400000*n - 7785443427588636312170304000000
```