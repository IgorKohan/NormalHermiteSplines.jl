# Selecting a good value of the scaling parameter
 
Approximating properties of the normal spline are getting better with the smaller value of the scaling ('shape') parameter ``\varepsilon``. In a case when value of this parameter is small enough the normal spline become similar to Duchon's ``D^m-spline`` [2]. Details are described in
[Comparison with Polyharmonic Splines](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Relation-to-Polyharmonic-Splines/).

However with decreasing value of the parameter ``\varepsilon`` the condition number of the corresponding problem Gram matrix is increasing and the problem becomes numerically unstable (see "uncertainty principle" of Schaback [4]). The Gram matrix of the interpolating problem even can lost its positive definiteness property if ``\varepsilon`` is small. Also, as it was pointed out in [3] the RBF interpolation with small value of the "shape" parameter may cause the Runge phenomenon i.e. undesirable interpolant oscillations which are most likely observed at the domain border. 

Therefore, when choosing the value of the ``\varepsilon``, a compromise is needed. In practice, it is necessary to choose such value of the scaling parameter that estimation of the Gram matrix condition number is a relatively small number and the value of interpolation result significant digits estimation is a good enough.  As a rule, the heuristic algorithm implemented within the interpolation procedure produces a good estimation of the scaling parameter value (this algorithm applies if the value of the scaling parameter was not provided explicitly in creation of the reproducing kernel object).

The following API functions are useful for selecting a suitable value of the scaling parameter

- ```estimate_cond```
- ```get_epsilon```
- ```estimate_accuracy```
- ```estimate_epsilon```  

As example let's consider interpolation of function ``f (x,y)`` (function ``N7`` from [5])

```math
f (x,y) = \frac{2}{3}cos(10x)sin(10y) + \frac{1}{3}sin(10xy)
```
```@raw html
<img src="../images/parameter-choice/p-cf.png" width="256"/>
```  ```@raw html
<img src="../images/parameter-choice/p-grid.png" width="197"/>
```
sampled on set ``\Chi`` of 100 pseudo-random nodes uniformly distributed on unit square ``\Omega = [0,1]^2``.

```
    using NormalHermiteSplines
    ....
    spline = interpolate(nodes, u, rk)
    ε = get_epsilon(spline)
    κ_1 = estimate_cond(spline)
    significant_digits = estimate_accuracy(spline)
    σ = evaluate(spline, grid)
    ....
```
the complete example code can be found in [Example usage](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Usage/#D-interpolation-case-2/).

The normal splines ``\sigma (x,y)`` are evaluated on a regular grid ``E_S= \{(x_i, y_i)\, , \ i=1, \dots, S \}`` of ``S = 101 × 101`` points. The interpolation error is measured by means of the Root Mean Square Error (``RMSE``)

```math
  RMSE = \sqrt { \frac{1}{S} \sum_{i=1}^S |f (x_i,y_i) - \sigma (x_i,y_i) |^2 } \ ,
```
and the Maximum Absolute Arror (``MAE``)
```math
  MAE = \max_{1 \le i \le S} |f (x_i,y_i) - \sigma (x_i,y_i) | \ .
```
Estimation of the number of the interpolant 'significant digits' (``SD``) is calculated as
```math
SD = -[\log_{10}(\max_{(x_k, y_k) \in \Chi} |f (x_k,y_k) - \sigma (x_k,y_k) |)] - 1 \ ,
```
here ``(x_k, y_k) \in \Chi`` are interplation nodes and ``[\cdot]`` denotes an integer part of number.

Estimation of the 1-norm condition number ``\kappa_1`` of the interpolation problem Gram matrix is obtained by function ```estimate_cond``` that implements the procedure described in [1]. It requires ``O(N^2)`` operations.

The results of this function interpolation with reproducing kernel ```RK_H0``` are displayed in Table I, results of interpolation with reproducing kernel ```RK_H1``` are displayed in Table II and results of interpolation received with reproducing kernel ```RK_H2``` – in Table III. The first row of each table corresponds to results obtained with automatically selected value of the scaling parameter ``\varepsilon``. 

Notice that a two-dimensional normal spline built with ```RK_H0``` reproducing kernel is an element of the Bessel potential space ``H^{3/2}_\varepsilon (R^2)``, an element of that space can be treated as a bounded continuous function. Correspondingly, a two-dimensional normal spline built with ```RK_H1``` reproducing kernel is an element of the Bessel potential space ``H^{5/2}_\varepsilon (R^2)``, an element of that space can be treated as a bounded with its derivatives continuously differentiable function. Further, a two-dimensional normal spline built with ```RK_H2``` reproducing kernel is an element of the Bessel potential space ``H^{7/2}_\varepsilon (R^2)``, an element of that space can be treated as a bounded with all its derivatives twice continuously differentiable function.  

Table I: Spline constructed with ```RK_H0``` reproducing kernel

|     ε      |``\kappa_1``|      SD      |     RMSE     |     MAE     |
|:----------:|:----------:|:------------:|:------------:|:-----------:|
|**1.7e+00** | **1.0e+05**|      **14**  |  **9.2e-02** |  **6.3e-01**|
|  1.0e+01   |  1.0e+03   |      15      |   1.1e-01    |   4.7e-01   |
|  1.0e+00   |  1.0e+05   |      14      |   9.2e-02    |   6.5e-01   |
|  1.0e-01   |  1.0e+07   |      13      |   9.3e-02    |   6.8e-01   |
|  1.0e-02   |  1.0e+08   |      12      |   9.3e-02    |   6.9e-01   |
|  1.0e-03   |  1.0e+09   |      10      |   9.3e-02    |   6.9e-01   |
|  1.0e-04   |  1.0e+10   |       9      |   9.3e-02    |   6.9e-01   |
|  1.0e-05   |  1.0e+11   |       8      |   9.3e-02    |   6.9e-01   |
|  1.0e-06   |  1.0e+12   |       8      |   9.3e-02    |   6.9e-01   |
|  1.0e-07   |  1.0e+13   |       6      |   9.3e-02    |   6.9e-01   |
|  1.0e-08   |  1.0e+14   |       5      |   9.3e-02    |   6.9e-01   |
|  1.0e-09   |  1.0e+15   |       5      |   9.3e-02    |   6.9e-01   |
|  1.0e-10   |  1.0e+16   |       4      |   9.3e-02    |   6.9e-01   |
|  1.0e-11   |  1.0e+17   |       3      |   9.3e-02    |   6.9e-01   |
|  1.0e-12   |  1.0e+18   |       2      |   9.3e-02    |   6.9e-01   |
|  1.0e-13   |  1.0e+19   |       1      |   9.6e-02    |   7.1e-01   |
|  1.0e-14   |  —         |       —      |   —          |   —         |

These plots show normal splines constructed with ```RK_H0()``` (``\varepsilon = 1.7``),
```RK_H0(1.0e-12)``` and ```RK_H0(1.0e-13)``` reproducing kernels
```@raw html
<img src="../images/parameter-choice/p-cs,0.0,-.png" width="256"/>
``` ```@raw html
<img src="../images/parameter-choice/p-cs,1.0e-12,-.png" width="256"/>
```  ```@raw html
<img src="../images/parameter-choice/p-cs,1.0e-13,-.png" width="256"/>
```
 

Table II: Spline constructed with ```RK_H1``` reproducing kernel

|     ε      |``\kappa_1``|      SD      |     RMSE     |     MAE     |
|:----------:|:----------:|:------------:|:------------:|:-----------:|
|**2.6e+00** | **1.0e+08**|      **13**  |  **4.5e-02** |  **4.1e-01**|
|  1.0e+01   |  1.0e+06   |      14      |   5.4e-02    |   4.7e-01   |
|  1.0e+00   |  1.0e+09   |      11      |   4.4e-02    |   3.7e-01   |
|  1.0e-01   |  1.0e+12   |       8      |   4.4e-02    |   3.7e-01   |
|  1.0e-02   |  1.0e+15   |       5      |   4.4e-02    |   3.7e-01   |
|  1.0e-03   |  1.0e+19   |       2      |   4.4e-02    |   3.6e-01   |
|  1.0e-04   |  —         |       —      |   —          |   —         |

 

Table III: Spline constructed with ```RK_H2``` reproducing kernel

|     ε      |``\kappa_1``|      SD      |     RMSE     |     MAE     |
|:----------:|:----------:|:------------:|:------------:|:-----------:|
|**3.5e+00** | **1.0e+10**|      **12**  |  **2.2e-02** |  **1.8e-01**|
|  1.0e+01   |  1.0e+07   |      13      |   3.4e-02    |   3.2e-01   |
|  1.0e+00   |  1.0e+13   |       9      |   2.3e-02    |   2.7e-01   |
|  1.0e-01   |  1.0e+18   |       4      |   2.5e-02    |   3.4e-01   |
|  1.0e-02   |  —         |       —      |   —          |   —         |


For comparison, in Table III-D we presented the results of interpolation received by using extended precision arithmetic (extended precision ``Double64`` float type from the [DoubleFloats.jl](https://github.com/JuliaMath/DoubleFloats.jl) package has 31 significant decimal digits).  
    
Table III-D: Spline constructed with ```RK_H2``` reproducing kernel using extended Double64 arithmetic

|     ε      |``\kappa_1``|      SD      |     RMSE     |     MAE     |
|:----------:|:----------:|:------------:|:------------:|:-----------:|
|**3.5e+00** | **1.0e+10**|      **28**  |  **2.2e-02** |  **1.8e-01**|
|  1.0e+01   |  1.0e+07   |      29      |   3.4e-02    |   3.2e-01   |
|  1.0e+00   |  1.0e+13   |      25      |   2.3e-02    |   2.7e-01   |
|  1.0e-01   |  1.0e+18   |      20      |   2.5e-02    |   3.4e-01   |
|  1.0e-02   |  1.0e+23   |      15      |   2.5e-02    |   3.4e-01   |
|  1.0e-03   |  1.0e+28   |      10      |   2.5e-02    |   3.4e-01   |
|  1.0e-04   |  1.0e+33   |       5      |   2.5e-02    |   3.4e-01   |
|  1.0e-05   |  —         |       —      |   —          |   —         |
 

As can be seen from these tables, for small values of the scaling parameter ``\varepsilon`` interpolation results have the similar quality, and it make sense to chose such ``\varepsilon`` value that provides "good" values of the result significant digits estimation (at least 7 – 10 digits in most cases) and reasonable estimation of the problem Gram matrix condition number.

**References**

[1] C. Brás, W. Hager, J. Júdice, An investigation of feasible descent algorithms for estimating the condition number of a matrix. TOP 20, 2012.

[2] J. Duchon, Splines minimizing rotation-invariant semi-norms in Sobolev spaces, Lect. Notes in Math., Springer, Berlin, Vol. 571, 1977.

[3] B. Fornberg, J. Zuev, The Runge phenomenon and spatially variable shape parameters in RBF interpolation,
Comput. Math. Appl., Vol.54, No.3, 2007.

[4] R. Schaback, Error estimates and condition numbers for radial basis functions interpolation, Adv. in Comput. Math. 3, 1995.

[5] R. Renka, R. Brown, Algorithm 792: accuracy test of ACM algorithms for interpolation of scattered data in the plane, ACM Transactions on Mathematical Software (TOMS), Vol.25, No.1, 1999.
