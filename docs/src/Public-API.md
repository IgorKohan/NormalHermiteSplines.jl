# Public API

## API Summary

| Function                            | Description                                                                                        |
|:----------------------------------- |:-------------------------------------------------------------------------------------------------- |
|```prepare```                        |Prepare the spline by constructing and factoring a Gram matrix of the interpolation problem.        |
|```construct```                      |Construct the spline by calculating its coefficients.                                               |
|```interpolate```                    |Prepare and construct the spline.                                                                   |
|```evaluate```                       |Evaluate the spline value at the required locations                                                 |
|```evaluate_one```                   |Evaluate the spline value at the required location                                                  |
|```evaluate_gradient```              |Evaluate gradient of the spline at the required location.                                           |
|```evaluate_derivative```            |Evaluate the 1D spline derivative at the required location.                                         |
|```estimate_accuracy```              |Estimate accuracy of the function interpolation result.                                             |
|```estimate_cond```                  |Estimate the Gram matrix 1-norm condition number.                                                   |
|```estimate_epsilon```               |Estimate the 'scaling parameter' of Bessel potential space the spline being built in.               |
|```get_epsilon```                    |Get the 'scaling parameter' of Bessel potential space the normal spline was built in.               |
|```get_cond```                       |Get the Gram matrix spectral condition number.                                                      |

## Functions
```@docs
prepare
construct
interpolate
evaluate
evaluate_one
evaluate_gradient
evaluate_derivative
estimate_accuracy
estimate_cond
get_epsilon
estimate_epsilon
get_cond
```

## Types

### Bessel potential Space Reproducing Kernels

```@docs
RK_H0
RK_H1
RK_H2
```

### NormalSpline structure

```@docs
NormalSpline
```

## Index
```@index
Order = [:function, :type]
```
