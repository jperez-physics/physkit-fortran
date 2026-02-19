# physkit-fortran

A modular computational physics toolkit written in modern Fortran.

## Description

physkit is a personal and academic project aimed at building a reusable computational physics toolkit throughout my academic formation as a physicist.

This repository is intended to grow progressively as new numerical methods, physical models, and computational techniques are learned and implemented.

## Objectives

- Develop reusable scientific computing tools in Fortran
- Implement numerical methods used in physics
- Build simulation-ready modules for future academic work
- Transition from isolated scripts to structured scientific software

## Current Features and roadmap
## 1. Mathematical and physical constants ('physkit_constants')
- There are some currently, but the idea is to add more

## 2. Linear Algebra ('physkit_linalg')

### Vectors
- Dot product 🟢 
- Norm of a vector 🟢
- Vector normalization 🟢
- Cross product 🟢

### Matrices
- Matrix-vector multiplication 🟢
- Matrix-matrix multiplication 🟢
- Hadamard product 🟢
- Zero matrix 🟢
- Identity matrix 🟢
- Determinant
- Trace 🟢
- Matrix inverse
- Pseudoinverse

### Advanced
- Eigenvalues and eigenvectors
- LU decomposition 🟢
- QR decomposition
- Cholesky decomposition


## Numerical Calculus ('physkit_numerical')

### Nonlinear Equations
- Bisection method 🟢
- Newton-Raphson method 🟢
- Secant method 🟢

### Numerical Differentiation
- Forward difference 🟢
- Backward difference 🟢
- Central difference 🟢
- Central second difference 🟢
- Higher-order derivatives
- Partial derivatives

### Numerical Integration
- Rectangular rule 🟢
- Trapezoidal rule 🟢
- Simpson’s rule 🟢
- Composite Simpson’s rule 🟢
- Adaptive Simpson's 🟢

### Interpolation and Fitting
- Linear interpolation
- Polynomial interpolation
- Cubic splines
- Least squares fitting

### Series and Approximations
- Discrete summations
- Taylor series expansion

### Error and Stability
- Truncation error estimation
- Stability analysis
- Convergence checks


## ODE Methods ('physkit_ode')
- Euler method 🟢
- Runge-Kutta 2 (RK2) 🟢
- Runge-Kutta 4 (RK4) 🟢
- Velocity Verlet for 2nd-order 🟢
- Adaptive Runge-Kutta


## Special Functions ('physkit_special')
- Factorial 🟢
- Combinations and permutations 🟢
- Gamma function 🟢
- Beta function 🟢
- Bessel functions
- Legendre polynomials
- Hypergeometric functions


## Fourier Analysis ('physkit_fourier')
- Discrete Fourier Transform (DFT)
- Fast Fourier Transform (FFT)
- Inverse transform
- Signal filtering and smoothing

(This may change)

## Project philosophy

Instead of writing independent programs for each assignment or simulation, physkit focuses on building a structured toolkit that evolves over time.

This project represents the transition from learning programming to developing scientific computing tools for physics.

## Status

Early development.

This project will expand continuously as part of my academic formation in physics and computational science.

## License

MIT License


