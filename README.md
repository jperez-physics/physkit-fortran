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

## Examples

![Damped Oscillator Plots](example/damped_oscillator/damped_oscillator_plots.png)

## Build

This project uses [fpm](https://fpm.fortran-lang.org/) (Fortran Package Manager).

```bash
fpm build       # build the library
fpm test        # run all tests in test/
fpm run --example <name>   # run an example (damped_oscillator, gamma_function)
```


## Latest Features

- **Linear Algebra Expansion**: Integrated QR Decomposition and Eigenvalue solvers in `physkit_linalg`.
- **Special Functions Visualization**: Added Python-based plotting tools for the Gamma function.
- **Advanced Physics Examples**: Implemented solvers for the **Damped Oscillator**.
- **Enhanced Documentation**: Applied Doxygen-compatible comments across all modules for better maintainability.

## Current Features and Roadmap

Here is a list of the currently implemented modules and the planned roadmap for future additions. Features marked with 🟢 are functional and available.

| Module | Feature | Status |
| :--- | :--- | :--- |
| `physkit_constants` | Physical and mathematical constants list | 🟡 WIP |
| `physkit_linalg` | Vectors: Dot product | 🟢 Done |
| `physkit_linalg` | Vectors: Norm of a vector | 🟢 Done |
| `physkit_linalg` | Vectors: Vector normalization | 🟢 Done |
| `physkit_linalg` | Vectors: Cross product | 🟢 Done |
| `physkit_linalg` | Matrices: Matrix-vector multiplication | 🟢 Done |
| `physkit_linalg` | Matrices: Matrix-matrix multiplication | 🟢 Done |
| `physkit_linalg` | Matrices: Hadamard product | 🟢 Done |
| `physkit_linalg` | Matrices: Zero matrix | 🟢 Done |
| `physkit_linalg` | Matrices: Identity matrix | 🟢 Done |
| `physkit_linalg` | Matrices: Determinant | 🟢 Done |
| `physkit_linalg` | Matrices: Trace | 🟢 Done |
| `physkit_linalg` | Matrices: Gram Schmidt orthogonalization | 🟢 Done |
| `physkit_linalg` | Matrices: Matrix inverse | ⚪ Planned |
| `physkit_linalg` | Matrices: Pseudoinverse | ⚪ Planned |
| `physkit_linalg` | Advanced: LU decomposition | 🟢 Done |
| `physkit_linalg` | Advanced: Eigenvalues and eigenvectors | 🟢 Done |
| `physkit_linalg` | Advanced: QR decomposition | 🟢 Done |
| `physkit_linalg` | Advanced: Cholesky decomposition | ⚪ Planned |
| `physkit_numerical` | Nonlinear: Bisection method | 🟢 Done |
| `physkit_numerical` | Nonlinear: Newton-Raphson method | 🟢 Done |
| `physkit_numerical` | Nonlinear: Secant method | 🟢 Done |
| `physkit_numerical` | Derivative: Forward difference | 🟢 Done |
| `physkit_numerical` | Derivative: Backward difference | 🟢 Done |
| `physkit_numerical` | Derivative: Central difference | 🟢 Done |
| `physkit_numerical` | Derivative: Central second difference | 🟢 Done |
| `physkit_numerical` | Derivative: Higher-order derivatives | ⚪ Planned |
| `physkit_numerical` | Derivative: Partial derivatives | ⚪ Planned |
| `physkit_numerical` | Integral: Rectangular rule | 🟢 Done |
| `physkit_numerical` | Integral: Trapezoidal rule | 🟢 Done |
| `physkit_numerical` | Integral: Simpson’s rule | 🟢 Done |
| `physkit_numerical` | Integral: Composite Simpson’s rule | 🟢 Done |
| `physkit_numerical` | Integral: Adaptive Simpson's | 🟢 Done |
| `physkit_numerical` | Interpolation: Linear interpolation | ⚪ Planned |
| `physkit_numerical` | Interpolation: Polynomial interpolation | ⚪ Planned |
| `physkit_numerical` | Interpolation: Cubic splines | ⚪ Planned |
| `physkit_numerical` | Interpolation: Least squares fitting | ⚪ Planned |
| `physkit_numerical` | Series: Discrete summations | 🟢 Done |
| `physkit_numerical` | Series: Taylor series expansion | ⚪ Planned |
| `physkit_numerical` | Error/Stability: Truncation error estimation | ⚪ Planned |
| `physkit_numerical` | Error/Stability: Stability analysis | ⚪ Planned |
| `physkit_numerical` | Error/Stability: Convergence checks | ⚪ Planned |
| `physkit_ode` | Euler method | 🟢 Done |
| `physkit_ode` | Runge-Kutta 2 (RK2) | 🟢 Done |
| `physkit_ode` | Runge-Kutta 4 (RK4) | 🟢 Done |
| `physkit_ode` | Velocity Verlet for 2nd-order | 🟢 Done |
| `physkit_ode` | Adaptive Runge-Kutta | ⚪ Planned |
| `physkit_special` | Factorial | 🟢 Done |
| `physkit_special` | Combinations and permutations | 🟢 Done |
| `physkit_special` | Gamma function | 🟢 Done |
| `physkit_special` | Beta function | 🟢 Done |
| `physkit_special` | Bessel functions | ⚪ Planned |
| `physkit_special` | Legendre polynomials | ⚪ Planned |
| `physkit_special` | Hypergeometric functions | ⚪ Planned |
| `physkit_newtonian_mechanics` | Displacement | 🟢 Done |
| `physkit_newtonian_mechanics` | Distance | 🟢 Done |
| `physkit_newtonian_mechanics` | Kinetic energy | 🟢 Done |
| `physkit_newtonian_mechanics` | Linear momentum | 🟢 Done |
| `physkit_newtonian_mechanics` | Angular momentum | 🟢 Done |
| `physkit_newtonian_mechanics` | Work | 🟢 Done |
| `physkit_fourier` | Discrete Fourier Transform (DFT) | ⚪ Planned |
| `physkit_fourier` | Fast Fourier Transform (FFT) | ⚪ Planned |
| `physkit_fourier` | Inverse transform | ⚪ Planned |
| `physkit_fourier` | Signal filtering and smoothing | ⚪ Planned |

*(Note: The roadmap and structure may change as the project evolves.)*

## Project philosophy

Instead of writing independent programs for each assignment or simulation, physkit focuses on building a structured toolkit that evolves over time.

This project represents the transition from learning programming to developing scientific computing tools for physics.

## Status

Early development.

This project will expand continuously as part of my academic formation in physics and computational science.

## License

MIT License


