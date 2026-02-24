!> @author Jaime (and contributors)
!> @brief Module for numerical methods including differentiation, integration, and root-finding.
!> @details Implementation of common numerical algorithms for real and complex functions,
!>          providing tools for calculus and solving non-linear equations.
module physkit_numerical
    use physkit_constants, only: dp
    implicit none
    private

    ! Public methods
    public :: pk_sum, &
              pk_forward_difference, pk_backward_difference, pk_central_difference, pk_second_central_difference, &
              pk_rectangular_rule, pk_trapezoidal_rule, pk_composite_simpson, &
              pk_adaptative_simpson, &
              pk_bisection_method, pk_newton_raphson, pk_secant_method

    interface pk_sum
        module procedure pk_sum_real
    end interface

    interface pk_forward_difference
        module procedure pk_forward_difference_scalar
        module procedure pk_forward_difference_vector
        module procedure pk_forward_difference_multi
    end interface

    interface pk_backward_difference
        module procedure pk_backward_difference_scalar
        module procedure pk_backward_difference_vector
        module procedure pk_backward_difference_multi
    end interface

    interface pk_central_difference
        module procedure pk_central_difference_scalar
        module procedure pk_central_difference_vector
        module procedure pk_central_difference_multi
    end interface

    interface pk_second_central_difference
        module procedure pk_second_central_difference_scalar
        module procedure pk_second_central_difference_vector
        module procedure pk_second_central_difference_multi
    end interface

    interface pk_composite_simpson
        module procedure pk_composite_simpson_real
        module procedure pk_composite_simpson_complex
    end interface

    interface pk_adaptative_simpson
        module procedure pk_adaptative_simpson_real
        module procedure pk_adaptative_simpson_complex
    end interface

    ! Abstract interfaces for user-defined functions
    abstract interface
        !> @brief Interface for real scalar functions f(x)
        function f(x)
            use physkit_constants, only: dp
            real(dp) :: f
            real(dp), intent(in) :: x
        end function f

        !
        function f_vector(x)
            use physkit_constants, only: dp
            real(dp), dimension(:), allocatable :: f_vector
            real(dp), intent(in) :: x
        end function f_vector

        !> @brief Interface for multi-variable scalar functions f(x1, x2, ..., xn)
        function f_multi(x)
            use physkit_constants, only: dp
            real(dp) :: f_multi
            real(dp), dimension(:), intent(in) :: x
        end function f_multi

        function f_multi_vector(x)
            use physkit_constants, only: dp
            real(dp), dimension(:), allocatable :: f_multi_vector
            real(dp), dimension(:), intent(in) :: x
        end function f_multi_vector

        !> @brief Interface for complex scalar functions f(z)
        function f_complex(x)
            use physkit_constants, only: dp
            complex(dp) :: f_complex
            complex(dp), intent(in) :: x
        end function f_complex
    end interface

contains

    !##################################################
    ! Mathematical Series and Summation Methods
    !##################################################

    !=================================================
    !> @brief Classical mathematical summation from i to n for a real array.
    !> @param arr Array of real numbers.
    !> @param i Starting index.
    !> @param n Ending index.
    !> @return Sum of arr(k) from k=i to n.
    !=================================================
    function pk_sum_real(arr, i, n) result(total_sum)
        real(dp), dimension(:), intent(in) :: arr
        integer, intent(in) :: i, n
        real(dp) :: total_sum
        integer :: k

        total_sum = 0.0_dp
        do k = i, n
            total_sum = total_sum + arr(k)
        end do
    end function pk_sum_real

    !=================================================
    !> @brief Computes the first derivative using forward difference.
    !> @param x Point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param y Real function f(x).
    !> @return Numerical derivative f'(x).
    !=================================================
    function pk_forward_difference_scalar(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp) :: yp
        procedure(f) :: y
        
        yp = (y(x+dx) - y(x))/ (dx)
    end function pk_forward_difference_scalar


    !=================================================
    !> @brief Computes the first derivative using forward difference for vector-valued functions.
    !> @param x Point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param y Vector function f(x).
    !> @return Numerical derivative vector f'(x).
    !=================================================
    function pk_forward_difference_vector(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp), dimension(:), allocatable :: yp
        procedure(f_vector) :: y
        
        yp = (y(x+dx) - y(x))/ (dx)
    end function pk_forward_difference_vector


    !=================================================
    !> @brief Computes the partial derivative using forward difference for multi-variable functions.
    !> @param x Vector of point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param i Index of the variable to differentiate with respect to.
    !> @param y Multi-variable real function f(x).
    !> @return Numerical partial derivative df/dxi.
    !=================================================
    function pk_forward_difference_multi(x, dx, i, y) result(yp)
        real(dp), dimension(:), intent(in) :: x
        real(dp), intent(in) :: dx
        integer, intent(in) :: i
        procedure(f_multi) :: y
        real(dp) :: yp
        real(dp), dimension(size(x)) :: x_plus_dx

        x_plus_dx = x
        x_plus_dx(i) = x_plus_dx(i) + dx
        
        yp = (y(x_plus_dx) - y(x)) / dx
    end function pk_forward_difference_multi

    !=================================================
    !> @brief Computes the first derivative using backward difference.
    !> @param x Point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param y Real function f(x).
    !> @return Numerical derivative f'(x).
    !=================================================
    function pk_backward_difference_scalar(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp) :: yp
        procedure(f) :: y
        
        yp = (y(x) - y(x-dx))/ (dx)
    end function pk_backward_difference_scalar


    !=================================================
    !> @brief Computes the first derivative using backward difference for vector-valued functions.
    !> @param x Point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param y Vector function f(x).
    !> @return Numerical derivative vector f'(x).
    !=================================================
    function pk_backward_difference_vector(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp), dimension(:), allocatable :: yp
        procedure(f_vector) :: y
        
        yp = (y(x) - y(x-dx))/ (dx)
    end function pk_backward_difference_vector


    !=================================================
    !> @brief Computes the partial derivative using backward difference for multi-variable functions.
    !> @param x Vector of point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param i Index of the variable to differentiate with respect to.
    !> @param y Multi-variable real function f(x).
    !> @return Numerical partial derivative df/dxi.
    !=================================================
    function pk_backward_difference_multi(x, dx, i, y) result(yp)
        real(dp), dimension(:), intent(in) :: x
        real(dp), intent(in) :: dx
        integer, intent(in) :: i
        procedure(f_multi) :: y
        real(dp) :: yp
        real(dp), dimension(size(x)) :: x_minus_dx

        x_minus_dx = x
        x_minus_dx(i) = x_minus_dx(i) - dx
        
        yp = (y(x) - y(x_minus_dx)) / dx
    end function pk_backward_difference_multi

    !=================================================
    !> @brief Computes the first derivative using central difference (higher accuracy).
    !> @param x Point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param y Real function f(x).
    !> @return Numerical derivative f'(x).
    !=================================================
    function pk_central_difference_scalar(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp) :: yp
        procedure(f) :: y
        
        yp = (y(x + dx) - y(x - dx))/ (2*dx)
    end function pk_central_difference_scalar

    
    !=================================================
    !> @brief Computes the first derivative using central difference (higher accuracy) for vector-valued functions.
    !> @param x Point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param y Vector function f(x).
    !> @return Numerical derivative vector f'(x).
    !=================================================
    function pk_central_difference_vector(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        procedure(f_vector) :: y
        real(dp), dimension(:), allocatable :: yp
        
        yp = (y(x + dx) - y(x - dx))/ (2.0_dp * dx)
    end function pk_central_difference_vector


    !=================================================
    !> @brief Computes the partial derivative using central difference for multi-variable functions.
    !> @param x Vector of point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param i Index of the variable to differentiate with respect to.
    !> @param y Multi-variable real function f(x).
    !> @return Numerical partial derivative df/dxi.
    !=================================================
    function pk_central_difference_multi(x, dx, i, y) result(yp)
        real(dp), dimension(:), intent(in) :: x
        real(dp), intent(in) :: dx
        integer, intent(in) :: i
        procedure(f_multi) :: y
        real(dp) :: yp
        real(dp), dimension(size(x)) :: x_plus_dx, x_minus_dx

        x_plus_dx = x
        x_plus_dx(i) = x_plus_dx(i) + dx

        x_minus_dx = x
        x_minus_dx(i) = x_minus_dx(i) - dx
        
        yp = (y(x_plus_dx) - y(x_minus_dx)) / (2.0_dp * dx)
    end function pk_central_difference_multi

    !=================================================
    !> @brief Computes the second derivative using central difference.
    !> @param x Point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param y Real function f(x).
    !> @return Numerical second derivative f''(x).
    !=================================================

    function pk_second_central_difference_scalar(x, dx, y) result(ypp)
        real(dp), intent(in) :: x, dx
        real(dp) :: ypp
        procedure(f) :: y
        
        ypp = (y(x+dx) - 2.0_dp*y(x) + y(x-dx))/ (dx**2)
    end function pk_second_central_difference_scalar


    !=================================================
    !> @brief Computes the second derivative using central difference for vector-valued functions.
    !> @param x Point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param y Vector function f(x).
    !> @return Numerical second derivative vector f''(x).
    !=================================================
    function pk_second_central_difference_vector(x, dx, y) result(ypp)
        real(dp), intent(in) :: x, dx
        real(dp), dimension(:), allocatable :: ypp
        procedure(f_vector) :: y
        
        ypp = (y(x+dx) - 2.0_dp*y(x) + y(x-dx))/ (dx**2)
    end function pk_second_central_difference_vector


    !=================================================
    !> @brief Computes the second partial derivative using central difference for multi-variable functions.
    !> @param x Vector of point at which to evaluate the derivative.
    !> @param dx Step size.
    !> @param i Index of the variable to differentiate with respect to.
    !> @param y Multi-variable real function f(x).
    !> @return Numerical second partial derivative d^2f/dxi^2.
    !=================================================
    function pk_second_central_difference_multi(x, dx, i, y) result(ypp)
        real(dp), dimension(:), intent(in) :: x
        real(dp), intent(in) :: dx
        integer, intent(in) :: i
        procedure(f_multi) :: y
        real(dp) :: ypp
        real(dp), dimension(size(x)) :: x_plus_dx, x_minus_dx

        x_plus_dx = x
        x_plus_dx(i) = x_plus_dx(i) + dx

        x_minus_dx = x
        x_minus_dx(i) = x_minus_dx(i) - dx
        
        ypp = (y(x_plus_dx) - 2.0_dp*y(x) + y(x_minus_dx)) / (dx**2)
    end function pk_second_central_difference_multi

    !#################################################
    ! Numerical Integration Methods
    !#################################################

    !=================================================
    !> @brief Simple rectangular rule for numerical integration.
    !> @param x0 Lower integration limit.
    !> @param x1 Upper integration limit.
    !> @param N Number of sub-intervals.
    !> @param y Real function f(x).
    !> @return Numerical integral of f from x0 to x1.
    !=================================================
    function pk_rectangular_rule(x0, x1, N, y) result(Integral)
        real(dp), intent(in) :: x0, x1
        real(dp) :: dx, Integral
        integer, intent(in) :: N
        integer :: i
        procedure(f) :: y
        
        dx = (x1 - x0)/N
        Integral = 0.0_dp

        do i = 0, N - 1
            Integral = Integral + y(x0 + i*dx)*dx
        end do

    end function pk_rectangular_rule

    !=================================================
    !> @brief Trapezoidal rule for numerical integration.
    !> @param x0 Lower integration limit.
    !> @param x1 Upper integration limit.
    !> @param N Number of sub-intervals.
    !> @param y Real function f(x).
    !> @return Numerical integral of f from x0 to x1.
    !=================================================
    function pk_trapezoidal_rule(x0, x1, N, y) result(Integral)
        real(dp), intent(in) :: x0, x1
        real(dp) :: dx, Integral
        integer, intent(in) :: N
        integer :: i
        procedure(f) :: y
        
        dx = (x1 - x0)/N
        Integral = 0.5_dp*(y(x0) + y(x1))

        do i = 1, N - 1
            Integral = Integral + y(x0 + i*dx)
        end do

        Integral = Integral*dx

    end function pk_trapezoidal_rule

    !=================================================
    !> @brief Composite Simpson's 1/3 rule.
    !> @param x0 Lower integration limit.
    !> @param x1 Upper integration limit.
    !> @param N Number of sub-intervals (must be even).
    !> @param y Real function f(x).
    !> @return Numerical integral of f from x0 to x1.
    !=================================================
    function pk_composite_simpson_real(x0, x1, N, y) result(Integral)
        real(dp), intent(in) :: x0, x1
        real(dp) :: dx, Integral
        integer, intent(in) :: N
        integer :: i
        procedure(f) :: y
        
        if (N <= 0) stop "Error: N must be positive"
        if (mod(N, 2) /= 0) stop "Error: N must be even for Simpson's rule"

        dx = (x1 - x0)/N
        Integral = y(x0) + y(x1)

        do i = 1, N - 1
            if (mod(i, 2) == 0) then
                Integral = Integral + 2.0_dp*y(x0 + i*dx)
            else
                Integral = Integral + 4.0_dp*y(x0 + i*dx)
            end if
        end do

        Integral = Integral*dx/3.0_dp

    end function pk_composite_simpson_real

    !=================================================
    !> @brief Composite Simpson's 1/3 rule for complex functions.
    !> @param x0 Lower integration limit.
    !> @param x1 Upper integration limit.
    !> @param N Number of sub-intervals (must be even).
    !> @param y Complex function f(z).
    !> @return Numerical integral of f from x0 to x1.
    !=================================================
    function pk_composite_simpson_complex(x0, x1, N, y) result(Integral)
        complex(dp), intent(in) :: x0, x1
        complex(dp) :: dx, Integral
        integer, intent(in) :: N
        integer :: i
        procedure(f_complex) :: y
        
        if (N <= 0) stop "Error: N must be positive"
        if (mod(N, 2) /= 0) stop "Error: N must be even for Simpson's rule"

        dx = (x1 - x0)/N
        Integral = y(x0) + y(x1)

        do i = 1, N - 1
            if (mod(i, 2) == 0) then
                Integral = Integral + 2.0_dp*y(x0 + i*dx)
            else
                Integral = Integral + 4.0_dp*y(x0 + i*dx)
            end if
        end do

        Integral = Integral*dx/3.0_dp

    end function pk_composite_simpson_complex

    !=================================================
    !> @brief Adaptive Simpson's method for automatic step size control.
    !> @param x0 Lower integration limit.
    !> @param x1 Upper integration limit.
    !> @param tol Numerical tolerance.
    !> @param i Current recursion depth (start with 0).
    !> @param imax Maximum recursion depth.
    !> @param y Real function f(x).
    !> @return Numerical integral result.
    !=================================================
    recursive function pk_adaptative_simpson_real(x0, x1, tol, i, imax, y) result(Integral)
        real(dp), intent(in) :: x0, x1, tol
        real(dp) :: Integral, m, S, S1, S2
        integer, intent(in) :: i, imax
        procedure(f) :: y

        m = (x0 + x1)/2

        S = pk_simpson(x0, x1, y)
        S1 = pk_simpson(x0, m, y)
        S2 = pk_simpson(m, x1, y)

        if (abs(S1 + S2 - S) < 15.0_dp*tol .or. i >= imax) then
            Integral = S1 + S2 + ((S1 + S2 - S)/15.0_dp)
        else
            Integral = pk_adaptative_simpson(x0, m, tol/2, i+1, imax, y) + pk_adaptative_simpson(m, x1, tol/2, i+1, imax, y)
        end if

    contains

        function pk_simpson(x0, x1, y) result(Integral)
            real(dp), intent(in) :: x0, x1
            real(dp) :: Integral
            procedure(f) :: y
            
            Integral = (x1 - x0)/6.0_dp*(y(x0) + 4.0_dp*y((x0 + x1)/2.0_dp) + y(x1))

        end function pk_simpson

    end function pk_adaptative_simpson_real


    !=================================================
    !> @brief Adaptive Simpson's method for complex functions.
    !> @param x0 Lower integration limit.
    !> @param x1 Upper integration limit.
    !> @param tol Numerical tolerance.
    !> @param i Current recursion depth.
    !> @param imax Maximum recursion depth.
    !> @param y Complex function f(z).
    !> @return Numerical integral result.
    !=================================================
    recursive function pk_adaptative_simpson_complex(x0, x1, tol, i, imax, y) result(Integral)
        complex(dp), intent(in) :: x0, x1
        real(dp), intent(in) :: tol
        integer, intent(in) :: i, imax
        complex(dp) :: Integral, m, S, S1, S2
        procedure(f_complex) :: y 

        m = (x0 + x1)/2.0_dp

        S  = pk_simpson_complex(x0, x1, y)
        S1 = pk_simpson_complex(x0, m, y)
        S2 = pk_simpson_complex(m, x1, y)

        if (abs(S1 + S2 - S) < 15.0_dp*tol .or. i >= imax) then
            Integral = S1 + S2 + (S1 + S2 - S)/15.0_dp
        else
            Integral = pk_adaptative_simpson_complex(x0, m, tol/2.0_dp, i+1, imax, y) + &
                    pk_adaptative_simpson_complex(m, x1, tol/2.0_dp, i+1, imax, y)
        end if
    
    contains

        function pk_simpson_complex(x0, x1, y) result(Integral)
            complex(dp), intent(in) :: x0, x1
            complex(dp) :: Integral
            procedure(f_complex) :: y
            
            Integral = (x1 - x0)/6.0_dp*(y(x0) + 4.0_dp*y((x0 + x1)/2.0_dp) + y(x1))

        end function pk_simpson_complex

    end function pk_adaptative_simpson_complex


    !###################################################
    ! Non-linear Equation Solvers (Root-Finding)
    !###################################################

    !=================================================
    !> @brief Bisection method for finding roots of f(x) = 0.
    !> @param a Lower bound of the interval.
    !> @param b Upper bound of the interval.
    !> @param tol Convergence tolerance.
    !> @param imax Maximum number of iterations.
    !> @param y Real function f(x).
    !> @return Estimated root.
    !=================================================
    function pk_bisection_method(a, b, tol, imax, y) result(root)
        real(dp), intent(in) :: a, b, tol
        real(dp) :: a_local, b_local
        integer, intent(in) :: imax
        real(dp) :: root
        integer :: i
        procedure(f) :: y

        if (y(a)*y(b) >= 0.0_dp) then
            print *, "Error: y(a) and y(b) must have opposite signs"
            root = 0.0_dp
            return
        end if

        a_local = a
        b_local = b
        root = (a_local + b_local)/2.0_dp

        do i = 1, imax
            if (abs(y(root)) < tol) then
                return
            else if (y(root)*y(a_local) < 0.0_dp) then
                b_local = root
            else
                a_local = root
            end if
            root = (a_local + b_local)/2.0_dp
        end do

        print *, "Warning: maximum number of iterations reached, result may not be accurate"
    end function pk_bisection_method

    !=================================================
    !> @brief Newton-Raphson method for root-finding.
    !> @details Uses numerical derivative (central difference) if not provided.
    !> @param x0 Initial guess.
    !> @param tol Convergence tolerance.
    !> @param imax Maximum number of iterations.
    !> @param y Real function f(x).
    !> @return Estimated root.
    !=================================================
    function pk_newton_raphson(x0, tol, imax, y) result(root)
        real(dp), intent(in) :: x0, tol
        real(dp) :: root, x, fx, dfx
        integer, intent(in) :: imax
        integer :: i
        procedure(f) :: y

        x = x0

        do i = 1, imax
            fx = y(x)
            dfx = pk_central_difference(x, 1.0e-6_dp, y)
            if (abs(fx) < tol) then
                root = x
                return
            else if (dfx == 0.0_dp) then
                print *, "Error: derivative is zero"
                root = x
                return
            end if
            x = x - fx/dfx
        end do

        root = x
        print *, "Warning: maximum number of iterations reached, result may not be accurate"
    end function pk_newton_raphson

    !=================================================
    !> @brief Secant method for root-finding.
    !> @details Alternative to Newton-Raphson that doesn't require an explicit derivative.
    !> @param x0 First initial guess.
    !> @param x1 Second initial guess.
    !> @param tol Convergence tolerance.
    !> @param imax Maximum number of iterations.
    !> @param y Real function f(x).
    !> @return Estimated root.
    !=================================================
    function pk_secant_method(x0, x1, tol, imax, y) result(root)
    real(dp), intent(in) :: x0, x1, tol
    real(dp) :: root, x_new, x0_local, x1_local, fx0, fx1
    integer, intent(in) :: imax
    integer :: i
    procedure(f) :: y
        
        fx0 = y(x0)
        fx1 = y(x1)
        x0_local = x0
        x1_local = x1

        do i = 1, imax
            if (abs(fx1 - fx0) < 1.0e-14_dp) then
                print*, "Error: division by zero in secant method"
                root = x1_local
                return
            end if
            x_new = x1_local - ((x1_local - x0_local)/(fx1 - fx0)) * fx1
        
            if (abs(x_new - x1_local) < tol) then
                root = x_new
                return
            end if

            x0_local = x1_local
            x1_local = x_new
            fx0 = fx1
            fx1 = y(x1_local)
        end do

        root = x1_local
        print*, "Warning: maximum number of iterations reached, result may not be accurate"
    end function pk_secant_method

end module physkit_numerical