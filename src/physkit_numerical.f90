module physkit_numerical
    use physkit_constants, only: dp
    implicit none
    private
    public :: pk_bisection_method, pk_forward_difference, pk_backward_difference, pk_central_difference, pk_second_central_difference, pk_rectangular_rule, pk_trapezoidal_rule, pk_composite_simpson, pk_adaptative_simpson

    interface
        function f(x)
            use physkit_constants, only: dp
            real(dp) :: f
            real(dp), intent(in) :: x
        end function f
    end interface

contains

    !##################################################
    !
    ! Numerical differentiation methods
    !
    !##################################################

    !=================================================
    ! Foward difference approximation
    !=================================================
    ! x: point at which to evaluate the derivative
    ! dx: step size
    ! y: function
    ! yp: numerical derivative of f at x
    !=================================================
    function pk_forward_difference(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp) :: yp
        procedure(f) :: y
        
        yp = (y(x+dx) - y(x))/ (dx)
    end function pk_forward_difference

    !=================================================
    ! Backward difference approximation
    !=================================================
    ! x: point at which to evaluate the derivative
    ! dx: step size
    ! y: function
    ! yp: numerical derivative of f at x
    !=================================================
    function pk_backward_difference(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp) :: yp
        procedure(f) :: y
        
        yp = (y(x) - y(x-dx))/ (dx)
    end function pk_backward_difference

    !=================================================
    ! Central difference approximation
    !=================================================
    ! x: point at which to evaluate the derivative
    ! dx: step size
    ! y: function
    ! yp: numerical derivative of f at x
    !=================================================
    function pk_central_difference(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp) :: yp
        procedure(f) :: y
        
        yp = (y(x + dx) - y(x - dx))/ (2*dx)
    end function pk_central_difference

    !=================================================
    ! Central second difference approximation
    !=================================================
    ! x: point at which to evaluate the second derivative
    ! dx: step size
    ! y: function
    ! ypp: numerical second derivative of f at x
    !=================================================
    function pk_second_central_difference(x, dx, y) result(ypp)
        real(dp), intent(in) :: x, dx
        real(dp) :: ypp
        procedure(f) :: y
        
        ypp = (y(x+dx) - 2.0_dp*y(x) + y(x-dx))/ (dx**2)
    end function pk_second_central_difference

    !#################################################
    !
    ! Numerical integration methods
    !
    !################################################

    !=================================================
    ! Rectangular rule
    !=================================================
    ! x0: down limit
    ! x1: upper limit
    ! N: number of subdivisions
    ! dx: step size
    ! y: function
    ! Integral: numerical integral of f from x0 to x1
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
    ! Trapezoidal rule
    !=================================================
    ! x0: down limit
    ! x1: upper limit
    ! N: number of subdivisions
    ! dx: step size
    ! y: function
    ! Integral: numerical integral of f from x0 to x1
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
    ! Simpson's rule
    !=================================================
    ! x0: down limit
    ! x1: upper limit
    ! y: function
    ! Integral: numerical integral of f from x0 to x1
    !=================================================
    function pk_simpson(x0, x1, y) result(Integral)
        real(dp), intent(in) :: x0, x1
        real(dp) :: Integral
        procedure(f) :: y
        
        Integral = (x1 - x0)/6.0_dp*(y(x0) + 4.0_dp*y((x0 + x1)/2.0_dp) + y(x1))

    end function pk_simpson

    !=================================================
    ! Composite Simpson's rule
    !=================================================
    ! x0: down limit
    ! x1: upper limit
    ! N: number of subdivisions (must be even)
    ! dx: step size
    ! y: function
    ! Integral: numerical integral of f from x0 to x1
    !=================================================
    function pk_composite_simpson(x0, x1, N, y) result(Integral)
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

    end function pk_composite_simpson

    !=================================================
    ! Adaptative Simpson's
    !=================================================
    ! x0: down limit
    ! x1: upper limit
    ! tol: tolerance
    ! imax: maximum iterations
    ! y: function
    ! Integral: numerical integral of f from x0 to x1
    !=================================================
    recursive function pk_adaptative_simpson(x0, x1, tol, i, imax, y) result(Integral)
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

    end function pk_adaptative_simpson

    !###################################################
    !
    ! Non-linear equations solvers
    !
    !###################################################

    !=================================================
    ! Bisection method
    !=================================================
    ! a: lower bound
    ! b: upper bound
    ! tol: tolerance
    ! imax: maximum iterations
    ! root: output root of f
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

    end function pk_bisection_method

    !------------------------------------------------
    ! Newton-Raphson method
    !------------------------------------------------
    ! x0: initial guess
    ! tol: tolerance
    ! imax: maximum iterations
    ! root: output root of f
    !------------------------------------------------
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

    end function pk_newton_raphson

end module physkit_numerical