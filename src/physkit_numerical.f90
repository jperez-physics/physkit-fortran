module physkit_numerical
    use physkit_constants, only: dp
    implicit none
    private
    public :: pk_fdaprox, pk_bdaprox, pk_cdaprox, pk_cddaprox, pk_rectangular, pk_trapezoidal, pk_csimpson, pk_asimpson

    interface
        function f(x)
            use physkit_constants, only: dp
            real(dp) :: f
            real(dp), intent(in) :: x
        end function f
    end interface

contains

    !=================================================
    ! Foward difference approximation
    !=================================================
    ! x: point at which to evaluate the derivative
    ! dx: step size
    ! y: function
    ! yp: numerical derivative of f at x
    !=================================================
    function pk_fdaprox(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp) :: yp
        procedure(f) :: y
        
        yp = (y(x+dx) - y(x))/ (dx)
    end function pk_fdaprox

    !=================================================
    ! Backward difference approximation
    !=================================================
    ! x: point at which to evaluate the derivative
    ! dx: step size
    ! y: function
    ! yp: numerical derivative of f at x
    !=================================================
    function pk_bdaprox(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp) :: yp
        procedure(f) :: y
        
        yp = (y(x) - y(x-dx))/ (dx)
    end function pk_bdaprox

    !=================================================
    ! Central difference approximation
    !=================================================
    ! x: point at which to evaluate the derivative
    ! dx: step size
    ! y: function
    ! yp: numerical derivative of f at x
    !=================================================
    function pk_cdaprox(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp) :: yp
        procedure(f) :: y
        
        yp = (y(x + dx) - y(x - dx))/ (2*dx)
    end function pk_cdaprox

    !=================================================
    ! Central second difference approximation
    !=================================================
    ! x: point at which to evaluate the second derivative
    ! dx: step size
    ! y: function
    ! ypp: numerical second derivative of f at x
    !=================================================
    function pk_cddaprox(x, dx, y) result(ypp)
        real(dp), intent(in) :: x, dx
        real(dp) :: ypp
        procedure(f) :: y
        
        ypp = (y(x+dx) - 2.0_dp*y(x) + y(x-dx))/ (dx**2)
    end function pk_cddaprox

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
    function pk_rectangular(x0, x1, N, y) result(Integral)
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

    end function pk_rectangular

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
    function pk_trapezoidal(x0, x1, N, y) result(Integral)
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

    end function pk_trapezoidal

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
    function pk_csimpson(x0, x1, N, y) result(Integral)
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

    end function pk_csimpson

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
    recursive function pk_asimpson(x0, x1, tol, i, imax, y) result(Integral)
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
        Integral = pk_asimpson(x0, m, tol/2, i+1, imax, y) + pk_asimpson(m, x1, tol/2, i+1, imax, y)
    end if

    end function pk_asimpson

end module physkit_numerical