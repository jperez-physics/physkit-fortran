module physkit_numerical
    use physkit_constants, only: dp
    implicit none
    private
    public :: pk_fdaprox, pk_bdaprox, pk_cdaprox, pk_cddaprox

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
    ! h: step size
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
    ! h: step size
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
    ! h: step size
    ! y: function
    ! yp: numerical derivative of f at x
    !=================================================
    function pk_cdaprox(x, dx, y) result(yp)
        real(dp), intent(in) :: x, dx
        real(dp) :: yp
        procedure(f) :: y
        
        yp = (y(x+dx) - y(x - dx))/ (2*dx)
    end function pk_cdaprox

    !=================================================
    ! Central second difference approximation
    !=================================================
    ! x: point at which to evaluate the second derivative
    ! h: step size
    ! y: function
    ! ypp: numerical second derivative of f at x
    !=================================================
    function pk_cddaprox(x, dx, y) result(ypp)
        real(dp), intent(in) :: x, dx
        real(dp) :: ypp
        procedure(f) :: y
        
        ypp = (y(x+dx) - 2.0_dp*y(x) + y(x-dx))/ (dx**2)
    end function pk_cddaprox

end module physkit_numerical