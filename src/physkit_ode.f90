module physkit_ode
    use physkit_constants, only: dp
    implicit none
    private
    public :: pk_euler_method, pk_rk2, pk_rk4, pk_velocity_verlet

    interface
        function f(x, y)
            use physkit_constants, only: dp
            real(dp) :: f
            real(dp), intent(in) :: x, y
        end function f
    end interface

contains

    !=================================================
    ! Euler Method for ODE Integration
    !=================================================
    ! dx: step size
    ! x0: initial x value
    ! x1: final x value
    ! yinitial: initial y value at x0
    ! ynext: output y value at x1
    ! yfunc: function that computes dy/dx = f(x, y)
    !=================================================
    subroutine pk_euler_method(dx, x0, x1, yinitial, ynext, yfunc)
        real(dp), intent(in)  :: dx, yinitial, x0, x1
        real(dp), intent(out) :: ynext
        procedure(f) :: yfunc

        real(dp) :: y, x
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        x = x0

        do j = 1, N
            y = y + dx * yfunc(x, y)
            x = x + dx
        end do

        ynext = y

    end subroutine pk_euler_method

    !=================================================
    ! RK2 Method for ODE Integration
    !=================================================
    ! dx: step size
    ! x0: initial x value
    ! x1: final x value
    ! yinitial: initial y value at x0
    ! ynext: output y value at x1
    ! yfunc: function that computes dy/dx = f(x, y)
    !=================================================
    subroutine pk_rk2(dx, x0, x1, yinitial, ynext, yfunc)
        real(dp), intent(in)  :: dx, yinitial, x0, x1
        real(dp), intent(out) :: ynext
        procedure(f) :: yfunc

        real(dp) :: y, x, k1, k2
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        x = x0

        do j = 1, N
            k1 = yfunc(x, y)
            k2 = yfunc(x + dx/2.0_dp, y + dx*k1/2.0_dp)
            y = y + dx * k2
            x = x + dx
        end do

        ynext = y

    end subroutine pk_rk2

    !=================================================
    ! RK4 Method for ODE Integration
    !=================================================
    ! dx: step size
    ! x0: initial x value
    ! x1: final x value
    ! yinitial: initial y value at x0
    ! ynext: output y value at x1
    ! yfunc: function that computes dy/dx = f(x, y)
    !=================================================
    subroutine pk_rk4(dx, x0, x1, yinitial, ynext, yfunc)
        real(dp), intent(in)  :: dx, yinitial, x0, x1
        real(dp), intent(out) :: ynext
        procedure(f) :: yfunc

        real(dp) :: y, x, k1, k2, k3, k4
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        x = x0

        do j = 1, N
            k1 = yfunc(x, y)
            k2 = yfunc(x + dx/2.0_dp, y + dx*k1/2.0_dp)
            k3 = yfunc(x + dx/2.0_dp, y + dx*k2/2.0_dp)
            k4 = yfunc(x + dx,        y + dx*k3)
            y = y + dx*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
            x = x + dx
        end do

        ynext = y

    end subroutine pk_rk4

    !=================================================
    ! Velocity Verlet Method for 2nd-order ODE
    !=================================================
    ! dx: step size
    ! x0: initial independent variable (time) value
    ! x1: final independent variable (time) value
    ! yinitial: initial position at x0
    ! vinitial: initial velocity at x0
    ! ynext: output position at x1
    ! vnext: output velocity at x1
    ! yfunc: function that computes acceleration a(x, y) = d2y/dx2
    ! NOTE: Here 'yfunc(x,y)' must return the acceleration (second order)
    !=================================================
    subroutine pk_velocity_verlet(dx, x0, x1, yinitial, vinitial, ynext, vnext, yfunc)
        real(dp), intent(in)  :: dx, yinitial, vinitial, x0, x1
        real(dp), intent(out) :: ynext, vnext
        procedure(f) :: yfunc

        real(dp) :: y, v, x, a, anext
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        v = vinitial
        x = x0

        do j = 1, N
            a = yfunc(x, y)
            y = y + v*dx + 0.5_dp*a*dx*dx
            anext = yfunc(x + dx, y)
            v = v + 0.5_dp*(a + anext)*dx
            x = x + dx
        end do

        ynext = y
        vnext = v

    end subroutine pk_velocity_verlet

end module physkit_ode
