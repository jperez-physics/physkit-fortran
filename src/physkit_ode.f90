module physkit_ode
    implicit none
    private
    public :: euler, rk2, rk4, verlet

    interface
        function f(x, y)
            real :: f
            real, intent(in) :: x, y
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
    subroutine euler(dx, x0, x1, yinitial, ynext, yfunc)
        real, intent(in)  :: dx, yinitial, x0, x1
        real, intent(out) :: ynext
        procedure(f) :: yfunc

        real :: y, x
        integer :: j, N

        if (dx == 0.0) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        x = x0

        do j = 1, N
            y = y + dx * yfunc(x, y)
            x = x + dx
        end do

        ynext = y

    end subroutine euler

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
    subroutine rk2(dx, x0, x1, yinitial, ynext, yfunc)
        real, intent(in)  :: dx, yinitial, x0, x1
        real, intent(out) :: ynext
        procedure(f) :: yfunc

        real :: y, x, k1, k2
        integer :: j, N

        if (dx == 0.0) stop "Error: dx cannot ser cero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        x = x0

        do j = 1, N
            k1 = yfunc(x, y)
            k2 = yfunc(x + dx/2.0, y + dx*k1/2.0)
            y = y + dx * k2
            x = x + dx
        end do

        ynext = y

    end subroutine rk2

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
    subroutine rk4(dx, x0, x1, yinitial, ynext, yfunc)
        real, intent(in)  :: dx, yinitial, x0, x1
        real, intent(out) :: ynext
        procedure(f) :: yfunc

        real :: y, x, k1, k2, k3, k4
        integer :: j, N

        if (dx == 0.0) stop "Error: dx cannot ser cero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        x = x0

        do j = 1, N
            k1 = yfunc(x, y)
            k2 = yfunc(x + dx/2.0, y + dx*k1/2.0)
            k3 = yfunc(x + dx/2.0, y + dx*k2/2.0)
            k4 = yfunc(x + dx,     y + dx*k3)
            y = y + dx*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
            x = x + dx
        end do

        ynext = y

    end subroutine rk4

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
    ! NOTE: Here 'yfunc(x,y)' must return la acceleration (second order)
    !=================================================
    subroutine verlet(dx, x0, x1, yinitial, vinitial, ynext, vnext, yfunc)
        real, intent(in)  :: dx, yinitial, vinitial, x0, x1
        real, intent(out) :: ynext, vnext
        procedure(f) :: yfunc

        real :: y, v, x, a, anext
        integer :: j, N

        if (dx == 0.0) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        v = vinitial
        x = x0

        do j = 1, N
            a = yfunc(x, y)
            y = y + v*dx + 0.5*a*dx*dx
            anext = yfunc(x + dx, y)
            v = v + 0.5*(a + anext)*dx
            x = x + dx
        end do

        ynext = y
        vnext = v

    end subroutine verlet

end module physkit_ode
