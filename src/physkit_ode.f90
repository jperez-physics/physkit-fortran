module physkit_ode
    implicit none
    private
    public :: euler

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

end module physkit_ode
