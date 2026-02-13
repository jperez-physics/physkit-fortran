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
    ! yinitial: initial value of y at x0
    ! x0: initial value of x
    ! N: number of steps
    ! ynext: output value of y after N steps
    ! yfunc: function that computes dy/dx = f(x, y)
    !=================================================
    subroutine euler(dx, yinitial, x0, N, ynext, yfunc)
        real, intent(in)  :: dx, yinitial, x0
        integer, intent(in) :: N
        real, intent(out) :: ynext
        procedure(f) :: yfunc

        real :: y, x
        integer :: j

        y = yinitial
        x = x0

        do j = 1, N
            y = y + dx * yfunc(x, y)
            x = x + dx
        end do

        ynext = y

    end subroutine euler

end module physkit_ode
