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
    ! Euler Method
    !=================================================
    ! dt: time step
    ! yinitial: initial value of y
    ! t0: initial time
    ! N: number of steps
    ! ynext: output variable for the next value of y
    ! yfunc: function that computes the derivative dy/dt = f(x, y)
    !=================================================
    subroutine euler(dt, yinitial, t0, N, yfunc, ynext)
        real, intent(in)  :: dt, yinitial, t0
        integer, intent(in) :: N
        real, intent(out) :: ynext
        procedure(f) :: yfunc

        real :: y, x
        integer :: j

        y = yinitial
        t = t0

        do j = 1, N
            y = y + dt * yfunc(x, y)
            x = x + dt
        end do

        ynext = y

    end subroutine euler

end module physkit_ode
