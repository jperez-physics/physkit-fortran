program ode_example
    use physkit_ode
    implicit none

    real :: y

    y = 0.0

    call euler(0.1, 0.0, 10.0, 1.0, y, f)
    print *, "Euler result:", y

    call rk2(0.1, 0.0, 10.0, 1.0, y, f)
    print *, "RK2 result:", y

    call rk4(0.1, 0.0, 10.0, 1.0, y, f)
    print *, "RK4 result:", y

contains

    function f(x, y)
        real :: f
        real, intent(in) :: x, y
        f = exp(x)
    end function f

end program ode_example
