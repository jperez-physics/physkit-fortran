program test
    use physkit_ode
    implicit none

    real :: y

    y = 0.0

    call euler(0.01, 1.0, 0.0, 100, y, f)

    print *, y

contains

    function f(x, y)
        real :: f
        real, intent(in) :: x, y
        f = exp(x)
    end function f

end program test
