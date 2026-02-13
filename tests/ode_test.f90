program ode_test
    use physkit_constants, only: dp
    use physkit_ode
    implicit none

    real(dp) :: y

    y = 0.0_dp

    call euler(0.1_dp, 0.0_dp, 10.0_dp, 1.0_dp, y, f)
    print *, "Euler result:", y

    call rk2(0.1_dp, 0.0_dp, 10.0_dp, 1.0_dp, y, f)
    print *, "RK2 result:", y

    call rk4(0.1_dp, 0.0_dp, 10.0_dp, 1.0_dp, y, f)
    print *, "RK4 result:", y

contains

    function f(x, y)
        real(dp) :: f
        real(dp), intent(in) :: x, y
        f = exp(x)
    end function f

end program ode_test