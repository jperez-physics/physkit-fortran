program gamma_function
    use physkit_constants, only: dp
    use physkit_special, only: pk_gamma
    implicit none

    integer :: i
    real(dp) :: x, y
    real(dp), parameter :: x_start = -4.5_dp
    real(dp), parameter :: x_end = 4.5_dp
    integer, parameter :: n_points = 2000
    real(dp) :: dx

    dx = (x_end - x_start) / real(n_points - 1, dp)

    open(unit=10, file='gamma_function.dat', status='replace')
    write(10, '(A15, A20)') '# x', 'Gamma(x)'

    do i = 0, n_points - 1
        x = x_start + real(i, dp) * dx
        
        y = pk_gamma(x)
        
        if (abs(y) > 20.0_dp) then
            write(10, '(F15.6, A20)') x, 'NaN'
        else
            write(10, '(F15.6, F20.6)') x, y
        end if
    end do

    close(10)

    print *, "Gamma function data saved to gamma_function.dat"

end program gamma_function