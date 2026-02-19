program gamma_graph
    use physkit_constants, only: dp
    use physkit_special, only: pk_gamma_real
    implicit none

    real(dp) :: x, y
    integer :: i, min, max, N

    min = -5
    max = 5
    N = 400

    open(unit=10, file='gamma_graph.dat')
    do i = min * N, max * N
        x = real(i, dp)/real(N, dp)
        y = pk_gamma_real(x)
        
        if (y > 5e5_dp .or. y < -5e5_dp) then
            write(10,*) x, "NaN"
        else
            write(10,*) x, y
        end if
    end do
    close(10)

    call execute_command_line("python plot.py")
    
end program gamma_graph