program gamma_graph
    use physkit_constants, only: dp
    use physkit_special, only: pk_gamma
    implicit none

    complex(dp) :: z, zaxis
    real(dp) :: xaxis, yaxis
    integer :: i, j, min, max, N

    min = -2
    max = 2
    N = 10

    open(unit=10, file='gamma_graph.dat')
    do i = N*min, N*max
        xaxis = real(i, dp) / real(N, dp)

        do j = N*min, N*max
            yaxis = real(j, dp) / real(N, dp)

            z = cmplx(xaxis, yaxis, kind=dp)

            zaxis = pk_gamma(z)

            if (abs(real(zaxis)) > 5e5_dp) then
                write(10,*) xaxis, yaxis, 'NaN'
            else
                write(10,*) xaxis, yaxis, real(zaxis)
            end if

        end do
    end do

    close(10)

    call execute_command_line("python plot.py")
    
end program gamma_graph