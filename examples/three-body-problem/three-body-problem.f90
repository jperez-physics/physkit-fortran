program three_body_problem
    use physkit_ode
    use physkit_constants, only: dp
    implicit none

    real(dp), parameter :: G = 1.0_dp
    real(dp), dimension(12) :: y0, y_final
    real(dp) :: t0, t1, dt
    integer :: i, n_steps
    character(len=100) :: filename = "three-body-problem.dat"

    ! Masses are assumed equal to 1.0
    y0(1:2)   = [ 0.97000436_dp, -0.24308753_dp] ! r1
    y0(3:4)   = [ 0.46620368_dp,  0.43236573_dp] ! v1

    y0(5:6)   = [-0.97000436_dp,  0.24308753_dp] ! r2
    y0(7:8)   = [ 0.46620368_dp,  0.43236573_dp] ! v2

    y0(9:10)  = [ 0.0_dp,         0.0_dp]         ! r3
    y0(11:12) = [-0.93240737_dp, -0.86473146_dp] ! v3

    t0 = 0.0_dp
    t1 = 6.3259_dp ! Approx period
    dt = 0.01_dp
    n_steps = nint((t1 - t0) / dt)

    print *, "Starting Three-Body Problem simulation (Figure-8 orbit)..."
    
    open(10, file=filename, status='replace')
    write(10, '(A)') "# t x1 y1 x2 y2 x3 y3"
    
    y_final = y0
    do i = 0, n_steps
        t1 = t0 + i * dt
        write(10, '(7F12.6)') t1, y_final(1), y_final(2), y_final(5), y_final(6), y_final(9), y_final(10)
        
        call pk_rk4(dt, t1, t1 + dt, y_final, y_final, three_body_derivs)
    end do
    close(10)

    print *, "Simulation finished. Trajectories saved to ", trim(filename)
    call execute_command_line("python plot.py")

contains

    function three_body_derivs(t, y) result(dydt)
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:), allocatable :: dydt
        
        real(dp), dimension(2) :: r1, r2, r3, v1, v2, v3
        real(dp), dimension(2) :: a1, a2, a3
        real(dp) :: r12, r13, r23
        
        allocate(dydt(size(y)))

        r1 = y(1:2)
        v1 = y(3:4)
        r2 = y(5:6)
        v2 = y(7:8)
        r3 = y(9:10)
        v3 = y(11:12)

        r12 = sqrt(sum((r1 - r2)**2))
        r13 = sqrt(sum((r1 - r3)**2))
        r23 = sqrt(sum((r2 - r3)**2))

        ! Accelerations (assuming m1=m2=m3=1)
        a1 = -G * ( (r1 - r2)/r12**3 + (r1 - r3)/r13**3 )
        a2 = -G * ( (r2 - r1)/r12**3 + (r2 - r3)/r23**3 )
        a3 = -G * ( (r3 - r1)/r13**3 + (r3 - r2)/r23**3 )

        dydt(1:2)   = v1
        dydt(3:4)   = a1
        dydt(5:6)   = v2
        dydt(7:8)   = a2
        dydt(9:10)  = v3
        dydt(11:12) = a3

    end function three_body_derivs

end program three_body_problem