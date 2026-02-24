program damped_oscillator
    use physkit_constants, only: dp
    use physkit_ode
    implicit none

    real(dp) :: t, dt, t_end
    real(dp), dimension(2) :: y
    real(dp), dimension(2) :: y_out
    real(dp) :: b, k, m, acc
    integer :: i, n_steps

    ! Physical Parameters
    m = 1.0_dp
    k = 10.0_dp
    b = 0.5_dp
    
    ! Simulation Parameters
    dt = 0.01_dp
    t_end = 10.0_dp
    n_steps = nint(t_end / dt)

    ! Initial conditions [position, velocity]
    y = [1.0_dp, 0.0_dp]
    t = 0.0_dp

    ! Open file for recording data
    open(10, file='damped_oscillator.dat', status='replace')
    write(10, '(A12, A15, A15, A15)') "# Time", "Position", "Velocity", "Acceleration"

    ! Initial state record
    acc = -(k*y(1) + b*y(2))/m
    write(10, '(4F15.6)') t, y(1), y(2), acc

    print *, "Starting simulation of a damped harmonic oscillator..."
    print *, "Saving data to 'damped_oscillator.dat'..."

    do i = 1, n_steps
        ! Integrate one step (dt) at a time using RK4 system solver
        call pk_rk4(dt, t, t + dt, y, y_out, oscillator_system)
        
        ! Update state
        y = y_out
        t = t + dt
        
        ! Calculate current acceleration
        acc = -(k*y(1) + b*y(2))/m
        
        ! Write to file: Time, Position, Velocity, Acceleration
        write(10, '(4F15.6)') t, y(1), y(2), acc
    end do

    close(10)
    print *, "Simulation finished successfully."

contains

    ! y(1) is position, y(2) is velocity
    function oscillator_system(t, y) result(dy)
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:), allocatable :: dy
        
        allocate(dy(size(y)))
        
        ! dy(1) = v = y(2)
        ! dy(2) = a = -k/m*x - b/m*v
        dy(1) = y(2)
        dy(2) = -(k*y(1) + b*y(2))/m
        
    end function oscillator_system

end program damped_oscillator