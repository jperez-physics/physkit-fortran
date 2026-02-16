program ode_test
    use physkit_ode
    use physkit_constants, only: dp, pi
    implicit none

    real(dp) :: x0, x1, dx, y_init, y_out, expected, v_init, v_out
    real(dp) :: tol_euler, tol_rk2, tol_rk4, tol_verlet

    ! Tolerances
    tol_euler = 5.0e-3_dp
    tol_rk2 = 1.0e-6_dp
    tol_rk4 = 1.0e-12_dp
    tol_verlet = 1.0e-4_dp

    x0 = 0.0_dp
    x1 = 1.0_dp
    dx = 0.001_dp
    y_init = 1.0_dp
    
    print *, "Running functionality tests for physkit_ode..."

    ! 1. Euler Method
    ! dy/dx = y, y(0)=1 => y(1) = e^1
    call pk_euler_method(dx, x0, x1, y_init, y_out, f_exp)
    expected = exp(x1)
    print *, "pk_euler_method (y'=y, x=[0,1], dx=0.01): ", y_out, " Expected: ", expected
    if (abs(y_out - expected) < tol_euler) then
        print *, "  PASSED"
    else
        print *, "  FAILED: Error = ", abs(y_out - expected)
    end if

    ! 2. RK2 Method
    call pk_rk2(dx, x0, x1, y_init, y_out, f_exp)
    print *, "pk_rk2 (y'=y, x=[0,1], dx=0.01): ", y_out, " Expected: ", expected
    if (abs(y_out - expected) < tol_rk2) then
        print *, "  PASSED"
    else
        print *, "  FAILED: Error = ", abs(y_out - expected)
    end if

    ! 3. RK4 Method
    call pk_rk4(dx, x0, x1, y_init, y_out, f_exp)
    print *, "pk_rk4 (y'=y, x=[0,1], dx=0.01): ", y_out, " Expected: ", expected
    if (abs(y_out - expected) < tol_rk4) then
        print *, "  PASSED"
    else
        print *, "  FAILED: Error = ", abs(y_out - expected)
    end if

    ! 4. Velocity Verlet
    ! d2y/dt2 = -y (Harmonic oscillator)
    ! y(0) = 1, v(0) = 0
    ! y(t) = cos(t), v(t) = -sin(t)
    x0 = 0.0_dp
    x1 = pi ! Half period
    y_init = 1.0_dp
    v_init = 0.0_dp
    
    call pk_velocity_verlet(dx, x0, x1, y_init, v_init, y_out, v_out, f_harm)
    expected = cos(x1) ! -1.0
    print *, "pk_velocity_verlet (y''=-y, x=[0,pi], dx=0.01) Position: ", y_out, " Expected: ", expected
    if (abs(y_out - expected) < tol_verlet) then
        print *, "  PASSED"
    else
        print *, "  FAILED: Error = ", abs(y_out - expected)
    end if
    
    expected = -sin(x1) ! 0.0
    print *, "pk_velocity_verlet (y''=-y, x=[0,pi], dx=0.01) Velocity: ", v_out, " Expected: ", expected
    if (abs(v_out - expected) < tol_verlet) then
        print *, "  PASSED"
    else
        print *, "  FAILED: Error = ", abs(v_out - expected)
    end if

    print *, "Tests finished."

contains

    function f_exp(x, y) result(dydx)
        real(dp), intent(in) :: x, y
        real(dp) :: dydx
        dydx = y
    end function f_exp

    function f_harm(x, y) result(d2ydx2)
        real(dp), intent(in) :: x, y
        real(dp) :: d2ydx2
        d2ydx2 = -y
    end function f_harm

end program ode_test