program ode_test
    use physkit_ode
    use physkit_constants, only: dp, pi
    implicit none

    real(dp) :: x0, x1, dx, y_init, y_out, expected, v_init, v_out
    real(dp), dimension(2) :: y_init_sol, y_out_sol, expected_sol, v_init_sol, v_out_sol
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

    ! Scalar Tests
    print *, ""
    print *, "--- SCALAR TESTS ---"
    ! 1. Euler Method
    call pk_euler_method(dx, x0, x1, y_init, y_out, f_exp)
    expected = exp(x1)
    call print_result("pk_euler_method", dx, y_out, expected, tol_euler)

    ! 2. RK2 Method
    call pk_rk2(dx, x0, x1, y_init, y_out, f_exp)
    call print_result("pk_rk2", dx, y_out, expected, tol_rk2)

    ! 3. RK4 Method
    call pk_rk4(dx, x0, x1, y_init, y_out, f_exp)
    call print_result("pk_rk4", dx, y_out, expected, tol_rk4)

    ! System (Vector) Tests
    print *, ""
    print *, "--- SYSTEM TESTS (Harmonic Oscillator) ---"
    x0 = 0.0_dp
    x1 = 1.0_dp
    y_init_sol = [1.0_dp, 0.0_dp]
    expected_sol = [cos(x1), -sin(x1)]

    ! 1. Euler System
    call pk_euler_method(dx, x0, x1, y_init_sol, y_out_sol, f_system_harm)
    call print_result_system("pk_euler_method system", dx, y_out_sol, expected_sol, tol_euler)

    ! 2. RK2 System
    call pk_rk2(dx, x0, x1, y_init_sol, y_out_sol, f_system_harm)
    call print_result_system("pk_rk2 system", dx, y_out_sol, expected_sol, tol_rk2)

    ! 3. RK4 System
    call pk_rk4(dx, x0, x1, y_init_sol, y_out_sol, f_system_harm)
    call print_result_system("pk_rk4 system", dx, y_out_sol, expected_sol, tol_rk4)

    ! 4. Velocity Verlet (Special case for 2nd order)
    print *, ""
    print *, "--- VELOCITY VERLET TEST ---"
    x0 = 0.0_dp
    x1 = 1.0_dp
    y_init = 1.0_dp
    v_init = 0.0_dp
    
    call pk_velocity_verlet(dx, x0, x1, y_init, v_init, y_out, v_out, f_harm)
    print *, "Scalar Verlet:"
    call print_result("  Position", dx, y_out, cos(x1), tol_verlet)
    call print_result("  Velocity", dx, v_out, -sin(x1), tol_verlet)

    ! Velocity Verlet System
    y_init_sol = [1.0_dp, 0.5_dp]
    v_init_sol = [0.0_dp, -0.1_dp]
    call pk_velocity_verlet(dx, x0, x1, y_init_sol, v_init_sol, y_out_sol, v_out_sol, f_system_acc_harm)
    print *, "System Verlet:"
    ! Analytical solution for y'' = -y: y(t) = y(0)*cos(t) + v(0)*sin(t)
    expected_sol = [y_init_sol(1)*cos(x1) + v_init_sol(1)*sin(x1), &
                    y_init_sol(2)*cos(x1) + v_init_sol(2)*sin(x1)]
    call print_result_system("  Position", dx, y_out_sol, expected_sol, tol_verlet)
    
    ! v(t) = -y(0)*sin(t) + v(0)*cos(t)
    expected_sol = [-y_init_sol(1)*sin(x1) + v_init_sol(1)*cos(x1), &
                    -y_init_sol(2)*sin(x1) + v_init_sol(2)*cos(x1)]
    call print_result_system("  Velocity", dx, v_out_sol, expected_sol, tol_verlet)

    print *, ""
    print *, "Tests finished."

contains

    subroutine print_result(label, dx, res, exp_val, tol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: dx, res, exp_val, tol
        real(dp) :: rel_err
        character(len=6) :: status
        
        rel_err = abs((res - exp_val) / exp_val) * 100.0_dp
        if (abs(res - exp_val) < tol) then
            status = "PASSED"
        else
            status = "FAILED"
        end if
        
        write(*, "(A25, ' - ', F6.4, ' - ', A6, ' - ', F12.8, ' - ', F12.8, ' - ', F10.6, ' %')") &
            label, dx, status, res, exp_val, rel_err
    end subroutine print_result

    subroutine print_result_system(label, dx, res, exp_val, tol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: dx, tol
        real(dp), dimension(:), intent(in) :: res, exp_val
        real(dp) :: rel_err
        character(len=6) :: status
        
        ! Mean relative error for system
        rel_err = sum(abs((res - exp_val) / exp_val)) / size(res) * 100.0_dp
        if (all(abs(res - exp_val) < tol)) then
            status = "PASSED"
        else
            status = "FAILED"
        end if

        write(*, "(A25, ' - ', F6.4, ' - ', A6, ' - [', F12.8, ',', F12.8, '] - [', F12.8, ',', F12.8, '] - ', F10.6, ' %')") &
            label, dx, status, res(1), res(2), exp_val(1), exp_val(2), rel_err
    end subroutine print_result_system

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

    function f_system_harm(x, y) result(dydx_system)
        real(dp), intent(in) :: x
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:), allocatable :: dydx_system
        
        allocate(dydx_system(size(y)))
        dydx_system(1) = y(2)
        dydx_system(2) = -y(1)
    end function f_system_harm

    function f_system_acc_harm(x, y) result(d2ydx2_system)
        real(dp), intent(in) :: x
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:), allocatable :: d2ydx2_system
        
        allocate(d2ydx2_system(size(y)))
        d2ydx2_system(1) = -y(1)
        d2ydx2_system(2) = -y(2)
    end function f_system_acc_harm

end program ode_test