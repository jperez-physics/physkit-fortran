program numerical_test
    use physkit_numerical
    use physkit_constants, only: dp, pi
    implicit none

    real(dp) :: x, dx, result_real, expected_real, tol
    real(dp) :: a, b
    integer :: N
    complex(dp) :: z0, z1, result_complex, expected_complex
    
    tol = 1.0e-5_dp
    dx = 1.0e-5_dp

    print *, "Running functionality tests for physkit_numerical..."

    !##################################################
    ! Differentiation Tests
    !##################################################
    print *, "--- Differentiation Tests ---"
    
    ! Test f(x) = x^2, f'(x) = 2x at x = 2.0 -> 4.0
    x = 2.0_dp
    expected_real = 4.0_dp

    ! Forward Difference
    result_real = pk_forward_difference(x, dx, f_poly)
    print *, "pk_forward_difference (x^2, x=2): ", result_real, " Expected: ", expected_real
    if (abs(result_real - expected_real) < 1.0e-4_dp) then ! Lower precision for forward diff
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Backward Difference
    result_real = pk_backward_difference(x, dx, f_poly)
    print *, "pk_backward_difference (x^2, x=2): ", result_real, " Expected: ", expected_real
    if (abs(result_real - expected_real) < 1.0e-4_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Central Difference
    result_real = pk_central_difference(x, dx, f_poly)
    print *, "pk_central_difference (x^2, x=2): ", result_real, " Expected: ", expected_real
    if (abs(result_real - expected_real) < 1.0e-8_dp) then ! Higher precision for central
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Second Central Difference
    ! f(x) = x^2, f''(x) = 2.0
    expected_real = 2.0_dp
    result_real = pk_second_central_difference(x, dx, f_poly)
    print *, "pk_second_central_difference (x^2, x=2): ", result_real, " Expected: ", expected_real
    if (abs(result_real - expected_real) < 1.0e-6_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    !##################################################
    ! Integration Tests
    !##################################################
    print *, "--- Integration Tests ---"
    
    ! Integrate f(x) = x^2 from 0 to 1 -> 1/3 = 0.333333...
    a = 0.0_dp
    b = 1.0_dp
    expected_real = 1.0_dp / 3.0_dp
    N = 100

    ! Rectangular Rule
    result_real = pk_rectangular_rule(a, b, N, f_poly)
    print *, "pk_rectangular_rule (x^2, 0-1, N=100): ", result_real
    if (abs(result_real - expected_real) < 1.0e-2_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Trapezoidal Rule
    result_real = pk_trapezoidal_rule(a, b, N, f_poly)
    print *, "pk_trapezoidal_rule (x^2, 0-1, N=100): ", result_real
    if (abs(result_real - expected_real) < 1.0e-4_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Simpson's Rule (Single interval)
    result_real = pk_simpson(a, b, f_poly)
    print *, "pk_simpson (x^2, 0-1): ", result_real
    if (abs(result_real - expected_real) < 1.0e-8_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Composite Simpson's Rule
    result_real = pk_composite_simpson(a, b, N, f_poly)
    print *, "pk_composite_simpson (x^2, 0-1, N=100): ", result_real
    if (abs(result_real - expected_real) < 1.0e-8_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Adaptive Simpson's Rule
    result_real = pk_adaptative_simpson(a, b, 1.0e-6_dp, 0, 20, f_poly)
    print *, "pk_adaptative_simpson (x^2, 0-1): ", result_real
    if (abs(result_real - expected_real) < 1.0e-6_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Complex Simpson's Rule
    ! Integrate z from 0 to 1+i -> z^2/2 | 0 to 1+i = (1+i)^2/2 = (1 + 2i - 1)/2 = i
    ! Using f_linear_complex(z) = z
    z0 = (0.0_dp, 0.0_dp)
    z1 = (1.0_dp, 1.0_dp)
    expected_complex = (0.0_dp, 1.0_dp) 
    
    result_complex = pk_simpson_complex(z0, z1, f_linear_complex)
    print *, "pk_simpson_complex (z, 0 to 1+i): ", result_complex
    if (abs(result_complex - expected_complex) < 1.0e-6_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED: Expected ", expected_complex, " Got ", result_complex
    end if

    ! Adaptative Simpson's Complex
    result_complex = pk_adaptative_simpson_complex(z0, z1, 1.0e-6_dp, 0, 20, f_linear_complex)
    print *, "pk_adaptative_simpson_complex (z, 0 to 1+i): ", result_complex
    if (abs(result_complex - expected_complex) < 1.0e-6_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED: Expected ", expected_complex, " Got ", result_complex
    end if


    !##################################################
    ! Root Finding Tests
    !##################################################
    print *, "--- Root Finding Tests ---"
    
    ! Find root of x^2 - 4 = 0 in [0, 3] -> x = 2
    a = 0.0_dp
    b = 3.0_dp
    expected_real = 2.0_dp
    
    ! Bisection
    result_real = pk_bisection_method(a, b, 1.0e-6_dp, 100, f_root)
    print *, "pk_bisection_method (x^2-4, [0,3]): ", result_real
    if (abs(result_real - expected_real) < 1.0e-5_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Newton-Raphson
    result_real = pk_newton_raphson(1.5_dp, 1.0e-6_dp, 100, f_root)
    print *, "pk_newton_raphson (x^2-4, x0=1.5): ", result_real
    if (abs(result_real - expected_real) < 1.0e-5_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Secant Method
    result_real = pk_secant_method(1.0_dp, 3.0_dp, 1.0e-6_dp, 100, f_root)
    print *, "pk_secant_method (x^2-4, x0=1, x1=3): ", result_real
    if (abs(result_real - expected_real) < 1.0e-5_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    print *, "Tests finished."

contains

    !-------------------------------------------------------------------------
    ! Real Functions
    !-------------------------------------------------------------------------
    function f_poly(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = x**2
    end function f_poly

    function f_root(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = x**2 - 4.0_dp
    end function f_root

    !-------------------------------------------------------------------------
    ! Complex Functions
    !-------------------------------------------------------------------------
    function f_linear_complex(z) result(y)
        complex(dp), intent(in) :: z
        complex(dp) :: y
        y = z
    end function f_linear_complex

end program numerical_test