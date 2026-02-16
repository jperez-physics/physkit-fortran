program special_test
    use physkit_special
    use physkit_constants, only: dp, pi
    implicit none

    real(dp) :: result_real, expected_real
    complex(dp) :: z, result_complex, expected_complex
    integer :: n, r, result_int, expected_int
    real(dp) :: tol = 1.0e-5_dp

    print *, "Running functionality tests for physkit_special..."

    !##################################################
    ! Gamma Function Tests
    !##################################################
    print *, "--- Gamma Function Tests ---"

    ! Test 1: pk_gamma_real(1.0) = 1.0
    result_real = pk_gamma_real(1.0_dp)
    expected_real = 1.0_dp
    print *, "pk_gamma_real(1.0): ", result_real, " Expected: ", expected_real
    if (abs(result_real - expected_real) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Test 2: pk_gamma_real(2.0) = 1.0
    result_real = pk_gamma_real(2.0_dp)
    expected_real = 1.0_dp
    print *, "pk_gamma_real(2.0): ", result_real, " Expected: ", expected_real
    if (abs(result_real - expected_real) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Test 3: pk_gamma_real(0.5) = sqrt(pi)
    result_real = pk_gamma_real(0.5_dp)
    expected_real = sqrt(pi)
    print *, "pk_gamma_real(0.5): ", result_real, " Expected: ", expected_real
    if (abs(result_real - expected_real) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Test 4: pk_gamma(1.0 + 0.0i) = 1.0
    z = cmplx(1.0_dp, 0.0_dp, kind=dp)
    result_complex = pk_gamma(z)
    expected_complex = cmplx(1.0_dp, 0.0_dp, kind=dp)
    print *, "pk_gamma((1.0, 0.0)): ", result_complex
    if (abs(result_complex - expected_complex) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    !##################################################
    ! Factorial Tests
    !##################################################
    print *, "--- Factorial Tests ---"

    ! n = 0 -> 1
    n = 0
    result_int = pk_factorial(n)
    expected_int = 1
    print *, "pk_factorial(0): ", result_int, " Expected: ", expected_int
    if (result_int == expected_int) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! n = 5 -> 120
    n = 5
    result_int = pk_factorial(n)
    expected_int = 120
    print *, "pk_factorial(5): ", result_int, " Expected: ", expected_int
    if (result_int == expected_int) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    !##################################################
    ! Permutation Tests
    !##################################################
    print *, "--- Permutation Tests ---"

    ! P(5, 2) = 5! / (5-2)! = 120 / 6 = 20
    n = 5
    r = 2
    result_int = pk_permutation(n, r)
    expected_int = 20
    print *, "pk_permutation(5, 2): ", result_int, " Expected: ", expected_int
    if (result_int == expected_int) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    !##################################################
    ! Combination Tests
    !##################################################
    print *, "--- Combination Tests ---"

    ! C(5, 2) = 5! / (2! * 3!) = 120 / (2 * 6) = 10
    n = 5
    r = 2
    result_int = pk_combination(n, r)
    expected_int = 10
    print *, "pk_combination(5, 2): ", result_int, " Expected: ", expected_int
    if (result_int == expected_int) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    print *, "Tests finished."

end program special_test
