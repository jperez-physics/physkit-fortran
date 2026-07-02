!> @author Jaime (and contributors)
!> @brief Module for special mathematical functions.
!> @details Implementation of common special functions like Factorial, Gamma (real and complex),
!>          Beta, Permutations, and Combinations.
module physkit_special
    use physkit_numerical
    use physkit_ode
    use physkit_linalg
    use physkit_constants, only: dp, pi, pk_success, pk_err_invalid_argument
    implicit none
    private
    public :: pk_factorial, pk_gamma, pk_beta, pk_permutation, pk_combination

    interface pk_gamma
        module procedure pk_gamma_real
        module procedure pk_gamma_complex
    end interface

contains

    !=================================================
    !> @brief Computes the factorial of a non-negative integer n!
    !> @param n Non-negative integer.
    !> @param[out] ierr Optional error code (pk_success or pk_err_invalid_argument).
    !> @return Factorial of n (-1 on error).
    !=================================================
    function pk_factorial(n, ierr) result(factorial)
    integer, intent(in) :: n
    integer, intent(out), optional :: ierr
    integer :: factorial, i

    if (present(ierr)) ierr = pk_success

    if (n < 0) then
        if (present(ierr)) then
            ierr = pk_err_invalid_argument
            factorial = -1
            return
        else
            error stop "physkit_special: pk_factorial: n must be positive or zero"
        end if
    else if (n == 0) then
        factorial = 1
    else
        factorial = 1
        do i = 1, n
        factorial = factorial * i
        end do
    end if

    end function

    !=================================================
    !> @brief Gamma function for real arguments using Lanczos approximation.
    !> @param z Real input value.
    !> @return Gamma(z).
    !=================================================
    recursive function pk_gamma_real(z) result(gamma)
        real(dp), intent(in) :: z
        real(dp) :: gamma, A, t, zz
        integer :: k

        real(dp), parameter :: g = 7.0_dp
        real(dp), dimension(0:8), parameter :: c = (/ &
                                                    0.99999999999980993_dp, &
                                                    676.5203681218851_dp, &
                                                    -1259.1392167224028_dp, &
                                                    771.32342877765313_dp, &
                                                    -176.61502916214059_dp, &
                                                    12.507343278686905_dp, &
                                                    -0.13857109526572012_dp, &
                                                    9.9843695780195716e-6_dp, &
                                                    1.5056327351493116e-7_dp /)

        if (z < 0.5_dp) then
            gamma = pi / (sin(pi*z) * pk_gamma_real(1.0_dp - z))
            return
        end if

        zz = z - 1.0_dp

        A = c(0)
        do k = 1, 8
            A = A + c(k)/(zz + real(k, dp))
        end do

        t = zz + g + 0.5_dp

        gamma = sqrt(2.0_dp*pi) * t**(zz + 0.5_dp) * exp(-t) * A

    end function pk_gamma_real

    !=================================================
    !> @brief Gamma function for complex arguments using Lanczos approximation.
    !> @param z Complex input value.
    !> @return Gamma(z).
    !=================================================
    recursive function pk_gamma_complex(z) result(gamma)
        complex(dp), intent(in) :: z
        complex(dp) :: zz, gamma, A, t
        integer :: k
        real(dp), parameter :: g = 7.0_dp
        real(dp), dimension(0:8), parameter :: c = (/ &
                                                    0.99999999999980993_dp, &
                                                    676.5203681218851_dp, &
                                                    -1259.1392167224028_dp, &
                                                    771.32342877765313_dp, &
                                                    -176.61502916214059_dp, &
                                                    12.507343278686905_dp, &
                                                    -0.13857109526572012_dp, &
                                                    9.9843695780195716e-6_dp, &
                                                    1.5056327351493116e-7_dp /)

        if (real(z) < 0.5_dp) then
            gamma = pi / (sin(pi*z) * pk_gamma(1.0_dp - z))
            return
        end if

        zz = z - 1.0_dp

        A = c(0)
        do k = 1, 8
            A = A + c(k)/(zz + real(k, dp))
        end do

        t = zz + g + 0.5_dp

        gamma = sqrt(2.0_dp*pi) * t**(zz + 0.5_dp) * exp(-t) * A

    end function pk_gamma_complex

    !=================================================
    !> @brief Beta function B(x, y).
    !> @param x First real input value.
    !> @param y Second real input value.
    !> @return Value of the Beta function.
    !=================================================
    function pk_beta(x, y) result(beta)
        real(dp), intent(in) :: x, y
        real(dp) :: beta

        beta = pk_gamma_real(x) * pk_gamma_real(y) / pk_gamma_real(x + y)

    end function pk_beta

    !=================================================
    !> @brief Calculates the number of permutations P(n, r).
    !> @param n Total number of items.
    !> @param r Number of items to choose.
    !> @param[out] ierr Optional error code (pk_success or pk_err_invalid_argument).
    !> @return Number of permutations (-1 on invalid input).
    !=================================================
    function pk_permutation(n, r, ierr) result(permutation)
        integer, intent(in) :: n, r
        integer, intent(out), optional :: ierr
        integer :: permutation

        if (present(ierr)) ierr = pk_success

        if (n < 0 .or. r < 0 .or. r > n) then
            if (present(ierr)) then
                ierr = pk_err_invalid_argument
                permutation = -1
                return
            else
                error stop "physkit_special: pk_permutation: n and r must be non-negative and r <= n"
            end if
        end if

        permutation = pk_factorial(n) / pk_factorial(n - r)

    end function pk_permutation

    !=================================================
    !> @brief Calculates the number of combinations C(n, r).
    !> @param n Total number of items.
    !> @param r Number of items to choose.
    !> @param[out] ierr Optional error code (pk_success or pk_err_invalid_argument).
    !> @return Number of combinations (-1 on invalid input).
    !=================================================
    function pk_combination(n, r, ierr) result(combination)
        integer, intent(in) :: n, r
        integer, intent(out), optional :: ierr
        integer :: combination

        if (present(ierr)) ierr = pk_success

        if (n < 0 .or. r < 0 .or. r > n) then
            if (present(ierr)) then
                ierr = pk_err_invalid_argument
                combination = -1
                return
            else
                error stop "physkit_special: pk_combination: n and r must be non-negative and r <= n"
            end if
        end if

        combination = pk_factorial(n) / (pk_factorial(r) * pk_factorial(n - r))

    end function pk_combination

    ! Bessel function


end module physkit_special


! Obsolete functions

    !=================================================
    ! Gamma function
    !=================================================
    ! z: input value
    ! gamma: output value
    !=================================================
    !recursive function pk_gamma_real(z) result(gamma)
    !    real(dp), intent(in) :: z
    !    real(dp) :: gamma
!
    !    if (z == 0.0_dp) then
    !        gamma = 1.0_dp
    !    else if (z < 1.0_dp) then
    !        gamma = pk_gamma_real(z + 1.0_dp) / z
    !    else
    !        ! Integrate from epsilon to 1-epsilon to avoid singularities at endpoints
    !        gamma = pk_adaptive_simpson(1.0e-12_dp, 1.0_dp - 1.0e-12_dp, 1.0e-8_dp, 0, 20, integrating)
    !    end if
!
    !contains
!
    !    function integrating(t) result(f)
    !        real(dp), intent(in) :: t
    !        real(dp) :: f
    !        
    !        f = ((t**(z - 1.0_dp)) * exp(-t/(1.0_dp - t))) / ((1.0_dp - t)**(z + 1.0_dp))
    !    end function integrating
!
    !end function pk_gamma_real