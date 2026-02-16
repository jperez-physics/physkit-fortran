module physkit_special
    use physkit_numerical
    use physkit_ode
    use physkit_linalg
    use physkit_constants, only: dp, pi
    implicit none
    private
    public :: pk_factorial, pk_gamma_real, pk_gamma, pk_permutation, pk_combination

contains

    !=================================================
    ! Factorial
    !=================================================
    ! n: input integer
    ! factorial: output factorial of n
    !=================================================
    function pk_factorial(n) result(factorial)
    integer, intent(in) :: n
    integer :: factorial, i

    if (n < 0) then
        print*, 'Error: n must be postive or zero'
        factorial = -1
        return
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
    ! Gamma function
    !=================================================
    ! z: input value
    ! gamma: output value
    !=================================================
    recursive function pk_gamma_real(z) result(gamma)
        real(dp), intent(in) :: z
        real(dp) :: gamma

        if (z < 1.0_dp) then
            gamma = pk_gamma_real(z + 1.0_dp) / z
        else
            ! Integrate from epsilon to 1-epsilon to avoid singularities at endpoints
            gamma = pk_adaptative_simpson(1.0e-12_dp, 1.0_dp - 1.0e-12_dp, 1.0e-8_dp, 0, 20, integrating)
        end if

    contains

        function integrating(t) result(f)
            real(dp), intent(in) :: t
            real(dp) :: f
            
            f = ((t**(z - 1.0_dp)) * exp(-t/(1.0_dp - t))) / ((1.0_dp - t)**(z + 1.0_dp))
        end function integrating

    end function pk_gamma_real

    !=================================================
    ! Gamma function for complex numbers
    !=================================================
    ! z: input value
    ! gamma: output value
    !=================================================
    recursive function pk_gamma(z) result(gamma)
        complex(dp), intent(in) :: z
        complex(dp) :: gamma

        if (real(z) < 1.0_dp) then
            gamma = pk_gamma(z + 1.0_dp) / z
        else
            ! Integrate from epsilon to 1-epsilon along real axis
             gamma = pk_adaptative_simpson_complex(cmplx(1.0e-12_dp, 0.0_dp, kind=dp), &
                                                   cmplx(1.0_dp - 1.0e-12_dp, 0.0_dp, kind=dp), &
                                                   1.0e-8_dp, 0, 20, integrating)
        end if

    contains

        function integrating(t) result(f)
            complex(dp), intent(in) :: t
            complex(dp) :: f
            ! Integral_0^1 (log(1/t))^(z-1) dt
            f = (log(1.0_dp/t))**(z - 1.0_dp)
        end function integrating

    end function pk_gamma

    !=================================================
    ! Beta function
    !=================================================
    ! x: input value
    ! beta: output value
    !=================================================
    function pk_beta(x, y) result(beta)
        real(dp), intent(in) :: x, y
        real(dp) :: beta

        beta = pk_gamma_real(x) * pk_gamma_real(y) / pk_gamma_real(x + y)

    end function pk_beta

    !=================================================
    ! Permutation
    !=================================================
    ! n: total number of items
    ! r: number of items to choose
    ! permutation: output number of permutations
    !=================================================
    function pk_permutation(n, r) result(permutation)
        integer, intent(in) :: n, r
        integer :: permutation

        if (n < 0 .or. r < 0 .or. r > n) then
            print*, 'Error: n and r must be non-negative and r must be less than or equal to n'
            permutation = -1
            return
        end if

        permutation = pk_factorial(n) / pk_factorial(n - r)

    end function pk_permutation

    !=================================================
    ! Combination
    !=================================================
    ! n: total number of items
    ! r: number of items to choose
    ! combination: output number of combinations
    !=================================================
    function pk_combination(n, r) result(combination)
        integer, intent(in) :: n, r
        integer :: combination

        if (n < 0 .or. r < 0 .or. r > n) then
            print*, 'Error: n and r must be non-negative and r must be less than or equal to n'
            combination = -1
            return
        end if

        combination = pk_factorial(n) / (pk_factorial(r) * pk_factorial(n - r))

    end function pk_combination


end module physkit_special