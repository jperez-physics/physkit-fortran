
program numerical_test
    use physkit_numerical
    use physkit_constants, only: dp
    implicit none

    real(dp) :: x, dx, fnum, fann, dfnum, dfan, d2num, d2an
    real(dp) :: a, b, I_rect, I_trap, I_comp, I_adapt, I_exact
    real(dp) :: tol
    integer :: N_rect, N_comp, imax

    ! Parameters for derivative tests
    x = 2.0_dp
    dx = 1.0e-6_dp
    tol = 1.0e-6_dp

    ! Forward / Backward / Central derivatives
    dfnum = pk_forward_difference(x, dx, f)
    dfan = cos(x)
    print *, "pk_forward_difference: numerical =", dfnum, "analytical =", dfan, "err =", abs(dfnum-dfan)

    dfnum = pk_backward_difference(x, dx, f)
    print *, "pk_backward_difference: numerical =", dfnum, "analytical =", dfan, "err =", abs(dfnum-dfan)

    dfnum = pk_central_difference(x, dx, f)
    print *, "pk_central_difference: numerical =", dfnum, "analytical =", dfan, "err =", abs(dfnum-dfan)

    ! Second derivative
    d2num = pk_second_central_difference(x, dx, f)
    d2an = -sin(x)
    print *, "pk_second_central_difference: numerical =", d2num, "analytical =", d2an, "err =", abs(d2num-d2an)

    print *, "==============================================================================="
    ! Parameters for integral tests
    a = 0.0_dp
    b = acos(-1.0_dp)
    I_exact = 2.0_dp

    N_rect = 100
    N_comp = 100  ! must be even for composite Simpson
    imax = 20

    I_rect = pk_rectangular_rule(a, b, N_rect, f)
    print *, "pk_rectangular_rule: N=", N_rect, "I=", I_rect, "err=", abs(I_rect - I_exact)

    I_trap = pk_trapezoidal_rule(a, b, N_rect, f)
    print *, "pk_trapezoidal_rule: N=", N_rect, "I=", I_trap, "err=", abs(I_trap - I_exact)

    I_comp = pk_composite_simpson(a, b, N_comp, f)
    print *, "pk_composite_simpson: N=", N_comp, "I=", I_comp, "err=", abs(I_comp - I_exact)

    I_adapt = pk_adaptative_simpson(a, b, 1.0e-8_dp, 0, imax, f)
    print *, "pk_adaptative_simpson: tol=1e-8 I=", I_adapt, "err=", abs(I_adapt - I_exact)

contains

    real(dp) function f(x)
        real(dp), intent(in) :: x
        f = sin(x)
    end function f

end program numerical_test