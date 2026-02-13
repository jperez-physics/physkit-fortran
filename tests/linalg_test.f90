program linalg_test
    use physkit_linalg
    use physkit_constants, only: dp
    implicit none

    real(dp), allocatable :: a(:), b(:)
    real(dp) :: res, expect, tol
    integer :: n

    n = 3
    allocate(a(n), b(n))
    a = [1.0_dp, 2.0_dp, 3.0_dp]
    b = [4.0_dp, 3.0_dp, 2.0_dp]

    res = dot(a, b)

    print *, "Dot product result:", res

end program linalg_test
