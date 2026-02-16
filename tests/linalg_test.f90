program linalg_test
    use physkit_linalg
    use physkit_constants, only: dp
    implicit none

    real(dp), allocatable :: v1(:), v2(:), v3(:), v_out(:)
    real(dp), allocatable :: m1(:,:), m2(:,:), m3(:,:), m_out(:,:)
    real(dp) :: s, tol
    integer :: i, j

    tol = 1.0e-6_dp
    print *, "Running functionality tests for physkit_linalg..."

    !##################################################
    ! Vector Tests
    !##################################################
    print *, "--- Vector Tests ---"
    
    allocate(v1(3), v2(3), v3(3))
    v1 = [1.0_dp, 2.0_dp, 3.0_dp]
    v2 = [4.0_dp, 5.0_dp, 6.0_dp]

    ! Dot Product
    s = pk_dot_product(v1, v2)
    ! 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    print *, "pk_dot_product ([1,2,3], [4,5,6]): ", s
    if (abs(s - 32.0_dp) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED: Expected 32.0, Got ", s
    end if

    ! Vector Norm
    s = pk_vector_norm(v1)
    ! sqrt(1+4+9) = sqrt(14) approx 3.741657
    print *, "pk_vector_norm ([1,2,3]): ", s
    if (abs(s - sqrt(14.0_dp)) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Cross Product
    v_out = pk_cross_product(v1, v2)
    ! [1,2,3] x [4,5,6] = [-3, 6, -3]
    ! x = 2*6 - 3*5 = 12 - 15 = -3
    ! y = 3*4 - 1*6 = 12 - 6 = 6
    ! z = 1*5 - 2*4 = 5 - 8 = -3
    print *, "pk_cross_product ([1,2,3], [4,5,6]): ", v_out
    if (abs(v_out(1) - (-3.0_dp)) < tol .and. &
        abs(v_out(2) - ( 6.0_dp)) < tol .and. &
        abs(v_out(3) - (-3.0_dp)) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Vector Normalize
    v_out = pk_vector_normalize(v1)
    s = pk_vector_norm(v_out)
    print *, "pk_vector_normalize ([1,2,3]) norm: ", s
    if (abs(s - 1.0_dp) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    !##################################################
    ! Matrix Tests
    !##################################################
    print *, "--- Matrix Tests ---"
    
    allocate(m1(2,2), m2(2,2))
    m1(1,:) = [1.0_dp, 2.0_dp]
    m1(2,:) = [3.0_dp, 4.0_dp]
    
    m2(1,:) = [2.0_dp, 0.0_dp]
    m2(2,:) = [1.0_dp, 2.0_dp]

    ! Trace
    s = pk_trace(m1)
    ! 1 + 4 = 5
    print *, "pk_trace ([[1,2],[3,4]]): ", s
    if (abs(s - 5.0_dp) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Identity Matrix
    allocate(m3(3,3))
    call pk_identity_matrix(3, m3)
    print *, "pk_identity_matrix (3x3):"
    do i=1,3
        print *, m3(i,:)
    end do
    s = pk_trace(m3)
    if (abs(s - 3.0_dp) < tol .and. m3(1,2) == 0.0_dp) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Zero Matrix
    call pk_zero_matrix(3, 3, m3)
    print *, "pk_zero_matrix (3x3) norm: ", pk_vector_norm(reshape(m3, [9])) ! treat as vector to check all zeros
    if (pk_vector_norm(reshape(m3, [9])) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Matrix-Vector Product
    ! m1 = [[1,2],[3,4]], v = [1,1] (using first 2 of v1 which is treated as 3D above, let's realloc)
    deallocate(v_out, v1)
    allocate(v1(2), v_out(2))
    v1 = [1.0_dp, 1.0_dp]
    call pk_matrix_vector_product(m1, v1, v_out)
    ! [1*1 + 2*1, 3*1 + 4*1] = [3, 7]
    print *, "pk_matrix_vector_product ([[1,2],[3,4]], [1,1]): ", v_out
    if (abs(v_out(1) - 3.0_dp) < tol .and. abs(v_out(2) - 7.0_dp) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    ! Matrix-Matrix Product
    ! m1 = [[1,2],[3,4]]
    ! m2 = [[2,0],[1,2]]
    ! C = m1 * m2
    ! c11 = 1*2 + 2*1 = 4
    ! c12 = 1*0 + 2*2 = 4
    ! c21 = 3*2 + 4*1 = 10
    ! c22 = 3*0 + 4*2 = 8
    allocate(m_out(2,2))
    call pk_matrix_matrix_product(m1, m2, m_out)
    print *, "pk_matrix_matrix_product:"
    do i=1,2
        print *, m_out(i,:)
    end do
    if (abs(m_out(1,1) - 4.0_dp) < tol .and. &
        abs(m_out(1,2) - 4.0_dp) < tol .and. &
        abs(m_out(2,1) - 10.0_dp) < tol .and. &
        abs(m_out(2,2) - 8.0_dp) < tol) then
        print *, "  PASSED"
    else
        print *, "  FAILED"
    end if

    print *, "Tests finished."

end program linalg_test
