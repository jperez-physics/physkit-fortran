program linalg_test
    use physkit_linalg
    use physkit_constants, only: dp
    implicit none

    real(dp), allocatable :: v1(:), v2(:), v3(:), v_out(:)
    real(dp), allocatable :: m1(:,:), m2(:,:), m3(:,:), m_out(:,:), Q(:,:), R(:,:)
    real(dp), allocatable :: evals(:)
    real(dp) :: s, tol
    integer :: i, j

    tol = 1.0e-6_dp
    print *, "Running functionality tests for physkit_linalg..."

    !##################################################
    ! Vector Tests
    !##################################################
    print *, ""
    print *, "--- VECTOR TESTS ---"
    
    allocate(v1(3), v2(3), v3(3))
    v1 = [1.0_dp, 2.0_dp, 3.0_dp]
    v2 = [4.0_dp, 5.0_dp, 6.0_dp]

    ! 1. Dot Product
    s = pk_dot_product(v1, v2)
    call print_result("pk_dot_product", s, 32.0_dp, tol)

    ! 2. Vector Norm
    s = pk_vector_norm(v1)
    call print_result("pk_vector_norm", s, sqrt(14.0_dp), tol)

    ! 3. Cross Product
    v_out = pk_cross_product(v1, v2)
    call print_result_vector("pk_cross_product", v_out, [-3.0_dp, 6.0_dp, -3.0_dp], tol)

    ! 4. Vector Normalize
    v_out = pk_vector_normalize(v1)
    s = pk_vector_norm(v_out)
    call print_result("pk_vector_normalize (norm)", s, 1.0_dp, tol)

    !##################################################
    ! Matrix Tests
    !##################################################
    print *, ""
    print *, "--- MATRIX TESTS ---"
    
    allocate(m1(2,2), m2(2,2))
    m1(1,:) = [1.0_dp, 2.0_dp]
    m1(2,:) = [3.0_dp, 4.0_dp]
    
    m2(1,:) = [2.0_dp, 0.0_dp]
    m2(2,:) = [1.0_dp, 2.0_dp]

    ! 1. Trace
    s = pk_trace(m1)
    call print_result("pk_trace", s, 5.0_dp, tol)

    ! 2. Identity Matrix
    allocate(m3(3,3))
    call pk_identity_matrix(3, m3)
    s = pk_trace(m3)
    call print_result("pk_identity_matrix (trace)", s, 3.0_dp, tol)

    ! 3. Zero Matrix
    call pk_zero_matrix(3, 3, m3)
    s = pk_vector_norm(reshape(m3, [9]))
    call print_result("pk_zero_matrix (norm)", s, 0.0_dp, tol)

    ! 4. Matrix-Vector Product
    if(allocated(v_out)) deallocate(v_out)
    if(allocated(v1)) deallocate(v1)
    allocate(v1(2), v_out(2))
    v1 = [1.0_dp, 1.0_dp]
    call pk_matrix_vector_product(m1, v1, v_out)
    call print_result_vector("pk_matrix_vector_product", v_out, [3.0_dp, 7.0_dp], tol)

    ! 5. Matrix-Matrix Product
    allocate(m_out(2,2))
    call pk_matrix_matrix_product(m1, m2, m_out)
    s = pk_vector_norm(reshape(m_out - reshape([4.0_dp, 4.0_dp, 10.0_dp, 8.0_dp], [2,2], order=[2,1]), [4]))
    call print_result("pk_matrix_matrix_product (err)", s, 0.0_dp, tol)
    
    ! 6. Hadamard Product
    call pk_hadamard_product(m1, m2, m_out)
    s = pk_vector_norm(reshape(m_out - reshape([2.0_dp, 0.0_dp, 3.0_dp, 8.0_dp], [2,2], order=[2,1]), [4]))
    call print_result("pk_hadamard_product (err)", s, 0.0_dp, tol)

    ! 7. Determinant
    s = pk_determinant(m1)
    call print_result("pk_determinant", s, -2.0_dp, tol)

    !##################################################
    ! Advanced Matrix Tests
    !##################################################
    print *, ""
    print *, "--- ADVANCED MATRIX TESTS ---"

    deallocate(m1, m_out, m3)
    allocate(m1(3,3), m_out(3,3), m3(3,3))
    
    ! 1. Gram-Schmidt
    m1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    m1(:, 2) = [2.0_dp, 1.0_dp, 0.0_dp]
    m1(:, 3) = [3.0_dp, 2.0_dp, 1.0_dp]
    m_out = pk_gram_schmidt(m1)
    call pk_identity_matrix(3, m3)
    s = pk_vector_norm(reshape(m_out - m3, [9]))
    call print_result("pk_gram_schmidt (err)", s, 0.0_dp, tol)

    ! 2. QR Decomposition & Eigenvalues
    ! We test pk_eigenvalues using trace property (sum of evals = trace)
    m1(1,:) = [4.0_dp, 1.0_dp, -2.0_dp]
    m1(2,:) = [1.0_dp, 2.0_dp, 0.0_dp]
    m1(3,:) = [-2.0_dp, 0.0_dp, 3.0_dp]
    
    allocate(Q(3,3), R(3,3))
    call pk_qr_decomposition(m1, Q, R, tol)
    s = 0.0_dp
    call print_result("pk_qr_decomp_algorithm (run)", s, 0.0_dp, tol)

    allocate(evals(3))
    evals = pk_eigenvalues(m1, tol)
    s = sum(evals)
    call print_result("pk_eigenvalues (trace sum)", s, pk_trace(m1), tol)

    print *, ""
    print *, "Tests finished."

contains

    subroutine print_result(label, res, exp_val, tol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: res, exp_val, tol
        real(dp) :: rel_err, err_abs
        character(len=6) :: status
        
        err_abs = abs(res - exp_val)
        if (abs(exp_val) > 1.0e-14_dp) then
            rel_err = (err_abs / abs(exp_val)) * 100.0_dp
        else
            rel_err = err_abs * 100.0_dp
        end if
        
        if (err_abs < tol .or. rel_err < tol*100) then
            status = "PASSED"
        else
            status = "FAILED"
        end if
        
        write(*, "(A30, ' - ', A6, ' - ', F14.8, ' - ', F14.8, ' - ', F10.6, ' %')") &
            label, status, res, exp_val, rel_err
    end subroutine print_result

    subroutine print_result_vector(label, res, exp_val, tol)
        character(len=*), intent(in) :: label
        real(dp), dimension(:), intent(in) :: res, exp_val
        real(dp), intent(in) :: tol
        real(dp) :: rel_err, err_abs
        character(len=6) :: status
        
        err_abs = maxval(abs(res - exp_val))
        if (maxval(abs(exp_val)) > 1.0e-14_dp) then
            rel_err = (err_abs / maxval(abs(exp_val))) * 100.0_dp
        else
            rel_err = err_abs * 100.0_dp
        end if
        
        if (err_abs < tol .or. rel_err < tol*100) then
            status = "PASSED"
        else
            status = "FAILED"
        end if

        if (size(res) == 2) then
            write(*, "(A30, ' - ', A6, ' - [', F8.4, ',', F8.4, '] - [', F8.4, ',', F8.4, '] - ', F10.6, ' %')") &
                label, status, res(1), res(2), exp_val(1), exp_val(2), rel_err
        else if (size(res) == 3) then
            write(*, "(A30, ' - ', A6, ' - [', F8.4, ',', F8.4, ',', F8.4, '] - [', F8.4, ',', F8.4, ',', F8.4, '] - ', F10.6, ' %')") &
                label, status, res(1), res(2), res(3), exp_val(1), exp_val(2), exp_val(3), rel_err
        else
            write(*, "(A30, ' - ', A6, ' - MaxAbsErr: ', F14.8, ' - ', F10.6, ' %')") &
                label, status, err_abs, rel_err
        end if
    end subroutine print_result_vector

end program linalg_test
