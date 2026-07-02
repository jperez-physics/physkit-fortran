!> @author Jaime (and contributors)
!> @brief Module for linear algebra operations.
!> @details Provides basic and advanced linear algebra tools including vector operations,
!>          matrix multiplications, products, and decompositions.
!>
!>          Error handling convention: procedures that can fail accept an
!>          optional `integer, intent(out) :: ierr` argument (last in the
!>          list). If the caller passes `ierr`, the procedure reports errors
!>          through it (see pk_success / pk_err_* in physkit_constants) and
!>          returns a well-defined "safe" result instead of leaving output
!>          arguments undefined. If the caller omits `ierr`, the procedure
!>          falls back to `error stop` with a descriptive message, matching
!>          the library's original fail-fast behavior.
module physkit_linalg
    use physkit_constants, only: dp, pk_success, pk_err_dimension_mismatch, &
                                  pk_err_invalid_argument, pk_err_singular, pk_err_no_convergence
    use physkit_numerical, only: pk_sum
    implicit none
    private

    ! Public methods
    public :: pk_dot_product, pk_vector_norm, pk_cross_product, pk_vector_normalize, &
              pk_matrix_vector_product, pk_matrix_matrix_product, pk_hadamard_product, &
              pk_zero_matrix, pk_identity_matrix, pk_trace, pk_determinant, pk_gram_schmidt, &
              pk_qr_decomposition, pk_eigenvalues

contains

    !=================================================
    !> @brief Computes the dot product of two vectors.
    !> @param a First vector.
    !> @param b Second vector.
    !> @param[out] ierr Optional error code (pk_success or pk_err_dimension_mismatch).
    !> @return Scaler dot product result (0.0 on error).
    !=================================================
    function pk_dot_product(a, b, ierr) result(dotprod)
        real(dp), intent(in) :: a(:), b(:)
        integer, intent(out), optional :: ierr
        real(dp) :: dotprod
        integer :: i, N

        if (present(ierr)) ierr = pk_success

        if (size(a) /= size(b)) then
            if (present(ierr)) then
                ierr = pk_err_dimension_mismatch
                dotprod = 0.0_dp
                return
            else
                error stop "physkit_linalg: pk_dot_product: vectors must be of the same size"
            end if
        end if

        N = size(a)
        dotprod = 0.0_dp

        do i = 1, N
            dotprod = dotprod + a(i) * b(i)
        end do

    end function pk_dot_product

    !=================================================
    !> @brief Computes the Euclidean norm (L2 norm) of a vector.
    !> @param a Input vector.
    !> @param[out] ierr Optional error code (pk_success or pk_err_invalid_argument).
    !> @return Euclidean norm of the vector (0.0 on error).
    !=================================================
    function pk_vector_norm(a, ierr) result(norm)
        real(dp), intent(in) :: a(:)
        integer, intent(out), optional :: ierr
        real(dp) :: norm
        integer :: i, N

        if (present(ierr)) ierr = pk_success

        if (size(a) == 0) then
            if (present(ierr)) then
                ierr = pk_err_invalid_argument
                norm = 0.0_dp
                return
            else
                error stop "physkit_linalg: pk_vector_norm: vector must be non-empty"
            end if
        end if

        N = size(a)
        norm = 0.0_dp

        do i = 1, N
            norm = norm + a(i) * a(i)
        end do

        norm = sqrt(norm)

    end function pk_vector_norm

    !=================================================
    !> @brief Computes the cross product of two 3-dimensional vectors.
    !> @param a First 3D vector.
    !> @param b Second 3D vector.
    !> @param[out] ierr Optional error code (pk_success or pk_err_dimension_mismatch).
    !> @return Resulting 3D vector perpendicular to both input vectors (zero vector on error).
    !=================================================
    function pk_cross_product(a, b, ierr) result(crossprod)
        real(dp), intent(in) :: a(:), b(:)
        integer, intent(out), optional :: ierr
        real(dp) :: crossprod(3)

        if (present(ierr)) ierr = pk_success

        if (size(a) /= 3 .or. size(b) /= 3) then
            if (present(ierr)) then
                ierr = pk_err_dimension_mismatch
                crossprod = 0.0_dp
                return
            else
                error stop "physkit_linalg: pk_cross_product: vectors must be 3D"
            end if
        end if

        crossprod(1) = a(2)*b(3) - a(3)*b(2)
        crossprod(2) = a(3)*b(1) - a(1)*b(3)
        crossprod(3) = a(1)*b(2) - a(2)*b(1)

    end function pk_cross_product

    !=================================================
    !> @brief Normalizes a vector to have unit length.
    !> @param a Input vector to be normalized.
    !> @param[out] ierr Optional error code (pk_success, pk_err_invalid_argument or pk_err_singular).
    !> @return Unit vector in the same direction as a (zero vector on error).
    !=================================================
    function pk_vector_normalize(a, ierr) result(normalization)
        real(dp), intent(in) :: a(:)
        integer, intent(out), optional :: ierr
        real(dp) :: normalization(size(a))
        real(dp) :: norm
        integer :: i, N

        if (present(ierr)) ierr = pk_success

        if (size(a) == 0) then
            if (present(ierr)) then
                ierr = pk_err_invalid_argument
                return
            else
                error stop "physkit_linalg: pk_vector_normalize: vector must be non-empty"
            end if
        end if

        N = size(a)
        norm = pk_vector_norm(a)

        if (norm == 0.0_dp) then
            if (present(ierr)) then
                ierr = pk_err_singular
                normalization = 0.0_dp
                return
            else
                error stop "physkit_linalg: pk_vector_normalize: cannot normalize zero vector"
            end if
        end if

        do i = 1, N
            normalization(i) = a(i) / norm
        end do

    end function pk_vector_normalize

    !=================================================
    !> @brief Performs matrix-vector multiplication y = A * x.
    !> @param A Input matrix of size (M, N).
    !> @param x Input vector of size (N).
    !> @param y Output vector of size (M).
    !> @param[out] ierr Optional error code (pk_success or pk_err_dimension_mismatch).
    !=================================================
    subroutine pk_matrix_vector_product(A, x, y, ierr)
        real(dp), intent(in) :: A(:, :), x(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out), optional :: ierr
        integer :: i, j, N, M

        if (present(ierr)) ierr = pk_success

        N = size(x)
        M = size(A, 1)

        if (size(A, 2) /= N) then
            if (present(ierr)) then
                ierr = pk_err_dimension_mismatch
                y = 0.0_dp
                return
            else
                error stop "physkit_linalg: pk_matrix_vector_product: incompatible matrix and vector dimensions"
            end if
        end if

        do i = 1, M
            y(i) = 0.0_dp
            do j = 1, N
                y(i) = y(i) + A(i, j) * x(j)
            end do
        end do

    end subroutine pk_matrix_vector_product

    !=================================================
    !> @brief Performs matrix-matrix multiplication C = A * B.
    !> @param A Input matrix of size (N, M).
    !> @param B Input matrix of size (M, P).
    !> @param C Output matrix of size (N, P).
    !> @param[out] ierr Optional error code (pk_success or pk_err_dimension_mismatch).
    !=================================================
    subroutine pk_matrix_matrix_product(A, B, C, ierr)
        real(dp), intent(in) :: A(:, :), B(:, :)
        real(dp), intent(out) :: C(:, :)
        integer, intent(out), optional :: ierr
        integer :: i, j, k, N, M, P

        if (present(ierr)) ierr = pk_success

        N = size(A, 1)
        M = size(A, 2)
        P = size(B, 2)

        if (size(B, 1) /= M) then
            if (present(ierr)) then
                ierr = pk_err_dimension_mismatch
                C = 0.0_dp
                return
            else
                error stop "physkit_linalg: pk_matrix_matrix_product: incompatible matrix dimensions"
            end if
        end if

        do i = 1, N
            do j = 1, P
                C(i, j) = 0.0_dp
                do k = 1, M
                    C(i, j) = C(i, j) + A(i, k) * B(k, j)
                end do
            end do
        end do

    end subroutine pk_matrix_matrix_product

    !=================================================
    !> @brief Computes the Hadamard product (element-wise multiplication) of two matrices C = A ⊙ B.
    !> @param A First matrix.
    !> @param B Second matrix.
    !> @param C Output matrix containing element-wise products.
    !> @param[out] ierr Optional error code (pk_success or pk_err_dimension_mismatch).
    !=================================================
    subroutine pk_hadamard_product(A, B, C, ierr)
        real(dp), intent(in) :: A(:, :), B(:, :)
        real(dp), intent(out) :: C(:, :)
        integer, intent(out), optional :: ierr
        integer :: i, j, N, M, P, Q

        if (present(ierr)) ierr = pk_success

        N = size(A, 1)
        M = size(A, 2)
        P = size(B, 1)
        Q = size(B, 2)

        if (P /= N .or. Q /= M) then
            if (present(ierr)) then
                ierr = pk_err_dimension_mismatch
                C = 0.0_dp
                return
            else
                error stop "physkit_linalg: pk_hadamard_product: incompatible matrix dimensions"
            end if
        end if

        do i = 1, N
            do j = 1, M
                C(i, j) = A(i, j) * B(i, j)
            end do
        end do

    end subroutine pk_hadamard_product

    !=================================================
    !> @brief Initializes a matrix with zeros.
    !> @param N Number of rows.
    !> @param M Number of columns.
    !> @param Z Output matrix of size (N, M).
    !=================================================
    subroutine pk_zero_matrix(N, M, Z)
        integer, intent(in) :: N, M
        real(dp), intent(out) :: Z(:, :)
        integer :: i, j

        do i = 1, N
            do j = 1, M
                Z(i, j) = 0.0_dp
            end do
        end do

    end subroutine pk_zero_matrix

    !=================================================
    !> @brief Generates an identity matrix of size N x N.
    !> @param N Dimension of the identity matrix.
    !> @param Iout Output identity matrix.
    !=================================================
    subroutine pk_identity_matrix(N, Iout)
        integer, intent(in) :: N
        real(dp), intent(out) :: Iout(:, :)
        integer :: i, j

        do i = 1, N
            do j = 1, N
                if (i == j) then
                    Iout(i, j) = 1.0_dp
                else
                    Iout(i, j) = 0.0_dp
                end if
            end do
        end do

    end subroutine pk_identity_matrix

    !=================================================
    !> @brief Computes the trace of a square matrix (sum of diagonal elements).
    !> @param A Square matrix.
    !> @param[out] ierr Optional error code (pk_success or pk_err_dimension_mismatch).
    !> @return Sum of elements on the main diagonal (0.0 on error).
    !=================================================
    function pk_trace(A, ierr) result(trace)
        real(dp), intent(in) :: A(:, :)
        integer, intent(out), optional :: ierr
        real(dp) :: trace
        integer :: i, N, M

        if (present(ierr)) ierr = pk_success

        N = size(A, 1)
        M = size(A, 2)
        trace = 0.0_dp

        if (N /= M) then
            if (present(ierr)) then
                ierr = pk_err_dimension_mismatch
                return
            else
                error stop "physkit_linalg: pk_trace: matrix must be square"
            end if
        end if

        do i = 1, N
            trace = trace + A(i, i)
        end do

    end function pk_trace

    !=================================================
    !> @brief Performs LU decomposition of a square matrix A = L * U.
    !> @param A Input square matrix.
    !> @param L Output lower triangular matrix with unit diagonal.
    !> @param U Output upper triangular matrix.
    !> @param[out] ierr Optional error code (pk_success or pk_err_dimension_mismatch).
    !=================================================
    subroutine pk_lu_decomposition(A, L, U, ierr)
        real(dp), intent(in) :: A(:, :)
        real(dp), intent(out) :: L(:, :), U(:, :)
        integer, intent(out), optional :: ierr
        real(dp) :: sum
        integer :: i, j, k, N, M

        if (present(ierr)) ierr = pk_success

        N = size(A, 1)
        M = size(A, 2)

        if (N /= M) then
            if (present(ierr)) then
                ierr = pk_err_dimension_mismatch
                L = 0.0_dp
                U = 0.0_dp
                return
            else
                error stop "physkit_linalg: pk_lu_decomposition: matrix must be square"
            end if
        end if

        call pk_identity_matrix(N, L)
        call pk_zero_matrix(N, M, U)

        do i = 1, N
            do j = i, N
                sum = 0.0_dp
                do k = 1, i - 1
                    sum = sum + L(i, k) * U(k, j)
                end do
                U(i, j) = A(i, j) - sum
            end do

            do j = i + 1, N
                sum = 0.0_dp
                do k = 1, i - 1
                    sum = sum + L(j, k) * U(k, i)
                end do
                L(j, i) = (A(j, i) - sum) / U(i, i)
            end do
        end do

    end subroutine pk_lu_decomposition

    !=================================================
    !> @brief Computes the determinant of a square matrix using LU decomposition.
    !> @param A Input square matrix.
    !> @param[out] ierr Optional error code (pk_success or pk_err_dimension_mismatch).
    !> @return Determinant of A (0.0 on error).
    !=================================================
    function pk_determinant(A, ierr) result(detA)
        real(dp), intent(in) :: A(:, :)
        integer, intent(out), optional :: ierr
        real(dp) :: detA
        real(dp), allocatable :: L(:, :), U(:, :)
        real(dp) :: detU
        integer :: i, n

        if (present(ierr)) ierr = pk_success

        n = size(A, 1)

        if (n /= size(A, 2)) then
            if (present(ierr)) then
                ierr = pk_err_dimension_mismatch
                detA = 0.0_dp
                return
            else
                error stop "physkit_linalg: pk_determinant: matrix must be square"
            end if
        end if

        allocate(L(n, n), U(n, n))
        call pk_lu_decomposition(A, L, U)

        detU = 1.0_dp
        do i = 1, n
            detU = detU * U(i, i)
        end do
        detA = detU

        deallocate(L, U)

    end function pk_determinant

    !=================================================
    !> @brief Performs Gram-Schmidt orthogonalization on a set of vectors.
    !> @param v Matrix whose columns are the input vectors.
    !> @return Matrix whose columns are the orthogonalized vectors.
    !=================================================
    function pk_gram_schmidt(v) result(u)
        real(dp), intent(in) :: v(:, :)
        real(dp) :: u(size(v, 1), size(v, 2))
        real(dp) :: proy(size(v, 1))
        real(dp) :: scalar
        integer :: j, k, num_vecs

        num_vecs = size(v, 2)

        if (num_vecs > 0) then
            u(:, 1) = v(:, 1)
        end if

        do k = 2, num_vecs
            u(:, k) = v(:, k)
            do j = 1, k - 1
                scalar = pk_dot_product(v(:, k), u(:, j)) / pk_dot_product(u(:, j), u(:, j))
                proy = scalar * u(:, j)
                u(:, k) = u(:, k) - proy
            end do
        end do

    end function pk_gram_schmidt

    !=================================================
    !> @brief Computes the QR decomposition of a matrix using the
    !>        modified Gram-Schmidt process, in a single pass.
    !> @details Decomposes A into Q (orthonormal columns) and R (upper
    !>          triangular) such that A = Q*R. This is a single-shot
    !>          decomposition, not an iterative algorithm; to compute
    !>          eigenvalues via the QR algorithm see pk_eigenvalues.
    !> @param[in]  A Input square matrix to decompose.
    !> @param[out] Q Orthonormal matrix.
    !> @param[out] R Upper triangular matrix.
    !=================================================
    subroutine pk_qr_decomposition(A, Q, R)
        real(dp), intent(in) :: A(:, :)
        real(dp), intent(out) :: Q(:, :), R(:, :)
        real(dp) :: q_temp(size(A, 1), size(A, 2))
        real(dp) :: Qt(size(A, 2), size(A, 1))
        integer :: j, M

        M = size(A, 2)
        q_temp = pk_gram_schmidt(A)
        do j = 1, M
            Q(:, j) = q_temp(:, j) / pk_vector_norm(q_temp(:, j))
        end do

        Qt = transpose(Q)
        call pk_matrix_matrix_product(Qt, A, R)

    end subroutine pk_qr_decomposition

    !=================================================
    !> @brief Estimates the eigenvalues of a square matrix using the
    !>        (unshifted) iterative QR algorithm.
    !> @details Repeatedly applies A_k+1 = R_k * Q_k, where Q_k, R_k are
    !>          the QR decomposition of A_k, until the sub-diagonal entries
    !>          vanish (convergence) or max_iter is reached. This method
    !>          reliably converges only for matrices with real eigenvalues
    !>          of distinct magnitude (e.g. symmetric matrices); matrices
    !>          with complex-conjugate eigenvalue pairs will not converge
    !>          to a diagonal form and max_iter guards against looping
    !>          forever in that case.
    !> @param[in] A Input square matrix.
    !> @param[in] tol Convergence tolerance for the sub-diagonal entries.
    !> @param[in] max_iter Optional cap on QR iterations (default: 1000).
    !> @param[out] ierr Optional error code (pk_success or pk_err_no_convergence).
    !>             If omitted, non-convergence is reported as a warning
    !>             (not fatal) and the best estimate is still returned.
    !> @return lambda Vector of estimated eigenvalues (real part only).
    !=================================================
    function pk_eigenvalues(A, tol, max_iter, ierr) result(lambda)
        real(dp), intent(in) :: A(:, :)
        real(dp), intent(in) :: tol
        integer, intent(in), optional :: max_iter
        integer, intent(out), optional :: ierr
        real(dp) :: Q(size(A, 1), size(A, 2)), R(size(A, 1), size(A, 2))
        real(dp) :: A_k(size(A, 1), size(A, 2))
        real(dp) :: sub_diag(size(A, 1))
        real(dp) :: lambda(size(A, 1))
        integer :: i, iter, iter_limit, N

        if (present(ierr)) ierr = pk_success

        N = size(A, 1)
        iter_limit = 1000
        if (present(max_iter)) iter_limit = max_iter

        A_k = A
        sub_diag = tol + 1.0_dp ! To start the loop
        iter = 0

        do while (maxval(abs(sub_diag)) > tol .and. iter < iter_limit)
            call pk_qr_decomposition(A_k, Q, R)
            call pk_matrix_matrix_product(R, Q, A_k)

            sub_diag = 0.0_dp
            do i = 1, N - 1
                sub_diag(i) = A_k(i+1, i)
            end do
            iter = iter + 1
        end do

        if (maxval(abs(sub_diag)) > tol) then
            if (present(ierr)) then
                ierr = pk_err_no_convergence
            else
                print *, "Warning: pk_eigenvalues did not converge within ", iter_limit, &
                          " iterations (matrix may have complex eigenvalues). Returning best estimate."
            end if
        end if

        do i = 1, N
            lambda(i) = A_k(i, i)
        end do
    end function pk_eigenvalues


end module physkit_linalg
