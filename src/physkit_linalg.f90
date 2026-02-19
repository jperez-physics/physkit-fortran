!> @author Jaime (and contributors)
!> @brief Module for linear algebra operations.
!> @details Provides basic and advanced linear algebra tools including vector operations,
!>          matrix multiplications, products, and decompositions.
module physkit_linalg
    use physkit_constants, only: dp
    implicit none
    private

    ! Public methods
    public :: pk_dot_product, pk_vector_norm, pk_cross_product, pk_vector_normalize, &
              pk_matrix_vector_product, pk_matrix_matrix_product, pk_hadamard_product, &
              pk_zero_matrix, pk_identity_matrix, pk_trace

contains

    !=================================================
    !> @brief Computes the dot product of two vectors.
    !> @param a First vector.
    !> @param b Second vector.
    !> @return Scaler dot product result.
    !=================================================
    function pk_dot_product(a, b) result(dotprod)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: dotprod
        integer :: i, N

        if (size(a) /= size(b)) stop "Error: vectors must be of the same size"
        N = size(a)
        dotprod = 0.0_dp

        do i = 1, N
            dotprod = dotprod + a(i) * b(i)
        end do

    end function pk_dot_product

    !=================================================
    !> @brief Computes the Euclidean norm (L2 norm) of a vector.
    !> @param a Input vector.
    !> @return Euclidean norm of the vector.
    !=================================================
    function pk_vector_norm(a) result(norm)
        real(dp), intent(in) :: a(:)
        real(dp) :: norm
        integer :: i, N

        if (size(a) == 0) stop "Error: vector must be non-empty"
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
    !> @return Resulting 3D vector perpendicular to both input vectors.
    !=================================================
    function pk_cross_product(a, b) result(crossprod)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: crossprod(3)
        integer :: i, N

        if (size(a) /= 3 .or. size(b) /= 3) stop "Error: vectors must be 3D"
        crossprod(1) = a(2)*b(3) - a(3)*b(2)
        crossprod(2) = a(3)*b(1) - a(1)*b(3)
        crossprod(3) = a(1)*b(2) - a(2)*b(1)

    end function pk_cross_product

    !=================================================
    !> @brief Normalizes a vector to have unit length.
    !> @param a Input vector to be normalized.
    !> @return Unit vector in the same direction as a.
    !=================================================
    function pk_vector_normalize(a) result(normalization)
        real(dp), intent(in) :: a(:)
        real(dp) :: normalization(size(a))
        real(dp) :: norm
        integer :: i, N

        if (size(a) == 0) stop "Error: vector must be non-empty"
        N = size(a)
        norm = pk_vector_norm(a)

        if (norm == 0.0_dp) stop "Error: cannot normalize zero vector"

        do i = 1, N
            normalization(i) = a(i) / norm
        end do

    end function pk_vector_normalize

    !=================================================
    !> @brief Performs matrix-vector multiplication y = A * x.
    !> @param A Input matrix of size (M, N).
    !> @param x Input vector of size (N).
    !> @param y Output vector of size (M).
    !=================================================
    subroutine pk_matrix_vector_product(A, x, y)
        real(dp), intent(in) :: A(:, :), x(:)
        real(dp), intent(out) :: y(:)
        integer :: i, j, N, M

        N = size(x)
        M = size(A, 1)

        if (size(A, 2) /= N) stop "Error: incompatible matrix and vector dimensions"

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
    !=================================================
    subroutine pk_matrix_matrix_product(A, B, C)
        real(dp), intent(in) :: A(:, :), B(:, :)
        real(dp), intent(out) :: C(:, :)
        integer :: i, j, k, N, M, P

        N = size(A, 1)
        M = size(A, 2)
        P = size(B, 2)

        if (size(B, 1) /= M) then
            print*, "Error: incompatible matrix dimensions"
            return
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
    !=================================================
    subroutine pk_hadamard_product(A, B, C)
        real(dp), intent(in) :: A(:, :), B(:, :)
        real(dp), intent(out) :: C(:, :)
        integer :: i, j, N, M, P, Q

        N = size(A, 1)
        M = size(A, 2)
        P = size(B, 1)
        Q = size(B, 2)

        if (P /= N .or. Q /= M) then
            print*, "Error: incompatible matrix dimensions"
            return
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
    !> @return Sum of elements on the main diagonal. Returns -1 on non-square matrices.
    !=================================================
    function pk_trace(A) result(trace)
        real(dp), intent(in) :: A(:, :)
        real(dp) :: trace
        integer :: i, N, M

        N = size(A, 1)
        M = size(A, 2)
        trace = 0.0_dp
        
        if (N /= M) then 
            print*, "Error: incompatible matrix dimensions"
            trace = -1.0_dp
            return
        else
            do i = 1, N
                trace = trace + A(i, i)
            end do
        end if

    end function pk_trace

    !=================================================
    !> @brief Performs LU decomposition of a square matrix A = L * U.
    !> @param A Input square matrix.
    !> @param L Output lower triangular matrix with unit diagonal.
    !> @param U Output upper triangular matrix.
    !=================================================
    subroutine pk_lu_decomposition(A, L, U)
        real(dp), intent(in) :: A(:, :)
        real(dp), intent(out) :: L(:, :), U(:, :)
        integer :: i, j, k, N, M

        N = size(A, 1)
        M = size(A, 2)

        if (N /= M) then
            print*, "Error: incompatible matrix dimensions"
            return
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

end module physkit_linalg