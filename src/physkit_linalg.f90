module physkit_linalg
    use physkit_constants, only: dp
    implicit none
    private
    public :: pk_dot, pk_norm, pk_cross, pk_normalize, pk_prodmatvec, pk_prodmatmat, pk_zeromatrix, pk_identitymatrix

contains

    !=================================================
    ! Dot product of two vectors
    !=================================================
    ! a: first vector (array)
    ! b: second vector (array)
    ! dotprod: dot product of a and b
    !=================================================
    function pk_dot(a, b) result(dotprod)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: dotprod
        integer :: i, N

        if (size(a) /= size(b)) stop "Error: vectors must be of the same size"
        N = size(a)
        dotprod = 0.0_dp

        do i = 1, N
            dotprod = dotprod + a(i) * b(i)
        end do

    end function pk_dot

    !=================================================
    ! Norm of a vector
    !=================================================
    ! a: vector (array)
    ! norm: norm of a vector a
    !=================================================
    function pk_norm(a) result(norm)
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

    end function pk_norm

    !=================================================
    ! Cross product of two 3D vectors
    !=================================================
    ! a: first vector (array)
    ! b: second vector (array)
    ! crossprod: cross product of a and b
    !=================================================
    function pk_cross(a, b) result(crossprod)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: crossprod(3)
        integer :: i, N

        if (size(a) /= 3 .or. size(b) /= 3) stop "Error: vectors must be 3D"
        crossprod(1) = a(2)*b(3) - a(3)*b(2)
        crossprod(2) = a(3)*b(1) - a(1)*b(3)
        crossprod(3) = a(1)*b(2) - a(2)*b(1)

    end function pk_cross

    !=================================================
    ! Normalization of a vector
    !=================================================
    ! a: vector (array)
    ! normalization: normalization of a vector a
    !=================================================
    function pk_normalize(a) result(normalization)
        real(dp), intent(in) :: a(:)
        real(dp) :: normalization(size(a))
        real(dp) :: norm
        integer :: i, N

        if (size(a) == 0) stop "Error: vector must be non-empty"
        N = size(a)
        norm = pk_norm(a)

        if (norm == 0.0_dp) stop "Error: cannot normalize zero vector"

        do i = 1, N
            normalization(i) = a(i) / norm
        end do

    end function pk_normalize

    !=================================================
    ! Matrix-vector multiplication
    !=================================================
    ! A: matrix (2D array)
    ! x: vector (array)
    ! y: output vector (array)
    !=================================================
    subroutine pk_prodmatvec(A, x, y)
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

    end subroutine pk_prodmatvec

    !=================================================
    ! Matrix-matrix multiplication
    !=================================================
    ! A: first matrix (2D array)
    ! B: second matrix (2D array)
    ! C: output matrix (2D array)
    !=================================================
    subroutine pk_prodmatmat(A, B, C)
        real(dp), intent(in) :: A(:, :), B(:, :)
        real(dp), intent(out) :: C(:, :)
        integer :: i, j, k, N, M, P

        N = size(A, 1)
        M = size(A, 2)
        P = size(B, 2)

        if (size(B, 1) /= M) stop "Error: incompatible matrix dimensions"

        do i = 1, N
            do j = 1, P
                C(i, j) = 0.0_dp
                do k = 1, M
                    C(i, j) = C(i, j) + A(i, k) * B(k, j)
                end do
            end do
        end do

    end subroutine pk_prodmatmat

    !=================================================
    ! Zero Matrix
    !=================================================
    ! N: number of rows
    ! M: number of columns
    ! Z: output zero matrix (2D array)
    !=================================================
    subroutine pk_zeromatrix(N, M, Z)
        integer, intent(in) :: N, M
        real(dp), intent(out) :: Z(:, :)
        integer :: i, j

        do i = 1, N
            do j = 1, M
                Z(i, j) = 0.0_dp
            end do
        end do

    end subroutine pk_zeromatrix

    !=================================================
    ! Identity Matrix
    !=================================================
    ! N: size of the identity matrix (N x N)
    ! I: output identity matrix (2D array)
    !=================================================
    subroutine pk_identitymatrix(N, I)
        integer, intent(in) :: N
        real(dp), intent(out) :: I(:, :)
        integer :: i, j

        do i = 1, N
            do j = 1, N
                if (i == j) then
                    I(i, j) = 1.0_dp
                else
                    I(i, j) = 0.0_dp
                end if
            end do
        end do

    end subroutine pk_identitymatrix

end module physkit_linalg