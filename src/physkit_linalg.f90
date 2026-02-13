module physkit_linalg
    use physkit_constants, only: dp
    implicit none
    private
    public :: pk_dot, pk_norm, pk_cross, pk_normalize

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

end module physkit_linalg