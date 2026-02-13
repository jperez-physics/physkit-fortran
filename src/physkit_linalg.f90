module physkit_linalg
    use physkit_constants, only: dp
    implicit none
    private
    public :: dot

contains

    !=================================================
    ! Dot product of two vectors
    !=================================================
    ! a: first vector (array)
    ! b: second vector (array)
    ! Returns: dot product sum_i a(i)*b(i)
    ! Implemented using arrays and a DO loop
    !=================================================
    function dot(a, b) result(dotprod)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: dotprod
        integer :: i, N

        if (size(a) /= size(b)) stop "Error: vectors must be of the same size"
        N = size(a)
        dotprod = 0.0_dp

        do i = 1, N
            dotprod = dotprod + a(i) * b(i)
        end do

    end function dot

end module physkit_linalg