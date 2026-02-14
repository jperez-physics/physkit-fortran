program linalg_test
    use physkit_linalg
    use physkit_constants, only: dp
    implicit none

    integer, parameter :: n3 = 3
    integer, parameter :: n_rows = 2, n_cols = 3

    real(dp), allocatable :: a(:), b(:), anorm(:), cross(:)
    real(dp), allocatable :: x(:), y(:)
    real(dp), allocatable :: matA(:,:), matB(:,:), matC(:,:), matZ(:,:), matI(:,:)
    real(dp) :: dotres, normres
    integer :: i, j

    allocate(a(n3), b(n3))
    a = [1.0_dp, 2.0_dp, 3.0_dp]
    b = [4.0_dp, 3.0_dp, 2.0_dp]

    ! Dot product
    dotres = pk_dot_product(a, b)
    print *, "pk_dot_product(a,b) =", dotres

    ! Vector norm
    normres = pk_vector_norm(a)
    print *, "pk_vector_norm(a) =", normres

    ! Cross product (3D)
    cross = pk_cross_product(a, b)
    print *, "pk_cross_product(a,b) ="
    do i = 1, n3
        print *, '  element(', i, ') =', cross(i)
    end do

    ! Vector normalization
    anorm = pk_vector_normalize(a)
    print *, "pk_vector_normalize(a) ="
    do i = 1, n3
        print *, '  element(', i, ') =', anorm(i)
    end do

    ! Matrix-vector product
    allocate(matA(n_rows, n_cols), x(n_cols), y(n_rows))
    matA = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], shape(matA))
    x = [1.0_dp, 0.0_dp, -1.0_dp]
    call pk_matrix_vector_product(matA, x, y)
    print *, "pk_matrix_vector_product(matA,x) ="
    do i = 1, n_rows
        print *, '  y(', i, ') =', y(i)
    end do

    ! Matrix-matrix product
    allocate(matB(n_cols, n_rows), matC(n_rows, n_rows))
    matB = reshape([7.0_dp, 8.0_dp, 9.0_dp, 10.0_dp, 11.0_dp, 12.0_dp], shape(matB))
    call pk_matrix_matrix_product(matA, matB, matC)
    print *, "pk_matrix_matrix_product(matA,matB) ="
    do i = 1, n_rows
        write(*,'(100(f10.4,1x))') (matC(i,j), j=1,n_rows)
    end do

    ! Zero matrix
    allocate(matZ(n_rows, n_cols))
    call pk_zero_matrix(n_rows, n_cols, matZ)
    print *, "pk_zero_matrix(2,3) ="
    do i = 1, n_rows
        write(*,'(100(f10.4,1x))') (matZ(i,j), j=1,n_cols)
    end do

    ! Identity matrix
    allocate(matI(n3, n3))
    call pk_identity_matrix(n3, matI)
    print *, "pk_identity_matrix(3) ="
    do i = 1, n3
        write(*,'(100(f10.4,1x))') (matI(i,j), j=1,n3)
    end do

end program linalg_test
