!> @author Jaime (and contributors)
!> @brief Module for numerical integration of Ordinary Differential Equations (ODEs).
!> @details Provides implementations of classic methods such as Euler, 2nd and 4th order Runge-Kutta,
!>          and the Velocity Verlet method for second-order equations.
module physkit_ode
    use physkit_constants, only: dp
    implicit none
    private

    ! Public module functions
    public :: pk_euler_method, pk_rk2, pk_rk4, pk_velocity_verlet

    ! Generic interface for the 4th order Runge-Kutta method
    interface pk_rk4
        module procedure pk_rk4_scalar
        module procedure pk_rk4_system
    end interface

    ! Abstract interface for scalar ODE function f(x, y) = dy/dx
    interface
        function f(x, y)
            use physkit_constants, only: dp
            real(dp) :: f
            real(dp), intent(in) :: x, y
        end function f

    ! Abstract interface for system of ODEs function f_system(x, y) = dy/dx
        function f_system(x, y)
            use physkit_constants, only: dp
            real(dp), dimension(:), allocatable :: f_system
            real(dp), intent(in) :: x
            real(dp), dimension(:), intent(in) :: y
        end function f_system
    end interface

contains

    !=================================================
    !> @brief Euler's method for ODE integration (1st order).
    !> @param dx Integration step size.
    !> @param x0 Initial value of the independent variable (e.g., initial time).
    !> @param x1 Final value of the independent variable (e.g., final time).
    !> @param yinitial Initial condition y(x0).
    !> @param ynext Integration result y(x1).
    !> @param yfunc Function defining the derivative dy/dx = f(x, y).
    !=================================================
    subroutine pk_euler_method(dx, x0, x1, yinitial, ynext, yfunc)
        real(dp), intent(in)  :: dx, yinitial, x0, x1
        real(dp), intent(out) :: ynext
        procedure(f) :: yfunc

        real(dp) :: y, x
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        x = x0

        do j = 1, N
            y = y + dx * yfunc(x, y)
            x = x + dx
        end do

        ynext = y

    end subroutine pk_euler_method

    !=================================================
    !> @brief 2nd order Runge-Kutta (RK2) method for ODE integration.
    !> @param dx Integration step size.
    !> @param x0 Initial value of the independent variable.
    !> @param x1 Final value of the independent variable.
    !> @param yinitial Initial condition y(x0).
    !> @param ynext Integration result y(x1).
    !> @param yfunc Function defining the derivative dy/dx = f(x, y).
    !=================================================
    subroutine pk_rk2(dx, x0, x1, yinitial, ynext, yfunc)
        real(dp), intent(in)  :: dx, yinitial, x0, x1
        real(dp), intent(out) :: ynext
        procedure(f) :: yfunc

        real(dp) :: y, x, k1, k2
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        x = x0

        do j = 1, N
            k1 = yfunc(x, y)
            k2 = yfunc(x + dx/2.0_dp, y + dx*k1/2.0_dp)
            y = y + dx * k2
            x = x + dx
        end do

        ynext = y

    end subroutine pk_rk2

    !=================================================
    !> @brief 4th order Runge-Kutta (RK4) method for scalar ODEs.
    !> @param dx Integration step size.
    !> @param x0 Initial value of the independent variable.
    !> @param x1 Final value of the independent variable.
    !> @param yinitial Initial condition y(x0).
    !> @param ynext Integration result y(x1).
    !> @param yfunc Function defining the derivative dy/dx = f(x, y).
    !=================================================
    subroutine pk_rk4_scalar(dx, x0, x1, yinitial, ynext, yfunc)
        real(dp), intent(in)  :: dx, yinitial, x0, x1
        real(dp), intent(out) :: ynext
        procedure(f) :: yfunc

        real(dp) :: y, x, k1, k2, k3, k4
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        x = x0

        do j = 1, N
            k1 = yfunc(x, y)
            k2 = yfunc(x + dx/2.0_dp, y + dx*k1/2.0_dp)
            k3 = yfunc(x + dx/2.0_dp, y + dx*k2/2.0_dp)
            k4 = yfunc(x + dx,        y + dx*k3)
            y = y + dx*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
            x = x + dx
        end do

        ynext = y

    end subroutine pk_rk4_scalar

    !=================================================
    !> @brief 4th order Runge-Kutta (RK4) method for systems of ODEs.
    !> @param dx Integration step size.
    !> @param x0 Initial value of the independent variable.
    !> @param x1 Final value of the independent variable.
    !> @param yinitial Initial condition vector y(x0).
    !> @param ynext Integration result vector y(x1).
    !> @param yfunc Function defining the derivative vector dy/dx = f(x, y).
    !=================================================
    subroutine pk_rk4_system(dx, x0, x1, yinitial, ynext, yfunc)
        real(dp), intent(in)  :: dx, x0, x1
        real(dp), dimension(:), intent(in) :: yinitial
        real(dp), dimension(:), intent(out) :: ynext
        procedure(f_system) :: yfunc

        real(dp), dimension(size(yinitial)) :: y, k1, k2, k3, k4
        real(dp) :: x
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        x = x0

        do j = 1, N
            k1 = yfunc(x, y)
            k2 = yfunc(x + dx/2.0_dp, y + dx*k1/2.0_dp)
            k3 = yfunc(x + dx/2.0_dp, y + dx*k2/2.0_dp)
            k4 = yfunc(x + dx,        y + dx*k3)
            y = y + dx*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
            x = x + dx
        end do

        ynext = y

    end subroutine pk_rk4_system

    !=================================================
    !> @brief Velocity Verlet method for 2nd order ODEs.
    !> @details Ideal for dynamics problems where acceleration depends on position.
    !> @param dx Integration step size.
    !> @param x0 Initial value of the independent variable (e.g., time t0).
    !> @param x1 Final value of the independent variable (e.g., time t1).
    !> @param yinitial Initial position y(x0).
    !> @param vinitial Initial velocity v(x0) = dy/dx(x0).
    !> @param ynext Final position y(x1).
    !> @param vnext Final velocity v(x1).
    !> @param yfunc Function that calculates the acceleration a(x, y) = d2y/dx2.
    ! NOTE: Here 'yfunc(x,y)' must return the acceleration (second order)
    !=================================================
    subroutine pk_velocity_verlet(dx, x0, x1, yinitial, vinitial, ynext, vnext, yfunc)
        real(dp), intent(in)  :: dx, yinitial, vinitial, x0, x1
        real(dp), intent(out) :: ynext, vnext
        procedure(f) :: yfunc

        real(dp) :: y, v, x, a, anext
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = yinitial
        v = vinitial
        x = x0

        do j = 1, N
            a = yfunc(x, y)
            y = y + v*dx + 0.5_dp*a*dx*dx
            anext = yfunc(x + dx, y)
            v = v + 0.5_dp*(a + anext)*dx
            x = x + dx
        end do

        ynext = y
        vnext = v

    end subroutine pk_velocity_verlet

end module physkit_ode