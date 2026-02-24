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

    ! Generic interfface for Euler method
    interface pk_euler_method
        module procedure pk_euler_method_scalar
        module procedure pk_euler_method_system
    end interface

    ! Generic interface for the 2th order Runge-Kutta method
    interface pk_rk2
        module procedure pk_rk2_scalar
        module procedure pk_rk2_system
    end interface

    ! Generic interface for the 4th order Runge-Kutta method
    interface pk_rk4
        module procedure pk_rk4_scalar
        module procedure pk_rk4_system
    end interface

    !
    interface pk_velocity_verlet
        module procedure pk_velocity_verlet_scalar
        module procedure pk_velocity_verlet_system
    end interface

    ! Abstract interface for scalar ODE function f(x, y) = dy/dx
    abstract interface
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
    !> @param y_init Initial condition y(x0).
    !> @param y_out Integration result y(x1).
    !> @param f_ode Function defining the derivative dy/dx = f(x, y).
    !=================================================
    subroutine pk_euler_method_scalar(dx, x0, x1, y_init, y_out, f_ode)
        real(dp), intent(in)  :: dx, y_init, x0, x1
        real(dp), intent(out) :: y_out
        procedure(f) :: f_ode

        real(dp) :: y, x
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = y_init
        x = x0

        do j = 1, N
            y = y + dx * f_ode(x, y)
            x = x + dx
        end do

        y_out = y

    end subroutine pk_euler_method_scalar

    !
    subroutine pk_euler_method_system(dx, x0, x1, y_init, y_out, f_ode)
        real(dp), intent(in)  :: dx, x0, x1
        real(dp), dimension(:), intent(in) :: y_init
        real(dp), dimension(:), intent(out) :: y_out
        procedure(f_system) :: f_ode

        real(dp), dimension(size(y_init)) :: y
        real(dp) :: x
        integer :: j, N
        
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = y_init
        x = x0

        do j = 1, N
            y = y + dx * f_ode(x, y)
            x = x + dx
        end do

        y_out = y

    end subroutine pk_euler_method_system

    !=================================================
    !> @brief 2nd order Runge-Kutta (RK2) method for ODE integration.
    !> @param dx Integration step size.
    !> @param x0 Initial value of the independent variable.
    !> @param x1 Final value of the independent variable.
    !> @param y_init Initial condition y(x0).
    !> @param y_out Integration result y(x1).
    !> @param f_ode Function defining the derivative dy/dx = f(x, y).
    !=================================================
    subroutine pk_rk2_scalar(dx, x0, x1, y_init, y_out, f_ode)
        real(dp), intent(in)  :: dx, y_init, x0, x1
        real(dp), intent(out) :: y_out
        procedure(f) :: f_ode

        real(dp) :: y, x, k1, k2
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = y_init
        x = x0

        do j = 1, N
            k1 = f_ode(x, y)
            k2 = f_ode(x + dx/2.0_dp, y + dx*k1/2.0_dp)
            y = y + dx * k2
            x = x + dx
        end do

        y_out = y

    end subroutine pk_rk2_scalar

    !
    subroutine pk_rk2_system(dx, x0, x1, y_init, y_out, f_ode)
        real(dp), intent(in) :: dx, x0, x1
        real(dp), dimension(:), intent(in) :: y_init
        real(dp), dimension(:), intent(out) :: y_out
        procedure(f_system) :: f_ode
        
        real(dp), dimension(size(y_init)) :: y, k1, k2
        real(dp) :: x
        integer :: j, N

        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = y_init
        x = x0

        do j = 1, N
            k1 = f_ode(x, y)
            k2 = f_ode(x + dx/2.0_dp, y + dx*k1/2.0_dp)
            y = y + dx * k2
            x = x + dx
        end do

        y_out = y

    end subroutine pk_rk2_system

    !=================================================
    !> @brief 4th order Runge-Kutta (RK4) method for scalar ODEs.
    !> @param dx Integration step size.
    !> @param x0 Initial value of the independent variable.
    !> @param x1 Final value of the independent variable.
    !> @param y_init Initial condition y(x0).
    !> @param y_out Integration result y(x1).
    !> @param f_ode Function defining the derivative dy/dx = f(x, y).
    !=================================================
    subroutine pk_rk4_scalar(dx, x0, x1, y_init, y_out, f_ode)
        real(dp), intent(in)  :: dx, y_init, x0, x1
        real(dp), intent(out) :: y_out
        procedure(f) :: f_ode

        real(dp) :: y, x, k1, k2, k3, k4
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = y_init
        x = x0

        do j = 1, N
            k1 = f_ode(x, y)
            k2 = f_ode(x + dx/2.0_dp, y + dx*k1/2.0_dp)
            k3 = f_ode(x + dx/2.0_dp, y + dx*k2/2.0_dp)
            k4 = f_ode(x + dx,        y + dx*k3)
            y = y + dx*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
            x = x + dx
        end do

        y_out = y

    end subroutine pk_rk4_scalar

    !
    subroutine pk_rk4_system(dx, x0, x1, y_init, y_out, f_ode)
        real(dp), intent(in)  :: dx, x0, x1
        real(dp), dimension(:), intent(in) :: y_init
        real(dp), dimension(:), intent(out) :: y_out
        procedure(f_system) :: f_ode

        real(dp), dimension(size(y_init)) :: y, k1, k2, k3, k4
        real(dp) :: x
        integer :: j, N

        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = y_init
        x = x0

        do j = 1, N
            k1 = f_ode(x, y)
            k2 = f_ode(x + dx/2.0_dp, y + dx*k1/2.0_dp)
            k3 = f_ode(x + dx/2.0_dp, y + dx*k2/2.0_dp)
            k4 = f_ode(x + dx,        y + dx*k3)
            y = y + dx*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
            x = x + dx
        end do

        y_out = y

    end subroutine pk_rk4_system

    !=================================================
    !> @brief Velocity Verlet method for 2nd order ODEs.
    !> @details Ideal for dynamics problems where acceleration depends on position.
    !> @param dx Integration step size.
    !> @param x0 Initial value of the independent variable (e.g., time t0).
    !> @param x1 Final value of the independent variable (e.g., time t1).
    !> @param y_init Initial position y(x0).
    !> @param v_init Initial velocity v(x0) = dy/dx(x0).
    !> @param y_out Final position y(x1).
    !> @param v_out Final velocity v(x1).
    !> @param f_ode Function that calculates the acceleration a(x, y) = d2y/dx2.
    ! NOTE: Here 'f_ode(x,y)' must return the acceleration (second order)
    !=================================================
    subroutine pk_velocity_verlet_scalar(dx, x0, x1, y_init, v_init, y_out, v_out, f_ode)
        real(dp), intent(in)  :: dx, y_init, v_init, x0, x1
        real(dp), intent(out) :: y_out, v_out
        procedure(f) :: f_ode

        real(dp) :: y, v, x, a_now, a_new
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = y_init
        v = v_init
        x = x0

        do j = 1, N
            a_now = f_ode(x, y)
            y = y + v*dx + 0.5_dp*a_now*dx*dx
            a_new = f_ode(x + dx, y)
            v = v + 0.5_dp*(a_now + a_new)*dx
            x = x + dx
        end do

        y_out = y
        v_out = v

    end subroutine pk_velocity_verlet_scalar

    !
    subroutine pk_velocity_verlet_system(dx, x0, x1, y_init, v_init, y_out, v_out, f_ode)
        real(dp), intent(in)  :: dx, x0, x1
        real(dp), dimension(:), intent(in) :: y_init, v_init
        real(dp), dimension(:), intent(out) :: y_out, v_out
        procedure(f_system) :: f_ode

        real(dp), dimension(size(y_init)) :: y, v, a_now, a_new
        real(dp) :: x
        integer :: j, N
        if (dx == 0.0_dp) stop "Error: dx cannot be zero"
        N = nint((x1 - x0) / dx)
        y = y_init
        v = v_init
        x = x0

        do j = 1, N
            a_now = f_ode(x, y)
            y = y + v*dx + 0.5_dp*a_now*dx*dx
            a_new = f_ode(x + dx, y)
            v = v + 0.5_dp*(a_now + a_new)*dx
            x = x + dx
        end do

        y_out = y
        v_out = v

    end subroutine pk_velocity_verlet_system

end module physkit_ode