!> @author Jaime (and contributors)
!> @brief Module for Newtonian Mechanics calculations.
!> @details Provides basic functions for point particle dynamics, including distance,
!>          displacement, energy, and momentum calculations.
module physkit_newtonian_mechanics
    use physkit_constants, only: dp, pk_success, pk_err_dimension_mismatch
    use physkit_ode
    use physkit_linalg
    use physkit_numerical

    implicit none
    private

    ! Public methods
    public :: pk_distance, pk_displacement, pk_kinetic_energy, &
              pk_linear_momentum, pk_angular_momentum, pk_work

    ! Abstract interface for vector functions (e.g., position as function of time)
    abstract interface
        function pk_vector_func(t) result(r)
            use physkit_constants, only: dp
            real(dp), intent(in) :: t
            real(dp) :: r(3)
        end function pk_vector_func
    end interface

contains

    !=================================================
    !> @brief Calculates the displacement vector between two points.
    !> @param r1 Initial position vector.
    !> @param r2 Final position vector.
    !> @param[out] ierr Optional error code (pk_success or pk_err_dimension_mismatch).
    !> @return Displacement vector (r2 - r1). Zero vector on error.
    !=================================================
    function pk_displacement(r1, r2, ierr) result(displacement)
        real(dp), intent(in) :: r1(:), r2(:)
        integer, intent(out), optional :: ierr
        real(dp) :: displacement(size(r1))

        if (present(ierr)) ierr = pk_success

        if (size(r1) /= size(r2)) then
            if (present(ierr)) then
                ierr = pk_err_dimension_mismatch
                displacement = 0.0_dp
                return
            else
                error stop "physkit_newtonian_mechanics: pk_displacement: vector sizes must be equal"
            end if
        end if

        displacement = r2 - r1

    end function pk_displacement

    !=================================================
    !> @brief Calculates the Euclidean distance between two points in space.
    !> @param r1 Coordinates of the first point.
    !> @param r2 Coordinates of the second point.
    !> @param[out] ierr Optional error code (pk_success or pk_err_dimension_mismatch).
    !> @return Scalar distance between r1 and r2 (0.0 on error).
    !=================================================
    function pk_distance(r1, r2, ierr) result(distance)
        real(dp), intent(in) :: r1(:), r2(:)
        integer, intent(out), optional :: ierr
        real(dp) :: distance

        if (present(ierr)) ierr = pk_success

        if (size(r1) /= size(r2)) then
            if (present(ierr)) then
                ierr = pk_err_dimension_mismatch
                distance = 0.0_dp
                return
            else
                error stop "physkit_newtonian_mechanics: pk_distance: vector sizes must be equal"
            end if
        end if

        distance = sqrt(sum((r2 - r1)**2))

    end function pk_distance

    !=================================================
    !> @brief Computes the kinetic energy of a point mass.
    !> @param m Mass of the particle.
    !> @param v Velocity vector (vx, vy, vz).
    !> @return Scalar kinetic energy T = 1/2 * m * v^2.
    !=================================================
    function pk_kinetic_energy(m, v) result(T)
        real(dp), intent(in) :: m, v(3)
        real(dp) :: T

        T = 0.5_dp * m * pk_dot_product(v, v)

    end function pk_kinetic_energy

    !=================================================
    !> @brief Computes the linear momentum of a point mass.
    !> @param m Mass of the particle.
    !> @param v Velocity vector.
    !> @return Linear momentum vector p = m * v.
    !=================================================
    function pk_linear_momentum(m, v) result(linear_momentum)
        real(dp), intent(in) :: m, v(3)
        real(dp) :: linear_momentum(3)

        linear_momentum = m * v

    end function pk_linear_momentum

    !=================================================
    !> @brief Computes the angular momentum of a particle.
    !> @param r Position vector relative to an origin.
    !> @param P Linear momentum vector.
    !> @return Angular momentum vector L = r x P.
    !=================================================
    function pk_angular_momentum(r, P) result(angular_momentum)
        real(dp), intent(in) :: r(3), P(3)
        real(dp) :: angular_momentum(3)

        angular_momentum = pk_cross_product(r, P)

    end function pk_angular_momentum

    !=================================================
    !> @brief Calculates the work done by a constant force over a displacement.
    !> @param F Constant force vector.
    !> @param r Displacement vector.
    !> @return Scalar work W = F . r.
    !=================================================
    function pk_work(F, r) result(work)
        real(dp), intent(in) :: F(3), r(3)
        real(dp) :: work

        work = pk_dot_product(F, r)

    end function pk_work

end module physkit_newtonian_mechanics
