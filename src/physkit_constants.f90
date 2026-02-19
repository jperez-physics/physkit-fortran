!> @author Jaime (and contributors)
!> @brief Common mathematical and physical constants for Physkit.
!> @details This module provides a central location for fundamental constants and
!>          numerical precision parameters used across the library.
module physkit_constants
    implicit none

    !=================================================
    !> @name Precision Parameters
    !> @{
    !=================================================
    !> @brief Double precision (approx. 15 decimal digits).
    integer, parameter :: dp = selected_real_kind(15, 307)
    !> @}

    !=================================================
    !> @name Mathematical Constants
    !> @{
    !=================================================
    real(dp), parameter :: pi = 3.14159265358979323846_dp        !< Pi number.
    real(dp), parameter :: euler = 2.71828182845904523536_dp     !< Euler's number (base of natural logarithm).
    real(dp), parameter :: tau = 2.0_dp * pi                     !< Tau number (2 * Pi).
    real(dp), parameter :: golden_ratio = 1.61803398874989484820_dp !< Golden ratio (phi).
    real(dp), parameter :: euler_mascheroni = 0.57721566490153286060_dp !< Euler-Mascheroni constant (gamma).
    real(dp), parameter :: catalan = 0.91596559417721901505_dp   !< Catalan's constant.
    real(dp), parameter :: apery_constant = 1.20205690315959428540_dp !< Apery's constant (zeta(3)).
    !> @}

    !=================================================
    !> @name Physical Constants (SI Units)
    !> @{
    !=================================================
    real(dp), parameter :: speed_of_light = 2.99792458e8_dp         !< Speed of light in vacuum [m/s].
    real(dp), parameter :: gravitational_constant = 6.67430e-11_dp  !< Newtonian constant of gravitation [m^3 kg^-1 s^-2].
    real(dp), parameter :: planck_constant = 6.62607015e-34_dp      !< Planck constant [J s].
    real(dp), parameter :: hbar = planck_constant / (2.0_dp * pi)  !< Reduced Planck constant [J s].
    real(dp), parameter :: elementary_charge = 1.602176634e-19_dp   !< Elementary charge [C].
    real(dp), parameter :: boltzmann_constant = 1.380649e-23_dp     !< Boltzmann constant [J/K].
    real(dp), parameter :: avogadro_number = 6.02214076e23_dp       !< Avogadro number [mol^-1].
    real(dp), parameter :: vacuum_permittivity = 8.8541878188e-12_dp !< Vacuum permittivity (epsilon_0) [F/m].
    real(dp), parameter :: vacuum_permeability = 1.25663706127e-6_dp !< Vacuum permeability (mu_0) [N/A^2].
    real(dp), parameter :: fine_structure_constant = 7.2973525643e-3_dp !< Fine-structure constant (alpha) [dimensionless].
    real(dp), parameter :: electron_mass = 9.1093837139e-31_dp      !< Electron mass [kg].
    real(dp), parameter :: proton_mass = 1.67262192595e-27_dp        !< Proton mass [kg].
    real(dp), parameter :: neutron_mass = 1.67492750056e-27_dp       !< Neutron mass [kg].
    real(dp), parameter :: atomic_unit_mass = 1.66053906660e-27_dp   !< Atomic mass constant (u) [kg].
    real(dp), parameter :: electron_volt = 1.602176634e-19_dp       !< Electron volt [J].
    real(dp), parameter :: gas_constant = 8.314462618_dp            !< Molar gas constant [J mol^-1 K^-1].
    real(dp), parameter :: stefan_boltzmann = 5.670374419e-8_dp     !< Stefan-Boltzmann constant [W m^-2 K^-4].
    real(dp), parameter :: rydberg_constant = 10973731.568157_dp    !< Rydberg constant [m^-1].
    real(dp), parameter :: standard_gravity = 9.80665_dp            !< Standard acceleration of gravity [m/s^2].
    !> @}

end module physkit_constants
