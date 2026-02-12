module physkit_constants
  implicit none

  !=================================================
  ! Precision
  !=================================================
  integer, parameter :: dp = selected_real_kind(15, 307)

  !=================================================
  ! Mathematical constants
  !=================================================
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: euler = 2.71828182845904523536_dp
  real(dp), parameter :: tau = 2.0_dp * pi
  real(dp), parameter :: golden_ratio = 1.61803398874989484820_dp
  real(dp), parameter :: sqrt2 = 1.41421356237309504880_dp
  real(dp), parameter :: sqrt3 = 1.73205080756887729353_dp
  real(dp), parameter :: euler_mascheroni = 0.57721566490153286060_dp
  real(dp), parameter :: catalan = 0.91596559417721901505_dp
  real(dp), parameter :: apery_constant = 1.20205690315959428540_dp

  !=================================================
  ! Physical constants (SI units)
  !=================================================

  ! Speed of light in vacuum (m/s)
  real(dp), parameter :: speed_of_light = 2.99792458e8_dp

  ! Gravitational constant (m^3 kg^-1 s^-2)
  real(dp), parameter :: gravitational_constant = 6.67430e-11_dp

  ! Planck constant (J s)
  real(dp), parameter :: planck_constant = 6.62607015e-34_dp

  ! Reduced Planck constant (J s)
  real(dp), parameter :: hbar = planck_constant / (2.0_dp * pi)

  ! Elementary charge (C)
  real(dp), parameter :: elementary_charge = 1.602176634e-19_dp

  ! Boltzmann constant (J/K)
  real(dp), parameter :: boltzmann_constant = 1.380649e-23_dp

  ! Avogadro number (mol^-1)
  real(dp), parameter :: avogadro_number = 6.02214076e23_dp

  ! Vacuum permittivity (F/m)
  real(dp), parameter :: vacuum_permittivity = 8.8541878188e-12_dp

  ! Vacuum permeability (N/A^2)
  real(dp), parameter :: vacuum_permeability = 1.25663706127e-6_dp

  ! Fine-structure constant (dimensionless)
  real(dp), parameter :: fine_structure_constant = 7.2973525643e-3_dp

  ! Electron mass (kg)
  real(dp), parameter :: electron_mass = 9.1093837139e-31_dp

  ! Proton mass (kg)
  real(dp), parameter :: proton_mass = 1.67262192595e-27_dp

  ! Neutron mass (kg)
  real(dp), parameter :: neutron_mass = 1.67492750056e-27_dp

  ! Atomic unit of mass (kg)
  real(dp), parameter :: atomic_unit_mass = 9.1093837139e-31_dp

  ! Electron volt (J)
  real(dp), parameter :: electron_volt = 1.602176634e-19_dp

  ! Molar gas constant (J mol^-1 K^-1)
  real(dp), parameter :: gas_constant = 8.314462618_dp

  ! Stefan-Boltzmann constant (W m^-2 K^-4)
  real(dp), parameter :: stefan_boltzmann = 5.670374419e-8_dp

  ! Rydberg constant (m^-1)
  real(dp), parameter :: rydberg_constant = 10973731.568157_dp

  ! Standard acceleration of gravity (m/s^2)
  real(dp), parameter :: standard_gravity = 9.80665_dp
end module physkit_constants
