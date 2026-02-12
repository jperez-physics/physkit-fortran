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

end module physkit_constants
