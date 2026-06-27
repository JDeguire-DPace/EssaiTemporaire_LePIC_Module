module mod_constants
  implicit none
  private
  public :: qe, eps0, c, pi, amu, eps0_si

  integer, parameter :: rk = 8
  real(rk), parameter :: qe=1.60217646d-19 ! Coulombs
  real(rk), parameter :: c=2.99792458d8 ! m/s
  real(rk), parameter :: eps0_si = 8.854187817d-12
  real(rk) :: eps0 = eps0_si
  real(rk), parameter :: pi=4.d0*datan(1.d0)
  real(rk), parameter :: amu=1.66053886d-27

end module mod_constants