module mod_legacy_particle_globals
  use iso_fortran_env, only: real64, int32
  use mod_config, only: npart
  use mod_config, only: Config
  implicit none
  public

  ! These are the “globals” load_part_OMP expects
  integer(int32) :: np_cell = 0
  integer(int32) :: n_cell  = 0

  real(real64) :: x_load = 0.0_real64
  real(real64) :: ymax   = 0.0_real64
  real(real64) :: zmax   = 0.0_real64

  real(real64) :: Pabs = 0.0_real64

  integer(int32) :: ixl_pow = 1, ixr_pow = 1
  integer(int32) :: iyl_pow = 1, iyr_pow = 1
  integer(int32) :: izl_pow = 1, izr_pow = 1

  real(real64) :: vt0(npart) = 0.0_real64
  real(real64) :: Nm (npart) = 0.0_real64

  

end module mod_legacy_particle_globals