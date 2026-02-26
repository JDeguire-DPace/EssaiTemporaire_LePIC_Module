module mod_legacy_particle_globals
  use iso_fortran_env, only: real64, int32
  implicit none

  ! Globals expected by legacy load_part
  integer(int32) :: np_cell = 0
  integer(int32) :: n_cell  = 0

  real(real64) :: x_load = 0.0_real64
  real(real64) :: ymax   = 0.0_real64
  real(real64) :: zmax   = 0.0_real64

  real(real64) :: Pabs   = 0.0_real64

  integer(int32) :: ixl_pow = 1
  integer(int32) :: ixr_pow = 1
  integer(int32) :: iyl_pow = 1
  integer(int32) :: iyr_pow = 1
  integer(int32) :: izl_pow = 1
  integer(int32) :: izr_pow = 1

  ! Per-species arrays (allocate at runtime once ntype is known)
  real(real64), allocatable :: vt0(:)
  real(real64), allocatable :: Nm(:)

contains

  subroutine legacy_globals_ensure_species(ntype)
    integer(int32), intent(in) :: ntype
    if (.not. allocated(vt0)) allocate(vt0(ntype))
    if (.not. allocated(Nm))  allocate(Nm(ntype))
    if (size(vt0) /= ntype) then
      deallocate(vt0); allocate(vt0(ntype))
    end if
    if (size(Nm) /= ntype) then
      deallocate(Nm); allocate(Nm(ntype))
    end if
    vt0 = 0.0_real64
    Nm  = 1.0_real64
  end subroutine legacy_globals_ensure_species

end module mod_legacy_particle_globals