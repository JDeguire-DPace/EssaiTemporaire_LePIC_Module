! module mod_legacy_particle_globals
!   use iso_fortran_env, only: real64, int32
!   implicit none
!   public

!   ! -----------------------------
!   ! Legacy "part_info" globals that control loading
!   ! (these used to live in COMMON blocks)
!   ! -----------------------------
!   integer(int32) :: flag_restart = 0
!   integer(int32) :: flag_inj     = 1
!   integer(int32) :: opt_inj      = 1
!   integer(int32) :: flag_heat    = 0
!   integer(int32) :: flag_pdf     = 0
!   integer(int32) :: ptype_pdf    = 0

!   real(real64)   :: n0   = 0.0_real64
!   real(real64)   :: ngas = 0.0_real64

!   ! -----------------------------
!   ! Counters / geometry
!   ! -----------------------------
!   integer(int32) :: np_cell = 0
!   integer(int32) :: n_cell  = 0

!   real(real64)   :: x_load = 0.0_real64
!   real(real64)   :: ymax   = 0.0_real64
!   real(real64)   :: zmax   = 0.0_real64

!   ! heating window indices (legacy uses these)
!   integer(int32) :: ixl_pow = 1, ixr_pow = 1
!   integer(int32) :: iyl_pow = 1, iyr_pow = 1
!   integer(int32) :: izl_pow = 1, izr_pow = 1

!   ! Power absorption (legacy heating)
!   real(real64)   :: Pabs   = 0.0_real64

!   ! Per-species arrays used by legacy load_part
!   real(real64), allocatable :: vt0(:)
!   real(real64), allocatable :: Nm(:)

! contains

!   subroutine legacy_globals_init(ntype)
!     integer(int32), intent(in) :: ntype
!     if (allocated(vt0)) deallocate(vt0)
!     if (allocated(Nm))  deallocate(Nm)
!     allocate(vt0(ntype))
!     allocate(Nm(ntype))
!     vt0 = 0.0_real64
!     Nm  = 1.0_real64
!   end subroutine legacy_globals_init

! end module mod_legacy_particle_globals