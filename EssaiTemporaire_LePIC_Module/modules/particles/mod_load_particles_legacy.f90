! module mod_load_particles_legacy
!   use iso_fortran_env, only: real64, int32
!   use omp_lib,         only: omp_get_max_threads
!   use mod_config,      only: Config
!   use mod_particles,   only: ParticleSet
!   use mod_legacy_particle_globals, only: &
!       np_cell, n_cell, x_load, ymax, zmax, Pabs, &
!       ixl_pow, ixr_pow, iyl_pow, iyr_pow, izl_pow, izr_pow, &
!       flag_restart, flag_inj, opt_inj, flag_heat, flag_pdf, ptype_pdf, &
!       n0, ngas, legacy_globals_init
!   implicit none
!   private
!   public :: load_particles_legacy

!   interface
!     subroutine load_part(n,h,bcnd,np,vxp,ntype,nmax,kq,ni0,np_tot,nproc,iseed,sum_dEk,Nh,mpi_rank,nproc_mpi)
!       use iso_fortran_env, only: real64
!       implicit none
!       integer :: mpi_rank, nproc_mpi
!       integer :: ntype, nmax, n(3), nproc
!       integer :: Nh(nproc)
!       real(real64) :: h(3)
!       integer :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
!       integer :: np_tot(ntype,nproc)
!       real(real64) :: vxp(6,nmax,ntype,nproc)
!       real(real64) :: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)
!       real(real64) :: ni0(*)
!       real(real64) :: sum_dEk(nproc)
!       real(real64) :: kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)
!       integer :: iseed(nproc)
!     end subroutine load_part
!   end interface

! contains

!   subroutine load_particles_legacy(cfg, mpi_rank, mpi_size, n, h, bcnd, kq, ntype_all, ntype_trk, iseed, part)
!     type(Config),   intent(in)    :: cfg
!     integer,        intent(in)    :: mpi_rank, mpi_size
!     integer(int32), intent(in)    :: n(3)
!     real(real64),   intent(in)    :: h(3)
!     integer,        intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
!     real(real64),   intent(in)    :: kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)
!     integer(int32), intent(in)    :: ntype_all
!     integer(int32), intent(in)    :: ntype_trk
!     integer(int32), intent(inout) :: iseed(:)
!     type(ParticleSet), intent(inout) :: part(:,:)

!     integer(int32) :: nproc
!     integer(int32) :: nmax, jmax_est
!     real(real64)   :: tmp
!     integer        :: iproc, it

!     integer(int32) :: ntype_call
!     integer(int32) :: ix_load

!     real(real64), allocatable :: vxp(:,:,:,:)
!     real(real64), allocatable :: np_arr(:,:,:,:,:)
!     integer,      allocatable :: np_tot(:,:)
!     integer,      allocatable :: Nh(:)
!     real(real64), allocatable :: sum_dEk(:)
!     real(real64), allocatable :: ni0_legacy(:)

!     ! runtime OpenMP thread count
!     nproc = int(max(1, omp_get_max_threads()), int32)

!     ! call legacy loader with charged species only
!     ntype_call = ntype_trk

!     ! sanity checks
!     if (size(part,1) /= ntype_trk) then
!       error stop "load_particles_legacy: part first dim must equal ntype_trk"
!     end if
!     if (size(part,2) < nproc) then
!       error stop "load_particles_legacy: part second dim < omp_get_max_threads()"
!     end if
!     if (size(iseed) < nproc) then
!       error stop "load_particles_legacy: iseed(:) smaller than omp_get_max_threads()"
!     end if

!     ! ---- Config -> legacy globals ----
!     np_cell      = int(cfg%np_cell, int32)
!     n0           = cfg%n0
!     ngas         = cfg%ngas
!     Pabs         = cfg%Pabs
!     flag_restart = int(cfg%flag_restart, int32)
!     flag_heat    = int(cfg%flag_heat, int32)
!     flag_pdf     = int(cfg%flag_pdf, int32)
!     ptype_pdf    = int(cfg%ptype_pdf, int32)

!     flag_inj = int(cfg%flag_inj, int32)
!     if (flag_inj == 0_int32) flag_inj = 1_int32

!     opt_inj = int(cfg%opt_inj, int32)
!     if (opt_inj == 0_int32) opt_inj = 1_int32

!     ! ---- Legacy-style particle loading region: x <= x_load ----
!     x_load = cfg%x_load
!     if (x_load <= 0.0_real64) then
!       x_load = real(n(1), real64) * h(1)
!     end if

!     ! ix_load = floor(x_load/dx) + 1, clamped to [1, n(1)+1]
!     ix_load = int(floor(x_load / h(1)), int32) + 1_int32
!     if (ix_load < 1_int32) ix_load = 1_int32
!     if (ix_load > n(1)+1_int32) ix_load = n(1)+1_int32

!     ! number of x-cells in load region is (ix_load-1)
!     n_cell = (ix_load - 1_int32) * n(2) * n(3)

!     ymax = real(n(2), real64) * h(2)
!     zmax = real(n(3), real64) * h(3)

!     ixl_pow = clamp_index(cfg%xl_pow, h(1), n(1))
!     ixr_pow = clamp_index(cfg%xr_pow, h(1), n(1))
!     iyl_pow = clamp_index(cfg%yl_pow, h(2), n(2))
!     iyr_pow = clamp_index(cfg%yr_pow, h(2), n(2))
!     izl_pow = clamp_index(cfg%zl_pow, h(3), n(3))
!     izr_pow = clamp_index(cfg%zr_pow, h(3), n(3))

!     call legacy_globals_init(ntype_call)

!     ! ---- acceptance probabilities (NOT normalized) ----
!     allocate(ni0_legacy(max(10_int32, ntype_call)))
!     ni0_legacy = 0.0_real64
!     do it = 1, ntype_call
!       ni0_legacy(it) = 1.0_real64
!     end do

!     ! ---- sizing ----
!     tmp = real(np_cell,real64) * real(n_cell,real64) / max(1.0_real64, real(mpi_size*nproc,real64))
!     jmax_est = max(1_int32, int(nint(tmp), int32))
!     nmax = int(1.3_real64*real(jmax_est,real64), int32) + 32_int32

!     allocate(vxp(6, nmax, ntype_call, nproc))
!     allocate(np_arr(0:n(1)+2,0:n(2)+2,0:n(3)+2, ntype_call, nproc))
!     allocate(np_tot(ntype_call, nproc))
!     allocate(Nh(nproc))
!     allocate(sum_dEk(nproc))

!     vxp     = 0.0_real64
!     np_arr  = 0.0_real64
!     np_tot  = 0
!     Nh      = 0
!     sum_dEk = 0.0_real64

!     if (mpi_rank == 0) then
!       write(*,*) " "
!       write(*,'(a)') "Loading particles..."
!     end if

!     call load_part( n=int(n), h=h, bcnd=bcnd, np=np_arr, vxp=vxp, &
!                     ntype=int(ntype_call), nmax=int(nmax), kq=kq, ni0=ni0_legacy, np_tot=np_tot, &
!                     nproc=int(nproc), iseed=iseed, sum_dEk=sum_dEk, Nh=Nh, &
!                     mpi_rank=mpi_rank, nproc_mpi=mpi_size )

!     ! Copy into ParticleSet (only tracked)
!     do iproc = 1, nproc
!       do it = 1, ntype_call
!         call part(it,iproc)%from_vxp(vxp(:,:,it,iproc), int(np_tot(it,iproc),int32), int(it,int32))
!       end do
!     end do

!     if (mpi_rank == 0) then
!       do it = 1, ntype_call
!         write(*,'(a,i0,a,i0)') "species ", it, ": ", sum(np_tot(it,1:nproc))
!       end do
!     end if

!     deallocate(vxp, np_arr, np_tot, Nh, sum_dEk, ni0_legacy)

!   end subroutine load_particles_legacy


!   pure integer(int32) function clamp_index(x, dx, nx) result(ix)
!     real(real64), intent(in) :: x, dx
!     integer(int32), intent(in) :: nx
!     integer(int32) :: raw

!     if (dx <= 0.0_real64) then
!       ix = 1_int32
!       return
!     end if

!     raw = int(floor(x/dx), int32) + 1_int32
!     if (raw < 1_int32) raw = 1_int32
!     if (raw > nx)      raw = nx
!     ix = raw
!   end function clamp_index

! end module mod_load_particles_legacy