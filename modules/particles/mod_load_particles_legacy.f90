module mod_load_particles_legacy
  use iso_fortran_env, only: real64, int32
  use mod_config,   only: Config
  use mod_particles, only: ParticleSet
  use mod_legacy_particle_globals, only: np_cell, n_cell, x_load, ymax, zmax, Pabs, &
                                        ixl_pow, ixr_pow, iyl_pow, iyr_pow, izl_pow, izr_pow, &
                                        legacy_globals_init
  implicit none
  private
  public :: load_particles_legacy

  interface
    subroutine load_part(n,h,bcnd,np,vxp,ntype,nmax,kq,ni0,np_tot,nproc,iseed,sum_dEk,Nh,mpi_rank,nproc_mpi)
      use iso_fortran_env, only: real64
      implicit none
      integer :: mpi_rank, nproc_mpi
      integer :: ntype, nmax, n(3), nproc
      integer :: Nh(nproc)
      real(real64) :: h(3)
      integer :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
      integer :: np_tot(ntype,nproc)
      real(real64) :: vxp(6,nmax,ntype,nproc)
      real(real64) :: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)
      real(real64) :: ni0(*)      ! legacy header uses npart; we pass size-1
      real(real64) :: sum_dEk(nproc)
      real(real64) :: kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)
      integer :: iseed(nproc)
    end subroutine load_part
  end interface

contains

  subroutine load_particles_legacy(cfg, mpi_rank, mpi_size, n, h, bcnd, kq, ntype, part)
    type(Config),  intent(in)    :: cfg
    integer,       intent(in)    :: mpi_rank, mpi_size
    integer(int32),intent(in)    :: n(3)
    real(real64),  intent(in)    :: h(3)
    integer,       intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),  intent(in)    :: kq  (0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32),intent(in)    :: ntype
    type(ParticleSet), intent(inout) :: part(ntype, cfg%omp_rank_max)

    ! locals
    integer(int32) :: nproc
    integer(int32) :: nmax, jmax_est
    real(real64)   :: tmp
    integer        :: ip, it

    real(real64), allocatable :: vxp(:,:,:,:)
    real(real64), allocatable :: np_arr(:,:,:,:,:)
    integer,      allocatable :: np_tot(:,:)
    integer,      allocatable :: Nh(:)
    integer,      allocatable :: iseed(:)
    real(real64), allocatable :: sum_dEk(:)

    real(real64) :: ni0(1)

    nproc = int(max(1, cfg%omp_rank_max), int32)

    ! -------------------------------------------------------
    ! 1) set legacy globals expected by load_part
    ! -------------------------------------------------------
    np_cell = int(cfg%np_cell, int32)

    ! count physical cells (using bcnd>=1 as "inside")
    n_cell  = int(count(bcnd(1:n(1),1:n(2),1:n(3)) >= 1), int32)

    x_load  = cfg%x_load
    if (x_load <= 0.0_real64) x_load = real(n(1),real64)*h(1)

    ymax    = real(n(2),real64)*h(2)
    zmax    = real(n(3),real64)*h(3)

    Pabs    = cfg%Pabs

    ! power window indices: convert physical x -> cell index in [1..n]
    ixl_pow = clamp_index(cfg%xl_pow, h(1), n(1))
    ixr_pow = clamp_index(cfg%xr_pow, h(1), n(1))
    iyl_pow = clamp_index(cfg%yl_pow, h(2), n(2))
    iyr_pow = clamp_index(cfg%yr_pow, h(2), n(2))
    izl_pow = clamp_index(cfg%zl_pow, h(3), n(3))
    izr_pow = clamp_index(cfg%zr_pow, h(3), n(3))

    ! allocate vt0(:) and Nm(:) arrays for this ntype
    call legacy_globals_init(ntype)

    ! -------------------------------------------------------
    ! 2) choose nmax (must be >= particles/species/thread)
    ! legacy uses: jmax = NINT(np_cell*n_cell/(mpi_size*nproc))
    ! -------------------------------------------------------
    tmp = real(np_cell,real64) * real(n_cell,real64) / max(1.0_real64, real(mpi_size*nproc,real64))
    jmax_est = max(1_int32, int(nint(tmp), int32))
    nmax = int(1.3_real64*real(jmax_est,real64), int32) + 32_int32

    ! -------------------------------------------------------
    ! 3) allocate legacy arrays for the loader
    ! -------------------------------------------------------
    allocate(vxp(6, nmax, ntype, nproc))
    allocate(np_arr(0:n(1)+2,0:n(2)+2,0:n(3)+2, ntype, nproc))
    allocate(np_tot(ntype, nproc))
    allocate(Nh(nproc), iseed(nproc))
    allocate(sum_dEk(nproc))

    vxp     = 0.0_real64
    np_arr  = 0.0_real64
    np_tot  = 0
    Nh      = 0
    iseed   = 0
    sum_dEk = 0.0_real64

    do ip = 1, nproc
      iseed(ip) = 13579 + 97*ip + 1000*max(0, mpi_rank)
    end do

    ni0 = 0.0_real64

    ! -------------------------------------------------------
    ! 4) call legacy loader
    ! -------------------------------------------------------
    call load_part( n=int(n), h=h, bcnd=bcnd, np=np_arr, vxp=vxp, &
                    ntype=int(ntype), nmax=int(nmax), kq=kq, ni0=ni0, np_tot=np_tot, &
                    nproc=int(nproc), iseed=iseed, sum_dEk=sum_dEk, Nh=Nh, &
                    mpi_rank=mpi_rank, nproc_mpi=mpi_size )

    ! -------------------------------------------------------
    ! 5) copy vxp -> modern ParticleSet
    ! -------------------------------------------------------
    do ip = 1, nproc
      do it = 1, ntype
        call part(it,ip)%from_vxp(vxp(:,:,it,ip), int(np_tot(it,ip),int32), int(it,int32))
      end do
    end do

    if (mpi_rank == 0) then
      write(*,'(a,i0,a,i0,a,i0)') "load_particles_legacy: ntype=", ntype, &
                                  " nproc=", nproc, " nmax=", nmax
    end if

    deallocate(vxp, np_arr, np_tot, Nh, iseed, sum_dEk)

  end subroutine load_particles_legacy


  pure integer(int32) function clamp_index(x, dx, nx) result(ix)
    real(real64), intent(in) :: x, dx
    integer(int32), intent(in) :: nx
    integer(int32) :: raw
    if (dx <= 0.0_real64) then
      ix = 1_int32
      return
    end if
    raw = int(floor(x/dx), int32) + 1_int32
    if (raw < 1_int32) raw = 1_int32
    if (raw > nx) raw = nx
    ix = raw
  end function clamp_index

end module mod_load_particles_legacy