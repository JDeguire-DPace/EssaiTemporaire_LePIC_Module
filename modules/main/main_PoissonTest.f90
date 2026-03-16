program main_PoissonTest
  use iso_fortran_env, only: int32, real64
  use mpi

  use mod_config,              only: Config
  use mod_readConditions,      only: read_input
  use mod_constants,           only: eps0
  use mod_domain,              only: Domain
  use mod_fields,              only: Fields
  use mod_boundary,            only: build_boundary
  use mod_density,             only: build_rho_from_np
  use mod_poisson_test_driver, only: run_poisson_actual
  use mod_output_2d,           only: write_scalar_planes

  implicit none

  integer :: ierr, rank, comm
  type(Config) :: cfg
  type(Domain) :: dom
  type(Fields) :: fld

  integer(int32) :: ix_plane, iy_plane, iz_plane
  integer(int32) :: ntype_test
  real(real64), allocatable :: np_test(:,:,:,:)

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)

  call read_input(cfg, rank)

  eps0 = cfg%k_eps0 * 8.854187817d-12
  if (rank == 0 .and. cfg%k_eps0 /= 1.0_real64) then
    write(*,*) 'eps0 HAS BEEN RE-SCALED: k_eps0 = ', cfg%k_eps0
  end if

  call dom%init_from_config(cfg)
  call dom%allocate_masks_domain()

  call fld%allocate_from_domain(dom)
  call fld%zero()

  call build_boundary(dom, cfg, fld, rank)

  ! ------------------------------------------------------------
  ! Simple test rho from a single synthetic charged species
  ! ------------------------------------------------------------
  ntype_test = 1_int32
  call fld%allocate_species_density(dom, ntype_test)

  allocate(np_test(0:dom%n(1)+2, 0:dom%n(2)+2, 0:dom%n(3)+2, ntype_test))
  np_test = 0.0_real64

  ! Example: put a compact density blob at the center node
  np_test(dom%n(1)/2+1, dom%n(2)/2+1, dom%n(3)/2+1, 1) = 1.0e14_real64

  fld%np = np_test

  call build_rho_from_np( &
       n      = int(dom%n, int32), &
       np_red = fld%np, &
       charge = [ -1.602176634d-19 ], &
       ntype  = 1, &
       rho    = fld%rho )

  if (rank == 0) then
    write(*,*) ' '
    write(*,*) 'DEBUG TEST rho'
    write(*,'(a,es16.8)') 'sum      = ', sum(fld%rho(0:dom%n(1)+1,0:dom%n(2)+1,0:dom%n(3)+1))
    write(*,'(a,es16.8)') 'sum(abs) = ', sum(abs(fld%rho(0:dom%n(1)+1,0:dom%n(2)+1,0:dom%n(3)+1)))
    write(*,'(a,es16.8)') 'max      = ', maxval(fld%rho(0:dom%n(1)+1,0:dom%n(2)+1,0:dom%n(3)+1))
    write(*,'(a,es16.8)') 'min      = ', minval(fld%rho(0:dom%n(1)+1,0:dom%n(2)+1,0:dom%n(3)+1))
  end if

  call run_poisson_actual( &
       comm        = comm, &
       n_in        = int(dom%n, int32), &
       h           = dom%h, &
       bcnd_global = dom%bcnd, &
       rhs_global  = fld%rho(0:dom%n(1)+1,0:dom%n(2)+1,0:dom%n(3)+1), &
       phi_global  = fld%phi, &
       ncycl       = 100, &
       eps         = cfg%eps, &
       omega       = cfg%omega, &
       ng          = cfg%ng )

  ix_plane = int(dom%n(1)/2 + 1, int32)
  iy_plane = int(dom%n(2)/2,     int32)   ! legacy phi_xz convention
  iz_plane = int(dom%n(3)/2 + 1, int32)

  if (rank == 0) then
    call write_scalar_planes( &
         f        = fld%phi, &
         n        = int(dom%n, int32), &
         ix_plane = ix_plane, &
         iy_plane = iy_plane, &
         iz_plane = iz_plane, &
         every    = 1_int32, &
         prefix   = '../Output/Output_2D/phi_test' )

    write(*,*) ' '
    write(*,*) 'DEBUG TEST phi'
    write(*,'(a,es16.8)') 'sum      = ', sum(fld%phi)
    write(*,'(a,es16.8)') 'sum(abs) = ', sum(abs(fld%phi))
    write(*,'(a,es16.8)') 'max      = ', maxval(fld%phi)
    write(*,'(a,es16.8)') 'min      = ', minval(fld%phi)
  end if

  if (allocated(np_test)) deallocate(np_test)

  call fld%destroy()

  call MPI_Finalize(ierr)
end program main_PoissonTest
