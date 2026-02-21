program main_poisson_stub
  use iso_fortran_env, only: int32, real64
  use mpi
  use mod_config,         only: Config
  use mod_readConditions, only: read_input
  use mod_domain,         only: Domain
  use mod_fields,         only: Fields
  use mod_poisson_stub_driver, only: run_poisson_stub
  use mod_boundary,       only: build_boundary
  implicit none

  integer :: ierr, rank
  integer :: comm
  type(Config) :: cfg
  type(Domain) :: dom
  type(Fields) :: fld
  integer :: izl, ix, iy, iz

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)

  ! Read conditions.inp -> cfg
  call read_input(cfg, rank)

  ! Domain + fields
  call dom%init_from_config(cfg)
  call dom%allocate_masks_domain()
  call fld%allocate_from_domain(dom)
  call fld%zero()

  ! Build boundary (fills dom%bcnd, dom%h, and sets wall phi in fld%phi)
  call build_boundary(dom, cfg, fld, rank)

  ! Run stub Poisson solve (scatter/build RHS/optional legacy pdesolver/gather)
  call run_poisson_stub(comm, dom%n, dom%h, dom%bcnd, fld%phi)
  
  if (rank == 0) then
  
        izl = dom%n(3)/2 + 1
        open(12,file='../Output/Output_2D/phi1_xy.mco')
        write(12,*) dom%n(1), dom%n(2)
    102   format(800(ES12.4,1x))
        do iy=dom%n(2)+1,1,-1
          write(12,102) (fld%phi(ix,iy,izl), ix=1,dom%n(1)+1,1)
          
        end do
        close(12)

        open(12,file='../Output/Output_2D/phi1_xz.mco')
        write(12,*) dom%n(1), dom%n(3)
        do iz=dom%n(3)+1,1,-1
          
          write(12,102) (fld%phi(ix,dom%n(2)/2+1,iz), ix=1,dom%n(1)+1,1)
        end do
        close(12)

        ! ix = n(1)/2 + 1
        ! if (inp%flag_grd == 1) ix = outp%ixg
        ! open(12,file='../Output/Output_2D/phi_yz.mco')
        ! write(12,*) n(2)/outp%every, n(3)/outp%every
        ! do iz=n(3)+1,1,-1*outp%every
        !   write(12,102) (u(ix,iy,iz), iy=1,n(2)+1,outp%every)
        ! end do
        ! close(12)
  end if

  call MPI_Finalize(ierr)
end program main_poisson_stub
