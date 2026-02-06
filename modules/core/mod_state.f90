module mod_state
  use iso_fortran_env, only: int32, real64
  use mpi

  use mod_config,         only: Config
  use mod_domain,         only: Domain
  use mod_boundary,       only: build_boundary
  use mod_readConditions, only: read_input

  use mod_chemistryState, only: ChemistryState
  use mod_reactionsDB,    only: ReactionsDB
  use mod_part_info,      only: npart

  use mod_fields,         only: Fields
  use mod_poisson_decomp, only: PoissonDecomp
  use mod_charge_weights, only: build_kq

  use mod_debug_checks,   only: checkpoint_poisson_decomp, checkpoint_kq

  implicit none
  private
  public :: State

  type :: State
    type(Config)        :: cfg
    type(Domain)        :: dom
    type(Fields)        :: fld
    type(PoissonDecomp) :: pdec

    type(ChemistryState) :: chem
    type(ReactionsDB)    :: rxn

    integer :: mpi_rank = -1
    integer :: mpi_size = -1
    integer :: comm     = MPI_COMM_NULL

  contains
    procedure :: init
    procedure :: build_boundary_only
    procedure :: init_chemistry
    procedure :: finalize
  end type State

contains

  subroutine init(self, comm_in)
    class(State), intent(inout) :: self
    integer,      intent(in)    :: comm_in
    integer :: ierr

    self%comm = comm_in
    call MPI_Comm_rank(self%comm, self%mpi_rank, ierr)
    call MPI_Comm_size(self%comm, self%mpi_size, ierr)

    ! Read config (rank-aware)
    call read_input(self%cfg, self%mpi_rank)

    ! Domain
    call self%dom%init_from_config(self%cfg)

    ! If you still store masks/geometry arrays on domain:
    call self%dom%allocate_masks_domain()

    ! Fields
    call self%fld%allocate_from_domain(self%dom)
    call self%fld%zero()

    ! Chemistry / reactions
    call self%init_chemistry()
  end subroutine init


  subroutine build_boundary_only(self)
    class(State), intent(inout) :: self
    integer :: ierr
    integer(int32), allocatable :: bcnd_i32(:,:,:)
    real(real64), allocatable   :: phi_before(:,:,:)
    real(real64) :: maxdiff

    if (self%mpi_rank == 0) write(*,*) "Initialization complete. Building boundary conditions..."
    if (self%mpi_rank == 0) write(*,*) "Building boundary..."

    ! build_boundary is expected to fill:
    !   - self%dom%bcnd (and other geometry flags)
    !   - possibly self%fld%phi (or at least it must exist)
    call build_boundary(self%dom, self%cfg, self%fld, self%mpi_rank)

    ! Initialize Poisson Z-slab decomposition
    call self%pdec%init( int(self%dom%n(1), int32), &
                         int(self%dom%n(2), int32), &
                         int(self%dom%n(3), int32), self%comm )

    ! Save phi for a round-trip identity check later
    allocate(phi_before(lbound(self%fld%phi,1):ubound(self%fld%phi,1), &
                        lbound(self%fld%phi,2):ubound(self%fld%phi,2), &
                        lbound(self%fld%phi,3):ubound(self%fld%phi,3)))
    phi_before = self%fld%phi

    ! IMPORTANT: avoid int(dom%bcnd) as a massive temporary on some compilers.
    allocate(bcnd_i32(lbound(self%dom%bcnd,1):ubound(self%dom%bcnd,1), &
                      lbound(self%dom%bcnd,2):ubound(self%dom%bcnd,2), &
                      lbound(self%dom%bcnd,3):ubound(self%dom%bcnd,3)))
    bcnd_i32 = int(self%dom%bcnd, int32)

    ! Scatter global -> local (phi + bcnd)
    call self%pdec%scatter_from_global(self%fld%phi, bcnd_i32)

    ! Tiny checkpoint to ensure slice bounds are valid and values look sane
    call checkpoint_poisson_decomp( self%mpi_rank, self%comm, &
                                    int(self%dom%n(3), int32), self%pdec%k0, self%pdec%m, &
                                    self%pdec%phi_dom, self%pdec%bcnd_dom )

    ! Build kq charge weights (global), legacy logic: 1 in plasma, 2 elsewhere
    if (self%mpi_rank == 0) write(*,*) "Building kq charge weights for Poisson solver..."
    call build_kq(bcnd_i32, self%fld%kq)
    call checkpoint_kq(self%fld%kq, self%comm, self%mpi_rank, "after build_kq")
    if (self%mpi_rank == 0) write(*,*) "Building kq charge weights for Poisson solver... DONE"

    ! Round-trip test: gather local phi_dom back to global phi and compare
    call self%pdec%gather_to_global(self%fld%phi)

    if (self%mpi_rank == 0) then
      maxdiff = maxval(abs(self%fld%phi - phi_before))
      write(*,'(a,1p,e12.3)') "PHI round-trip max|diff| = ", maxdiff
    end if

    deallocate(phi_before)
    deallocate(bcnd_i32)

    if (self%mpi_rank == 0) then
      write(*,'(a,3(i0,1x))')       "n = ", self%dom%n(1), self%dom%n(2), self%dom%n(3)
      write(*,'(a,3(1p,e12.4,1x))') "h(m) = ", self%dom%h(1), self%dom%h(2), self%dom%h(3)
      write(*,'(a,3(1p,e12.4,1x))') "box(m) = ", self%dom%xmax, self%dom%ymax, self%dom%zmax
      write(*,'(a,1p,e12.4)')       "Sg(m^2) = ", self%dom%Sg
      write(*,'(a,4(i0,1x))')       "flags(pbc,pbcz,nmn,die) = ", &
           self%dom%flag_pbc, self%dom%flag_pbcz, self%dom%flag_nmn, self%dom%flag_die
    end if

    call MPI_Barrier(self%comm, ierr)
  end subroutine build_boundary_only


  subroutine init_chemistry(self)
    class(State), intent(inout) :: self

    if (self%mpi_rank == 0) write(*,*) "Initializing chemistry/reactions..."

    call self%chem%init(npart, self%rxn%ncol_mx, self%rxn%npt_mx)
    self%chem%ngas = self%cfg%ngas
    call self%rxn%load(self%chem, trim(self%cfg%rname), self%mpi_rank)

    if (self%mpi_rank == 0) then
      write(*,'(a,i0)') "rxn%ntype = ", self%rxn%ntype
      write(*,'(a,i0)') "rxn%n_neu = ", self%rxn%n_neu
      write(*,'(a,i0)') "chem%ncol = ", self%chem%ncol
      write(*,'(a,i0)') "chem%sig_npt_mx = ", self%chem%sig_npt_mx
    end if
  end subroutine init_chemistry


  subroutine finalize(self)
    class(State), intent(inout) :: self
    call self%pdec%destroy()
    call self%fld%destroy()
    call self%rxn%destroy()
  end subroutine finalize

end module mod_state
