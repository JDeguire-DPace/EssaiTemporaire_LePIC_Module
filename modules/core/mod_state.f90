module mod_state
  use iso_fortran_env, only: int32, real64
  use mpi

  use mod_config,         only: Config
  use mod_readConditions, only: read_input

  use mod_domain,         only: Domain
  use mod_boundary,       only: build_boundary

  use mod_fields,         only: Fields
  use mod_poisson_decomp, only: PoissonDecomp
  use mod_charge_weights, only: build_kq
  use mod_debug_checks,   only: checkpoint_poisson_decomp, checkpoint_kq

  use mod_chemistryState, only: ChemistryState
  use mod_reactionsDB,    only: ReactionsDB
  use mod_part_info,      only: npart

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

    ! Read config
    call read_input(self%cfg, self%mpi_rank)

    ! Domain setup
    call self%dom%init_from_config(self%cfg)
    call self%dom%allocate_masks_domain()   ! your current domain allocation

    ! Fields setup
    call self%fld%allocate_from_domain(self%dom)
    call self%fld%zero()

    ! Chemistry / reactions
    call self%init_chemistry()
  end subroutine init


  subroutine build_boundary_only(self)
    class(State), intent(inout) :: self
    integer :: ierr
    real(real64), allocatable :: phi_rt(:,:,:)
    real(real64) :: lmax, gmax

    if (self%mpi_rank == 0) write(*,*) "Building boundary..."
    call build_boundary(self%dom, self%cfg, self%fld, self%mpi_rank)

    ! ------------------------------------------------------------
    ! Poisson Z-slab decomposition (legacy structure)
    ! ------------------------------------------------------------
    write(*,*) "Initializing Poisson decomposition..."
    call self%pdec%init(int(self%dom%n(1),int32), int(self%dom%n(2),int32), int(self%dom%n(3),int32), self%comm)
    
    write(*,*) "Scattering phi and bcnd to local Poisson decomposition arrays..." 
    call self%pdec%scatter_from_global(self%fld%phi, self%dom%bcnd)


    write(*,*) "Checking Poisson decomposition and local arrays..."
    call checkpoint_poisson_decomp(self%mpi_rank, self%comm, &
                                   int(self%dom%n(3),int32), self%pdec%k0, self%pdec%m, &
                                   self%pdec%phi_dom, self%pdec%bcnd_dom)

    ! ------------------------------------------------------------
    ! Charge weights kq (global, exact legacy logic)
    ! ------------------------------------------------------------
    write(*,*) "Building kq charge weights for Poisson solver..."
    call build_kq(self%dom%bcnd, self%fld%kq)
    write(*,*) "Checking kq charge weights for Poisson solver..."
    call checkpoint_kq(self%fld%kq, self%comm, self%mpi_rank, "after build_kq", self%pdec)
    write(*,*) "Building kq charge weights for Poisson solver... DONE"

    ! ------------------------------------------------------------
    ! Tiny checkpoint: phi round-trip (scatter -> gather)
    ! ------------------------------------------------------------
    allocate(phi_rt(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2))
    phi_rt = -999.0_real64

    call self%pdec%gather_phi_to_global(phi_rt)

    lmax = maxval(abs(phi_rt - self%fld%phi))
    call MPI_Allreduce(lmax, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, self%comm, ierr)

    if (self%mpi_rank == 0) write(*,'(a,1p,e12.3)') "PHI round-trip max|diff| = ", gmax
    deallocate(phi_rt)

    ! ------------------------------------------------------------
    ! Print domain summary (rank 0)
    ! ------------------------------------------------------------
    if (self%mpi_rank == 0) then
      write(*,'(a,3(i0,1x))')         "n = ", self%dom%n(1), self%dom%n(2), self%dom%n(3)
      write(*,'(a,3(1p,e12.4,1x))')   "h(m) = ", self%dom%h(1), self%dom%h(2), self%dom%h(3)
      write(*,'(a,3(1p,e12.4,1x))')   "box(m) = ", self%dom%xmax, self%dom%ymax, self%dom%zmax
      write(*,'(a,1p,e12.4)')         "Sg(m^2) = ", self%dom%Sg
      write(*,'(a,4(i0,1x))')         "flags(pbc,pbcz,nmn,die) = ", &
           self%dom%flag_pbc, self%dom%flag_pbcz, self%dom%flag_nmn, self%dom%flag_die
    end if
  end subroutine build_boundary_only


  subroutine init_chemistry(self)
    class(State), intent(inout) :: self

    if (self%mpi_rank == 0) write(*,*) "Initializing chemistry/reactions..."

    ! 1) init chemistry container (species properties)
    call self%chem%init(npart, self%rxn%ncol_mx, self%rxn%npt_mx)

    ! legacy reader expects ngas available via ChemistryState
    self%chem%ngas = self%cfg%ngas

    ! 2) load reaction tables into ReactionsDB
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
    ! later: destroy domain arrays / close IO if needed
  end subroutine finalize

end module mod_state
