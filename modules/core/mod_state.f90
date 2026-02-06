module mod_state
  use iso_fortran_env, only: real64
  use mpi
  use mod_config,         only: Config
  use mod_domain,         only: Domain
  use mod_boundary,       only: build_boundary
  use mod_readConditions, only: read_input
  use mod_chemistryState, only: ChemistryState
  use mod_reactionsDB,   only: ReactionsDB
  use mod_part_info,      only: npart
  use mod_fields
  
  implicit none
  private
  public :: State

  type :: State
    type(Config) :: cfg
    type(Domain) :: dom
    type(Fields) :: fld

    type(ChemistryState) :: chem
    type(ReactionsDB)    :: rxn   ! <--- moved out of State

    ! --- MPI ---
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
    call self%dom%allocate_masks_domain()
    
    ! Fields setup
    call self%fld%allocate_from_domain(self%dom)
    call self%fld%zero()

    ! Chemistry / reactions
    call self%init_chemistry()
  end subroutine init


  subroutine build_boundary_only(self)
    class(State), intent(inout) :: self

    if (self%mpi_rank == 0) write(*,*) "Building boundary..."
    call build_boundary(self%dom, self%cfg, self%fld ,self%mpi_rank)

    if (self%mpi_rank == 0) then
      write(*,'(a,3(i0,1x))') "n = ", self%dom%n(1), self%dom%n(2), self%dom%n(3)
      write(*,'(a,3(1p,e12.4,1x))') "h(m) = ", self%dom%h(1), self%dom%h(2), self%dom%h(3)
      write(*,'(a,3(1p,e12.4,1x))') "box(m) = ", self%dom%xmax, self%dom%ymax, self%dom%zmax
      write(*,'(a,1p,e12.4)') "Sg(m^2) = ", self%dom%Sg
      write(*,'(a,4(i0,1x))') "flags(pbc,pbcz,nmn,die) = ", &
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
    call self%fld%destroy()
    call self%rxn%destroy()
    ! later: also destroy domain arrays / close IO if needed
  end subroutine finalize

end module mod_state
