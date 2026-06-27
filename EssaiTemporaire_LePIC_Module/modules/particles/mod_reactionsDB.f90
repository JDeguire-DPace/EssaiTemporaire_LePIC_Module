module mod_reactionsDB
  use iso_fortran_env, only: real64
  use mod_part_info,      only: npart                 ! legacy PARAM sizing (for now)
  use mod_chemistryState, only: ChemistryState
  use mod_reactions,      only: read_reactions
  implicit none
  private
  public :: ReactionsDB

  type :: ReactionsDB
    ! Keep the legacy “max sizes” here (NOT in State)
    integer :: ncol_mx = 200
    integer :: npt_mx  = 10000

    ! Cross sections / reaction bookkeeping
    real(real64), allocatable :: sig(:,:)       ! (npt_mx, ncol_mx)
    real(real64), allocatable :: sig_Er(:)      ! (npt_mx)
    integer,      allocatable :: sig_list(:,:)  ! (npart, ncol_mx)
    real(real64), allocatable :: sig_Eex(:,:)   ! (ncol_mx, 2)
    integer,      allocatable :: sig_type(:)    ! (ncol_mx)
    integer,      allocatable :: col_info(:,:)  ! (ncol_mx, 10)

    ! Surface reactions
    integer,      allocatable :: scol_rank(:,:) ! (npart, ncol_mx)
    real(real64), allocatable :: scol_info(:,:) ! (ncol_mx, 4)

    ! Null-collision max (sig*v)_max
    real(real64), allocatable :: sigv_mx(:,:)   ! (npart, ncol_mx)

    ! Initial densities read from chemistry file
    real(real64), allocatable :: ni0(:)         ! (npart)

    ! Useful summary numbers (mirrors legacy)
    integer :: ntype = 0
    integer :: n_neu = 0

  contains
    procedure :: allocate_tables => rxn_allocate_tables
    procedure :: load            => rxn_load
    procedure :: destroy         => rxn_destroy
  end type ReactionsDB

contains

  subroutine rxn_allocate_tables(self)
    class(ReactionsDB), intent(inout) :: self

    if (.not. allocated(self%sig)) then
      allocate(self%sig(self%npt_mx, self%ncol_mx))
      allocate(self%sig_Er(self%npt_mx))
      allocate(self%sig_list(npart, self%ncol_mx))
      allocate(self%sig_Eex(self%ncol_mx, 2))
      allocate(self%sig_type(self%ncol_mx))
      allocate(self%col_info(self%ncol_mx, 10))

      allocate(self%scol_rank(npart, self%ncol_mx))
      allocate(self%scol_info(self%ncol_mx, 4))

      allocate(self%sigv_mx(npart, self%ncol_mx))
      allocate(self%ni0(npart))
    end if
  end subroutine rxn_allocate_tables


  subroutine rxn_load(self, chem, rname, mpi_rank)
    class(ReactionsDB),  intent(inout) :: self
    type(ChemistryState),intent(inout) :: chem
    character(len=*),    intent(in)    :: rname
    integer,             intent(in)    :: mpi_rank

    call self%allocate_tables()

    call read_reactions(chem, self%sig, self%sig_Er, self%sig_list, self%sig_Eex, &
                        self%ncol_mx, self%sig_type, self%npt_mx, trim(rname),     &
                        self%col_info, self%scol_rank, self%scol_info, self%sigv_mx,&
                        self%ni0, self%ntype, self%n_neu, mpi_rank)
  end subroutine rxn_load


  subroutine rxn_destroy(self)
    class(ReactionsDB), intent(inout) :: self

    if (allocated(self%sig))       deallocate(self%sig)
    if (allocated(self%sig_Er))    deallocate(self%sig_Er)
    if (allocated(self%sig_list))  deallocate(self%sig_list)
    if (allocated(self%sig_Eex))   deallocate(self%sig_Eex)
    if (allocated(self%sig_type))  deallocate(self%sig_type)
    if (allocated(self%col_info))  deallocate(self%col_info)
    if (allocated(self%scol_rank)) deallocate(self%scol_rank)
    if (allocated(self%scol_info)) deallocate(self%scol_info)
    if (allocated(self%sigv_mx))   deallocate(self%sigv_mx)
    if (allocated(self%ni0))       deallocate(self%ni0)

    self%ntype = 0
    self%n_neu = 0
  end subroutine rxn_destroy

end module mod_reactionsDB
