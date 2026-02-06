module mod_chemistryState
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: ChemistryState

  ! ---- indices used by legacy col_info / sig_Eex ----
  integer, parameter, public :: ind_nre = 1
  integer, parameter, public :: ind_nby = 2
  integer, parameter, public :: ind_Eth = 1
  integer, parameter, public :: ind_dE  = 2

  type :: ChemistryState
    ! sizes (legacy)
    integer :: npart   = 0
    integer :: ncol_mx = 0
    integer :: npt_mx  = 0

    ! gas / chemistry bookkeeping
    integer :: ntype = 0
    integer :: ncol  = 0
    integer :: nscol = 0
    integer :: n_neu = 0

    integer :: sig_npt_mx = 0
    integer, allocatable :: p_ncol(:)     ! (npart)
    integer, allocatable :: p_nscol(:)    ! (npart)

    integer :: tag_neg  = 0
    integer :: tag_neu  = 0
    integer :: tag_beam = 0

    ! species properties (size npart)
    character(len=6), allocatable :: pname(:)
    real(real64),      allocatable :: mass(:)
    real(real64),      allocatable :: charge(:)
    real(real64),      allocatable :: Ti(:)

    ! densities (size npart)
    real(real64), allocatable :: ni0(:)

    ! You used ngas in the legacy reader (global). Keep it here for now:
    real(real64) :: ngas = 0.0_real64
  contains
    procedure :: init => chem_init
  end type ChemistryState

contains

  subroutine chem_init(self, npart, ncol_mx, npt_mx)
    class(ChemistryState), intent(inout) :: self
    integer, intent(in) :: npart, ncol_mx, npt_mx

    self%npart   = npart
    self%ncol_mx = ncol_mx
    self%npt_mx  = npt_mx

    allocate(self%pname(npart), self%mass(npart), self%charge(npart), self%Ti(npart))
    allocate(self%ni0(npart))
    allocate(self%p_ncol(npart), self%p_nscol(npart))

    self%pname   = ''
    self%mass    = 0.0_real64
    self%charge  = 0.0_real64
    self%Ti      = 0.0_real64
    self%ni0     = 0.0_real64
    self%p_ncol  = 0
    self%p_nscol = 0

    self%ntype = 0
    self%ncol  = 0
    self%nscol = 0
    self%n_neu = 0
    self%sig_npt_mx = 0
    self%tag_neg  = 0
    self%tag_neu  = 0
    self%tag_beam = 0
    self%ngas = 0.0_real64
  end subroutine chem_init

end module mod_chemistryState
