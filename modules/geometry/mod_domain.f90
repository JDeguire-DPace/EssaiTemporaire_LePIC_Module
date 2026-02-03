module mod_domain
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: Domain

  type :: Domain
    integer :: n(3) = 0
    real(real64) :: h(3) = 0.0_real64

    ! geometry extents (filled by legacy generate_boundary for now)
    real(real64) :: xmax = 0.0_real64, ymax = 0.0_real64, zmax = 0.0_real64

    ! ---- NEW: boundary/geometry metadata produced by generate_boundary ----
    real(real64) :: Lgy = 0.0_real64, Lgz = 0.0_real64
    real(real64) :: Sg  = 0.0_real64
    real(real64) :: xg1 = 0.0_real64
    integer      :: ixg = 0

    integer :: flag_pbc  = 0
    integer :: flag_pbcz = 0
    integer :: flag_nmn  = 0
    integer :: flag_die  = 0

    ! wall material types per label (0:ngrid). Allocated in build_boundary.
    integer, allocatable :: dtype(:)

    ! core arrays (legacy layout: ghosted 0:n+2)
    integer,      allocatable :: bcnd(:,:,:)
    real(real64), allocatable :: u(:,:,:)

  contains
    procedure :: init_from_config
    procedure :: allocate_fields
  end type Domain

contains

  subroutine init_from_config(self, cfg)
    use mod_config, only: Config
    implicit none
    class(Domain), intent(inout) :: self
    type(Config),  intent(in)    :: cfg

    self%n = cfg%n

    ! Keep h/xmax/ymax/zmax unset here: legacy generate_boundary defines them
    self%h = 0.0_real64
    self%xmax = 0.0_real64
    self%ymax = 0.0_real64
    self%zmax = 0.0_real64
  end subroutine init_from_config

  subroutine allocate_fields(self)
    implicit none
    class(Domain), intent(inout) :: self

    if (any(self%n <= 0)) error stop "Domain%allocate_fields: n not initialized"

    allocate(self%bcnd(0:self%n(1)+2, 0:self%n(2)+2, 0:self%n(3)+2))
    allocate(self%u   (0:self%n(1)+2, 0:self%n(2)+2, 0:self%n(3)+2))

    self%bcnd = 0
    self%u    = 0.0_real64
  end subroutine allocate_fields

end module mod_domain
