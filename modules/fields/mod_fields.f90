module mod_fields
  use iso_fortran_env, only: real64, int32
  use mod_domain,      only: Domain
  implicit none
  private
  public :: Fields

  type :: Fields
    real(real64), allocatable :: phi(:,:,:)
    real(real64), allocatable :: E(:,:,:,:)
    real(real64), allocatable :: kq(:,:,:)
    real(real64), allocatable :: rho(:,:,:)
    real(real64), allocatable :: np(:,:,:,:)
  contains
    procedure :: allocate_from_domain
    procedure :: allocate_species_density
    procedure :: zero
    procedure :: destroy
  end type Fields

contains

  subroutine allocate_from_domain(self, dom)
    class(Fields), intent(inout) :: self
    type(Domain),  intent(in)    :: dom
    integer :: nx, ny, nz

    nx = dom%n(1)
    ny = dom%n(2)
    nz = dom%n(3)

    if (.not. allocated(self%phi)) allocate(self%phi(0:nx+2, 0:ny+2, 0:nz+2))
    if (.not. allocated(self%E))   allocate(self%E(3, 0:nx+2, 0:ny+2, 0:nz+2))
    if (.not. allocated(self%kq))  allocate(self%kq(0:nx+2, 0:ny+2, 0:nz+2))
    if (.not. allocated(self%rho)) allocate(self%rho(0:nx+2, 0:ny+2, 0:nz+2))
  end subroutine allocate_from_domain

  subroutine allocate_species_density(self, dom, ntype)
    class(Fields), intent(inout) :: self
    type(Domain),  intent(in)    :: dom
    integer(int32), intent(in)   :: ntype
    integer :: nx, ny, nz

    nx = dom%n(1)
    ny = dom%n(2)
    nz = dom%n(3)

    if (allocated(self%np)) deallocate(self%np)
    allocate(self%np(0:nx+2, 0:ny+2, 0:nz+2, ntype))
    self%np = 0.0_real64
  end subroutine allocate_species_density

  subroutine zero(self)
    class(Fields), intent(inout) :: self
    if (allocated(self%phi)) self%phi = 0.0_real64
    if (allocated(self%E))   self%E   = 0.0_real64
    if (allocated(self%kq))  self%kq  = 0.0_real64
    if (allocated(self%rho)) self%rho = 0.0_real64
    if (allocated(self%np))  self%np  = 0.0_real64
  end subroutine zero

  subroutine destroy(self)
    class(Fields), intent(inout) :: self
    if (allocated(self%phi)) deallocate(self%phi)
    if (allocated(self%E))   deallocate(self%E)
    if (allocated(self%kq))  deallocate(self%kq)
    if (allocated(self%rho)) deallocate(self%rho)
    if (allocated(self%np))  deallocate(self%np)
  end subroutine destroy

end module mod_fields