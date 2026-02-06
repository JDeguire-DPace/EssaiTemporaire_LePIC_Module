module mod_fields
  use iso_fortran_env, only: real64
  use mod_domain,      only: Domain
  implicit none
  private
  public :: Fields

  type :: Fields
    ! Keep it minimal: phi and E are the first things used everywhere
    real(real64), allocatable :: phi(:,:,:)   ! (0:nx+2,0:ny+2,0:nz+2)
    real(real64), allocatable :: E(:,:,:,:)   ! (3,0:nx+2,0:ny+2,0:nz+2)

  contains
    procedure :: allocate_from_domain
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

    if (.not. allocated(self%phi)) then
      allocate(self%phi(0:nx+2, 0:ny+2, 0:nz+2))
    end if

    if (.not. allocated(self%E)) then
      allocate(self%E(3, 0:nx+2, 0:ny+2, 0:nz+2))
    end if
  end subroutine allocate_from_domain


  subroutine zero(self)
    class(Fields), intent(inout) :: self
    if (allocated(self%phi)) self%phi = 0.0_real64
    if (allocated(self%E))   self%E   = 0.0_real64
  end subroutine zero


  subroutine destroy(self)
    class(Fields), intent(inout) :: self
    if (allocated(self%phi)) deallocate(self%phi)
    if (allocated(self%E))   deallocate(self%E)
  end subroutine destroy

end module mod_fields
