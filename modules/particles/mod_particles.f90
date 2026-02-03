module mod_particles
  use iso_fortran_env, only: real64, int32
  implicit none
  private
  public :: ParticleSet

  type :: ParticleSet
    integer(int32) :: n = 0
    integer(int32) :: nmax = 0
    real(real64), allocatable :: x(:), y(:), z(:)
    real(real64), allocatable :: vx(:), vy(:), vz(:)
    real(real64), allocatable :: w(:)
    integer(int32), allocatable :: sp(:)
  contains
    procedure :: allocate_pset
    procedure :: clear
  end type ParticleSet

contains
  subroutine allocate_pset(self, nmax)
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: nmax
    self%nmax = nmax
    self%n    = 0
    allocate(self%x(nmax), self%y(nmax), self%z(nmax))
    allocate(self%vx(nmax), self%vy(nmax), self%vz(nmax))
    allocate(self%w(nmax))
    allocate(self%sp(nmax))
    call self%clear()
  end subroutine

  subroutine clear(self)
    class(ParticleSet), intent(inout) :: self
    if (allocated(self%x))  self%x  = 0.0_real64
    if (allocated(self%vx)) self%vx = 0.0_real64
    if (allocated(self%w))  self%w  = 0.0_real64
    if (allocated(self%sp)) self%sp = 0
    self%n = 0
  end subroutine
end module mod_particles
