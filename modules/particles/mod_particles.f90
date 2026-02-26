module mod_particles
  use iso_fortran_env, only: real64, int32
  implicit none
  private
  public :: ParticleSet

  type :: ParticleSet
    integer(int32) :: n    = 0
    integer(int32) :: nmax = 0
    real(real64), allocatable :: x(:), y(:), z(:)
    real(real64), allocatable :: vx(:), vy(:), vz(:)
    real(real64), allocatable :: w(:)
    integer(int32), allocatable :: sp(:)
  contains
    procedure :: allocate_pset
    procedure :: ensure_capacity
    procedure :: from_vxp
    procedure :: clear
    procedure :: destroy
  end type ParticleSet

contains

  subroutine allocate_pset(self, nmax)
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: nmax

    call self%destroy()
    self%nmax = max(0_int32, nmax)
    self%n    = 0

    if (self%nmax > 0_int32) then
      allocate(self%x(self%nmax), self%y(self%nmax), self%z(self%nmax))
      allocate(self%vx(self%nmax), self%vy(self%nmax), self%vz(self%nmax))
      allocate(self%w(self%nmax))
      allocate(self%sp(self%nmax))
      call self%clear()
    end if
  end subroutine allocate_pset


  subroutine ensure_capacity(self, nneed)
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: nneed
    integer(int32) :: newmax

    if (nneed <= 0_int32) return

    if (.not. allocated(self%x)) then
      newmax = max(32_int32, nneed)
      call self%allocate_pset(newmax)
      return
    end if

    if (self%nmax < nneed) then
      ! growth factor to reduce realloc frequency
      newmax = max(nneed, int(1.3_real64*real(self%nmax,real64), int32) + 32_int32)
      call self%allocate_pset(newmax)
    end if
  end subroutine ensure_capacity


  subroutine from_vxp(self, vxp, np_this, ptype)
    class(ParticleSet), intent(inout) :: self
    real(real64),       intent(in)    :: vxp(:,:)
    integer(int32),     intent(in)    :: np_this
    integer(int32),     intent(in)    :: ptype
    integer(int32) :: i

    if (size(vxp,1) /= 6) then
      error stop "ParticleSet%from_vxp: vxp first dimension must be 6"
    end if

    if (np_this <= 0_int32) then
      self%n = 0_int32
      return
    end if

    call self%ensure_capacity(np_this)

    self%n = np_this
    do i = 1, np_this
      self%x(i)  = vxp(1,i)
      self%y(i)  = vxp(2,i)
      self%z(i)  = vxp(3,i)
      self%vx(i) = vxp(4,i)
      self%vy(i) = vxp(5,i)
      self%vz(i) = vxp(6,i)
      self%w(i)  = 1.0_real64
      self%sp(i) = ptype
    end do
  end subroutine from_vxp


  subroutine clear(self)
    class(ParticleSet), intent(inout) :: self
    if (allocated(self%x))  self%x  = 0.0_real64
    if (allocated(self%y))  self%y  = 0.0_real64
    if (allocated(self%z))  self%z  = 0.0_real64
    if (allocated(self%vx)) self%vx = 0.0_real64
    if (allocated(self%vy)) self%vy = 0.0_real64
    if (allocated(self%vz)) self%vz = 0.0_real64
    if (allocated(self%w))  self%w  = 0.0_real64
    if (allocated(self%sp)) self%sp = 0
    self%n = 0_int32
  end subroutine clear


  subroutine destroy(self)
    class(ParticleSet), intent(inout) :: self
    if (allocated(self%x))  deallocate(self%x)
    if (allocated(self%y))  deallocate(self%y)
    if (allocated(self%z))  deallocate(self%z)
    if (allocated(self%vx)) deallocate(self%vx)
    if (allocated(self%vy)) deallocate(self%vy)
    if (allocated(self%vz)) deallocate(self%vz)
    if (allocated(self%w))  deallocate(self%w)
    if (allocated(self%sp)) deallocate(self%sp)
    self%n = 0_int32
    self%nmax = 0_int32
  end subroutine destroy

end module mod_particles