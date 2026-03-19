module mod_particles
  use iso_fortran_env, only: real64, int32
  implicit none
  private
  public :: ParticleSet

  type :: ParticleSet
    integer(int32) :: n    = 0
    integer(int32) :: nmax = 0

    ! Particle data
    real(real64), allocatable :: x(:), y(:), z(:)
    real(real64), allocatable :: vx(:), vy(:), vz(:)
    real(real64), allocatable :: w(:)
    integer(int32), allocatable :: sp(:)

    ! Sorting / binning metadata
    integer(int32), allocatable :: cell_id(:)     ! size nmax
    integer(int32), allocatable :: cell_count(:)  ! size ncells
    integer(int32), allocatable :: cell_start(:)  ! size ncells
    integer(int32) :: ncells = 0

  contains
    procedure :: allocate_pset
    procedure :: ensure_capacity
    procedure :: ensure_cell_storage
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
    self%n    = 0_int32

    if (self%nmax > 0_int32) then
      allocate(self%x(self%nmax), self%y(self%nmax), self%z(self%nmax))
      allocate(self%vx(self%nmax), self%vy(self%nmax), self%vz(self%nmax))
      allocate(self%w(self%nmax))
      allocate(self%sp(self%nmax))
      allocate(self%cell_id(self%nmax))

      call self%clear()
    end if
  end subroutine allocate_pset


  subroutine ensure_capacity(self, nneed)
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: nneed

    integer(int32) :: newmax, ncopy
    real(real64), allocatable :: tx(:), ty(:), tz(:)
    real(real64), allocatable :: tvx(:), tvy(:), tvz(:), tw(:)
    integer(int32), allocatable :: tsp(:), tcell(:)

    if (nneed <= 0_int32) return

    if (.not. allocated(self%x)) then
      newmax = max(32_int32, nneed)
      call self%allocate_pset(newmax)
      return
    end if

    if (self%nmax >= nneed) return

    newmax = max(nneed, int(1.3_real64 * real(self%nmax, real64), int32) + 32_int32)
    ncopy  = self%n

    allocate(tx(newmax), ty(newmax), tz(newmax))
    allocate(tvx(newmax), tvy(newmax), tvz(newmax))
    allocate(tw(newmax))
    allocate(tsp(newmax))
    allocate(tcell(newmax))

    tx    = 0.0_real64 ; ty    = 0.0_real64 ; tz    = 0.0_real64
    tvx   = 0.0_real64 ; tvy   = 0.0_real64 ; tvz   = 0.0_real64
    tw    = 0.0_real64
    tsp   = 0_int32
    tcell = 0_int32

    if (ncopy > 0_int32) then
      tx(1:ncopy)    = self%x(1:ncopy)
      ty(1:ncopy)    = self%y(1:ncopy)
      tz(1:ncopy)    = self%z(1:ncopy)
      tvx(1:ncopy)   = self%vx(1:ncopy)
      tvy(1:ncopy)   = self%vy(1:ncopy)
      tvz(1:ncopy)   = self%vz(1:ncopy)
      tw(1:ncopy)    = self%w(1:ncopy)
      tsp(1:ncopy)   = self%sp(1:ncopy)
      if (allocated(self%cell_id)) tcell(1:ncopy) = self%cell_id(1:ncopy)
    end if

    call move_alloc(tx,    self%x)
    call move_alloc(ty,    self%y)
    call move_alloc(tz,    self%z)
    call move_alloc(tvx,   self%vx)
    call move_alloc(tvy,   self%vy)
    call move_alloc(tvz,   self%vz)
    call move_alloc(tw,    self%w)
    call move_alloc(tsp,   self%sp)
    call move_alloc(tcell, self%cell_id)

    self%nmax = newmax
  end subroutine ensure_capacity


  subroutine ensure_cell_storage(self, ncells_in)
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: ncells_in

    if (ncells_in <= 0_int32) return

    if (.not. allocated(self%cell_count)) then
      allocate(self%cell_count(ncells_in), self%cell_start(ncells_in))
      self%ncells = ncells_in
    else if (self%ncells /= ncells_in) then
      deallocate(self%cell_count, self%cell_start)
      allocate(self%cell_count(ncells_in), self%cell_start(ncells_in))
      self%ncells = ncells_in
    end if

    self%cell_count = 0_int32
    self%cell_start = 0_int32
  end subroutine ensure_cell_storage


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

    if (allocated(self%cell_id)) self%cell_id(1:self%n) = 0_int32
  end subroutine from_vxp


  subroutine clear(self)
    class(ParticleSet), intent(inout) :: self

    if (allocated(self%x))         self%x         = 0.0_real64
    if (allocated(self%y))         self%y         = 0.0_real64
    if (allocated(self%z))         self%z         = 0.0_real64
    if (allocated(self%vx))        self%vx        = 0.0_real64
    if (allocated(self%vy))        self%vy        = 0.0_real64
    if (allocated(self%vz))        self%vz        = 0.0_real64
    if (allocated(self%w))         self%w         = 0.0_real64
    if (allocated(self%sp))        self%sp        = 0_int32
    if (allocated(self%cell_id))   self%cell_id   = 0_int32
    if (allocated(self%cell_count)) self%cell_count = 0_int32
    if (allocated(self%cell_start)) self%cell_start = 0_int32

    self%n = 0_int32
  end subroutine clear


  subroutine destroy(self)
    class(ParticleSet), intent(inout) :: self

    if (allocated(self%x))          deallocate(self%x)
    if (allocated(self%y))          deallocate(self%y)
    if (allocated(self%z))          deallocate(self%z)
    if (allocated(self%vx))         deallocate(self%vx)
    if (allocated(self%vy))         deallocate(self%vy)
    if (allocated(self%vz))         deallocate(self%vz)
    if (allocated(self%w))          deallocate(self%w)
    if (allocated(self%sp))         deallocate(self%sp)
    if (allocated(self%cell_id))    deallocate(self%cell_id)
    if (allocated(self%cell_count)) deallocate(self%cell_count)
    if (allocated(self%cell_start)) deallocate(self%cell_start)

    self%n      = 0_int32
    self%nmax   = 0_int32
    self%ncells = 0_int32
  end subroutine destroy

end module mod_particles