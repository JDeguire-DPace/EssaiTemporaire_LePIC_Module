module mod_particles
  use iso_fortran_env, only: real64, int32, int8
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

    ! Collision-related per-particle flags
    integer(int8),  allocatable :: flag_dead(:)   ! size nmax
    integer(int32), allocatable :: flag_cex(:)    ! size nmax

  contains
    procedure :: allocate_pset
    procedure :: ensure_capacity
    procedure :: ensure_cell_storage
    procedure :: from_vxp
    procedure :: clear
    procedure :: destroy
  end type ParticleSet


contains

  subroutine allocate_pset(self, nmax_in)
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: nmax_in

    if (allocated(self%x)) call self%destroy()

    self%n    = 0_int32
    self%nmax = max(0_int32, nmax_in)
    self%ncells = 0_int32

    if (self%nmax <= 0_int32) return

    allocate(self%x(self%nmax), self%y(self%nmax), self%z(self%nmax))
    allocate(self%vx(self%nmax), self%vy(self%nmax), self%vz(self%nmax))
    allocate(self%w(self%nmax))
    allocate(self%sp(self%nmax))

    allocate(self%cell_id(self%nmax))
    allocate(self%flag_dead(self%nmax))
    allocate(self%flag_cex(self%nmax))

    self%x = 0.0_real64
    self%y = 0.0_real64
    self%z = 0.0_real64
    self%vx = 0.0_real64
    self%vy = 0.0_real64
    self%vz = 0.0_real64
    self%w  = 0.0_real64
    self%sp = 0_int32

    self%cell_id   = 0_int32
    self%flag_dead = 0_int8
    self%flag_cex  = 0_int32
  end subroutine allocate_pset

  subroutine ensure_capacity(self, needed)
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: needed

    integer(int32) :: new_nmax, old_nmax, ncopy
    real(real64), allocatable :: x_new(:), y_new(:), z_new(:)
    real(real64), allocatable :: vx_new(:), vy_new(:), vz_new(:)
    real(real64), allocatable :: w_new(:)
    integer(int32), allocatable :: sp_new(:), cell_id_new(:), flag_cex_new(:)
    integer(int8),  allocatable :: flag_dead_new(:)

    if (needed <= self%nmax) return

    old_nmax = self%nmax
    new_nmax = max(needed, max(1_int32, 2_int32*old_nmax))
    ncopy    = self%n

    allocate(x_new(new_nmax), y_new(new_nmax), z_new(new_nmax))
    allocate(vx_new(new_nmax), vy_new(new_nmax), vz_new(new_nmax))
    allocate(w_new(new_nmax))
    allocate(sp_new(new_nmax))
    allocate(cell_id_new(new_nmax))
    allocate(flag_dead_new(new_nmax))
    allocate(flag_cex_new(new_nmax))

    x_new = 0.0_real64
    y_new = 0.0_real64
    z_new = 0.0_real64
    vx_new = 0.0_real64
    vy_new = 0.0_real64
    vz_new = 0.0_real64
    w_new  = 0.0_real64
    sp_new = 0_int32
    cell_id_new   = 0_int32
    flag_dead_new = 0_int8
    flag_cex_new  = 0_int32

    if (old_nmax > 0_int32) then
      if (allocated(self%x))  x_new(1:ncopy) = self%x(1:ncopy)
      if (allocated(self%y))  y_new(1:ncopy) = self%y(1:ncopy)
      if (allocated(self%z))  z_new(1:ncopy) = self%z(1:ncopy)
      if (allocated(self%vx)) vx_new(1:ncopy) = self%vx(1:ncopy)
      if (allocated(self%vy)) vy_new(1:ncopy) = self%vy(1:ncopy)
      if (allocated(self%vz)) vz_new(1:ncopy) = self%vz(1:ncopy)
      if (allocated(self%w))  w_new(1:ncopy) = self%w(1:ncopy)
      if (allocated(self%sp)) sp_new(1:ncopy) = self%sp(1:ncopy)

      if (allocated(self%cell_id))   cell_id_new(1:ncopy)   = self%cell_id(1:ncopy)
      if (allocated(self%flag_dead)) flag_dead_new(1:ncopy) = self%flag_dead(1:ncopy)
      if (allocated(self%flag_cex))  flag_cex_new(1:ncopy)  = self%flag_cex(1:ncopy)

      if (allocated(self%x))         deallocate(self%x)
      if (allocated(self%y))         deallocate(self%y)
      if (allocated(self%z))         deallocate(self%z)
      if (allocated(self%vx))        deallocate(self%vx)
      if (allocated(self%vy))        deallocate(self%vy)
      if (allocated(self%vz))        deallocate(self%vz)
      if (allocated(self%w))         deallocate(self%w)
      if (allocated(self%sp))        deallocate(self%sp)
      if (allocated(self%cell_id))   deallocate(self%cell_id)
      if (allocated(self%flag_dead)) deallocate(self%flag_dead)
      if (allocated(self%flag_cex))  deallocate(self%flag_cex)
    end if

    call move_alloc(x_new, self%x)
    call move_alloc(y_new, self%y)
    call move_alloc(z_new, self%z)
    call move_alloc(vx_new, self%vx)
    call move_alloc(vy_new, self%vy)
    call move_alloc(vz_new, self%vz)
    call move_alloc(w_new, self%w)
    call move_alloc(sp_new, self%sp)
    call move_alloc(cell_id_new, self%cell_id)
    call move_alloc(flag_dead_new, self%flag_dead)
    call move_alloc(flag_cex_new, self%flag_cex)

    self%nmax = new_nmax
  end subroutine ensure_capacity


  subroutine ensure_cell_storage(self, ncells_in)
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: ncells_in

    if (self%ncells == ncells_in) return

    if (allocated(self%cell_count)) deallocate(self%cell_count)
    if (allocated(self%cell_start)) deallocate(self%cell_start)

    self%ncells = max(0_int32, ncells_in)

    if (self%ncells <= 0_int32) return

    allocate(self%cell_count(self%ncells))
    allocate(self%cell_start(self%ncells))

    self%cell_count = 0_int32
    self%cell_start = 0_int32
  end subroutine ensure_cell_storage


  subroutine from_vxp(self, vxp_in, npar, species_id)
    class(ParticleSet), intent(inout) :: self
    real(real64),       intent(in)    :: vxp_in(:,:)
    integer(int32),     intent(in)    :: npar
    integer(int32),     intent(in)    :: species_id

    integer(int32) :: i

    if (size(vxp_in,1) /= 6) then
      error stop 'from_vxp: first dimension of vxp_in must be 6'
    end if

    if (size(vxp_in,2) < npar) then
      error stop 'from_vxp: second dimension of vxp_in is smaller than npar'
    end if

    call self%ensure_capacity(npar)

    self%n = npar

    do i = 1, npar
      self%x(i)  = vxp_in(1,i)
      self%y(i)  = vxp_in(2,i)
      self%z(i)  = vxp_in(3,i)
      self%vx(i) = vxp_in(4,i)
      self%vy(i) = vxp_in(5,i)
      self%vz(i) = vxp_in(6,i)
      self%w(i)  = 1.0_real64
      self%sp(i) = species_id
    end do

    self%cell_id(1:npar)   = 0_int32
    self%flag_dead(1:npar) = 0_int8
    self%flag_cex(1:npar)  = 0_int32
  end subroutine from_vxp

  subroutine clear(self)
    class(ParticleSet), intent(inout) :: self

    self%n = 0_int32

    if (allocated(self%cell_id))   self%cell_id   = 0_int32
    if (allocated(self%flag_dead)) self%flag_dead = 0_int8
    if (allocated(self%flag_cex))  self%flag_cex  = 0_int32

    if (allocated(self%cell_count)) self%cell_count = 0_int32
    if (allocated(self%cell_start)) self%cell_start = 0_int32
  end subroutine clear


  subroutine destroy(self)
    class(ParticleSet), intent(inout) :: self

    if (allocated(self%x))         deallocate(self%x)
    if (allocated(self%y))         deallocate(self%y)
    if (allocated(self%z))         deallocate(self%z)
    if (allocated(self%vx))        deallocate(self%vx)
    if (allocated(self%vy))        deallocate(self%vy)
    if (allocated(self%vz))        deallocate(self%vz)
    if (allocated(self%w))         deallocate(self%w)
    if (allocated(self%sp))        deallocate(self%sp)

    if (allocated(self%cell_id))   deallocate(self%cell_id)
    if (allocated(self%cell_count)) deallocate(self%cell_count)
    if (allocated(self%cell_start)) deallocate(self%cell_start)

    if (allocated(self%flag_dead)) deallocate(self%flag_dead)
    if (allocated(self%flag_cex))  deallocate(self%flag_cex)

    self%n      = 0_int32
    self%nmax   = 0_int32
    self%ncells = 0_int32
  end subroutine destroy

end module mod_particles