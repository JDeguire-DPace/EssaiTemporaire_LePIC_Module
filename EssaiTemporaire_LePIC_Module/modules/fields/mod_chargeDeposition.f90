module mod_chargeDeposition
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  implicit none
  private

  public :: clear_np_thread
  public :: deposit_particle_set_to_np_thread

  ! Keep this .false. for production speed.
  ! Switch to .true. only while debugging bad particle positions.
  logical, parameter :: DEPOSITION_SAFETY_CHECKS = .false.

contains

  subroutine clear_np_thread(n, ntype, nproc, np_thread)
    integer(int32), intent(in)    :: n(3)
    integer(int32), intent(in)    :: ntype, nproc
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)

    np_thread = 0.0_real64
  end subroutine clear_np_thread


  subroutine deposit_particle_set_to_np_thread(part, n, h, kq, Nm_species, np_local)
    ! Fast production version of charge deposition.
    ! Same deposition convention as the previous modular code, but the
    ! expensive debug/error guards are compile-time disabled by default.
    class(ParticleSet), intent(in)    :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(in)    :: Nm_species
    real(real64),       intent(inout) :: np_local(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    integer(int32) :: i
    integer(int32) :: ix, iy, iz
    real(real64)   :: x, y, z
    real(real64)   :: px, py, pz
    real(real64)   :: wx2, wy2, wz2
    real(real64)   :: k1
    real(real64)   :: w1, w2, w3, w4, w5, w6, w7, w8
    real(real64)   :: xmax_loc, ymax_loc, zmax_loc
    logical        :: has_dead

    if (.not. allocated(part%x)) return
    if (part%n <= 0_int32) return

    xmax_loc = real(n(1), real64) * h(1)
    ymax_loc = real(n(2), real64) * h(2)
    zmax_loc = real(n(3), real64) * h(3)

    k1 = Nm_species / (h(1) * h(2) * h(3))
    has_dead = allocated(part%flag_dead)

    do i = 1, part%n
      if (has_dead) then
        if (part%flag_dead(i) /= 0) cycle
      end if

      x = part%x(i)
      y = part%y(i)
      z = part%z(i)

      ! Cheap guard kept for particles slightly outside due to roundoff.
      ! This matches the previous behavior where tiny excursions were skipped.
      if (x < -1.0e-12_real64 .or. x > xmax_loc + 1.0e-12_real64) cycle
      if (y < -1.0e-12_real64 .or. y > ymax_loc + 1.0e-12_real64) cycle
      if (z < -1.0e-12_real64 .or. z > zmax_loc + 1.0e-12_real64) cycle

      if (DEPOSITION_SAFETY_CHECKS) then
        if (.not. (x == x .and. y == y .and. z == z)) then
          write(*,*) 'NaN particle position in deposition'
          write(*,*) 'i = ', i
          write(*,*) 'x,y,z = ', x, y, z
          error stop 'deposit_particle_set_to_np_thread: NaN position'
        end if

        if (x < 0.0_real64 .or. x > xmax_loc .or. &
            y < 0.0_real64 .or. y > ymax_loc .or. &
            z < 0.0_real64 .or. z > zmax_loc) then
          write(*,*) 'Out-of-range particle position in deposition'
          write(*,*) 'i = ', i
          write(*,*) 'x,y,z = ', x, y, z
          write(*,*) 'xmax,ymax,zmax = ', xmax_loc, ymax_loc, zmax_loc
          error stop 'deposit_particle_set_to_np_thread: particle out of bounds'
        end if
      end if

      ix = int(x / h(1), int32) + 1_int32
      iy = int(y / h(2), int32) + 1_int32
      iz = int(z / h(3), int32) + 1_int32

      if (DEPOSITION_SAFETY_CHECKS) then
        if (ix < 1_int32 .or. ix > n(1)+1_int32 .or. &
            iy < 1_int32 .or. iy > n(2)+1_int32 .or. &
            iz < 1_int32 .or. iz > n(3)+1_int32) then
          write(*,*) 'Bad deposition cell index'
          write(*,*) 'i = ', i
          write(*,*) 'ix,iy,iz = ', ix, iy, iz
          write(*,*) 'x,y,z = ', x, y, z
          write(*,*) 'n = ', n
          error stop 'deposit_particle_set_to_np_thread: invalid cell index'
        end if
      end if

      px = (real(ix, real64) * h(1) - x) / h(1)
      py = (real(iy, real64) * h(2) - y) / h(2)
      pz = (real(iz, real64) * h(3) - z) / h(3)

      wx2 = 1.0_real64 - px
      wy2 = 1.0_real64 - py
      wz2 = 1.0_real64 - pz

      w1 = k1 * px  * py  * pz
      w2 = k1 * wx2 * py  * pz
      w3 = k1 * wx2 * wy2 * pz
      w4 = k1 * px  * wy2 * pz
      w5 = k1 * px  * py  * wz2
      w6 = k1 * wx2 * py  * wz2
      w7 = k1 * wx2 * wy2 * wz2
      w8 = k1 * px  * wy2 * wz2

      np_local(ix  ,iy  ,iz  ) = np_local(ix  ,iy  ,iz  ) + kq(ix  ,iy  ,iz  ) * w1
      np_local(ix+1,iy  ,iz  ) = np_local(ix+1,iy  ,iz  ) + kq(ix+1,iy  ,iz  ) * w2
      np_local(ix+1,iy+1,iz  ) = np_local(ix+1,iy+1,iz  ) + kq(ix+1,iy+1,iz  ) * w3
      np_local(ix  ,iy+1,iz  ) = np_local(ix  ,iy+1,iz  ) + kq(ix  ,iy+1,iz  ) * w4
      np_local(ix  ,iy  ,iz+1) = np_local(ix  ,iy  ,iz+1) + kq(ix  ,iy  ,iz+1) * w5
      np_local(ix+1,iy  ,iz+1) = np_local(ix+1,iy  ,iz+1) + kq(ix+1,iy  ,iz+1) * w6
      np_local(ix+1,iy+1,iz+1) = np_local(ix+1,iy+1,iz+1) + kq(ix+1,iy+1,iz+1) * w7
      np_local(ix  ,iy+1,iz+1) = np_local(ix  ,iy+1,iz+1) + kq(ix  ,iy+1,iz+1) * w8
    end do
  end subroutine deposit_particle_set_to_np_thread

end module mod_chargeDeposition
