module mod_chargeDeposition
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  implicit none
  private

  public :: clear_np_thread
  public :: deposit_particle_set_to_np_thread

contains
  !> hello this is a test
  !! @param[in] i am a test
  subroutine clear_np_thread(n, ntype, nproc, np_thread)
    integer(int32), intent(in)    :: n(3)
    integer(int32), intent(in)    :: ntype, nproc
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)

    np_thread = 0.0_real64
  end subroutine clear_np_thread


  subroutine deposit_particle_set_to_np_thread(part, n, h, kq, Nm_species, np_local)
    !===========================================================
    ! Modern version of legacy charge_deposition for one local
    ! particle bucket.
    !
    ! Legacy behavior preserved:
    !   - ix = INT(x/hx) + 1
    !   - trilinear weighting to the 8 surrounding nodes
    !   - scaling by Nm / (hx*hy*hz)
    !   - multiplication by kq(node)
    !   - accumulation into one local 3D density grid
    !
    ! This deposits species density, not total charge.
    !
    ! Extra safety checks added:
    !   - finite x,y,z
    !   - positions inside physical domain
    !   - valid ix,iy,iz before deposition
    !===========================================================
    class(ParticleSet), intent(in)    :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(in)    :: Nm_species
    real(real64),       intent(inout) :: np_local(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    integer(int32) :: i
    integer(int32) :: ix, iy, iz
    real(real64)   :: xp_new, yp_new, zp_new
    real(real64)   :: px, py, pz
    real(real64)   :: k1
    real(real64)   :: ki(8)
    real(real64)   :: xmax_loc, ymax_loc, zmax_loc

    if (.not. allocated(part%x)) return
    if (part%n <= 0_int32) return

    xmax_loc = real(n(1), real64) * h(1)
    ymax_loc = real(n(2), real64) * h(2)
    zmax_loc = real(n(3), real64) * h(3)

    k1 = Nm_species / (h(1) * h(2) * h(3))

    do i = 1, part%n

      if (allocated(part%flag_dead)) then
        if (part%flag_dead(i) /= 0) cycle
      end if

      xp_new = part%x(i)
      yp_new = part%y(i)
      zp_new = part%z(i)

      ! --------------------------------------------------------
      ! NaN / Inf guard
      ! --------------------------------------------------------
      if (.not. (xp_new == xp_new .and. yp_new == yp_new .and. zp_new == zp_new)) then
        write(*,*) 'NaN particle position in deposition'
        write(*,*) 'i = ', i
        write(*,*) 'x,y,z = ', xp_new, yp_new, zp_new
        if (allocated(part%sp)) write(*,*) 'species id = ', part%sp(i)
        error stop 'deposit_particle_set_to_np_thread: NaN position'
      end if

      if (xp_new < -1e-12_real64 .or. xp_new > xmax_loc + 1e-12_real64) cycle
      if (yp_new < -1e-12_real64 .or. yp_new > ymax_loc + 1e-12_real64) cycle
      if (zp_new < -1e-12_real64 .or. zp_new > zmax_loc + 1e-12_real64) cycle

      ! --------------------------------------------------------
      ! Physical-domain guard
      ! We allow x==xmax etc. because legacy can place particles on
      ! boundaries; INT(x/h)+1 then gives n+1, which is valid for
      ! deposition to nodes up to n+1.
      ! --------------------------------------------------------
      if (xp_new < 0.0_real64 .or. xp_new > xmax_loc .or. &
          yp_new < 0.0_real64 .or. yp_new > ymax_loc .or. &
          zp_new < 0.0_real64 .or. zp_new > zmax_loc) then
        write(*,*) 'Out-of-range particle position in deposition'
        write(*,*) 'i = ', i
        write(*,*) 'x,y,z = ', xp_new, yp_new, zp_new
        write(*,*) 'xmax,ymax,zmax = ', xmax_loc, ymax_loc, zmax_loc
        if (allocated(part%sp)) write(*,*) 'species id = ', part%sp(i)
        error stop 'deposit_particle_set_to_np_thread: particle out of bounds'
      end if

      ix = int(xp_new / h(1), int32) + 1_int32
      iy = int(yp_new / h(2), int32) + 1_int32
      iz = int(zp_new / h(3), int32) + 1_int32

      ! --------------------------------------------------------
      ! Index guard for node-based deposition
      ! Valid ix,iy,iz are 1..n+1 here, because we deposit to
      ! ix and ix+1, and the arrays exist up to n+2.
      ! --------------------------------------------------------
      if (ix < 1_int32 .or. ix > n(1)+1_int32 .or. &
          iy < 1_int32 .or. iy > n(2)+1_int32 .or. &
          iz < 1_int32 .or. iz > n(3)+1_int32) then
        write(*,*) 'Bad deposition cell index'
        write(*,*) 'i = ', i
        write(*,*) 'ix,iy,iz = ', ix, iy, iz
        write(*,*) 'x,y,z = ', xp_new, yp_new, zp_new
        write(*,*) 'n = ', n
        if (allocated(part%sp)) write(*,*) 'species id = ', part%sp(i)
        error stop 'deposit_particle_set_to_np_thread: invalid cell index'
      end if

      px = (real(ix, real64) * h(1) - xp_new) / h(1)
      py = (real(iy, real64) * h(2) - yp_new) / h(2)
      pz = (real(iz, real64) * h(3) - zp_new) / h(3)

      ki(1) = k1 * px * py * pz
      ki(2) = k1 * (1.0_real64 - px) * py * pz
      ki(3) = k1 * (1.0_real64 - px) * (1.0_real64 - py) * pz
      ki(4) = k1 * px * (1.0_real64 - py) * pz
      ki(5) = k1 * px * py * (1.0_real64 - pz)
      ki(6) = k1 * (1.0_real64 - px) * py * (1.0_real64 - pz)
      ki(7) = k1 * (1.0_real64 - px) * (1.0_real64 - py) * (1.0_real64 - pz)
      ki(8) = k1 * px * (1.0_real64 - py) * (1.0_real64 - pz)

      np_local(ix  ,iy  ,iz  ) = np_local(ix  ,iy  ,iz  ) + kq(ix  ,iy  ,iz  ) * ki(1)
      np_local(ix+1,iy  ,iz  ) = np_local(ix+1,iy  ,iz  ) + kq(ix+1,iy  ,iz  ) * ki(2)
      np_local(ix+1,iy+1,iz  ) = np_local(ix+1,iy+1,iz  ) + kq(ix+1,iy+1,iz  ) * ki(3)
      np_local(ix  ,iy+1,iz  ) = np_local(ix  ,iy+1,iz  ) + kq(ix  ,iy+1,iz  ) * ki(4)
      np_local(ix  ,iy  ,iz+1) = np_local(ix  ,iy  ,iz+1) + kq(ix  ,iy  ,iz+1) * ki(5)
      np_local(ix+1,iy  ,iz+1) = np_local(ix+1,iy  ,iz+1) + kq(ix+1,iy  ,iz+1) * ki(6)
      np_local(ix+1,iy+1,iz+1) = np_local(ix+1,iy+1,iz+1) + kq(ix+1,iy+1,iz+1) * ki(7)
      np_local(ix  ,iy+1,iz+1) = np_local(ix  ,iy+1,iz+1) + kq(ix  ,iy+1,iz+1) * ki(8)

    end do

  end subroutine deposit_particle_set_to_np_thread

end module mod_chargeDeposition