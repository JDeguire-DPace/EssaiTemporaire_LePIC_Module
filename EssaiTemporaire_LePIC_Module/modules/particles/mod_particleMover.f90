module mod_particleMover
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  implicit none
  private

  public :: move_particles_electrostatic
  public :: interpolate_E_trilinear

contains

  pure integer(int32) function clamp_index(i, ilo, ihi) result(ic)
    integer(int32), intent(in) :: i, ilo, ihi
    ic = max(ilo, min(ihi, i))
  end function clamp_index


  pure subroutine interpolate_E_trilinear(xp, yp, zp, n, h, E, Exp, Eyp, Ezp)
    real(real64),   intent(in)  :: xp, yp, zp
    integer(int32), intent(in)  :: n(3)
    real(real64),   intent(in)  :: h(3)
    real(real64),   intent(in)  :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(out) :: Exp, Eyp, Ezp

    integer(int32) :: ix, iy, iz
    real(real64)   :: px, py, pz
    real(real64)   :: wx2, wy2, wz2
    real(real64)   :: w1, w2, w3, w4, w5, w6, w7, w8

    ix = int(xp / h(1), int32) + 1_int32
    iy = int(yp / h(2), int32) + 1_int32
    iz = int(zp / h(3), int32) + 1_int32

    ix = clamp_index(ix, 0_int32, n(1)+1_int32)
    iy = clamp_index(iy, 0_int32, n(2)+1_int32)
    iz = clamp_index(iz, 0_int32, n(3)+1_int32)

    px = (real(ix, real64)*h(1) - xp) / h(1)
    py = (real(iy, real64)*h(2) - yp) / h(2)
    pz = (real(iz, real64)*h(3) - zp) / h(3)

    wx2 = 1.0_real64 - px
    wy2 = 1.0_real64 - py
    wz2 = 1.0_real64 - pz

    w1 = px  * py  * pz
    w2 = wx2 * py  * pz
    w3 = wx2 * wy2 * pz
    w4 = px  * wy2 * pz
    w5 = px  * py  * wz2
    w6 = wx2 * py  * wz2
    w7 = wx2 * wy2 * wz2
    w8 = px  * wy2 * wz2

    Exp = w1*E(1,ix  ,iy  ,iz  ) + w2*E(1,ix+1,iy  ,iz  ) + &
          w3*E(1,ix+1,iy+1,iz  ) + w4*E(1,ix  ,iy+1,iz  ) + &
          w5*E(1,ix  ,iy  ,iz+1) + w6*E(1,ix+1,iy  ,iz+1) + &
          w7*E(1,ix+1,iy+1,iz+1) + w8*E(1,ix  ,iy+1,iz+1)

    Eyp = w1*E(2,ix  ,iy  ,iz  ) + w2*E(2,ix+1,iy  ,iz  ) + &
          w3*E(2,ix+1,iy+1,iz  ) + w4*E(2,ix  ,iy+1,iz  ) + &
          w5*E(2,ix  ,iy  ,iz+1) + w6*E(2,ix+1,iy  ,iz+1) + &
          w7*E(2,ix+1,iy+1,iz+1) + w8*E(2,ix  ,iy+1,iz+1)

    Ezp = w1*E(3,ix  ,iy  ,iz  ) + w2*E(3,ix+1,iy  ,iz  ) + &
          w3*E(3,ix+1,iy+1,iz  ) + w4*E(3,ix  ,iy+1,iz  ) + &
          w5*E(3,ix  ,iy  ,iz+1) + w6*E(3,ix+1,iy  ,iz+1) + &
          w7*E(3,ix+1,iy+1,iz+1) + w8*E(3,ix  ,iy+1,iz+1)
  end subroutine interpolate_E_trilinear


  subroutine move_particles_electrostatic(part, n, h, E, q, m, dt)
    ! Fast version:
    ! - same electrostatic update as before
    ! - trilinear interpolation is inlined inside the particle loop
    ! - avoids one procedure call per particle
    class(ParticleSet), intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(in)    :: q
    real(real64),       intent(in)    :: m
    real(real64),       intent(in)    :: dt

    integer(int32) :: i
    integer(int32) :: ix, iy, iz
    real(real64)   :: qmdt
    real(real64)   :: xp, yp, zp
    real(real64)   :: px, py, pz
    real(real64)   :: wx2, wy2, wz2
    real(real64)   :: w1, w2, w3, w4, w5, w6, w7, w8
    real(real64)   :: Exp, Eyp, Ezp

    if (.not. allocated(part%x)) return
    if (part%n <= 0_int32) return

    qmdt = dt*q/m

    do i = 1, part%n
      xp = part%x(i)
      yp = part%y(i)
      zp = part%z(i)

      ix = int(xp / h(1), int32) + 1_int32
      iy = int(yp / h(2), int32) + 1_int32
      iz = int(zp / h(3), int32) + 1_int32

      ix = max(0_int32, min(n(1)+1_int32, ix))
      iy = max(0_int32, min(n(2)+1_int32, iy))
      iz = max(0_int32, min(n(3)+1_int32, iz))

      px = (real(ix, real64)*h(1) - xp) / h(1)
      py = (real(iy, real64)*h(2) - yp) / h(2)
      pz = (real(iz, real64)*h(3) - zp) / h(3)

      wx2 = 1.0_real64 - px
      wy2 = 1.0_real64 - py
      wz2 = 1.0_real64 - pz

      w1 = px  * py  * pz
      w2 = wx2 * py  * pz
      w3 = wx2 * wy2 * pz
      w4 = px  * wy2 * pz
      w5 = px  * py  * wz2
      w6 = wx2 * py  * wz2
      w7 = wx2 * wy2 * wz2
      w8 = px  * wy2 * wz2

      Exp = w1*E(1,ix  ,iy  ,iz  ) + w2*E(1,ix+1,iy  ,iz  ) + &
            w3*E(1,ix+1,iy+1,iz  ) + w4*E(1,ix  ,iy+1,iz  ) + &
            w5*E(1,ix  ,iy  ,iz+1) + w6*E(1,ix+1,iy  ,iz+1) + &
            w7*E(1,ix+1,iy+1,iz+1) + w8*E(1,ix  ,iy+1,iz+1)

      Eyp = w1*E(2,ix  ,iy  ,iz  ) + w2*E(2,ix+1,iy  ,iz  ) + &
            w3*E(2,ix+1,iy+1,iz  ) + w4*E(2,ix  ,iy+1,iz  ) + &
            w5*E(2,ix  ,iy  ,iz+1) + w6*E(2,ix+1,iy  ,iz+1) + &
            w7*E(2,ix+1,iy+1,iz+1) + w8*E(2,ix  ,iy+1,iz+1)

      Ezp = w1*E(3,ix  ,iy  ,iz  ) + w2*E(3,ix+1,iy  ,iz  ) + &
            w3*E(3,ix+1,iy+1,iz  ) + w4*E(3,ix  ,iy+1,iz  ) + &
            w5*E(3,ix  ,iy  ,iz+1) + w6*E(3,ix+1,iy  ,iz+1) + &
            w7*E(3,ix+1,iy+1,iz+1) + w8*E(3,ix  ,iy+1,iz+1)

      part%vx(i) = part%vx(i) + qmdt*Exp
      part%vy(i) = part%vy(i) + qmdt*Eyp
      part%vz(i) = part%vz(i) + qmdt*Ezp

      part%x(i)  = xp + dt*part%vx(i)
      part%y(i)  = yp + dt*part%vy(i)
      part%z(i)  = zp + dt*part%vz(i)
    end do

    if (allocated(part%cell_id))    part%cell_id(1:part%n) = 0_int32
    if (allocated(part%cell_count)) part%cell_count = 0_int32
    if (allocated(part%cell_start)) part%cell_start = 0_int32
  end subroutine move_particles_electrostatic

end module mod_particleMover
