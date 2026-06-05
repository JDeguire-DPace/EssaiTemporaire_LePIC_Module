module mod_particleMover
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  implicit none
  private

  public :: move_particles_electrostatic
  public :: move_particles_boris
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


  ! Boris push: half E-kick -> magnetic rotation -> half E-kick -> position advance.
  ! When n_B == [1,1,1] the field is uniform and Bi(1:3,1,1,1) is read directly;
  ! otherwise trilinear interpolation is performed on the B-field grid (h_B / n_B).
  subroutine move_particles_boris(part, n, h, E, n_B, h_B, Bi, q, m, dt)
    class(ParticleSet), intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32),     intent(in)    :: n_B(3)
    real(real64),       intent(in)    :: h_B(3)
    real(real64),       intent(in)    :: Bi(4,0:n_B(1)+2,0:n_B(2)+2,0:n_B(3)+2)
    real(real64),       intent(in)    :: q
    real(real64),       intent(in)    :: m
    real(real64),       intent(in)    :: dt

    integer(int32) :: i
    integer(int32) :: ix, iy, iz
    real(real64)   :: qm2dt
    real(real64)   :: xp, yp, zp
    real(real64)   :: px, py, pz
    real(real64)   :: wx2, wy2, wz2
    real(real64)   :: w1, w2, w3, w4, w5, w6, w7, w8
    real(real64)   :: Exp, Eyp, Ezp
    real(real64)   :: Bpx, Bpy, Bpz, Btot, k2
    real(real64)   :: v1x, v1y, v1z
    real(real64)   :: v3x, v3y, v3z
    real(real64)   :: v2x, v2y, v2z
    logical        :: uniform_B

    if (.not. allocated(part%x)) return
    if (part%n <= 0_int32) return

    qm2dt    = 0.5_real64 * dt * q / m
    uniform_B = (n_B(1) == 1)

    do i = 1, part%n
      xp = part%x(i)
      yp = part%y(i)
      zp = part%z(i)

      ! E-field interpolation
      ix = int(xp / h(1), int32) + 1_int32
      iy = int(yp / h(2), int32) + 1_int32
      iz = int(zp / h(3), int32) + 1_int32

      ix = max(0_int32, min(n(1)+1_int32, ix))
      iy = max(0_int32, min(n(2)+1_int32, iy))
      iz = max(0_int32, min(n(3)+1_int32, iz))

      px  = (real(ix, real64)*h(1) - xp) / h(1)
      py  = (real(iy, real64)*h(2) - yp) / h(2)
      pz  = (real(iz, real64)*h(3) - zp) / h(3)
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

      ! B-field at particle position
      if (uniform_B) then
        Bpx = Bi(1,1,1,1)
        Bpy = Bi(2,1,1,1)
        Bpz = Bi(3,1,1,1)
      else
        ix = int(xp / h_B(1), int32) + 1_int32
        iy = int(yp / h_B(2), int32) + 1_int32
        iz = int(zp / h_B(3), int32) + 1_int32

        ix = max(0_int32, min(n_B(1)+1_int32, ix))
        iy = max(0_int32, min(n_B(2)+1_int32, iy))
        iz = max(0_int32, min(n_B(3)+1_int32, iz))

        px  = (real(ix, real64)*h_B(1) - xp) / h_B(1)
        py  = (real(iy, real64)*h_B(2) - yp) / h_B(2)
        pz  = (real(iz, real64)*h_B(3) - zp) / h_B(3)
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

        Bpx = w1*Bi(1,ix  ,iy  ,iz  ) + w2*Bi(1,ix+1,iy  ,iz  ) + &
              w3*Bi(1,ix+1,iy+1,iz  ) + w4*Bi(1,ix  ,iy+1,iz  ) + &
              w5*Bi(1,ix  ,iy  ,iz+1) + w6*Bi(1,ix+1,iy  ,iz+1) + &
              w7*Bi(1,ix+1,iy+1,iz+1) + w8*Bi(1,ix  ,iy+1,iz+1)

        Bpy = w1*Bi(2,ix  ,iy  ,iz  ) + w2*Bi(2,ix+1,iy  ,iz  ) + &
              w3*Bi(2,ix+1,iy+1,iz  ) + w4*Bi(2,ix  ,iy+1,iz  ) + &
              w5*Bi(2,ix  ,iy  ,iz+1) + w6*Bi(2,ix+1,iy  ,iz+1) + &
              w7*Bi(2,ix+1,iy+1,iz+1) + w8*Bi(2,ix  ,iy+1,iz+1)

        Bpz = w1*Bi(3,ix  ,iy  ,iz  ) + w2*Bi(3,ix+1,iy  ,iz  ) + &
              w3*Bi(3,ix+1,iy+1,iz  ) + w4*Bi(3,ix  ,iy+1,iz  ) + &
              w5*Bi(3,ix  ,iy  ,iz+1) + w6*Bi(3,ix+1,iy  ,iz+1) + &
              w7*Bi(3,ix+1,iy+1,iz+1) + w8*Bi(3,ix  ,iy+1,iz+1)
      end if

      ! Boris rotation
      Btot = sqrt(Bpx*Bpx + Bpy*Bpy + Bpz*Bpz)
      k2   = qm2dt * 2.0_real64 / (1.0_real64 + (qm2dt*Btot)**2)

      v1x = part%vx(i) + qm2dt*Exp
      v1y = part%vy(i) + qm2dt*Eyp
      v1z = part%vz(i) + qm2dt*Ezp

      v3x = v1x + qm2dt*(v1y*Bpz - v1z*Bpy)
      v3y = v1y + qm2dt*(v1z*Bpx - v1x*Bpz)
      v3z = v1z + qm2dt*(v1x*Bpy - v1y*Bpx)

      v2x = v1x + k2*(v3y*Bpz - v3z*Bpy)
      v2y = v1y + k2*(v3z*Bpx - v3x*Bpz)
      v2z = v1z + k2*(v3x*Bpy - v3y*Bpx)

      part%vx(i) = v2x + qm2dt*Exp
      part%vy(i) = v2y + qm2dt*Eyp
      part%vz(i) = v2z + qm2dt*Ezp

      part%x(i) = xp + dt*part%vx(i)
      part%y(i) = yp + dt*part%vy(i)
      part%z(i) = zp + dt*part%vz(i)
    end do

    if (allocated(part%cell_id))    part%cell_id(1:part%n) = 0_int32
    if (allocated(part%cell_count)) part%cell_count = 0_int32
    if (allocated(part%cell_start)) part%cell_start = 0_int32
  end subroutine move_particles_boris

end module mod_particleMover
