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
    !=============================================================
    ! Trilinear interpolation of the electric field at particle
    ! position, adapted from the legacy part_expmover.f90 routine.
    !
    ! Legacy convention:
    !   ix = INT(x/hx) + 1
    !   px = (ix*hx - x)/hx
    !
    ! and similarly for y,z.
    !=============================================================
    real(real64),   intent(in)  :: xp, yp, zp
    integer(int32), intent(in)  :: n(3)
    real(real64),   intent(in)  :: h(3)
    real(real64),   intent(in)  :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(out) :: Exp, Eyp, Ezp

    integer(int32) :: ix, iy, iz
    real(real64)   :: px, py, pz
    real(real64)   :: ki(8)

    ix = int(xp / h(1), int32) + 1_int32
    iy = int(yp / h(2), int32) + 1_int32
    iz = int(zp / h(3), int32) + 1_int32

    ! Clamp so ix+1, iy+1, iz+1 remain valid
    ix = clamp_index(ix, 0_int32, n(1)+1_int32)
    iy = clamp_index(iy, 0_int32, n(2)+1_int32)
    iz = clamp_index(iz, 0_int32, n(3)+1_int32)

    px = (real(ix, real64)*h(1) - xp) / h(1)
    py = (real(iy, real64)*h(2) - yp) / h(2)
    pz = (real(iz, real64)*h(3) - zp) / h(3)

    ki(1) = px*py*pz
    ki(2) = (1.0_real64-px)*py*pz
    ki(3) = (1.0_real64-px)*(1.0_real64-py)*pz
    ki(4) = px*(1.0_real64-py)*pz
    ki(5) = px*py*(1.0_real64-pz)
    ki(6) = (1.0_real64-px)*py*(1.0_real64-pz)
    ki(7) = (1.0_real64-px)*(1.0_real64-py)*(1.0_real64-pz)
    ki(8) = px*(1.0_real64-py)*(1.0_real64-pz)

    Exp = ki(1)*E(1,ix  ,iy  ,iz  ) + &
          ki(2)*E(1,ix+1,iy  ,iz  ) + &
          ki(3)*E(1,ix+1,iy+1,iz  ) + &
          ki(4)*E(1,ix  ,iy+1,iz  ) + &
          ki(5)*E(1,ix  ,iy  ,iz+1) + &
          ki(6)*E(1,ix+1,iy  ,iz+1) + &
          ki(7)*E(1,ix+1,iy+1,iz+1) + &
          ki(8)*E(1,ix  ,iy+1,iz+1)

    Eyp = ki(1)*E(2,ix  ,iy  ,iz  ) + &
          ki(2)*E(2,ix+1,iy  ,iz  ) + &
          ki(3)*E(2,ix+1,iy+1,iz  ) + &
          ki(4)*E(2,ix  ,iy+1,iz  ) + &
          ki(5)*E(2,ix  ,iy  ,iz+1) + &
          ki(6)*E(2,ix+1,iy  ,iz+1) + &
          ki(7)*E(2,ix+1,iy+1,iz+1) + &
          ki(8)*E(2,ix  ,iy+1,iz+1)

    Ezp = ki(1)*E(3,ix  ,iy  ,iz  ) + &
          ki(2)*E(3,ix+1,iy  ,iz  ) + &
          ki(3)*E(3,ix+1,iy+1,iz  ) + &
          ki(4)*E(3,ix  ,iy+1,iz  ) + &
          ki(5)*E(3,ix  ,iy  ,iz+1) + &
          ki(6)*E(3,ix+1,iy  ,iz+1) + &
          ki(7)*E(3,ix+1,iy+1,iz+1) + &
          ki(8)*E(3,ix  ,iy+1,iz+1)
  end subroutine interpolate_E_trilinear


  subroutine move_particles_electrostatic(part, n, h, E, q, m, dt)
    !=============================================================
    ! First modern layer of the legacy part_mover:
    !
    !   1) interpolate E at particle position
    !   2) electrostatic leapfrog-style velocity update
    !   3) update particle position
    !
    ! Adapted from the legacy branch:
    !
    !   k1 = dt*q/(2*m)
    !   vx = vx + 2*k1*Ex
    !   vy = vy + 2*k1*Ey
    !   vz = vz + 2*k1*Ez
    !
    !   x = x + dt*vx
    !   y = y + dt*vy
    !   z = z + dt*vz
    !
    ! This routine does NOT yet apply:
    !   - wall / periodic BCs
    !   - particle deletion
    !   - magnetic field Boris push
    !   - MPI transfer
    !
    ! Sorting metadata is invalidated after motion.
    !=============================================================
    class(ParticleSet), intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(in)    :: q
    real(real64),       intent(in)    :: m
    real(real64),       intent(in)    :: dt

    integer(int32) :: i
    real(real64)   :: k1
    real(real64)   :: Exp, Eyp, Ezp

    if (.not. allocated(part%x)) return
    if (part%n <= 0_int32) return

    k1 = dt*q/(2.0_real64*m)

    do i = 1, part%n

      call interpolate_E_trilinear( &
           xp  = part%x(i), &
           yp  = part%y(i), &
           zp  = part%z(i), &
           n   = n, &
           h   = h, &
           E   = E, &
           Exp = Exp, &
           Eyp = Eyp, &
           Ezp = Ezp )

      ! Legacy electrostatic branch:
      part%vx(i) = part%vx(i) + 2.0_real64*k1*Exp
      part%vy(i) = part%vy(i) + 2.0_real64*k1*Eyp
      part%vz(i) = part%vz(i) + 2.0_real64*k1*Ezp

      part%x(i)  = part%x(i)  + dt*part%vx(i)
      part%y(i)  = part%y(i)  + dt*part%vy(i)
      part%z(i)  = part%z(i)  + dt*part%vz(i)

    end do

    ! Motion invalidates cell sorting metadata.
    if (allocated(part%cell_id))    part%cell_id(1:part%n) = 0_int32
    if (allocated(part%cell_count)) part%cell_count = 0_int32
    if (allocated(part%cell_start)) part%cell_start = 0_int32

  end subroutine move_particles_electrostatic

end module mod_particleMover