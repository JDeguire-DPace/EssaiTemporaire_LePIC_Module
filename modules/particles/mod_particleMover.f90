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
		integer(int32) :: ixB, iyB, izB

		real(real64) :: qm2dt
		real(real64) :: xp, yp, zp
		real(real64) :: px, py, pz
		real(real64) :: wx2, wy2, wz2
		real(real64) :: w1, w2, w3, w4, w5, w6, w7, w8

		real(real64) :: pxB, pyB, pzB
		real(real64) :: wx2B, wy2B, wz2B
		real(real64) :: b1, b2, b3, b4, b5, b6, b7, b8

		real(real64) :: Exp, Eyp, Ezp
		real(real64) :: Bpx, Bpy, Bpz

		real(real64) :: vminus_x, vminus_y, vminus_z
		real(real64) :: vprime_x, vprime_y, vprime_z
		real(real64) :: vplus_x, vplus_y, vplus_z

		real(real64) :: tx, ty, tz
		real(real64) :: sx, sy, sz
		real(real64) :: t2, inv_denom

		logical :: uniform_B, same_grid_B

		if (.not. allocated(part%x)) return
		if (part%n <= 0_int32) return

		qm2dt = 0.5_real64 * dt * q / m

		uniform_B = (n_B(1) == 1_int32 .and. n_B(2) == 1_int32 .and. n_B(3) == 1_int32)

		same_grid_B = .false.
		if (.not. uniform_B) then
			same_grid_B = all(n_B == n) .and. &
										abs(h_B(1)-h(1)) < 1.0e-14_real64 .and. &
										abs(h_B(2)-h(2)) < 1.0e-14_real64 .and. &
										abs(h_B(3)-h(3)) < 1.0e-14_real64
		end if

		do i = 1_int32, part%n

			xp = part%x(i)
			yp = part%y(i)
			zp = part%z(i)

			! ------------------------------------------------------------
			! E-field interpolation
			! ------------------------------------------------------------
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

			! ------------------------------------------------------------
			! B-field interpolation
			! ------------------------------------------------------------
			if (uniform_B) then

				Bpx = Bi(1,1,1,1)
				Bpy = Bi(2,1,1,1)
				Bpz = Bi(3,1,1,1)

			else if (same_grid_B) then

				! Reuse E-grid indices and weights.
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

			else

				ixB = int(xp / h_B(1), int32) + 1_int32
				iyB = int(yp / h_B(2), int32) + 1_int32
				izB = int(zp / h_B(3), int32) + 1_int32

				ixB = max(0_int32, min(n_B(1)+1_int32, ixB))
				iyB = max(0_int32, min(n_B(2)+1_int32, iyB))
				izB = max(0_int32, min(n_B(3)+1_int32, izB))

				pxB = (real(ixB, real64)*h_B(1) - xp) / h_B(1)
				pyB = (real(iyB, real64)*h_B(2) - yp) / h_B(2)
				pzB = (real(izB, real64)*h_B(3) - zp) / h_B(3)

				wx2B = 1.0_real64 - pxB
				wy2B = 1.0_real64 - pyB
				wz2B = 1.0_real64 - pzB

				b1 = pxB  * pyB  * pzB
				b2 = wx2B * pyB  * pzB
				b3 = wx2B * wy2B * pzB
				b4 = pxB  * wy2B * pzB
				b5 = pxB  * pyB  * wz2B
				b6 = wx2B * pyB  * wz2B
				b7 = wx2B * wy2B * wz2B
				b8 = pxB  * wy2B * wz2B

				Bpx = b1*Bi(1,ixB  ,iyB  ,izB  ) + b2*Bi(1,ixB+1,iyB  ,izB  ) + &
							b3*Bi(1,ixB+1,iyB+1,izB  ) + b4*Bi(1,ixB  ,iyB+1,izB  ) + &
							b5*Bi(1,ixB  ,iyB  ,izB+1) + b6*Bi(1,ixB+1,iyB  ,izB+1) + &
							b7*Bi(1,ixB+1,iyB+1,izB+1) + b8*Bi(1,ixB  ,iyB+1,izB+1)

				Bpy = b1*Bi(2,ixB  ,iyB  ,izB  ) + b2*Bi(2,ixB+1,iyB  ,izB  ) + &
							b3*Bi(2,ixB+1,iyB+1,izB  ) + b4*Bi(2,ixB  ,iyB+1,izB  ) + &
							b5*Bi(2,ixB  ,iyB  ,izB+1) + b6*Bi(2,ixB+1,iyB  ,izB+1) + &
							b7*Bi(2,ixB+1,iyB+1,izB+1) + b8*Bi(2,ixB  ,iyB+1,izB+1)

				Bpz = b1*Bi(3,ixB  ,iyB  ,izB  ) + b2*Bi(3,ixB+1,iyB  ,izB  ) + &
							b3*Bi(3,ixB+1,iyB+1,izB  ) + b4*Bi(3,ixB  ,iyB+1,izB  ) + &
							b5*Bi(3,ixB  ,iyB  ,izB+1) + b6*Bi(3,ixB+1,iyB  ,izB+1) + &
							b7*Bi(3,ixB+1,iyB+1,izB+1) + b8*Bi(3,ixB  ,iyB+1,izB+1)

			end if

			! ------------------------------------------------------------
			! Boris push without sqrt(|B|^2)
			! ------------------------------------------------------------

			vminus_x = part%vx(i) + qm2dt*Exp
			vminus_y = part%vy(i) + qm2dt*Eyp
			vminus_z = part%vz(i) + qm2dt*Ezp

			tx = qm2dt * Bpx
			ty = qm2dt * Bpy
			tz = qm2dt * Bpz

			t2 = tx*tx + ty*ty + tz*tz
			inv_denom = 1.0_real64 / (1.0_real64 + t2)

			sx = 2.0_real64 * tx * inv_denom
			sy = 2.0_real64 * ty * inv_denom
			sz = 2.0_real64 * tz * inv_denom

			vprime_x = vminus_x + (vminus_y*tz - vminus_z*ty)
			vprime_y = vminus_y + (vminus_z*tx - vminus_x*tz)
			vprime_z = vminus_z + (vminus_x*ty - vminus_y*tx)

			vplus_x = vminus_x + (vprime_y*sz - vprime_z*sy)
			vplus_y = vminus_y + (vprime_z*sx - vprime_x*sz)
			vplus_z = vminus_z + (vprime_x*sy - vprime_y*sx)

			part%vx(i) = vplus_x + qm2dt*Exp
			part%vy(i) = vplus_y + qm2dt*Eyp
			part%vz(i) = vplus_z + qm2dt*Ezp

			part%x(i) = xp + dt*part%vx(i)
			part%y(i) = yp + dt*part%vy(i)
			part%z(i) = zp + dt*part%vz(i)

		end do

		if (allocated(part%cell_id))    part%cell_id(1:part%n) = 0_int32
		if (allocated(part%cell_count)) part%cell_count = 0_int32
		if (allocated(part%cell_start)) part%cell_start = 0_int32

	end subroutine move_particles_boris

end module mod_particleMover
