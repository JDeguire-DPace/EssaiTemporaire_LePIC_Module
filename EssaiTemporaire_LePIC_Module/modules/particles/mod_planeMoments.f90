module mod_planeMoments
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  use mod_constants,   only: qe
  use mod_simParams,   only: SimParams
  implicit none
  private

  public :: compute_particle_plane_moments_species

contains

  subroutine compute_particle_plane_moments_species(part, nproc, n, h, ptype, mass_species, &
                                                    ix_plane, iy_plane, iz_plane,             &
                                                    data_xy, data_xz, data_yz, params)
    type(ParticleSet), intent(in)    :: part(:)
    integer(int32),    intent(in)    :: nproc
    integer(int32),    intent(in)    :: n(3)
    integer(int32),    intent(in)    :: ptype
    integer(int32),    intent(in)    :: ix_plane, iy_plane, iz_plane
    real(real64),      intent(in)    :: h(3)
    real(real64),      intent(in)    :: mass_species
    real(real64),      intent(inout) :: data_xy(5,0:n(1)+2,0:n(2)+2)
    real(real64),      intent(inout) :: data_xz(5,0:n(1)+2,0:n(3)+2)
    real(real64),      intent(inout) :: data_yz(5,0:n(2)+2,0:n(3)+2)
    type(SimParams),   intent(in)    :: params

    integer(int32), parameter :: np_avg=1, Tp_avg=2, u1_avg=3, u2_avg=4, u3_avg=5

    integer(int32) :: iproc, i
    integer(int32) :: ix, iy, iz
    real(real64)   :: x, y, z, vx, vy, vz
    real(real64)   :: px, py, pz
    real(real64)   :: pplane, v2
    real(real64)   :: wx1, wx2, wy1, wy2, wz1, wz2
    real(real64)   :: w11, w21, w22, w12
    real(real64)   :: u2, v2mean, thermal_v2

    real(real64), allocatable :: cnt_xy(:,:), cnt_xz(:,:), cnt_yz(:,:)
    real(real64), allocatable :: vx2_xy(:,:), vy2_xy(:,:), vz2_xy(:,:)
    real(real64), allocatable :: vx2_xz(:,:), vy2_xz(:,:), vz2_xz(:,:)
    real(real64), allocatable :: vx2_yz(:,:), vy2_yz(:,:), vz2_yz(:,:)

    allocate(cnt_xy(0:n(1)+2,0:n(2)+2), cnt_xz(0:n(1)+2,0:n(3)+2), cnt_yz(0:n(2)+2,0:n(3)+2))
    allocate(vx2_xy(0:n(1)+2,0:n(2)+2), vy2_xy(0:n(1)+2,0:n(2)+2), vz2_xy(0:n(1)+2,0:n(2)+2))
    allocate(vx2_xz(0:n(1)+2,0:n(3)+2), vy2_xz(0:n(1)+2,0:n(3)+2), vz2_xz(0:n(1)+2,0:n(3)+2))
    allocate(vx2_yz(0:n(2)+2,0:n(3)+2), vy2_yz(0:n(2)+2,0:n(3)+2), vz2_yz(0:n(2)+2,0:n(3)+2))

    cnt_xy = 0.0_real64
    cnt_xz = 0.0_real64
    cnt_yz = 0.0_real64

    vx2_xy = 0.0_real64 ; vy2_xy = 0.0_real64 ; vz2_xy = 0.0_real64
    vx2_xz = 0.0_real64 ; vy2_xz = 0.0_real64 ; vz2_xz = 0.0_real64
    vx2_yz = 0.0_real64 ; vy2_yz = 0.0_real64 ; vz2_yz = 0.0_real64

    do iproc = 1, nproc
      if (.not. allocated(part(iproc)%x)) cycle
      if (part(iproc)%n <= 0_int32) cycle

      do i = 1, part(iproc)%n
        x = part(iproc)%x(i) - part(iproc)%vx(i) * params%dt * 0.5_real64
        y = part(iproc)%y(i) - part(iproc)%vy(i) * params%dt * 0.5_real64
        z = part(iproc)%z(i) - part(iproc)%vz(i) * params%dt * 0.5_real64
        vx = part(iproc)%vx(i)
        vy = part(iproc)%vy(i)
        vz = part(iproc)%vz(i)
        v2 = vx*vx + vy*vy + vz*vz

        ix = int(x / h(1), int32) + 1_int32
        iy = int(y / h(2), int32) + 1_int32
        iz = int(z / h(3), int32) + 1_int32

        if (ix < 1_int32 .or. ix > n(1)) cycle
        if (iy < 1_int32 .or. iy > n(2)) cycle
        if (iz < 1_int32 .or. iz > n(3)) cycle

        px = (real(ix, real64)*h(1) - x) / h(1)
        py = (real(iy, real64)*h(2) - y) / h(2)
        pz = (real(iz, real64)*h(3) - z) / h(3)

        px = max(0.0_real64, min(1.0_real64, px))
        py = max(0.0_real64, min(1.0_real64, py))
        pz = max(0.0_real64, min(1.0_real64, pz))

        wx1 = px
        wx2 = 1.0_real64 - px
        wy1 = py
        wy2 = 1.0_real64 - py
        wz1 = pz
        wz2 = 1.0_real64 - pz

        ! --------------------------------------------------
        ! XY plane: interpolate in z to plane, then deposit
        ! bilinearly in x-y
        ! --------------------------------------------------
        if (iz == iz_plane .or. iz == iz_plane-1_int32) then
          if (iz == iz_plane) then
            pplane = pz
          else
            pplane = 1.0_real64 - pz
          end if

          w11 = pplane * wx1 * wy1
          w21 = pplane * wx2 * wy1
          w22 = pplane * wx2 * wy2
          w12 = pplane * wx1 * wy2

          cnt_xy(ix  ,iy  ) = cnt_xy(ix  ,iy  ) + w11
          cnt_xy(ix+1,iy  ) = cnt_xy(ix+1,iy  ) + w21
          cnt_xy(ix+1,iy+1) = cnt_xy(ix+1,iy+1) + w22
          cnt_xy(ix  ,iy+1) = cnt_xy(ix  ,iy+1) + w12

          data_xy(u1_avg,ix  ,iy  ) = data_xy(u1_avg,ix  ,iy  ) + w11*vx
          data_xy(u1_avg,ix+1,iy  ) = data_xy(u1_avg,ix+1,iy  ) + w21*vx
          data_xy(u1_avg,ix+1,iy+1) = data_xy(u1_avg,ix+1,iy+1) + w22*vx
          data_xy(u1_avg,ix  ,iy+1) = data_xy(u1_avg,ix  ,iy+1) + w12*vx

          data_xy(u2_avg,ix  ,iy  ) = data_xy(u2_avg,ix  ,iy  ) + w11*vy
          data_xy(u2_avg,ix+1,iy  ) = data_xy(u2_avg,ix+1,iy  ) + w21*vy
          data_xy(u2_avg,ix+1,iy+1) = data_xy(u2_avg,ix+1,iy+1) + w22*vy
          data_xy(u2_avg,ix  ,iy+1) = data_xy(u2_avg,ix  ,iy+1) + w12*vy

          data_xy(u3_avg,ix  ,iy  ) = data_xy(u3_avg,ix  ,iy  ) + w11*vz
          data_xy(u3_avg,ix+1,iy  ) = data_xy(u3_avg,ix+1,iy  ) + w21*vz
          data_xy(u3_avg,ix+1,iy+1) = data_xy(u3_avg,ix+1,iy+1) + w22*vz
          data_xy(u3_avg,ix  ,iy+1) = data_xy(u3_avg,ix  ,iy+1) + w12*vz

          vx2_xy(ix  ,iy  ) = vx2_xy(ix  ,iy  ) + w11*vx*vx
          vx2_xy(ix+1,iy  ) = vx2_xy(ix+1,iy  ) + w21*vx*vx
          vx2_xy(ix+1,iy+1) = vx2_xy(ix+1,iy+1) + w22*vx*vx
          vx2_xy(ix  ,iy+1) = vx2_xy(ix  ,iy+1) + w12*vx*vx

          vy2_xy(ix  ,iy  ) = vy2_xy(ix  ,iy  ) + w11*vy*vy
          vy2_xy(ix+1,iy  ) = vy2_xy(ix+1,iy  ) + w21*vy*vy
          vy2_xy(ix+1,iy+1) = vy2_xy(ix+1,iy+1) + w22*vy*vy
          vy2_xy(ix  ,iy+1) = vy2_xy(ix  ,iy+1) + w12*vy*vy

          vz2_xy(ix  ,iy  ) = vz2_xy(ix  ,iy  ) + w11*vz*vz
          vz2_xy(ix+1,iy  ) = vz2_xy(ix+1,iy  ) + w21*vz*vz
          vz2_xy(ix+1,iy+1) = vz2_xy(ix+1,iy+1) + w22*vz*vz
          vz2_xy(ix  ,iy+1) = vz2_xy(ix  ,iy+1) + w12*vz*vz
        end if

        ! --------------------------------------------------
        ! XZ plane: interpolate in y to plane, then deposit
        ! bilinearly in x-z
        ! --------------------------------------------------
        if (iy == iy_plane .or. iy == iy_plane-1_int32) then
          if (iy == iy_plane) then
            pplane = py
          else
            pplane = 1.0_real64 - py
          end if

          w11 = pplane * wx1 * wz1
          w21 = pplane * wx2 * wz1
          w22 = pplane * wx2 * wz2
          w12 = pplane * wx1 * wz2

          cnt_xz(ix  ,iz  ) = cnt_xz(ix  ,iz  ) + w11
          cnt_xz(ix+1,iz  ) = cnt_xz(ix+1,iz  ) + w21
          cnt_xz(ix+1,iz+1) = cnt_xz(ix+1,iz+1) + w22
          cnt_xz(ix  ,iz+1) = cnt_xz(ix  ,iz+1) + w12

          data_xz(u1_avg,ix  ,iz  ) = data_xz(u1_avg,ix  ,iz  ) + w11*vx
          data_xz(u1_avg,ix+1,iz  ) = data_xz(u1_avg,ix+1,iz  ) + w21*vx
          data_xz(u1_avg,ix+1,iz+1) = data_xz(u1_avg,ix+1,iz+1) + w22*vx
          data_xz(u1_avg,ix  ,iz+1) = data_xz(u1_avg,ix  ,iz+1) + w12*vx

          data_xz(u2_avg,ix  ,iz  ) = data_xz(u2_avg,ix  ,iz  ) + w11*vy
          data_xz(u2_avg,ix+1,iz  ) = data_xz(u2_avg,ix+1,iz  ) + w21*vy
          data_xz(u2_avg,ix+1,iz+1) = data_xz(u2_avg,ix+1,iz+1) + w22*vy
          data_xz(u2_avg,ix  ,iz+1) = data_xz(u2_avg,ix  ,iz+1) + w12*vy

          data_xz(u3_avg,ix  ,iz  ) = data_xz(u3_avg,ix  ,iz  ) + w11*vz
          data_xz(u3_avg,ix+1,iz  ) = data_xz(u3_avg,ix+1,iz  ) + w21*vz
          data_xz(u3_avg,ix+1,iz+1) = data_xz(u3_avg,ix+1,iz+1) + w22*vz
          data_xz(u3_avg,ix  ,iz+1) = data_xz(u3_avg,ix  ,iz+1) + w12*vz

          vx2_xz(ix  ,iz  ) = vx2_xz(ix  ,iz  ) + w11*vx*vx
          vx2_xz(ix+1,iz  ) = vx2_xz(ix+1,iz  ) + w21*vx*vx
          vx2_xz(ix+1,iz+1) = vx2_xz(ix+1,iz+1) + w22*vx*vx
          vx2_xz(ix  ,iz+1) = vx2_xz(ix  ,iz+1) + w12*vx*vx

          vy2_xz(ix  ,iz  ) = vy2_xz(ix  ,iz  ) + w11*vy*vy
          vy2_xz(ix+1,iz  ) = vy2_xz(ix+1,iz  ) + w21*vy*vy
          vy2_xz(ix+1,iz+1) = vy2_xz(ix+1,iz+1) + w22*vy*vy
          vy2_xz(ix  ,iz+1) = vy2_xz(ix  ,iz+1) + w12*vy*vy

          vz2_xz(ix  ,iz  ) = vz2_xz(ix  ,iz  ) + w11*vz*vz
          vz2_xz(ix+1,iz  ) = vz2_xz(ix+1,iz  ) + w21*vz*vz
          vz2_xz(ix+1,iz+1) = vz2_xz(ix+1,iz+1) + w22*vz*vz
          vz2_xz(ix  ,iz+1) = vz2_xz(ix  ,iz+1) + w12*vz*vz
        end if

        ! --------------------------------------------------
        ! YZ plane: interpolate in x to plane, then deposit
        ! bilinearly in y-z
        ! --------------------------------------------------
        if (ix == ix_plane .or. ix == ix_plane-1_int32) then
          if (ix == ix_plane) then
            pplane = px
          else
            pplane = 1.0_real64 - px
          end if

          w11 = pplane * wy1 * wz1
          w21 = pplane * wy2 * wz1
          w22 = pplane * wy2 * wz2
          w12 = pplane * wy1 * wz2

          cnt_yz(iy  ,iz  ) = cnt_yz(iy  ,iz  ) + w11
          cnt_yz(iy+1,iz  ) = cnt_yz(iy+1,iz  ) + w21
          cnt_yz(iy+1,iz+1) = cnt_yz(iy+1,iz+1) + w22
          cnt_yz(iy  ,iz+1) = cnt_yz(iy  ,iz+1) + w12

          data_yz(u1_avg,iy  ,iz  ) = data_yz(u1_avg,iy  ,iz  ) + w11*vx
          data_yz(u1_avg,iy+1,iz  ) = data_yz(u1_avg,iy+1,iz  ) + w21*vx
          data_yz(u1_avg,iy+1,iz+1) = data_yz(u1_avg,iy+1,iz+1) + w22*vx
          data_yz(u1_avg,iy  ,iz+1) = data_yz(u1_avg,iy  ,iz+1) + w12*vx

          data_yz(u2_avg,iy  ,iz  ) = data_yz(u2_avg,iy  ,iz  ) + w11*vy
          data_yz(u2_avg,iy+1,iz  ) = data_yz(u2_avg,iy+1,iz  ) + w21*vy
          data_yz(u2_avg,iy+1,iz+1) = data_yz(u2_avg,iy+1,iz+1) + w22*vy
          data_yz(u2_avg,iy  ,iz+1) = data_yz(u2_avg,iy  ,iz+1) + w12*vy

          data_yz(u3_avg,iy  ,iz  ) = data_yz(u3_avg,iy  ,iz  ) + w11*vz
          data_yz(u3_avg,iy+1,iz  ) = data_yz(u3_avg,iy+1,iz  ) + w21*vz
          data_yz(u3_avg,iy+1,iz+1) = data_yz(u3_avg,iy+1,iz+1) + w22*vz
          data_yz(u3_avg,iy  ,iz+1) = data_yz(u3_avg,iy  ,iz+1) + w12*vz

          vx2_yz(iy  ,iz  ) = vx2_yz(iy  ,iz  ) + w11*vx*vx
          vx2_yz(iy+1,iz  ) = vx2_yz(iy+1,iz  ) + w21*vx*vx
          vx2_yz(iy+1,iz+1) = vx2_yz(iy+1,iz+1) + w22*vx*vx
          vx2_yz(iy  ,iz+1) = vx2_yz(iy  ,iz+1) + w12*vx*vx

          vy2_yz(iy  ,iz  ) = vy2_yz(iy  ,iz  ) + w11*vy*vy
          vy2_yz(iy+1,iz  ) = vy2_yz(iy+1,iz  ) + w21*vy*vy
          vy2_yz(iy+1,iz+1) = vy2_yz(iy+1,iz+1) + w22*vy*vy
          vy2_yz(iy  ,iz+1) = vy2_yz(iy  ,iz+1) + w12*vy*vy

          vz2_yz(iy  ,iz  ) = vz2_yz(iy  ,iz  ) + w11*vz*vz
          vz2_yz(iy+1,iz  ) = vz2_yz(iy+1,iz  ) + w21*vz*vz
          vz2_yz(iy+1,iz+1) = vz2_yz(iy+1,iz+1) + w22*vz*vz
          vz2_yz(iy  ,iz+1) = vz2_yz(iy  ,iz+1) + w12*vz*vz
        end if

      end do
    end do

    do iy = 0, n(2)+2
      do ix = 0, n(1)+2
        if (cnt_xy(ix,iy) > 0.0_real64) then
          data_xy(np_avg,ix,iy) = cnt_xy(ix,iy)
          data_xy(u1_avg,ix,iy) = data_xy(u1_avg,ix,iy) / cnt_xy(ix,iy)
          data_xy(u2_avg,ix,iy) = data_xy(u2_avg,ix,iy) / cnt_xy(ix,iy)
          data_xy(u3_avg,ix,iy) = data_xy(u3_avg,ix,iy) / cnt_xy(ix,iy)

          u2 = data_xy(u1_avg,ix,iy)**2 + data_xy(u2_avg,ix,iy)**2 + data_xy(u3_avg,ix,iy)**2
          v2mean = (vx2_xy(ix,iy) + vy2_xy(ix,iy) + vz2_xy(ix,iy)) / cnt_xy(ix,iy)
          thermal_v2 = max(0.0_real64, v2mean - u2)
          data_xy(Tp_avg,ix,iy) = data_xy(Tp_avg,ix,iy) + &
                        mass_species * thermal_v2 / (3.0_real64 * qe)
        end if
      end do
    end do

    do iz = 0, n(3)+2
      do ix = 0, n(1)+2
        if (cnt_xz(ix,iz) > 0.0_real64) then
          data_xz(np_avg,ix,iz) = cnt_xz(ix,iz)
          data_xz(u1_avg,ix,iz) = data_xz(u1_avg,ix,iz) / cnt_xz(ix,iz)
          data_xz(u2_avg,ix,iz) = data_xz(u2_avg,ix,iz) / cnt_xz(ix,iz)
          data_xz(u3_avg,ix,iz) = data_xz(u3_avg,ix,iz) / cnt_xz(ix,iz)

          u2 = data_xz(u1_avg,ix,iz)**2 + data_xz(u2_avg,ix,iz)**2 + data_xz(u3_avg,ix,iz)**2
          v2mean = (vx2_xz(ix,iz) + vy2_xz(ix,iz) + vz2_xz(ix,iz)) / cnt_xz(ix,iz)
          thermal_v2 = max(0.0_real64, v2mean - u2)
          data_xz(Tp_avg,ix,iz) = data_xz(Tp_avg,ix,iz) + &
                      mass_species * thermal_v2 / (3.0_real64 * qe)
        end if
      end do
    end do

    do iz = 0, n(3)+2
      do iy = 0, n(2)+2
        if (cnt_yz(iy,iz) > 0.0_real64) then
          data_yz(np_avg,iy,iz) = cnt_yz(iy,iz)
          data_yz(u1_avg,iy,iz) = data_yz(u1_avg,iy,iz) / cnt_yz(iy,iz)
          data_yz(u2_avg,iy,iz) = data_yz(u2_avg,iy,iz) / cnt_yz(iy,iz)
          data_yz(u3_avg,iy,iz) = data_yz(u3_avg,iy,iz) / cnt_yz(iy,iz)

          u2 = data_yz(u1_avg,iy,iz)**2 + data_yz(u2_avg,iy,iz)**2 + data_yz(u3_avg,iy,iz)**2
          v2mean = (vx2_yz(iy,iz) + vy2_yz(iy,iz) + vz2_yz(iy,iz)) / cnt_yz(iy,iz)
          thermal_v2 = max(0.0_real64, v2mean - u2)
          data_yz(Tp_avg,iy,iz) = data_yz(Tp_avg,iy,iz) + &
                  mass_species * thermal_v2 / (3.0_real64 * qe)
        end if
      end do
    end do

    deallocate(cnt_xy, cnt_xz, cnt_yz)
    deallocate(vx2_xy, vy2_xy, vz2_xy)
    deallocate(vx2_xz, vy2_xz, vz2_xz)
    deallocate(vx2_yz, vy2_yz, vz2_yz)

  end subroutine compute_particle_plane_moments_species
end module mod_planeMoments