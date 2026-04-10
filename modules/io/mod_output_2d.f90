module mod_output_2d
  use iso_fortran_env, only: real64, int32
  use mod_particles,   only: ParticleSet
  use mod_constants,   only: qe
  implicit none
  private

  public :: write_plane_xy_scalar, write_plane_xz_scalar, write_plane_yz_scalar
  public :: write_plane_xy_scalar_2d, write_plane_xz_scalar_2d, write_plane_yz_scalar_2d
  public :: write_density_planes, write_scalar_planes
  public :: write_temperature_xz_particles
  public :: write_temperature_planes_particles

contains

  subroutine write_plane_xy_scalar(filename, f, n, iz_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: iz_plane
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer :: ix, iy, u

    u = 21
    open(u, file=filename, status='replace', action='write')
    write(u,*) n(1)/every, n(2)/every
    do iy = n(2)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(ix,iy,iz_plane), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xy_scalar


  subroutine write_plane_xz_scalar(filename, f, n, iy_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: iy_plane
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer :: ix, iz, u

    u = 22
    open(u, file=filename, status='replace', action='write')
    write(u,*) n(1)/every, n(3)/every
    do iz = n(3)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(ix,iy_plane,iz), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xz_scalar


  subroutine write_plane_yz_scalar(filename, f, n, ix_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: ix_plane
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer :: iy, iz, u

    u = 23
    open(u, file=filename, status='replace', action='write')
    write(u,*) n(2)/every, n(3)/every
    do iz = n(3)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(ix_plane,iy,iz), iy=1,n(2)+1,every )
    end do
    close(u)
  end subroutine write_plane_yz_scalar


  subroutine write_density_planes(np, n, ptype, ix_plane, iy_plane, iz_plane, every, prefix)
    character(len=*), intent(in) :: prefix
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: ptype
    integer(int32),   intent(in) :: ix_plane, iy_plane, iz_plane
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: np(:,:,:,:)

    real(real64), allocatable :: tmp(:,:,:)
    character(len=256) :: fxy, fxz, fyz

    allocate(tmp(0:n(1)+2,0:n(2)+2,0:n(3)+2))
    tmp = np(:,:,:,ptype)

    write(fxy,'(a,"_xy.mco")') trim(prefix)
    write(fxz,'(a,"_xz.mco")') trim(prefix)
    write(fyz,'(a,"_yz.mco")') trim(prefix)

    call write_plane_xy_scalar(fxy, tmp, n, iz_plane, every)
    call write_plane_xz_scalar(fxz, tmp, n, iy_plane, every)
    call write_plane_yz_scalar(fyz, tmp, n, ix_plane, every)

    deallocate(tmp)
  end subroutine write_density_planes


  subroutine write_scalar_planes(f, n, ix_plane, iy_plane, iz_plane, every, prefix)
    character(len=*), intent(in) :: prefix
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: ix_plane, iy_plane, iz_plane
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    character(len=256) :: fxy, fxz, fyz

    write(fxy,'(a,"_xy.mco")') trim(prefix)
    write(fxz,'(a,"_xz.mco")') trim(prefix)
    write(fyz,'(a,"_yz.mco")') trim(prefix)

    call write_plane_xy_scalar(fxy, f, n, iz_plane, every)
    call write_plane_xz_scalar(fxz, f, n, iy_plane, every)
    call write_plane_yz_scalar(fyz, f, n, ix_plane, every)
  end subroutine write_scalar_planes


  subroutine write_temperature_xz_particles(filename, part, nproc, h, n, iy_plane, every, mass_species)
    character(len=*), intent(in) :: filename
    type(ParticleSet), intent(in) :: part(:)
    integer(int32),    intent(in) :: nproc
    integer(int32),    intent(in) :: n(3)
    integer(int32),    intent(in) :: iy_plane
    integer(int32),    intent(in) :: every
    real(real64),      intent(in) :: h(3)
    real(real64),      intent(in) :: mass_species

    integer(int32) :: iproc, i
    integer(int32) :: ix, iy, iz
    integer :: u
    real(real64) :: x, y, z, vx, vy, vz
    real(real64) :: u2, v2mean, thermal_v2
    real(real64), allocatable :: cnt(:,:), sx(:,:), sy(:,:), sz(:,:), sv2(:,:), T(:,:)

    allocate(cnt(0:n(1)+2,0:n(3)+2))
    allocate(sx (0:n(1)+2,0:n(3)+2))
    allocate(sy (0:n(1)+2,0:n(3)+2))
    allocate(sz (0:n(1)+2,0:n(3)+2))
    allocate(sv2(0:n(1)+2,0:n(3)+2))
    allocate(T  (0:n(1)+2,0:n(3)+2))

    cnt = 0.0_real64
    sx  = 0.0_real64
    sy  = 0.0_real64
    sz  = 0.0_real64
    sv2 = 0.0_real64
    T   = 0.0_real64

    do iproc = 1, nproc
      if (.not. allocated(part(iproc)%x)) cycle
      if (part(iproc)%n <= 0_int32) cycle

      do i = 1, part(iproc)%n
        x = part(iproc)%x(i)
        y = part(iproc)%y(i)
        z = part(iproc)%z(i)

        ix = int(x / h(1), int32) + 1_int32
        iy = int(y / h(2), int32) + 1_int32
        iz = int(z / h(3), int32) + 1_int32

        if (ix < 1_int32 .or. ix > n(1)+1_int32) cycle
        if (iy /= iy_plane .and. iy /= iy_plane-1_int32) cycle
        if (iz < 1_int32 .or. iz > n(3)+1_int32) cycle

        vx = part(iproc)%vx(i)
        vy = part(iproc)%vy(i)
        vz = part(iproc)%vz(i)

        cnt(ix,iz) = cnt(ix,iz) + 1.0_real64
        sx(ix,iz)  = sx(ix,iz)  + vx
        sy(ix,iz)  = sy(ix,iz)  + vy
        sz(ix,iz)  = sz(ix,iz)  + vz
        sv2(ix,iz) = sv2(ix,iz) + (vx*vx + vy*vy + vz*vz)
      end do
    end do

    do iz = 1, n(3)+1
      do ix = 1, n(1)+1
        if (cnt(ix,iz) > 0.0_real64) then
          u2 = (sx(ix,iz)/cnt(ix,iz))**2 + &
               (sy(ix,iz)/cnt(ix,iz))**2 + &
               (sz(ix,iz)/cnt(ix,iz))**2

          v2mean = sv2(ix,iz) / cnt(ix,iz)
          thermal_v2 = max(0.0_real64, v2mean - u2)

          T(ix,iz) = mass_species * thermal_v2 / (3.0_real64 * qe)
        end if
      end do
    end do

    u = 24
    open(u, file=filename, status='replace', action='write')
    write(u,*) n(1)/every, n(3)/every
    do iz = n(3)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( T(ix,iz), ix=1,n(1)+1,every )
    end do
    close(u)

    deallocate(cnt, sx, sy, sz, sv2, T)
  end subroutine write_temperature_xz_particles


  subroutine write_temperature_planes_particles(part, nproc, n, h, ix_plane, iy_plane, iz_plane, every, mass_species, prefix)
    type(ParticleSet), intent(in) :: part(:)
    integer(int32),    intent(in) :: nproc
    integer(int32),    intent(in) :: n(3)
    integer(int32),    intent(in) :: ix_plane, iy_plane, iz_plane
    integer(int32),    intent(in) :: every
    real(real64),      intent(in) :: h(3)
    real(real64),      intent(in) :: mass_species
    character(len=*),  intent(in) :: prefix

    character(len=256) :: fxz

    ! First concrete step: XZ temperature plane
    write(fxz,'(a,"_xz.mco")') trim(prefix)
    call write_temperature_xz_particles(fxz, part, nproc, h, n, iy_plane, every, mass_species)

    ! XY and YZ can be added next with the same pattern.
  end subroutine write_temperature_planes_particles

  subroutine write_plane_xy_scalar_2d(filename, f, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2)
    integer :: ix, iy, u

    u = 31
    open(u, file=filename, status='replace', action='write')
    write(u,*) n(1)/every, n(2)/every
    do iy = n(2)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(ix,iy), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xy_scalar_2d


  subroutine write_plane_xz_scalar_2d(filename, f, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(3)+2)
    integer :: ix, iz, u

    u = 32
    open(u, file=filename, status='replace', action='write')
    write(u,*) n(1)/every, n(3)/every
    do iz = n(3)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(ix,iz), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xz_scalar_2d


  subroutine write_plane_yz_scalar_2d(filename, f, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(2)+2,0:n(3)+2)
    integer :: iy, iz, u

    u = 33
    open(u, file=filename, status='replace', action='write')
    write(u,*) n(2)/every, n(3)/every
    do iz = n(3)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(iy,iz), iy=1,n(2)+1,every )
    end do
    close(u)
  end subroutine write_plane_yz_scalar_2d

  



end module mod_output_2d