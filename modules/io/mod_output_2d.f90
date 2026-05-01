module mod_output_2d
  use iso_fortran_env, only: real64, int32
  use mod_particles,   only: ParticleSet
  use mod_constants,   only: qe

  implicit none
  private

  public :: write_plane_xy_scalar, write_plane_xz_scalar, write_plane_yz_scalar
  public :: write_plane_xy_scalar_2d, write_plane_xz_scalar_2d, write_plane_yz_scalar_2d
  public :: write_density_planes, write_scalar_planes

contains

  ! =========================
  ! 3D FIELD WRITERS (n+1)
  ! =========================

  subroutine write_plane_xy_scalar(filename, f, n, iz_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: iz_plane
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    integer :: ix, iy, u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) n(1)/every, n(2)/every

    do iy = n(2)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(ix,iy,iz_plane), ix=1,n(1)+1,every )
    end do

    close(u)
  end subroutine


  subroutine write_plane_xz_scalar(filename, f, n, iy_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: iy_plane
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    integer :: ix, iz, u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) n(1)/every, n(3)/every

    do iz = n(3)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(ix,iy_plane,iz), ix=1,n(1)+1,every )
    end do

    close(u)
  end subroutine


  subroutine write_plane_yz_scalar(filename, f, n, ix_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: ix_plane
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    integer :: iy, iz, u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) n(2)/every, n(3)/every

    do iz = n(3)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(ix_plane,iy,iz), iy=1,n(2)+1,every )
    end do

    close(u)
  end subroutine


  ! =========================
  ! WRAPPERS
  ! =========================

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
  end subroutine


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
  end subroutine


  ! =========================
  ! 2D DATA WRITERS (n+1)
  ! =========================

  subroutine write_plane_xy_scalar_2d(filename, f, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2)

    integer :: ix, iy, u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) n(1)/every, n(2)/every

    do iy = n(2)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(ix,iy), ix=1,n(1)+1,every )
    end do

    close(u)
  end subroutine


  subroutine write_plane_xz_scalar_2d(filename, f, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(3)+2)

    integer :: ix, iz, u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) n(1)/every, n(3)/every

    do iz = n(3)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(ix,iz), ix=1,n(1)+1,every )
    end do

    close(u)
  end subroutine


  subroutine write_plane_yz_scalar_2d(filename, f, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3)
    integer(int32),   intent(in) :: every
    real(real64),     intent(in) :: f(0:n(2)+2,0:n(3)+2)

    integer :: iy, iz, u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) n(2)/every, n(3)/every

    do iz = n(3)+1, 1, -every
      write(u,'(800(es18.10,1x))') ( f(iy,iz), iy=1,n(2)+1,every )
    end do

    close(u)
  end subroutine


end module mod_output_2d