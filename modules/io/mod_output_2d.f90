module mod_output_2d
  use iso_fortran_env, only: real64, int32
  implicit none
  private
  public :: write_plane_xy_scalar, write_plane_xz_scalar, write_plane_yz_scalar
  public :: write_density_planes

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

end module mod_output_2d