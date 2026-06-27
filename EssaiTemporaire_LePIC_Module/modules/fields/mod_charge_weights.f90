module mod_charge_weights
  use iso_fortran_env, only: real64, int32
  implicit none
  private
  public :: build_kq

contains

  subroutine build_kq(bcnd, kq)
    ! Accept whatever bcnd kind Domain uses (default integer),
    ! and cast inside if needed.
    integer,        intent(in)  :: bcnd(:,:,:)
    real(real64),   intent(out) :: kq(:,:,:)

    integer :: ix, iy, iz
    integer :: nx, ny, nz

    nx = ubound(kq,1)
    ny = ubound(kq,2)
    nz = ubound(kq,3)

    kq = 1.0_real64

    !$omp parallel do collapse(3) private(ix,iy,iz)
    do iz = lbound(kq,3), nz
      do iy = lbound(kq,2), ny
        do ix = lbound(kq,1), nx
          if (bcnd(ix,iy,iz) /= -1) kq(ix,iy,iz) = 2.0_real64
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine build_kq

end module mod_charge_weights
