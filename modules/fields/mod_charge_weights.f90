module mod_charge_weights
  use iso_fortran_env, only: real64, int32
  implicit none
  private
  public :: build_kq

contains

  subroutine build_kq(bcnd, kq)
    integer(int32), intent(in)  :: bcnd(:,:,:)
    real(real64),   intent(out) :: kq(:,:,:)
    integer :: ix, iy, iz
    integer :: ixu, iyu, izu

    if (any(lbound(bcnd) /= lbound(kq)) .or. any(ubound(bcnd) /= ubound(kq))) then
      error stop "build_kq: bcnd and kq bounds mismatch"
    end if

    ixu = ubound(kq,1); iyu = ubound(kq,2); izu = ubound(kq,3)

    kq = 1.0_real64
    !$omp parallel do collapse(3) private(ix,iy,iz)
    do iz = lbound(kq,3), izu
      do iy = lbound(kq,2), iyu
        do ix = lbound(kq,1), ixu
          if (bcnd(ix,iy,iz) /= -1_int32) kq(ix,iy,iz) = 2.0_real64
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine

end module mod_charge_weights
