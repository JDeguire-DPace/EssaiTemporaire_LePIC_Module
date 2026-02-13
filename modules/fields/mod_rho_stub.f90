module mod_rho_stub
  use iso_fortran_env, only: real64, int32
  implicit none
  private
  public :: build_rhs_gaussian

contains

  subroutine build_rhs_gaussian(dom_n, h, k0, m, rhs_loc)
    ! Build rhs_loc(0:nx+1,0:ny+1,0:m+1) from a Gaussian rho.
    !
    ! NOTE: legacy pdesolver expects b on (0:nx+1,0:ny+1,0:m+1).
    ! We build it in global coordinates using k0 offset.

    integer(int32), intent(in) :: dom_n(3)
    real(real64),   intent(in) :: h(3)
    integer(int32), intent(in) :: k0, m
    real(real64),   intent(out):: rhs_loc(0:dom_n(1)+1,0:dom_n(2)+1,0:m+1)

    integer :: i,j,kl
    integer(int32) :: kg
    real(real64) :: x,y,z, xc,yc,zc, sig, r2
    real(real64) :: amp

    rhs_loc = 0.0_real64

    ! Center (in cell-centered coordinates). Use physical coordinates.
    xc = 0.5_real64 * dom_n(1) * h(1)
    yc = 0.5_real64 * dom_n(2) * h(2)
    zc = 0.5_real64 * dom_n(3) * h(3)

    sig = 0.08_real64 * min(dom_n(1)*h(1), min(dom_n(2)*h(2), dom_n(3)*h(3)))
    amp = 1.0_real64   ! scale this; sign controls polarity

    do kl = 0, m+1
      kg = k0 + int(kl, int32)   ! global k index in [k0, k0+m+1]
      z  = real(kg, real64) * h(3)
      do j = 0, dom_n(2)+1
        y = real(j, real64) * h(2)
        do i = 0, dom_n(1)+1
          x = real(i, real64) * h(1)
          r2 = (x-xc)**2 + (y-yc)**2 + (z-zc)**2
          rhs_loc(i,j,kl) = amp * exp(-r2/(2.0_real64*sig*sig))
        end do
      end do
    end do
    rhs_loc = 0.0_real64  ! override with zero for now, to test homogeneous BCs and solver convergence without a source

    ! If you want "Poisson form" consistent with ∇²phi = -rho/eps0,
    ! you can choose rhs_loc = -rho/eps0 here. For now: just a source.
  end subroutine build_rhs_gaussian

end module mod_rho_stub
