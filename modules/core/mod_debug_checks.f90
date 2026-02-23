module mod_debug_checks
  use iso_fortran_env, only: real64, int32
  use mpi
  use mod_poisson_decomp, only: PoissonDecomp
  implicit none
  private
  public :: checkpoint_poisson_decomp, checkpoint_kq

contains

  subroutine checkpoint_poisson_decomp(rank, comm, nz, k0, m, phi_dom, bcnd_dom)
    integer, intent(in) :: rank, comm
    integer(int32), intent(in) :: nz, k0, m
    real(real64), intent(in) :: phi_dom(:,:,:)
    integer(int32), intent(in) :: bcnd_dom(:,:,:)

    integer :: ierr
    integer(int32) :: zlo, zhi
    real(real64) :: pmin, pmax
    integer(int32) :: bmin, bmax

    zlo = k0
    zhi = k0 + m + 2

    if (zlo < 0_int32 .or. zhi > nz+2_int32) then
      write(*,'(a,i0,a,i0,a,i0,a,i0)') "BAD SLICE rank=", rank, " k0=", k0, " m=", m, " z=[", zlo, ",", zhi, "]"
      call MPI_Abort(comm, 911, ierr)
    end if

    pmin = minval(phi_dom); pmax = maxval(phi_dom)
    bmin = minval(bcnd_dom); bmax = maxval(bcnd_dom)
  end subroutine checkpoint_poisson_decomp


  subroutine checkpoint_kq(kq, comm, rank, label, pdec)
    use iso_fortran_env, only: real64, int32
    use mpi
    use mod_poisson_decomp, only: PoissonDecomp
    implicit none

    real(real64), intent(in) :: kq(0:,0:,0:)   ! <--- force lb=0
    integer, intent(in) :: comm, rank
    character(*), intent(in) :: label
    type(PoissonDecomp), intent(in) :: pdec

    integer(int32) :: zlo, zhi
    integer :: ierr
    integer :: c1, c2, g1, g2
    integer :: nplanes, gplanes
    integer :: nz_kq_lo, nz_kq_hi

    nz_kq_lo = lbound(kq,3)
    nz_kq_hi = ubound(kq,3)

    ! Global z-slab owned by this rank in legacy gather sense:
    zlo = pdec%k0 + pdec%kl
    zhi = pdec%k0 + pdec%kr

    nplanes = zhi - zlo + 1

    call MPI_Reduce(nplanes, gplanes, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    c1 = count(kq(:,:,zlo:zhi) == 1.0_real64)
    c2 = count(kq(:,:,zlo:zhi) == 2.0_real64)

    call MPI_Reduce(c1, g1, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
    call MPI_Reduce(c2, g2, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  end subroutine checkpoint_kq

end module mod_debug_checks
