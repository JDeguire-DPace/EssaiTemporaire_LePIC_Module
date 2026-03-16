module mod_poisson_test_driver
  use iso_fortran_env, only: real64, int32
  use mpi
  use mod_poisson_decomp, only: PoissonDecomp
  implicit none
  private
  public :: run_poisson_actual, write_mco_real_xy

  interface
    subroutine pdesolver(u,b,bcnd,h,n,ncycl,eps,omega,k,ktot,res,ng,rank,nproc)
      use iso_fortran_env, only: real64
      implicit none
      integer, intent(in)    :: n(3), ng, ncycl, rank, nproc
      integer, intent(inout) :: k
      real(real64), intent(inout) :: h(3)
      real(real64), intent(inout) :: u(0:n(1)+2,0:n(2)+2,-1:n(3)/nproc+2)
      real(real64), intent(inout) :: b(0:n(1)+1,0:n(2)+1,0:n(3)/nproc+1)
      integer,      intent(inout) :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)/nproc+2)
      real(real64), intent(in)    :: eps, omega
      real(real64), intent(inout) :: ktot, res
    end subroutine pdesolver
  end interface

contains

  subroutine run_poisson_actual(comm, n_in, h, bcnd_global, rhs_global, phi_global, &
                                ncycl, eps, omega, ng)
    use ieee_arithmetic, only: ieee_is_nan
    implicit none

    integer, intent(in) :: comm
    integer(int32), intent(in) :: n_in(3)
    real(real64), intent(in) :: h(3)
    integer, intent(in) :: bcnd_global(0:,0:,0:)
    real(real64), intent(in) :: rhs_global(0:,0:,0:)
    real(real64), intent(inout) :: phi_global(0:,0:,0:)

    integer, intent(in) :: ncycl, ng
    real(real64), intent(in) :: eps, omega

    type(PoissonDecomp) :: pdec

    integer :: rank, nproc, ierr
    integer(int32) :: m
    integer :: k_it
    real(real64) :: ktot, res
    real(real64) :: h_mg(3)
    integer :: n_leg(3)

    real(real64), allocatable :: u_mg(:,:,:)
    real(real64), allocatable :: b_loc(:,:,:)
    integer,      allocatable :: bcnd_loc(:,:,:)

    real(real64) :: local_sum, global_sum
    integer :: kk, nan_k, flag_loc, flag_glob

    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, nproc, ierr)

    h_mg = h
    n_leg = int(n_in, kind=4)

    call pdec%init(n_in(1), n_in(2), n_in(3), comm)
    call pdec%scatter_from_global(phi_global, bcnd_global)
    call pdec%scatter_rhs_from_global(rhs_global)

    m = pdec%m

    allocate(u_mg(0:n_leg(1)+2, 0:n_leg(2)+2, -1:m+2))
    allocate(b_loc(0:n_leg(1)+1, 0:n_leg(2)+1,  0:m+1))
    allocate(bcnd_loc(0:n_leg(1)+2,0:n_leg(2)+2, 0:m+2))

    u_mg     = 0.0_real64
    b_loc    = 0.0_real64
    bcnd_loc = 0

    u_mg(:,:,0:m+2)     = pdec%phi_dom(:,:,0:m+2)
    bcnd_loc(:,:,0:m+2) = pdec%bcnd_dom(:,:,0:m+2)
    b_loc(:,:,:)        = pdec%rhs_dom(:,:,:)

    if (rank == 0) then
      write(*,*) ' '
      write(*,*) 'DEBUG POISSON TEST DRIVER rhs'
      write(*,'(a,es16.8)') 'sum      = ', sum(rhs_global)
      write(*,'(a,es16.8)') 'sum(abs) = ', sum(abs(rhs_global))
      write(*,'(a,es16.8)') 'max      = ', maxval(rhs_global)
      write(*,'(a,es16.8)') 'min      = ', minval(rhs_global)
    end if

    local_sum = sum(b_loc(1:n_leg(1)+1, 1:n_leg(2)+1, 1:m))
    call MPI_Allreduce(local_sum, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

    if (rank == 0) then
      write(*,'(a,es16.8)') 'Global RHS slab sum = ', global_sum
      write(*,'(a,es16.8,1x,es16.8)') 'pre: min/max u_mg = ', minval(u_mg), maxval(u_mg)
    end if

    k_it = 0
    ktot = 0.0_real64
    res  = 0.0_real64

    call pdesolver(u_mg, b_loc, bcnd_loc, h_mg, n_leg, ncycl, eps, omega, &
                   k_it, ktot, res, ng, rank, nproc)

    if (rank == 0) then
      write(*,'(a,es16.8,1x,es16.8)') 'post: min/max u_mg = ', minval(u_mg), maxval(u_mg)
    end if

    flag_loc = merge(1, 0, any(ieee_is_nan(u_mg)))
    call MPI_Allreduce(flag_loc, flag_glob, 1, MPI_INTEGER, MPI_MAX, comm, ierr)

    if (rank == 0) then
      write(*,'(a,l1)') 'NaN present after pdesolver? ', (flag_glob == 1)
    end if

    do kk = lbound(u_mg,3), ubound(u_mg,3)
      nan_k = count(ieee_is_nan(u_mg(:,:,kk)))
      if (nan_k > 0) then
        write(*,'(A,I3,A,I5,A,I10)') 'rank ', rank, ' NaNs in k=', kk, ' : ', nan_k
      end if
    end do

    pdec%phi_dom(:,:,0:m+2) = u_mg(:,:,0:m+2)
    call pdec%gather_phi_to_global(phi_global)

    call pdec%destroy()
    deallocate(u_mg, b_loc, bcnd_loc)
  end subroutine run_poisson_actual


  subroutine write_mco_real_xy(fname, a, nx, ny, kslice, rank)
    implicit none
    character(len=*), intent(in) :: fname
    real(real64),     intent(in) :: a(0:nx+2,0:ny+2,*)
    integer,          intent(in) :: nx, ny, kslice, rank
    integer :: i, j, u

    if (rank /= 0) return

    open(newunit=u, file=fname, status='replace', action='write', form='formatted')
    write(u,'(2I8)') nx+1, ny+1

    do j = ny+1, 1, -1
      do i = 1, nx+1
        write(u,'(ES24.16)', advance='no') a(i,j,kslice)
      end do
      write(u,*)
    end do

    close(u)
  end subroutine write_mco_real_xy

end module mod_poisson_test_driver
