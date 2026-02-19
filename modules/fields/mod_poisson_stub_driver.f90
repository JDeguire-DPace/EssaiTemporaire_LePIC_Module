module mod_poisson_stub_driver
  use iso_fortran_env, only: real64
  use mpi
  use mod_poisson_decomp, only: PoissonDecomp
  use mod_constants,    only: eps0, pi
  implicit none
  private
  public :: run_poisson_stub, write_mco_real_xy

  ! Give the compiler the interface of the legacy routine
  interface
    subroutine pdesolver(u,b,bcnd,h,n,ncycl,eps,omega,k,ktot,res,ng,rank,nproc)
      use iso_fortran_env, only: real64
      implicit none
      integer, intent(in) :: n(3), ng, ncycl, rank, nproc
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

  subroutine run_poisson_stub(comm, n, h, bcnd_global, phi_global)
    use iso_fortran_env,  only: real64
    use mpi
    use ieee_arithmetic,  only: ieee_is_nan
    use mod_poisson_decomp, only: PoissonDecomp
    use mod_part_info,      only: flag_pbc, flag_nmn
    implicit none

    integer, intent(in) :: comm
    integer, intent(in) :: n(3)
    real(real64), intent(in) :: h(3)
    integer, intent(in) :: bcnd_global(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64), intent(inout) :: phi_global(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    type(PoissonDecomp) :: pdec
    integer :: rank, nproc, ierr
    integer :: m
    real(real64), allocatable :: b_loc(:,:,:)
    integer,      allocatable :: bcnd_loc(:,:,:)
    real(real64), allocatable :: u_mg(:,:,:)

    integer :: ng, ncycl, k
    real(real64) :: eps, omega, ktot, res
    real(real64) :: h_mg(3)

    real(real64) :: local_sum, global_sum
    integer :: local_n, global_n
    integer :: kk, nan_k, ix, iy, iz
    integer :: flag_loc, flag_glob
    logical :: force_driver_pbc_ghosts

    ! ------------------------------------------------------------------
    ! IMPORTANT: legacy mg/sors uses MPI_COMM_WORLD internally.
    ! So for now treat comm as WORLD-consistent.
    ! ------------------------------------------------------------------
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    h_mg = h
    m    = n(3) / nproc

    allocate(u_mg(0:n(1)+2, 0:n(2)+2, -1:m+2))
    allocate(b_loc(0:n(1)+1, 0:n(2)+1,  0:m+1))
    allocate(bcnd_loc(0:n(1)+2,0:n(2)+2, 0:m+2))

    u_mg     = 0.0_real64
    b_loc    = 0.0_real64
    bcnd_loc = 0

    ! ------------------------------------------------------------------
    ! Legacy solver dependency: these flags live in mod_part_info
    ! Set them explicitly for this test.
    ! ------------------------------------------------------------------
    ! If your geometry implies periodic:
    flag_pbc = 1
    ! If you know flag_nmn should be 0/1, set it deterministically:
    ! flag_nmn = 0
    ! (leave as-is if you already manage it elsewhere)

    ! Local decomp + scatter
    call pdec%init(n(1), n(2), n(3), MPI_COMM_WORLD)
    call pdec%scatter_from_global(phi_global, bcnd_global)

    u_mg(:,:,0:m+2)     = pdec%phi_dom(:,:,0:m+2)
    bcnd_loc(:,:,0:m+2) = pdec%bcnd_dom(:,:,0:m+2)

    call build_rhs_sin(n, h, pdec%k0, m, b_loc, bcnd_loc)


    if (rank == 0) then
      if (rank == 0) then
        open(12,file='../DATA/DATA_2D/rho_xz.mco')
        102   format(800(ES12.4,1x))
          write(12,*) n(1), n(3)
          
          do iz=n(3)+1,1,-1
            
            write(12,102) (b_loc(ix,n(2)/2+1,iz), ix=1,n(1)+1,1)
          end do
          close(12)
      end if
    end if

    ! Diagnostics: RHS sum
    local_sum = sum(b_loc(1:n(1)+1, 1:n(2)+1, 1:m))
    call MPI_Allreduce(local_sum, global_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
      print *, "Global RHS sum = ", global_sum
      print *, "pre:  min/max u_mg= ", minval(u_mg), maxval(u_mg)
    end if

    ! ------------------------------------------------------------------
    ! Toggle: FORCE periodic ghost planes in the driver as a diagnostic.
    ! If this removes NaNs, the bug is in legacy periodic handling (flag/comm/gating).
    ! ------------------------------------------------------------------
    force_driver_pbc_ghosts = .true.
    if (force_driver_pbc_ghosts) then
      u_mg(:,:,0)   = u_mg(:,:,m)
      u_mg(:,:,m+1) = u_mg(:,:,1)
    end if


    ! Solver params
    ng    = 6
    ncycl = 100
    eps   = 1.0e-7_real64
    omega = 1.7_real64
    ktot  = 0.0_real64
    res   = 0.0_real64
    k     = 0
    
    call pdesolver(u_mg, b_loc, bcnd_loc, h_mg, n, ncycl, eps, omega, k, ktot, res, ng, rank, nproc)
    
    ! Gauge fix for periodic (remove mean)
    local_sum = sum(u_mg(1:n(1)+1, 1:n(2)+1, 1:m))
    local_n   = (n(1)+1)*(n(2)+1)*m
    call MPI_Allreduce(local_sum, global_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_n,   global_n,   1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! if (global_n > 0) then
    !   u_mg(1:n(1)+1, 1:n(2)+1, 1:m) = u_mg(1:n(1)+1, 1:n(2)+1, 1:m) - global_sum/real(global_n, real64)
    ! end if

    if (rank == 0) then
      print *, "post: min/max u_mg= ", minval(u_mg), maxval(u_mg)
    end if

    ! NaN check
    flag_loc = merge(1, 0, any(ieee_is_nan(u_mg)))
    call MPI_Allreduce(flag_loc, flag_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (rank == 0) print *, "NaN present after pdesolver? ", (flag_glob == 1)

    do kk = lbound(u_mg,3), ubound(u_mg,3)
      nan_k = count(ieee_is_nan(u_mg(:,:,kk)))
      if (nan_k > 0) then
        write(*,'(A,I3,A,I5,A,I10)') "rank ", rank, " NaNs in k=", kk, " : ", nan_k
      end if
    end do
    write(*,*) "rank",rank," total NaNs =", count(ieee_is_nan(u_mg))

    ! Back to decomp + gather
    pdec%phi_dom(:,:,lbound(pdec%phi_dom,3):ubound(pdec%phi_dom,3)) = &
        u_mg(:,:,lbound(pdec%phi_dom,3):ubound(pdec%phi_dom,3))
    call pdec%gather_phi_to_global(phi_global)

    call pdec%destroy()
    deallocate(u_mg, b_loc, bcnd_loc)
  end subroutine run_poisson_stub



  subroutine build_rhs_gaussian(n, h, k0, m, b_loc, bcnd_loc)
    integer, intent(in)          :: n(3)
    real(real64), intent(in)     :: h(3)
    integer, intent(in)          :: k0, m
    real(real64), intent(inout)  :: b_loc(0:n(1)+1,0:n(2)+1,0:m+1)
    integer, intent(in)          :: bcnd_loc(0:n(1)+2,0:n(2)+2,0:m+2)

    integer :: i,j,k
    real(real64) :: x,y,z
    real(real64) :: x0,y0,z0, sig2, amp

    x0 = 0.5_real64 * (n(1)*h(1))
    y0 = 0.5_real64 * (n(2)*h(2))
    z0 = 0.5_real64 * (n(3)*h(3))

    sig2 = (5.0e-2_real64/4)**2
    amp  = 1.0e-12_real64

    b_loc = 1.0e-6_real64! - 10.0_real64/(3.0_real64)   ! scale factor to get from a physical charge density (C/m^3) to the appropriate RHS units for the Poisson equation in V/m^2, assuming a grid spacing of order 1e-2 m. Adjust as needed for your specific grid and physical scenario.

    do k=0,m+1
      z = real(k0 + k, real64) * h(3)
      do j=0,n(2)+1
        y = real(j, real64) * h(2)
        do i=0,n(1)+1
          x = real(i, real64) * h(1)
          if (bcnd_loc(i,j,k) <= 0) then
            b_loc(i,j,k) = amp*exp(-((x-h(1))-x0)**2/(2*sig2)) * &
                    (x0**2-2.0_real64*x0*(x-h(1))-sig2+(x-h(1))**2)/(sqrt(2*3.141592*sig2**5))
          end if
        end do
      end do
    end do
    !b_loc = 0.0_real64  ! override with zero for now, to test homogeneous BCs and solver convergence without a source
    write(*,*) "Value of the sum", sum(b_loc)
  end subroutine build_rhs_gaussian

  subroutine build_rhs_sin(n, h, k0, m, b_loc, bcnd_loc)
    integer, intent(in)          :: n(3)
    real(real64), intent(in)     :: h(3)
    integer, intent(in)          :: k0, m
    real(real64), intent(inout)  :: b_loc(0:n(1)+1,0:n(2)+1,0:m+1)
    integer, intent(in)          :: bcnd_loc(0:n(1)+2,0:n(2)+2,0:m+2)

    integer :: i,j,k
    real(real64) :: x,y,z
    real(real64) :: x0,y0,z0, sig2, amp,L

    x0 = 0.5_real64 * (n(1)*h(1))
    y0 = 0.5_real64 * (n(2)*h(2))
    z0 = 0.5_real64 * (n(3)*h(3))

    
    amp  = 1.0e-6_real64
    L = 0.10_real64

    b_loc = 0.0_real64
    
    do k=0,m+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          x = real(i, real64) * h(1)
          if (bcnd_loc(i,j,k) <= 0) then
            b_loc(i,j,k) = -(4.0_real64*pi/L)**2 * sin(4.0_real64*pi*(x-h(1))/L)*eps0
          end if
        end do
      end do
    end do


  end subroutine build_rhs_sin

  subroutine build_rhs_constant(n, h, k0, m, b_loc, bcnd_loc)
    integer, intent(in)          :: n(3)
    real(real64), intent(in)     :: h(3)
    integer, intent(in)          :: k0, m
    real(real64), intent(inout)  :: b_loc(0:n(1)+1,0:n(2)+1,0:m+1)
    integer, intent(in)          :: bcnd_loc(0:n(1)+2,0:n(2)+2,0:m+2)

    integer :: i,j,k
    real(real64) :: x,y,z
    real(real64) :: x0,y0,z0, sig2, amp

    x0 = 0.5_real64 * (n(1)*h(1))
    y0 = 0.5_real64 * (n(2)*h(2))
    z0 = 0.5_real64 * (n(3)*h(3))

    b_loc = 1.0e-7_real64 
  end subroutine build_rhs_constant

  subroutine write_mco_real_xy(fname, a, nx, ny, kslice, rank)
    use iso_fortran_env, only: real64
    implicit none
    character(len=*), intent(in) :: fname
    real(real64),     intent(in) :: a(0:nx+2,0:ny+2,*)
    integer,          intent(in) :: nx, ny, kslice, rank
    integer :: i, j, u

    if (rank /= 0) return

    open(newunit=u, file=fname, status="replace", action="write", form="formatted")

    ! header consistent with your legacy-style 2D MCO slices (interior 1..nx+1, 1..ny+1)
    write(u,'(2I8)') nx+1, ny+1

    do j = ny+1, 1, -1
      do i = 1, nx+1
        write(u,'(ES24.16)', advance='no') a(i,j,kslice)
      end do
      write(u,*)
    end do

    close(u)
  end subroutine write_mco_real_xy

end module mod_poisson_stub_driver
