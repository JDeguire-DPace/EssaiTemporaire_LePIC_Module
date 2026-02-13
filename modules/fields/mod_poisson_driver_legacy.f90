module mod_poisson_driver_legacy
  use iso_fortran_env, only: real64, int32
  use mpi
  use mod_poisson_decomp, only: PoissonDecomp
  use mod_rho_stub,       only: build_rhs_gaussian
  implicit none
  private
  public :: solve_poisson_legacy_stub

  interface
    subroutine pdesolver(u,b,bcnd,h,n,ncycl,eps,omega,k,ktot,res,ng,rank,nproc)
      implicit none
      integer :: k, n(3), ng, ncycl, rank, nproc
      integer :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)/nproc+2)
      real(kind=8) :: h(3)
      real(kind=8) :: u(0:n(1)+2,0:n(2)+2,-1:n(3)/nproc+2)
      real(kind=8) :: b(0:n(1)+1,0:n(2)+1,0:n(3)/nproc+1)
      real(kind=8) :: eps, omega, ktot, res
    end subroutine pdesolver
  end interface

contains

  subroutine solve_poisson_legacy_stub(pdec, phi_global, bcnd_global, h, n, &
                                      ncycl, eps, omega, ng)
    type(PoissonDecomp), intent(inout) :: pdec
    real(real64), intent(inout)        :: phi_global(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)         :: bcnd_global(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64), intent(in)           :: h(3)
    integer(int32), intent(in)         :: n(3)
    integer, intent(in)                :: ncycl, ng
    real(real64), intent(in)           :: eps, omega

    integer :: k_it
    real(real64) :: ktot, res

    integer(int32) :: m
    real(real64), allocatable :: u_mg(:,:,:)
    real(real64), allocatable :: b_loc(:,:,:)
    integer,      allocatable :: bcnd_loc(:,:,:)

    ! Scatter global phi/bcnd into pdec local buffers
    call pdec%scatter_from_global(phi_global, bcnd_global)

    m = pdec%m

    ! Allocate legacy-shaped local arrays
    allocate(u_mg(0:n(1)+2, 0:n(2)+2, -1:m+2))
    allocate(b_loc(0:n(1)+1, 0:n(2)+1,  0:m+1))
    allocate(bcnd_loc(0:n(1)+2,0:n(2)+2, 0:m+2))

    ! Map decomposition arrays into MG arrays
    u_mg(:,:, -1)    = 0.0_real64
    u_mg(:,:, 0:m+2) = pdec%phi_dom(:,:,0:m+2)

    bcnd_loc(:,:,0:m+2) = pdec%bcnd_dom(:,:,0:m+2)

    ! Build stub RHS in local slab coordinates (global z offset = pdec%k0)
    call build_rhs_gaussian(n, h, pdec%k0, m, b_loc)

    ! Call legacy solver
    k_it = 0
    ktot = 0.0_real64
    res  = 0.0_real64

    call pdesolver(u_mg, b_loc, bcnd_loc, h, int(n,kind=4), ncycl, eps, omega, &
                   k_it, ktot, res, ng, pdec%rank, pdec%nproc)

    ! Copy solution back and gather to global
    pdec%phi_dom(:,:,0:m+2) = u_mg(:,:,0:m+2)
    call pdec%gather_phi_to_global(phi_global)

    deallocate(u_mg, b_loc, bcnd_loc)
  end subroutine solve_poisson_legacy_stub

end module mod_poisson_driver_legacy
