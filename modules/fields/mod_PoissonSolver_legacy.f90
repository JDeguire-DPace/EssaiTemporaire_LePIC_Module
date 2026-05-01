module mod_PoissonSolver_legacy
  use iso_fortran_env, only: real64, int32
  use mpi
  use mod_poisson_decomp, only: PoissonDecomp
  use mod_part_info, only: flag_pbc, flag_nmn
  implicit none
  private
  public :: solve_poisson_legacy

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

  subroutine solve_poisson_legacy(pdec, phi_global, bcnd_global, rhs_global, h, n_in, &
                                  ncycl, eps, omega, ng, flag_pbc_in, flag_nmn_in)
    type(PoissonDecomp), intent(inout) :: pdec
    real(real64), intent(inout)        :: phi_global(0:,0:,0:)
    integer(int32), intent(in)         :: bcnd_global(0:,0:,0:)
    real(real64), intent(in)           :: rhs_global(0:,0:,0:)
    real(real64), intent(in)           :: h(3)
    integer(int32), intent(in)         :: n_in(3)
    integer, intent(in)                :: ncycl, ng
    real(real64), intent(in)           :: eps, omega
    integer, intent(in)                :: flag_pbc_in, flag_nmn_in

    integer :: k_it
    real(real64) :: ktot, res
    integer(int32) :: m
    integer :: n_leg(3)
    real(real64) :: h_leg(3)

    real(real64), allocatable :: u_mg(:,:,:)
    real(real64), allocatable :: b_loc(:,:,:)
    integer,      allocatable :: bcnd_loc(:,:,:)

    flag_pbc = flag_pbc_in
    flag_nmn = flag_nmn_in

    call pdec%scatter_from_global(phi_global, bcnd_global)
    call pdec%scatter_rhs_from_global(rhs_global)

    m     = pdec%m
    n_leg = int(n_in, kind=4)
    h_leg = h

    allocate(u_mg(0:n_leg(1)+2, 0:n_leg(2)+2, -1:m+2))
    allocate(b_loc(0:n_leg(1)+1, 0:n_leg(2)+1,  0:m+1))
    allocate(bcnd_loc(0:n_leg(1)+2, 0:n_leg(2)+2, 0:m+2))

    u_mg = 0.0_real64
    b_loc = 0.0_real64
    bcnd_loc = 0

    u_mg(:,:,0:m+2)     = pdec%phi_dom(:,:,0:m+2)
    bcnd_loc(:,:,0:m+2) = pdec%bcnd_dom(:,:,0:m+2)
    b_loc(:,:,:)        = pdec%rhs_dom(:,:,:)

    ! if (pdec%rank == 0) then
    !   write(*,*) 'DEBUG LEGACY-DRIVER flags: flag_pbc, flag_nmn = ', flag_pbc, flag_nmn
    !   write(*,*) 'DEBUG LEGACY-DRIVER phi before pdesolver'
    !   write(*,'(a,es16.8)') 'sum      = ', sum(u_mg)
    !   write(*,'(a,es16.8)') 'sum(abs) = ', sum(abs(u_mg))
    !   write(*,'(a,es16.8)') 'max      = ', maxval(u_mg)
    !   write(*,'(a,es16.8)') 'min      = ', minval(u_mg)
    ! end if

    k_it = 0
    ktot = 0.0_real64
    res  = 0.0_real64



    call pdesolver(u_mg, b_loc, bcnd_loc, h_leg, n_leg, ncycl, eps, omega, &
                   k_it, ktot, res, ng, pdec%rank, pdec%nproc)

    ! if (pdec%rank == 0) then
    !   write(*,*) 'DEBUG LEGACY-DRIVER phi after pdesolver'
    !   write(*,'(a,es16.8)') 'sum      = ', sum(u_mg)
    !   write(*,'(a,es16.8)') 'sum(abs) = ', sum(abs(u_mg))
    !   write(*,'(a,es16.8)') 'max      = ', maxval(u_mg)
    !   write(*,'(a,es16.8)') 'min      = ', minval(u_mg)
    ! end if

    pdec%phi_dom(:,:,0:m+2) = u_mg(:,:,0:m+2)
    call pdec%gather_phi_to_global(phi_global)

    deallocate(u_mg, b_loc, bcnd_loc)
  end subroutine solve_poisson_legacy

end module mod_PoissonSolver_legacy
