module mod_collisions
  use iso_fortran_env,    only: int32, real64
  use mod_particles,      only: ParticleSet
  use mod_MCCcollisions,  only: perform_MCC_collisions
  use mod_BMCCcollisions, only: perform_BMCC_collisions

  implicit none
  private
  public :: perform_collisions_step

contains

  subroutine perform_collisions_step( &
      part, ntype_tracked, ntype_all, mass, Ti, Nm, p_ncol, sig_list, col_info, &
      sigv_mx, sig, sig_Er, sig_Eex, ni0, ns_coll, dt, nu_uplim, iseed, &
      mpi_rank, Pcoll, dom_volume)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in) :: ntype_tracked, ntype_all
    real(real64), intent(in) :: mass(:), Ti(:), Nm(:)
    integer(int32), intent(in) :: p_ncol(:)
    integer(int32), intent(in) :: sig_list(:,:), col_info(:,:)
    real(real64), intent(in) :: sigv_mx(:,:), sig(:,:), sig_Er(:), sig_Eex(:,:)
    real(real64), intent(in) :: ni0(:)
    integer(int32), intent(in) :: ns_coll
    real(real64), intent(in) :: dt, nu_uplim(:)
    integer(int32), intent(inout) :: iseed(:)
    integer(int32), intent(in) :: mpi_rank
    real(real64), intent(inout) :: Pcoll(:,:)
    real(real64), intent(in) :: dom_volume

    ! MCC: projectile (tracked) + neutral background target
    call perform_MCC_collisions( &
        part          = part, &
        ntype_tracked = ntype_tracked, &
        ntype_all     = ntype_all, &
        mass          = mass, &
        Ti            = Ti, &
        Nm            = Nm, &
        p_ncol        = p_ncol, &
        sig_list      = sig_list, &
        col_info      = col_info, &
        sigv_mx       = sigv_mx, &
        sig           = sig, &
        sig_Er        = sig_Er, &
        sig_Eex       = sig_Eex, &
        ni0           = ni0, &
        ns_coll       = ns_coll, &
        dt            = dt, &
        nu_uplim      = nu_uplim, &
        iseed         = iseed, &
        mpi_rank      = mpi_rank, &
        Pcoll         = Pcoll)

    ! BMCC: projectile (tracked) + charged tracked target (e.g. H+ + H-)
    call perform_BMCC_collisions( &
        part          = part, &
        ntype_tracked = ntype_tracked, &
        mass          = mass, &
        Nm            = Nm, &
        p_ncol        = p_ncol, &
        sig_list      = sig_list, &
        col_info      = col_info, &
        sigv_mx       = sigv_mx, &
        sig           = sig, &
        sig_Er        = sig_Er, &
        sig_Eex       = sig_Eex, &
        dom_volume    = dom_volume, &
        ns_coll       = ns_coll, &
        dt            = dt, &
        iseed         = iseed, &
        mpi_rank      = mpi_rank, &
        Pcoll         = Pcoll)

  end subroutine perform_collisions_step

end module mod_collisions
