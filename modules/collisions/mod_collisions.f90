module mod_collisions
  use iso_fortran_env,    only: int32, real64
  use mod_particles,      only: ParticleSet
  use mod_collisionDiagnostics, only: init_rxn_counts
  use mod_collisionsGwenael, only: perform_collisions_gwenael

  implicit none
  private
  public :: perform_collisions_step

contains

  subroutine perform_collisions_step( &
      part, n, h, ntype_tracked, ntype_all, mass, Ti, Nm, p_ncol, sig_list, col_info, &
      sigv_mx, sig, sig_Er, sig_Eex, ni0, ns_coll, dt, nu_uplim, iseed, &
      mpi_rank, Pcoll, dom_volume, np_red, bcnd)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in) :: n(3)
    real(real64),   intent(in) :: h(3)
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
    ! Deposited density grid (legacy np_red), needed by
    ! mod_collisionsGwenael to reproduce legacy's DSMC local-density
    ! estimate for charged (tracked) collision targets.
    real(real64), intent(in) :: np_red(0:,0:,0:,:)
    ! Boundary-condition flag grid (legacy bcnd), needed to restrict the
    ! np_mx running-max density estimate to interior cells, exactly like
    ! legacy calc_rho.f90.
    integer(int32), intent(in) :: bcnd(0:,0:,0:)

    call init_rxn_counts(int(size(sigv_mx,2), int32))

    ! The legacy collision algorithm (target search via cell-sorted
    ! neighbour lists, energy-dichotomy cross-section lookup, single
    ! null-collision draw across every reaction channel of a species,
    ! deferred particle kill/creation) is implemented in
    ! mod_collisionsGwenael. This routine no longer does anything else -
    ! see that module's header for what is/isn't reproduced exactly.
    call perform_collisions_gwenael( &
        part          = part, &
        n             = n, &
        h             = h, &
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
        Pcoll         = Pcoll, &
        dom_volume    = dom_volume, &
        np_red        = np_red, &
        bcnd          = bcnd)

  end subroutine perform_collisions_step

end module mod_collisions