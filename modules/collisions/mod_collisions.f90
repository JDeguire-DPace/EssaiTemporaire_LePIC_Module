module mod_collisions
  use iso_fortran_env,    only: int32, real64
  use mod_particles,      only: ParticleSet
  use mod_collisionDiagnostics, only: init_rxn_counts
  use mod_collisionsGwenael, only: perform_collisions_gwenael
  use mod_CoulombCollisions, only: perform_coulomb_collisions

  implicit none
  private
  public :: perform_collisions_step

contains

  subroutine perform_collisions_step( &
      part, n, h, ntype_tracked, ntype_all, mass, charge, Ti, Nm, p_ncol, sig_list, col_info, &
      sigv_mx, sig, sig_Er, sig_Eex, ni0, ns_coll, dt, nu_uplim, iseed, &
      mpi_rank, Pcoll, dom_volume, np_red, bcnd, flag_coulomb, n_e)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in) :: n(3)
    real(real64),   intent(in) :: h(3)
    integer(int32), intent(in) :: ntype_tracked, ntype_all
    real(real64), intent(in) :: mass(:), charge(:), Ti(:), Nm(:)
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
    real(real64), intent(in) :: np_red(0:,0:,0:,:)
    integer(int32), intent(in) :: bcnd(0:,0:,0:)
    integer(int32), intent(in) :: flag_coulomb
    real(real64),   intent(in) :: n_e   ! [m^-3] global average electron density

    integer(int32) :: nproc

    call init_rxn_counts(int(size(sigv_mx,2), int32))

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

    ! Coulomb small-angle scattering (Nanbu 2000), same cadence as MC collisions.
    ! Enabled via flag_coulomb == 1 in the Config / input file.
    if (flag_coulomb == 1_int32) then
      nproc = size(part, 2)
      call perform_coulomb_collisions( &
          part   = part, &
          ntype  = ntype_tracked, &
          nproc  = nproc, &
          mass   = mass(1:ntype_tracked), &
          charge = charge(1:ntype_tracked), &
          Nm     = Nm(1:ntype_tracked), &
          n_e    = n_e, &
          dt     = dt * real(ns_coll, real64), &
          iseed  = iseed)
    end if

  end subroutine perform_collisions_step

end module mod_collisions
