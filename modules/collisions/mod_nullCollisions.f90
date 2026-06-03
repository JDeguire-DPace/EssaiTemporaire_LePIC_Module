module mod_nullCollisions
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  use mod_RNG,         only: ran2
  use mod_collisionProducts, only: apply_collision_products

  implicit none
  private
  public :: perform_null_collisions_selection_test

  real(real64), parameter :: QE_ABS = 1.602176634e-19_real64
  real(real64), parameter :: PI_R8 = 3.1415926535897932384626433832795_real64
  
contains

  subroutine perform_null_collisions_selection_test( &
      part, ntype_tracked, ntype_all, mass, Ti, Nm, p_ncol, sig_list, col_info, &
      sigv_mx, sig, sig_Er, sig_Eex, ni0, ns_coll, dt, nu_uplim, iseed, &
      mpi_rank, Pcoll)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in) :: ntype_tracked, ntype_all
    real(real64), intent(in) :: mass(:), Nm(:)
    integer(int32), intent(in) :: p_ncol(:)
    integer(int32), intent(in) :: sig_list(:,:), col_info(:,:)
    real(real64), intent(in) :: sigv_mx(:,:), sig(:,:), sig_Er(:), sig_Eex(:,:)
    real(real64), intent(in) :: ni0(:)
    integer(int32), intent(in) :: ns_coll
    real(real64), intent(in) :: dt, nu_uplim(:)
    integer(int32), intent(inout) :: iseed(:)
    integer(int32), intent(in) :: mpi_rank
    real(real64), intent(inout) :: Pcoll(:,:)
    real(real64), intent(in) :: Ti(:)

    integer(int32) :: ptype, iproc, k, ip
    integer(int32) :: nproc, n_total, n_selected
    integer(int32) :: selected_total, accepted
    integer(int32) :: chosen_col, chosen_ttype
    integer(int32) :: npt_sig
    real(real64) :: nu_max, P_null, n_selected_real, rnd
    real(real64) :: sum_nu
    real(real64) :: target_vx, target_vy, target_vz


    ! mpi_rank is intentionally kept in the interface for compatibility.
    ! It is not used in this production/clean version.
    associate(dummy_rank => mpi_rank)
    end associate
    

    nproc = int(size(part,2), int32)
    npt_sig = count_valid_energy_points(sig_Er)

    if (npt_sig < 2_int32) then
      error stop "perform_null_collisions: invalid sig_Er table"
    end if

    do ptype = 1_int32, ntype_tracked

      if (p_ncol(ptype) <= 0_int32) cycle
      if (mass(ptype) <= 0.0_real64) cycle

      ! Match the legacy null-collision upper-frequency estimate:
      !   nu_max(ptype) = sum_t np_mx(t) * sigv_mx(ptype, reaction)
      ! Here ni0(ttype) is the modular equivalent of the target neutral density.
      ! Charged-target/DSMC reactions are skipped in this module and must be handled
      ! by a separate charged-target collision module.
      nu_max = estimate_legacy_nu_max( &
          ptype, ntype_tracked, ntype_all, p_ncol, sig_list, col_info, &
          sigv_mx, ni0)

      if (ptype <= size(nu_uplim)) then
        nu_max = min(nu_max, nu_uplim(ptype))
      end if

      if (nu_max <= 0.0_real64) cycle

      P_null = nu_max * real(ns_coll,real64) * dt
      if (P_null > 1.0_real64) then
        error stop "perform_null_collisions: P_null > 1"
      end if

      accepted = 0_int32
      selected_total = 0_int32

      !$omp parallel do default(shared) &
      !$omp private(iproc,k,ip,n_total,n_selected,n_selected_real,rnd,sum_nu,chosen_col,chosen_ttype,target_vx,target_vy,target_vz) &
      !$omp reduction(+:selected_total,accepted)
      do iproc = 1_int32, nproc

        if (.not. allocated(part(ptype,iproc)%x)) cycle

        n_total = part(ptype,iproc)%n
        if (n_total <= 0_int32) cycle

        n_selected_real = P_null * real(n_total, real64)
        n_selected = int(n_selected_real, int32)

        rnd = ran2(iseed(iproc))
        if (rnd < n_selected_real - real(n_selected,real64)) then
          n_selected = n_selected + 1_int32
        end if

        selected_total = selected_total + n_selected

        do k = 1_int32, n_selected

          if (n_total <= 0_int32) exit

          rnd = ran2(iseed(iproc))
          ip = int(real(n_total,real64) * rnd, int32) + 1_int32
          if (ip > n_total) ip = n_total

          call choose_neutral_collision_channel( &
              p              = part(ptype,iproc), &
              ip             = ip, &
              ptype          = ptype, &
              ntype_tracked  = ntype_tracked, &
              ntype_all      = ntype_all, &
              mass           = mass, &
              Ti             = Ti, &
              p_ncol         = p_ncol, &
              sig_list       = sig_list, &
              col_info       = col_info, &
              sig            = sig, &
              sig_Er         = sig_Er, &
              npt_sig        = npt_sig, &
              ni0            = ni0, &
              iseed          = iseed(iproc), &
              sum_nu         = sum_nu, &
              chosen_col     = chosen_col, &
              chosen_ttype   = chosen_ttype, &
              target_vx      = target_vx, &
              target_vy      = target_vy, &
              target_vz      = target_vz)

          rnd = ran2(iseed(iproc)) * nu_max

          if (rnd <= sum_nu .and. chosen_col > 0_int32) then

            accepted = accepted + 1_int32

            call apply_collision_products( &
                part        = part, &
                iproc       = iproc, &
                ptype       = ptype, &
                ip          = ip, &
                target_type = chosen_ttype, &
                target_vx   = target_vx, &
                target_vy   = target_vy, &
                target_vz   = target_vz, &
                c_ind       = chosen_col, &
                col_info    = col_info, &
                sig_Eex     = sig_Eex, &
                mass        = mass, &
                Nm          = Nm, &
                Pcoll       = Pcoll, &
                iseed       = iseed(iproc))

          end if

          ! Exclude already-selected particles from being selected again during
          ! the same collision sweep, matching the current modular selection logic.
          call swap_particle_in_set(part(ptype,iproc), ip, n_total)
          n_total = n_total - 1_int32

        end do

      end do
      !$omp end parallel do

    end do

  end subroutine perform_null_collisions_selection_test


  subroutine choose_neutral_collision_channel( &
      p, ip, ptype, ntype_tracked, ntype_all, mass, Ti, p_ncol, &
      sig_list, col_info, sig, sig_Er, npt_sig, ni0, iseed, &
      sum_nu, chosen_col, chosen_ttype, target_vx, target_vy, target_vz)

    type(ParticleSet), intent(in) :: p
    integer(int32), intent(in) :: ip, ptype, ntype_tracked, ntype_all
    real(real64), intent(in) :: mass(:), Ti(:)
    integer(int32), intent(in) :: p_ncol(:), npt_sig
    integer(int32), intent(in) :: sig_list(:,:), col_info(:,:)
    real(real64), intent(in) :: sig(:,:), sig_Er(:), ni0(:)
    integer(int32), intent(inout) :: iseed
    real(real64), intent(out) :: sum_nu
    integer(int32), intent(out) :: chosen_col, chosen_ttype
    real(real64), intent(out) :: target_vx, target_vy, target_vz

    integer(int32) :: icol, ind_col, n_re, ttype
    real(real64) :: vx1, vy1, vz1
    real(real64) :: vr, mu, Ekr, sig_p
    real(real64) :: nu_icol, target, cumulative
    real(real64) :: vx_t, vy_t, vz_t
    real(real64) :: nu_store(1000)
    real(real64) :: vx_store(1000), vy_store(1000), vz_store(1000)

    chosen_col = 0_int32
    chosen_ttype = 0_int32
    sum_nu = 0.0_real64

    target_vx = 0.0_real64
    target_vy = 0.0_real64
    target_vz = 0.0_real64

    nu_store = 0.0_real64
    vx_store = 0.0_real64
    vy_store = 0.0_real64
    vz_store = 0.0_real64

    vx1 = p%vx(ip)
    vy1 = p%vy(ip)
    vz1 = p%vz(ip)

    do icol = 1_int32, p_ncol(ptype)

      ind_col = sig_list(ptype,icol)
      n_re    = col_info(ind_col,1)
      ttype   = col_info(ind_col,2+n_re)

      if (ttype <= ntype_tracked) cycle
      if (ttype >  ntype_all)     cycle

      call sample_maxwellian_velocity(Ti(ttype), mass(ttype), iseed, vx_t, vy_t, vz_t)

      vr = sqrt((vx1 - vx_t)**2 + &
                (vy1 - vy_t)**2 + &
                (vz1 - vz_t)**2)

      mu = abs(mass(ptype))*abs(mass(ttype)) / &
          (abs(mass(ptype)) + abs(mass(ttype)))

      Ekr = 0.5_real64 * mu * vr * vr / QE_ABS
      sig_p = interpolate_sigma(Ekr, sig(:,ind_col), sig_Er, npt_sig)

      nu_icol = ni0(ttype) * sig_p * vr

      nu_store(icol) = nu_icol
      vx_store(icol) = vx_t
      vy_store(icol) = vy_t
      vz_store(icol) = vz_t

      sum_nu = sum_nu + nu_icol

    end do

    if (sum_nu <= 0.0_real64) return

    target = ran2(iseed) * sum_nu
    cumulative = 0.0_real64

    do icol = 1_int32, p_ncol(ptype)

      ind_col = sig_list(ptype,icol)
      n_re    = col_info(ind_col,1)
      ttype   = col_info(ind_col,2+n_re)

      if (ttype <= ntype_tracked) cycle
      if (ttype >  ntype_all)     cycle

      nu_icol = nu_store(icol)
      if (nu_icol <= 0.0_real64) cycle

      cumulative = cumulative + nu_icol

      if (target <= cumulative) then
        chosen_col = ind_col
        chosen_ttype = ttype

        target_vx = vx_store(icol)
        target_vy = vy_store(icol)
        target_vz = vz_store(icol)

        return
      end if

    end do

  end subroutine choose_neutral_collision_channel


  function estimate_legacy_nu_max( &
      ptype, ntype_tracked, ntype_all, p_ncol, sig_list, col_info, &
      sigv_mx, ni0) result(nu_max)

    integer(int32), intent(in) :: ptype, ntype_tracked, ntype_all
    integer(int32), intent(in) :: p_ncol(:)
    integer(int32), intent(in) :: sig_list(:,:), col_info(:,:)
    real(real64), intent(in) :: sigv_mx(:,:), ni0(:)

    real(real64) :: nu_max
    integer(int32) :: icol, ind_col, n_re, ttype

    nu_max = 0.0_real64

    do icol = 1_int32, p_ncol(ptype)

      ind_col = sig_list(ptype,icol)
      n_re    = col_info(ind_col,1)
      ttype   = col_info(ind_col,2+n_re)

      if (ttype <= ntype_tracked) cycle
      if (ttype >  ntype_all)     cycle

      nu_max = nu_max + ni0(ttype) * sigv_mx(ptype,ind_col)

    end do

  end function estimate_legacy_nu_max


  function interpolate_sigma(E, sig_col, sig_Er, npt) result(sig_p)

    real(real64), intent(in) :: E
    real(real64), intent(in) :: sig_col(:), sig_Er(:)
    integer(int32), intent(in) :: npt

    real(real64) :: sig_p
    integer(int32) :: lo, hi, mid

    if (npt < 2_int32) then
      sig_p = 0.0_real64
      return
    end if

    if (E <= sig_Er(1)) then
      sig_p = sig_col(1)
      return
    end if

    if (E >= sig_Er(npt)) then
      sig_p = sig_col(npt)
      return
    end if

    lo = 1_int32
    hi = npt

    do while (hi - lo > 1_int32)
      mid = (lo + hi) / 2_int32

      if (E <= sig_Er(mid)) then
        hi = mid
      else
        lo = mid
      end if
    end do

    sig_p = sig_col(lo) + &
        (E - sig_Er(lo)) * &
        (sig_col(hi) - sig_col(lo)) / &
        (sig_Er(hi) - sig_Er(lo))

  end function interpolate_sigma


  function count_valid_energy_points(sig_Er) result(npt)

    real(real64), intent(in) :: sig_Er(:)
    integer(int32) :: npt
    integer(int32) :: ipt

    npt = 0_int32

    do ipt = 1_int32, int(size(sig_Er), int32)
      if (sig_Er(ipt) > 0.0_real64) npt = ipt
    end do

  end function count_valid_energy_points


  subroutine swap_particle_in_set(p, i, j)

    type(ParticleSet), intent(inout) :: p
    integer(int32), intent(in) :: i, j

    if (i == j) return

    call swap_real(p%x(i),  p%x(j))
    call swap_real(p%y(i),  p%y(j))
    call swap_real(p%z(i),  p%z(j))
    call swap_real(p%vx(i), p%vx(j))
    call swap_real(p%vy(i), p%vy(j))
    call swap_real(p%vz(i), p%vz(j))

  end subroutine swap_particle_in_set


  pure subroutine swap_real(a,b)

    real(real64), intent(inout) :: a,b
    real(real64) :: tmp

    tmp = a
    a = b
    b = tmp

  end subroutine swap_real

  subroutine sample_maxwellian_velocity(T_eV, m, iseed, vx, vy, vz)
    real(real64), intent(in) :: T_eV, m
    integer(int32), intent(inout) :: iseed
    real(real64), intent(out) :: vx, vy, vz

    real(real64) :: sigma_v
    real(real64) :: dummy

    if (T_eV <= 0.0_real64 .or. m == 0.0_real64) then
      vx = 0.0_real64
      vy = 0.0_real64
      vz = 0.0_real64
      return
    end if

    sigma_v = sqrt(QE_ABS*T_eV/abs(m))

    call gaussian_pair(sigma_v, iseed, vx, vy)
    call gaussian_pair(sigma_v, iseed, vz, dummy)
  end subroutine sample_maxwellian_velocity

  subroutine gaussian_pair(sigma_v, iseed, g1, g2)
    real(real64), intent(in) :: sigma_v
    integer(int32), intent(inout) :: iseed
    real(real64), intent(out) :: g1, g2

    real(real64) :: u1, u2, r, theta

    u1 = max(ran2(iseed), 1.0e-300_real64)
    u2 = ran2(iseed)

    r = sigma_v * sqrt(-2.0_real64*log(u1))
    theta = 2.0_real64*PI_R8*u2

    g1 = r*cos(theta)
    g2 = r*sin(theta)
  end subroutine gaussian_pair

end module mod_nullCollisions