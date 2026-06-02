module mod_nullCollisions
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  use mod_RNG,         only: ran2
  use mod_collisionProducts, only: apply_collision_products

  implicit none
  private
  public :: perform_null_collisions_selection_test

  real(real64), parameter :: QE_ABS = 1.602176634e-19_real64
  real(real64), parameter :: NU_SAFETY_ION = 1.2_real64

contains

  subroutine perform_null_collisions_selection_test( &
      part, ntype_tracked, ntype_all, mass, Nm, p_ncol, sig_list, col_info, &
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

    integer(int32) :: ptype, iproc, k, ip
    integer(int32) :: nproc, n_total, n_selected
    integer(int32) :: selected_total, accepted
    integer(int32) :: chosen_col, chosen_ttype
    real(real64) :: nu_max, P_null, n_selected_real, rnd
    real(real64) :: sum_nu

    integer(int32) :: npt_sig

    nproc = int(size(part,2), int32)

    npt_sig = count_valid_energy_points(sig_Er)

    if (npt_sig < 2_int32) then
      error stop "perform_null_collisions: invalid sig_Er table"
    end if

    ! if (count_valid_energy_points(sig_Er) < 2_int32) then
    !   error stop "perform_null_collisions: invalid sig_Er table"
    ! end if

    do ptype = 1_int32, ntype_tracked

      if (p_ncol(ptype) <= 0_int32) cycle
      if (mass(ptype) <= 0.0_real64) cycle

      if (ptype == 1_int32) then

        nu_max = estimate_legacy_electron_nu_max( &
            ptype, ntype_tracked, ntype_all, p_ncol, sig_list, col_info, &
            sigv_mx, ni0)

        if (ptype <= size(nu_uplim)) then
          nu_max = min(nu_max, nu_uplim(ptype))
        end if

      else

        nu_max = estimate_dynamic_nu_max( &
            part, ptype, ntype_tracked, ntype_all, mass, p_ncol, &
            sig_list, col_info, sig, sig_Er, npt_sig, ni0)

        nu_max = NU_SAFETY_ION * nu_max

      end if

      if (nu_max <= 0.0_real64) cycle

      P_null = 1.0_real64 - exp(-nu_max * real(ns_coll,real64) * dt)

      accepted = 0_int32
      selected_total = 0_int32

      !$omp parallel do default(shared) &
      !$omp private(iproc,k,ip,n_total,n_selected,n_selected_real,rnd,sum_nu,chosen_col,chosen_ttype) &
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
              chosen_ttype   = chosen_ttype)

          rnd = ran2(iseed(iproc)) * nu_max

          if (rnd <= sum_nu .and. chosen_col > 0_int32) then

            accepted = accepted + 1_int32

            call apply_collision_products( &
                part        = part, &
                iproc       = iproc, &
                ptype       = ptype, &
                ip          = ip, &
                target_type = chosen_ttype, &
                target_vx   = 0.0_real64, &
                target_vy   = 0.0_real64, &
                target_vz   = 0.0_real64, &
                c_ind       = chosen_col, &
                col_info    = col_info, &
                sig_Eex     = sig_Eex, &
                mass        = mass, &
                Nm          = Nm, &
                Pcoll       = Pcoll, &
                iseed       = iseed(iproc))


          end if

          call swap_particle_in_set(part(ptype,iproc), ip, n_total)
          n_total = n_total - 1_int32

        end do

      end do
      !$omp end parallel do

    end do



  end subroutine perform_null_collisions_selection_test


  subroutine choose_neutral_collision_channel( &
      p, ip, ptype, ntype_tracked, ntype_all, mass, p_ncol, &
      sig_list, col_info, sig, sig_Er, npt_sig, ni0, iseed, &
      sum_nu, chosen_col, chosen_ttype)

    type(ParticleSet), intent(in) :: p
    integer(int32), intent(in) :: ip, ptype, ntype_tracked, ntype_all
    real(real64), intent(in) :: mass(:)
    integer(int32), intent(in) :: p_ncol(:), npt_sig
    integer(int32), intent(in) :: sig_list(:,:), col_info(:,:)
    real(real64), intent(in) :: sig(:,:), sig_Er(:), ni0(:)
    integer(int32), intent(inout) :: iseed
    real(real64), intent(out) :: sum_nu
    integer(int32), intent(out) :: chosen_col, chosen_ttype

    integer(int32) :: icol, ind_col, n_re, ttype
    real(real64) :: vx1, vy1, vz1
    real(real64) :: vr, mu, Ekr, sig_p
    real(real64) :: nu_icol, target, cumulative

    chosen_col = 0_int32
    chosen_ttype = 0_int32
    sum_nu = 0.0_real64

    vx1 = p%vx(ip)
    vy1 = p%vy(ip)
    vz1 = p%vz(ip)

    vr = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1)

    do icol = 1_int32, p_ncol(ptype)

      ind_col = sig_list(ptype,icol)
      n_re    = col_info(ind_col,1)
      ttype   = col_info(ind_col,2+n_re)

      if (ttype <= ntype_tracked) cycle
      if (ttype >  ntype_all)     cycle

      mu = abs(mass(ptype))*abs(mass(ttype)) / &
           (abs(mass(ptype)) + abs(mass(ttype)))

      Ekr = 0.5_real64 * mu * vr * vr / QE_ABS
      sig_p = interpolate_sigma(Ekr, sig(:,ind_col), sig_Er, npt_sig)

      sum_nu = sum_nu + ni0(ttype) * sig_p * vr

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

      mu = abs(mass(ptype))*abs(mass(ttype)) / &
           (abs(mass(ptype)) + abs(mass(ttype)))

      Ekr = 0.5_real64 * mu * vr * vr / QE_ABS
      sig_p = interpolate_sigma(Ekr, sig(:,ind_col), sig_Er, npt_sig)

      nu_icol = ni0(ttype) * sig_p * vr
      cumulative = cumulative + nu_icol

      if (target <= cumulative) then
        chosen_col = ind_col
        chosen_ttype = ttype
        return
      end if

    end do

  end subroutine choose_neutral_collision_channel


  function estimate_legacy_electron_nu_max( &
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

  end function estimate_legacy_electron_nu_max


  function estimate_dynamic_nu_max( &
      part, ptype, ntype_tracked, ntype_all, mass, p_ncol, &
      sig_list, col_info, sig, sig_Er, npt_sig, ni0) result(nu_max_dyn)

    type(ParticleSet), intent(in) :: part(:,:)
    integer(int32), intent(in) :: ptype, ntype_tracked, ntype_all
    real(real64), intent(in) :: mass(:)
    integer(int32), intent(in) :: p_ncol(:), npt_sig
    integer(int32), intent(in) :: sig_list(:,:), col_info(:,:)
    real(real64), intent(in) :: sig(:,:), sig_Er(:), ni0(:)

    real(real64) :: nu_max_dyn
    integer(int32) :: iproc, ip
    integer(int32) :: stride
    real(real64) :: nu_sum

    nu_max_dyn = 0.0_real64
    stride = 1000_int32

    do iproc = 1_int32, int(size(part,2), int32)

      if (.not. allocated(part(ptype,iproc)%x)) cycle

      do ip = 1_int32, part(ptype,iproc)%n, stride

        nu_sum = compute_particle_neutral_collision_frequency( &
            part(ptype,iproc), ip, ptype, ntype_tracked, ntype_all, &
            mass, p_ncol, sig_list, col_info, sig, sig_Er, npt_sig, ni0)

        if (nu_sum > nu_max_dyn) nu_max_dyn = nu_sum

      end do

    end do

  end function estimate_dynamic_nu_max


  function compute_particle_neutral_collision_frequency( &
      p, ip, ptype, ntype_tracked, ntype_all, mass, p_ncol, &
      sig_list, col_info, sig, sig_Er, npt_sig, ni0) result(sum_nu)

    type(ParticleSet), intent(in) :: p
    integer(int32), intent(in) :: ip, ptype, ntype_tracked, ntype_all
    real(real64), intent(in) :: mass(:)
    integer(int32), intent(in) :: p_ncol(:), npt_sig
    integer(int32), intent(in) :: sig_list(:,:), col_info(:,:)
    real(real64), intent(in) :: sig(:,:), sig_Er(:), ni0(:)

    real(real64) :: sum_nu
    integer(int32) :: icol, ind_col, n_re, ttype
    real(real64) :: vx1, vy1, vz1, vr, mu, Ekr, sig_p

    sum_nu = 0.0_real64

    vx1 = p%vx(ip)
    vy1 = p%vy(ip)
    vz1 = p%vz(ip)

    vr = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1)

    do icol = 1_int32, p_ncol(ptype)

      ind_col = sig_list(ptype,icol)
      n_re    = col_info(ind_col,1)
      ttype   = col_info(ind_col,2+n_re)

      if (ttype <= ntype_tracked) cycle
      if (ttype >  ntype_all)     cycle

      mu = abs(mass(ptype))*abs(mass(ttype)) / &
           (abs(mass(ptype)) + abs(mass(ttype)))

      Ekr = 0.5_real64 * mu * vr * vr / QE_ABS
      sig_p = interpolate_sigma(Ekr, sig(:,ind_col), sig_Er, npt_sig)

      sum_nu = sum_nu + ni0(ttype) * sig_p * vr

    end do

  end function compute_particle_neutral_collision_frequency


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

end module mod_nullCollisions