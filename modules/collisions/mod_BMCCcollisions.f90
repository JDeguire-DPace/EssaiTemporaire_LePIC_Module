module mod_BMCCcollisions
  ! Binary Monte Carlo Collision module.
  !
  ! Handles reactions where the collision TARGET is a tracked charged species
  ! (e.g. H+ + H- mutual neutralization, H+ + H- -> H2+ + e).
  ! The MCC module handles only neutral background targets; this module
  ! handles charged-charged pairs using actual particle velocities.
  !
  ! Algorithm follows the legacy collision_OMP for charged targets:
  !   - For each projectile species ptype, loop over reactions with charged ttype
  !   - Select a random target particle from part(ttype,iproc)
  !   - Null collision acceptance (nu_max computed from instantaneous target density)
  !   - Apply collision products for both particles

  use iso_fortran_env,  only: int32, int8, real64
  use mod_particles,    only: ParticleSet
  use mod_RNG,          only: ran2

  implicit none
  private
  public :: perform_BMCC_collisions

  real(real64), parameter :: QE_ABS = 1.602176634e-19_real64
  real(real64), parameter :: PI_R8  = 3.14159265358979323846_real64

contains

  subroutine perform_BMCC_collisions( &
      part, n, h, ntype_tracked, mass, Nm, p_ncol, sig_list, col_info, &
      sigv_mx, sig, sig_Er, sig_Eex, dom_volume, ns_coll, dt, &
      iseed, mpi_rank, Pcoll)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in) :: n(3)
    real(real64),   intent(in) :: h(3)
    integer(int32),    intent(in)    :: ntype_tracked
    real(real64),      intent(in)    :: mass(:), Nm(:)
    integer(int32),    intent(in)    :: p_ncol(:)
    integer(int32),    intent(in)    :: sig_list(:,:), col_info(:,:)
    real(real64),      intent(in)    :: sigv_mx(:,:)
    real(real64),      intent(in)    :: sig(:,:), sig_Er(:), sig_Eex(:,:)
    real(real64),      intent(in)    :: dom_volume
    integer(int32),    intent(in)    :: ns_coll
    real(real64),      intent(in)    :: dt
    integer(int32),    intent(inout) :: iseed(:)
    integer(int32),    intent(in)    :: mpi_rank
    real(real64),      intent(inout) :: Pcoll(:,:)

    integer(int32) :: ptype, iproc, nproc
    integer(int32) :: icol, ind_col, n_re, ttype
    integer(int32) :: ip, it, k
    integer(int32) :: n_total, n_selected
    real(real64)   :: nu_max, P_null, n_sel_real, rnd
    real(real64)   :: n_ttype_avg
    real(real64)   :: vx1, vy1, vz1, tvx, tvy, tvz
    real(real64)   :: vr, mu, Ekr, sig_p, nu_icol
    real(real64)   :: sum_nu, target_rnd, cumul
    integer(int32) :: chosen_col, chosen_ttype, chosen_it
    integer(int32) :: max_ncol

    ! Per-reaction storage (stack-allocated; size = max reactions per species)
    integer(int32) :: saved_it(size(sig_list,2))
    integer(int32) :: saved_ttype(size(sig_list,2))
    real(real64)   :: saved_tvx(size(sig_list,2))
    real(real64)   :: saved_tvy(size(sig_list,2))
    real(real64)   :: saved_tvz(size(sig_list,2))
    real(real64)   :: nu_store(size(sig_list,2))

    nproc    = int(size(part, 2), int32)
    max_ncol = int(size(sig_list, 2), int32)

    if (count_valid_pts(sig_Er) < 2_int32) return



    do ptype = 1_int32, ntype_tracked

      if (p_ncol(ptype) <= 0_int32)    cycle
      if (mass(ptype)   <= 0.0_real64) cycle

      ! Skip ptype entirely if it has no BMCC reactions (charged target)
      if (.not. has_bmcc_reactions(ptype, ntype_tracked, p_ncol, sig_list, col_info)) cycle

      !$omp parallel do default(shared) &
      !$omp private(iproc,ip,it,k,n_total,n_selected,n_sel_real,rnd, &
      !$omp         nu_max,P_null,n_ttype_avg, &
      !$omp         vx1,vy1,vz1,tvx,tvy,tvz,vr,mu,Ekr,sig_p,nu_icol, &
      !$omp         sum_nu,target_rnd,cumul, &
      !$omp         chosen_col,chosen_ttype,chosen_it, &
      !$omp         icol,ind_col,n_re,ttype) &
      !$omp schedule(static)
      do iproc = 1_int32, nproc
        block
          integer(int32) :: saved_it(max_ncol)
          integer(int32) :: saved_ttype(max_ncol)
          real(real64)   :: saved_tvx(max_ncol)
          real(real64)   :: saved_tvy(max_ncol)
          real(real64)   :: saved_tvz(max_ncol)
          real(real64)   :: nu_store(max_ncol)
          if (.not. allocated(part(ptype,iproc)%x)) cycle
          n_total = part(ptype,iproc)%n
          if (n_total <= 0_int32) cycle

          ! nu_max from instantaneous target density this iproc
          nu_max = 0.0_real64
          do icol = 1_int32, p_ncol(ptype)
            ind_col = sig_list(ptype, icol)
            n_re    = col_info(ind_col, 1)
            ttype   = col_info(ind_col, 2 + n_re)
            if (ttype < 1_int32 .or. ttype > ntype_tracked) cycle
            if (mass(ttype) == 0.0_real64) cycle
            if (.not. allocated(part(ttype,iproc)%x)) cycle
            if (part(ttype,iproc)%n <= 0_int32) cycle
            n_ttype_avg = real(part(ttype,iproc)%n, real64) * Nm(ttype) / dom_volume
            nu_max = nu_max + n_ttype_avg * sigv_mx(ptype, ind_col)
          end do
          if (nu_max <= 0.0_real64) cycle

          P_null = min(nu_max * real(ns_coll, real64) * dt, 1.0_real64)

          n_sel_real = P_null * real(n_total, real64)
          n_selected = int(n_sel_real, int32)
          rnd = ran2(iseed(iproc))
          if (rnd < n_sel_real - real(n_selected, real64)) n_selected = n_selected + 1_int32

          do k = 1_int32, n_selected

            ! Refresh n_total in case kills reduced the array
            n_total = part(ptype,iproc)%n
            if (n_total <= 0_int32) exit

            ! Random projectile
            rnd = ran2(iseed(iproc))
            ip  = int(real(n_total, real64) * rnd, int32) + 1_int32
            if (ip > n_total) ip = n_total

            vx1 = part(ptype,iproc)%vx(ip)
            vy1 = part(ptype,iproc)%vy(ip)
            vz1 = part(ptype,iproc)%vz(ip)

            sum_nu      = 0.0_real64
            nu_store    = 0.0_real64
            saved_it    = 0_int32
            saved_ttype = 0_int32

            ! Accumulate rates and save target particles
            do icol = 1_int32, p_ncol(ptype)
              ind_col = sig_list(ptype, icol)
              n_re    = col_info(ind_col, 1)
              ttype   = col_info(ind_col, 2 + n_re)

              if (ttype < 1_int32 .or. ttype > ntype_tracked) cycle
              if (mass(ttype) == 0.0_real64) cycle
              if (.not. allocated(part(ttype,iproc)%x)) cycle
              if (part(ttype,iproc)%n <= 0_int32) cycle

              ! Random target particle
              call select_local_target_particle( &
                  part(ttype,iproc), &
                  part(ptype,iproc)%x(ip), &
                  part(ptype,iproc)%y(ip), &
                  part(ptype,iproc)%z(ip), &
                  n, h, iseed(iproc), it)

              if (it <= 0_int32) cycle

              tvx = part(ttype,iproc)%vx(it)
              tvy = part(ttype,iproc)%vy(it)
              tvz = part(ttype,iproc)%vz(it)

              vr    = sqrt((vx1-tvx)**2 + (vy1-tvy)**2 + (vz1-tvz)**2)
              mu    = abs(mass(ptype))*abs(mass(ttype)) / (abs(mass(ptype))+abs(mass(ttype)))
              Ekr   = 0.5_real64 * mu * vr * vr / QE_ABS
              sig_p = interp_sigma(Ekr, sig, ind_col, sig_Er, count_valid_pts(sig_Er))

              n_ttype_avg = estimate_local_density(part(ttype,iproc), &
                  part(ptype,iproc)%x(ip), part(ptype,iproc)%y(ip), &
                  part(ptype,iproc)%z(ip), n, h, Nm(ttype))
              nu_icol     = n_ttype_avg * sig_p * vr

              nu_store(icol)   = nu_icol
              saved_it(icol)   = it
              saved_ttype(icol)= ttype
              saved_tvx(icol)  = tvx
              saved_tvy(icol)  = tvy
              saved_tvz(icol)  = tvz
              sum_nu = sum_nu + nu_icol
            end do

            ! Null collision test
            rnd = ran2(iseed(iproc)) * nu_max
            if (rnd > sum_nu) cycle   ! null collision — no swap needed, just skip

            ! Select reaction proportional to rates
            target_rnd   = ran2(iseed(iproc)) * sum_nu
            cumul        = 0.0_real64
            chosen_col   = 0_int32
            chosen_ttype = 0_int32
            chosen_it    = 0_int32
            tvx = 0.0_real64; tvy = 0.0_real64; tvz = 0.0_real64

            do icol = 1_int32, p_ncol(ptype)
              if (nu_store(icol) <= 0.0_real64) cycle
              cumul = cumul + nu_store(icol)
              if (target_rnd <= cumul) then
                chosen_col   = sig_list(ptype, icol)
                chosen_ttype = saved_ttype(icol)
                chosen_it    = saved_it(icol)
                tvx = saved_tvx(icol)
                tvy = saved_tvy(icol)
                tvz = saved_tvz(icol)
                exit
              end if
            end do

            if (chosen_col > 0_int32 .and. chosen_it > 0_int32) then
              call apply_bmcc_products( &
                  part, iproc, ptype, ip, chosen_ttype, chosen_it, &
                  tvx, tvy, tvz, chosen_col, col_info, sig_Eex, &
                  mass, Nm, Pcoll, iseed(iproc), ntype_tracked)
            end if

          end do  ! k
        end block
      end do  ! iproc
      !$omp end parallel do

    end do  ! ptype

  end subroutine perform_BMCC_collisions


  ! Apply collision products for a charged-charged reaction.
  ! Handles both reactants: projectile (ptype,ip) and target (ttype,it).
  ! Neutral products are simply discarded (not tracked).
  ! Charged products reuse the reactant slots or are appended as new particles.
  subroutine apply_bmcc_products( &
      part, iproc, ptype, ip, ttype, it, &
      tvx, tvy, tvz, c_ind, col_info, sig_Eex, &
      mass, Nm, Pcoll, iseed, ntype_tracked)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32),    intent(in)    :: iproc, ptype, ip, ttype, it, c_ind
    real(real64),      intent(in)    :: tvx, tvy, tvz
    integer(int32),    intent(in)    :: col_info(:,:)
    real(real64),      intent(in)    :: sig_Eex(:,:), mass(:), Nm(:)
    real(real64),      intent(inout) :: Pcoll(:,:)
    integer(int32),    intent(inout) :: iseed
    integer(int32),    intent(in)    :: ntype_tracked

    integer(int32) :: n_re, n_by, rt, i_by, btype, ib
    integer(int32) :: flag_proj_kept, flag_targ_kept
    real(real64)   :: vx1, vy1, vz1
    real(real64)   :: vx_cm, vy_cm, vz_cm
    real(real64)   :: vr, mu, Eth, dE_heavy, Erel_after
    real(real64)   :: sum_mass_inv, vp
    real(real64)   :: cos_th, sin_th, th, th_add, phi
    real(real64)   :: cos_th_s, phi_s
    real(real64)   :: ex, ey, ez, ex1, ey1, ez1
    real(real64)   :: rnd1, rnd2, v2old, v2new

    ! Guard: target particle index must still be valid
    if (it > part(ttype,iproc)%n) return

    n_re = col_info(c_ind, 1)
    n_by = col_info(c_ind, 2)
    rt   = col_info(c_ind, 2 + n_re + n_by + 1)

    vx1 = part(ptype,iproc)%vx(ip)
    vy1 = part(ptype,iproc)%vy(ip)
    vz1 = part(ptype,iproc)%vz(ip)

    mu = abs(mass(ptype))*abs(mass(ttype)) / (abs(mass(ptype))+abs(mass(ttype)))
    vr = sqrt((vx1-tvx)**2 + (vy1-tvy)**2 + (vz1-tvz)**2)

    vx_cm = (abs(mass(ptype))*vx1 + abs(mass(ttype))*tvx) / (abs(mass(ptype))+abs(mass(ttype)))
    vy_cm = (abs(mass(ptype))*vy1 + abs(mass(ttype))*tvy) / (abs(mass(ptype))+abs(mass(ttype)))
    vz_cm = (abs(mass(ptype))*vz1 + abs(mass(ttype))*tvz) / (abs(mass(ptype))+abs(mass(ttype)))

    Eth      = sig_Eex(c_ind, 1)
    dE_heavy = sig_Eex(c_ind, 2)

    Erel_after = 0.5_real64*mu*vr*vr - Eth*QE_ABS
    if (Erel_after < 0.0_real64) return

    ! Sum of 1/mass for tracked charged products only
    sum_mass_inv = 0.0_real64
    do i_by = 1_int32, n_by
      btype = col_info(c_ind, 2 + n_re + i_by)
      if (btype >= 1_int32 .and. btype <= ntype_tracked) then
        if (abs(mass(btype)) > 0.0_real64) &
            sum_mass_inv = sum_mass_inv + 1.0_real64 / abs(mass(btype))
      end if
    end do

    ! Scattering angles
    rnd1   = ran2(iseed); rnd2   = ran2(iseed)
    th_add = 2.0_real64*PI_R8 / real(max(n_by,1), real64)
    cos_th = 1.0_real64 - 2.0_real64*rnd1
    th     = acos(cos_th)
    phi    = 2.0_real64*PI_R8*rnd2
    rnd1   = ran2(iseed); rnd2   = ran2(iseed)
    cos_th_s = 1.0_real64 - 2.0_real64*rnd1
    phi_s    = 2.0_real64*PI_R8*rnd2

    flag_proj_kept = 0_int32
    flag_targ_kept = 0_int32

    do i_by = 1_int32, n_by
      btype = col_info(c_ind, 2 + n_re + i_by)

      ! Neutral products are not tracked — skip
      if (btype < 1_int32 .or. btype > ntype_tracked) then
        th = th + th_add
        cycle
      end if
      if (abs(mass(btype)) <= 0.0_real64) then
        th = th + th_add
        cycle
      end if

      cos_th = cos(th)
      sin_th = sin(th)
      th = th + th_add

      ex = cos_th
      ey = sin_th * sin(phi)
      ez = sin_th * cos(phi)
      call scatter_vec(ex1, ey1, ez1, ex, ey, ez, cos_th_s, phi_s)

      ! Assign product to existing slot or create new particle
      if (btype == ptype .and. flag_proj_kept == 0_int32) then
        ib = ip
        flag_proj_kept = 1_int32
      else if (btype == ttype .and. flag_targ_kept == 0_int32) then
        ib = it
        flag_targ_kept = 1_int32
      else
        call append_bmcc(part(ptype,iproc), part(btype,iproc), ip, ib)
      end if

      if (sum_mass_inv > 0.0_real64) then
        vp = sqrt(2.0_real64*(Erel_after/sum_mass_inv) / abs(mass(btype))**2)
      else
        vp = 0.0_real64
      end if

      if (dE_heavy > 0.0_real64 .and. abs(mass(btype)) > 10.0_real64*abs(mass(1))) then
        vp = sqrt(vp*vp + 2.0_real64*dE_heavy*QE_ABS/abs(mass(btype)))
      end if

      v2old = part(btype,iproc)%vx(ib)**2 + &
              part(btype,iproc)%vy(ib)**2 + &
              part(btype,iproc)%vz(ib)**2

      part(btype,iproc)%vx(ib) = vx_cm + vp*ex1
      part(btype,iproc)%vy(ib) = vy_cm + vp*ey1
      part(btype,iproc)%vz(ib) = vz_cm + vp*ez1

      v2new = part(btype,iproc)%vx(ib)**2 + &
              part(btype,iproc)%vy(ib)**2 + &
              part(btype,iproc)%vz(ib)**2

      Pcoll(btype,iproc) = Pcoll(btype,iproc) + &
          0.5_real64*Nm(btype)*abs(mass(btype))*(v2new - v2old)
    end do

    ! Kill reactants that did not appear as products
    if (flag_proj_kept == 0_int32) call kill_bmcc(part(ptype,iproc), ip)
    if (flag_targ_kept == 0_int32) call kill_bmcc(part(ttype,iproc), it)

  end subroutine apply_bmcc_products


  ! Append a new particle at the position of src(ip_src) into dst.
  subroutine append_bmcc(src, dst, ip_src, ip_new)
    type(ParticleSet), intent(in)    :: src
    type(ParticleSet), intent(inout) :: dst
    integer(int32),    intent(in)    :: ip_src
    integer(int32),    intent(out)   :: ip_new

    ip_new = dst%n + 1_int32
    if (ip_new > dst%nmax) then
      ip_new = dst%n   ! capacity exceeded: silently drop
      return
    end if
    dst%n           = ip_new
    dst%x(ip_new)   = src%x(ip_src)
    dst%y(ip_new)   = src%y(ip_src)
    dst%z(ip_new)   = src%z(ip_src)
    dst%vx(ip_new)  = 0.0_real64
    dst%vy(ip_new)  = 0.0_real64
    dst%vz(ip_new)  = 0.0_real64
    if (allocated(dst%flag_dead)) dst%flag_dead(ip_new) = 0_int8
    if (allocated(dst%flag_cex))  dst%flag_cex(ip_new)  = 0_int32
  end subroutine append_bmcc


  ! Remove particle ip by swap-with-last, decrement n.
  subroutine kill_bmcc(p, ip)
    type(ParticleSet), intent(inout) :: p
    integer(int32),    intent(in)    :: ip

    if (ip < 1_int32 .or. ip > p%n) return
    if (ip < p%n) then
      p%x(ip) = p%x(p%n);  p%y(ip) = p%y(p%n);  p%z(ip) = p%z(p%n)
      p%vx(ip)= p%vx(p%n); p%vy(ip)= p%vy(p%n); p%vz(ip)= p%vz(p%n)
      if (allocated(p%flag_dead)) p%flag_dead(ip) = p%flag_dead(p%n)
      if (allocated(p%flag_cex))  p%flag_cex(ip)  = p%flag_cex(p%n)
    end if
    p%n = p%n - 1_int32
  end subroutine kill_bmcc


  ! Rotate vector (vx,vy,vz) by polar angle costheta and azimuthal angle phi.
  subroutine scatter_vec(vx1, vy1, vz1, vx, vy, vz, costheta, phi)
    real(real64), intent(out) :: vx1, vy1, vz1
    real(real64), intent(in)  :: vx, vy, vz, costheta, phi

    real(real64) :: sintheta, cosphi, sinphi, v, vv

    sintheta = sqrt(max(1.0_real64 - costheta**2, 0.0_real64))
    cosphi   = cos(phi)
    sinphi   = sin(phi)
    v = sqrt(vx**2 + vy**2 + vz**2)

    if (v == 0.0_real64) then
      vx1 = 0.0_real64; vy1 = 0.0_real64; vz1 = 0.0_real64
      return
    end if

    if (abs(vy) > abs(vz)) then
      vv  = sqrt(vx**2 + vy**2)
      vx1 = vx*costheta + (vy*v*sinphi + vx*vz*cosphi)/vv*sintheta
      vy1 = vy*costheta + (-vx*v*sinphi + vy*vz*cosphi)/vv*sintheta
      vz1 = vz*costheta - vv*cosphi*sintheta
    else
      vv  = sqrt(vx**2 + vz**2)
      vx1 = vx*costheta + (vz*v*sinphi - vy*vx*cosphi)/vv*sintheta
      vy1 = vy*costheta + vv*cosphi*sintheta
      vz1 = vz*costheta - (vx*v*sinphi + vy*vz*cosphi)/vv*sintheta
    end if
  end subroutine scatter_vec


  function interp_sigma(E, sig, ind_col, sig_Er, npt) result(sig_p)
    real(real64),   intent(in) :: E, sig(:,:), sig_Er(:)
    integer(int32), intent(in) :: ind_col, npt
    real(real64) :: sig_p
    integer(int32) :: lo, hi, mid

    if (npt < 2_int32) then; sig_p = 0.0_real64; return; end if
    if (E <= sig_Er(1))   then; sig_p = sig(1,ind_col);   return; end if
    if (E >= sig_Er(npt)) then; sig_p = sig(npt,ind_col); return; end if

    lo = 1_int32; hi = npt
    do while (hi - lo > 1_int32)
      mid = (lo + hi) / 2_int32
      if (E <= sig_Er(mid)) then; hi = mid; else; lo = mid; end if
    end do
    sig_p = sig(lo,ind_col) + &
        (E - sig_Er(lo))*(sig(hi,ind_col) - sig(lo,ind_col)) / (sig_Er(hi) - sig_Er(lo))
  end function interp_sigma


  function count_valid_pts(sig_Er) result(npt)
    real(real64), intent(in) :: sig_Er(:)
    integer(int32) :: npt, ipt
    npt = 0_int32
    do ipt = 1_int32, int(size(sig_Er), int32)
      if (sig_Er(ipt) > 0.0_real64) npt = ipt
    end do
  end function count_valid_pts


  logical function has_bmcc_reactions(ptype, ntype_tracked, p_ncol, sig_list, col_info)
    integer(int32), intent(in) :: ptype, ntype_tracked
    integer(int32), intent(in) :: p_ncol(:), sig_list(:,:), col_info(:,:)
    integer(int32) :: icol, ind_col, n_re, ttype

    has_bmcc_reactions = .false.
    do icol = 1_int32, p_ncol(ptype)
      ind_col = sig_list(ptype, icol)
      n_re    = col_info(ind_col, 1)
      ttype   = col_info(ind_col, 2 + n_re)
      if (ttype >= 1_int32 .and. ttype <= ntype_tracked) then
        has_bmcc_reactions = .true.
        return
      end if
    end do
  end function has_bmcc_reactions

  subroutine select_local_target_particle(p, xp, yp, zp, n, h, iseed, it)
    type(ParticleSet), intent(in) :: p
    real(real64), intent(in) :: xp, yp, zp
    integer(int32), intent(in) :: n(3)
    real(real64), intent(in) :: h(3)
    integer(int32), intent(inout) :: iseed
    integer(int32), intent(out) :: it

    integer(int32) :: ix0, iy0, iz0
    integer(int32) :: ix, iy, iz
    integer(int32) :: j, trial, max_trial
    real(real64) :: rnd

    it = 0_int32

    if (p%n <= 0_int32) return

    ix0 = int(xp / h(1), int32) + 1_int32
    iy0 = int(yp / h(2), int32) + 1_int32
    iz0 = int(zp / h(3), int32) + 1_int32

    ix0 = max(1_int32, min(n(1), ix0))
    iy0 = max(1_int32, min(n(2), iy0))
    iz0 = max(1_int32, min(n(3), iz0))

    max_trial = min(2000_int32, max(100_int32, p%n))

    do trial = 1_int32, max_trial

      rnd = ran2(iseed)
      j = int(real(p%n, real64) * rnd, int32) + 1_int32
      if (j > p%n) j = p%n

      ix = int(p%x(j) / h(1), int32) + 1_int32
      iy = int(p%y(j) / h(2), int32) + 1_int32
      iz = int(p%z(j) / h(3), int32) + 1_int32

      if (abs(ix - ix0) <= 1_int32 .and. &
          abs(iy - iy0) <= 1_int32 .and. &
          abs(iz - iz0) <= 1_int32) then
        it = j
        return
      end if

    end do

  end subroutine select_local_target_particle

  function estimate_local_density(p, xp, yp, zp, n, h, Nm) result(ndens)
    type(ParticleSet), intent(in) :: p
    real(real64), intent(in) :: xp, yp, zp
    integer(int32), intent(in) :: n(3)
    real(real64), intent(in) :: h(3)
    real(real64), intent(in) :: Nm
    real(real64) :: ndens

    integer(int32) :: j, ix0, iy0, iz0, ix, iy, iz
    integer(int32) :: count_local
    real(real64) :: Vlocal

    ix0 = int(xp / h(1), int32) + 1_int32
    iy0 = int(yp / h(2), int32) + 1_int32
    iz0 = int(zp / h(3), int32) + 1_int32

    count_local = 0_int32

    do j = 1_int32, p%n
      ix = int(p%x(j) / h(1), int32) + 1_int32
      iy = int(p%y(j) / h(2), int32) + 1_int32
      iz = int(p%z(j) / h(3), int32) + 1_int32

      if (abs(ix - ix0) <= 1_int32 .and. &
          abs(iy - iy0) <= 1_int32 .and. &
          abs(iz - iz0) <= 1_int32) then
        count_local = count_local + 1_int32
      end if
    end do

    Vlocal = 27.0_real64 * h(1) * h(2) * h(3)

    ndens = real(count_local, real64) * Nm / Vlocal

  end function estimate_local_density

end module mod_BMCCcollisions
