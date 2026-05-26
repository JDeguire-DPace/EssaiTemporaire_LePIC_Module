module mod_collisions
  use iso_fortran_env, only: int32, real64
  use mpi

  use mod_particles,        only: ParticleSet
  use mod_particle_sorting, only: get_cell_particle_range
  use mod_utils,            only: stop_calculation
  use mod_RNG,              only: ran2

  implicit none
  private

  public :: CollisionWorkspace
  public :: scatter_vector
  public :: density_average_at_particle_cell
  public :: perform_collisions_step

  integer(int32), parameter :: ind_Eth = 1_int32
  real(real64),   parameter :: qe_si   = 1.60217646e-19_real64
  real(real64),   parameter :: pi      = acos(-1.0_real64)
  integer(int32), parameter :: MAX_COL = 128_int32

  logical, parameter :: DEBUG_COLLISIONS = .false.

  type :: CollisionWorkspace
    real(real64),    allocatable :: np_mx(:)
    real(real64),    allocatable :: np_mx_all(:)
    real(real64),    allocatable :: nu_max(:)
    integer(int32),  allocatable :: Nc(:,:)
    integer(int32),  allocatable :: np_add(:,:)
    integer(int32),  allocatable :: err_coll(:,:)
  contains
    procedure :: ensure_sizes => collision_workspace_ensure_sizes
    procedure :: clear        => collision_workspace_clear
  end type CollisionWorkspace

contains

  subroutine collision_workspace_ensure_sizes(self, ntype_tracked, ntype_total, nproc)
    class(CollisionWorkspace), intent(inout) :: self
    integer(int32), intent(in) :: ntype_tracked, ntype_total, nproc

    if (.not. allocated(self%np_mx))     allocate(self%np_mx(ntype_tracked))
    if (.not. allocated(self%np_mx_all)) allocate(self%np_mx_all(ntype_total))
    if (.not. allocated(self%nu_max))    allocate(self%nu_max(ntype_tracked))
    if (.not. allocated(self%Nc))        allocate(self%Nc(ntype_tracked,nproc))
    if (.not. allocated(self%np_add))    allocate(self%np_add(ntype_tracked,nproc))
    if (.not. allocated(self%err_coll))  allocate(self%err_coll(ntype_tracked,nproc))
  end subroutine collision_workspace_ensure_sizes


  subroutine collision_workspace_clear(self)
    class(CollisionWorkspace), intent(inout) :: self

    if (allocated(self%np_mx))     self%np_mx     = 0.0_real64
    if (allocated(self%np_mx_all)) self%np_mx_all = 0.0_real64
    if (allocated(self%nu_max))    self%nu_max    = 0.0_real64
    if (allocated(self%Nc))        self%Nc        = 0_int32
    if (allocated(self%np_add))    self%np_add    = 0_int32
    if (allocated(self%err_coll))  self%err_coll  = 0_int32
  end subroutine collision_workspace_clear


  subroutine scatter_vector(vx1, vy1, vz1, vx, vy, vz, costheta, phi)
    real(real64), intent(out) :: vx1, vy1, vz1
    real(real64), intent(in)  :: vx, vy, vz, costheta, phi

    real(real64) :: sintheta, cosphi, sinphi, v, vv

    sintheta = sqrt(max(1.0_real64 - costheta*costheta, 0.0_real64))
    sinphi   = sin(phi)
    cosphi   = cos(phi)
    v        = sqrt(vx*vx + vy*vy + vz*vz)

    vx1 = vx
    vy1 = vy
    vz1 = vz

    if (v == 0.0_real64) return

    if (abs(vy) > abs(vz)) then
      vv  = sqrt(vx*vx + vy*vy)
      vx1 = vx*costheta + (vy*v*sinphi + vx*vz*cosphi)/vv * sintheta
      vy1 = vy*costheta + (-vx*v*sinphi + vy*vz*cosphi)/vv * sintheta
      vz1 = vz*costheta - vv*cosphi*sintheta
    else
      vv  = sqrt(vx*vx + vz*vz)
      vx1 = vx*costheta + (vz*v*sinphi - vy*vx*cosphi)/vv * sintheta
      vy1 = vy*costheta + vv*cosphi*sintheta
      vz1 = vz*costheta - (vx*v*sinphi + vy*vz*cosphi)/vv * sintheta
    end if
  end subroutine scatter_vector


  pure real(real64) function density_average_at_particle_cell(np_red, ix, iy, iz, ttype) result(np_t)
    real(real64),   intent(in) :: np_red(:,:,:,:)
    integer(int32), intent(in) :: ix, iy, iz, ttype

    np_t = 0.125_real64 * ( &
         np_red(ix  ,iy  ,iz  ,ttype) + np_red(ix+1,iy  ,iz  ,ttype) + &
         np_red(ix  ,iy+1,iz  ,ttype) + np_red(ix+1,iy+1,iz  ,ttype) + &
         np_red(ix  ,iy  ,iz+1,ttype) + np_red(ix+1,iy  ,iz+1,ttype) + &
         np_red(ix  ,iy+1,iz+1,ttype) + np_red(ix+1,iy+1,iz+1,ttype) )
  end function density_average_at_particle_cell


  pure integer(int32) function legacy_linear_cell(ix, iy, iz, n) result(icell)
    integer(int32), intent(in) :: ix, iy, iz
    integer(int32), intent(in) :: n(3)

    icell = (ix - 1_int32) + n(1) * ((iy - 1_int32) + n(2) * (iz - 1_int32)) + 1_int32
  end function legacy_linear_cell


  subroutine pick_target_same_or_neighbor_cell(part, n, ici, iseed, itarget, found)
    class(ParticleSet), intent(in)    :: part
    integer(int32),     intent(in)    :: n(3)
    integer(int32),     intent(in)    :: ici
    integer(int32),     intent(inout) :: iseed
    integer(int32),     intent(out)   :: itarget
    logical,            intent(out)   :: found

    integer(int32) :: cells_to_try(7)
    integer(int32) :: k, icell, i0, i1, count
    real(real64)   :: rnd

    found   = .false.
    itarget = 0_int32

    if (.not. allocated(part%cell_count)) return
    if (.not. allocated(part%cell_start)) return

    cells_to_try(1) = ici
    cells_to_try(2) = ici + 1_int32
    cells_to_try(3) = ici - 1_int32
    cells_to_try(4) = ici + n(1)
    cells_to_try(5) = ici - n(1)
    cells_to_try(6) = ici + n(1)*n(2)
    cells_to_try(7) = ici - n(1)*n(2)

    do k = 1, 7
      icell = cells_to_try(k)
      if (icell < 1_int32) cycle
      if (icell > size(part%cell_count)) cycle

      call get_cell_particle_range(part, icell, i0, i1, count)
      if (count <= 0_int32) cycle

      rnd = ran2(iseed)
      itarget = i0 + int(rnd * real(count, real64), int32)
      if (itarget > i1) itarget = i1

      found = .true.
      return
    end do
  end subroutine pick_target_same_or_neighbor_cell


  pure real(real64) function target_density(ttype, np_red, np_mx_all, ix, iy, iz, ntype_tracked) result(np_t)
    integer(int32), intent(in) :: ttype, ix, iy, iz, ntype_tracked
    real(real64),   intent(in) :: np_red(:,:,:,:)
    real(real64),   intent(in) :: np_mx_all(:)

    if (ttype >= 1_int32 .and. ttype <= ntype_tracked) then
      np_t = density_average_at_particle_cell(np_red, ix, iy, iz, ttype)
    else if (ttype >= 1_int32 .and. ttype <= size(np_mx_all)) then
      np_t = np_mx_all(ttype)
    else
      np_t = 0.0_real64
    end if
  end function target_density


  subroutine neutral_velocity(vx, vy, vz, vt, iseed)
    real(real64),    intent(out)   :: vx, vy, vz
    real(real64),    intent(in)    :: vt
    integer(int32),  intent(inout) :: iseed

    real(real64) :: r1, r2, dummy

    r1 = ran2(iseed)
    r2 = ran2(iseed)
    call box_muller(vx, vy, vt, r1, r2)

    r1 = ran2(iseed)
    r2 = ran2(iseed)
    call box_muller(vz, dummy, vt, r1, r2)
  end subroutine neutral_velocity


  subroutine append_or_overwrite_product(part, btype, iproc, ib, x, y, z, vx, vy, vz)
    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32),    intent(in)    :: btype, iproc, ib
    real(real64),      intent(in)    :: x, y, z, vx, vy, vz

    call part(btype,iproc)%ensure_capacity(ib)

    if (ib > part(btype,iproc)%n) part(btype,iproc)%n = ib

    part(btype,iproc)%x(ib)  = x
    part(btype,iproc)%y(ib)  = y
    part(btype,iproc)%z(ib)  = z
    part(btype,iproc)%vx(ib) = vx
    part(btype,iproc)%vy(ib) = vy
    part(btype,iproc)%vz(ib) = vz

    if (allocated(part(btype,iproc)%flag_dead)) part(btype,iproc)%flag_dead(ib) = 0
    if (allocated(part(btype,iproc)%flag_cex))  part(btype,iproc)%flag_cex(ib)  = 0
    if (allocated(part(btype,iproc)%sp))        part(btype,iproc)%sp(ib)        = btype
    if (allocated(part(btype,iproc)%w))         part(btype,iproc)%w(ib)         = 1.0_real64
  end subroutine append_or_overwrite_product


  subroutine perform_collisions_step( &
        part, n, h, np_red, mass, charge, vt0, Nm, neutral_density, &
        p_ncol, sig, sig_Er, sig_list, sig_Eex, col_info, sigv_mx, &
        ns_coll, dt, nu_uplim, iseed, nproc_mpi, mpi_rank, workspace, Pcoll)

    type(ParticleSet),  intent(inout) :: part(:,:)
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: np_red(:,:,:,:)
    real(real64),       intent(in)    :: mass(:), charge(:), vt0(:), Nm(:)
    real(real64),       intent(in)    :: neutral_density(:)
    integer(int32),     intent(in)    :: p_ncol(:)
    real(real64),       intent(in)    :: sig(:,:), sig_Er(:), sig_Eex(:,:)
    integer(int32),     intent(in)    :: sig_list(:,:), col_info(:,:)
    real(real64),       intent(in)    :: sigv_mx(:,:)
    integer(int32),     intent(in)    :: ns_coll
    real(real64),       intent(in)    :: dt, nu_uplim(:)
    integer(int32),     intent(inout) :: iseed(:)
    integer(int32),     intent(in)    :: nproc_mpi, mpi_rank
    type(CollisionWorkspace), intent(inout) :: workspace
    real(real64),       intent(inout) :: Pcoll(:,:)

    integer(int32) :: ntype_tracked, ntype_total, nproc
    integer(int32) :: ptype, iproc, icol, ind_col, ttype
    integer(int32) :: n_re, n_by, sum_np_tot_local, sum_np_tot_global, ierr
    integer(int32) :: Nc_tmp
    real(real64)   :: Pmax, dNc, rnd, np_for_nu

    integer(int32) :: ic, ip, ix, iy, iz, ici
    integer(int32) :: itarget, c_ind, i_by, i_re, btype, ctype, ib
    integer(int32) :: chosen_icol, append_index
    logical        :: found_target, flag_coll
    !logical        :: target_is_tracked
    real(real64)   :: vx1, vy1, vz1, vx2, vy2, vz2
    real(real64)   :: mu, vr, Ekr, Eth, Ee, sum_mass_inv
    real(real64)   :: vx_cm, vy_cm, vz_cm
    real(real64)   :: np_t, sig_p, Ek_L, Ek_R, sig_L, sig_R
    real(real64)   :: nu(MAX_COL), sum_nu
    
    !!! MODIFICATION
    real(real64)   :: sort_arr(MAX_COL)
    integer(int32) :: indx(MAX_COL)
    logical        :: use_fast_electron_neutral
    integer(int32) :: e_neu_count, k
    integer(int32) :: e_neu_col(MAX_COL)
    integer(int32) :: e_neu_ttype(MAX_COL)

    integer(int32) :: saved_target_ip(size(mass))
    logical        :: target_selected(size(mass))
    logical        :: target_found_by_type(size(mass))
    real(real64)   :: saved_target_vx(size(mass))
    real(real64)   :: saved_target_vy(size(mass))
    real(real64)   :: saved_target_vz(size(mass))
    real(real64)   :: saved_vr(size(mass))
    real(real64)   :: saved_Ekr(size(mass))
    real(real64)   :: saved_mu(size(mass))
    real(real64)   :: saved_vx_cm(size(mass))
    real(real64)   :: saved_vy_cm(size(mass))
    real(real64)   :: saved_vz_cm(size(mass))

    integer(int32) :: ipt, ipt_L, ipt_R, ipt_M, npt
    real(real64)   :: th_add, costh, phi, ex, ey, ez, ex1, ey1, ez1
    real(real64)   :: costh_s, phi_s, vp, st
    real(real64)   :: x0, y0, z0
    real(real64)   :: nu_max_local
    real(real64)   :: v2old, v2new
    integer(int32) :: flag_add(10)
    integer(int32) :: accepted_total
    integer(int32) :: ineu

    ntype_tracked = int(size(part,1), int32)
    nproc         = int(size(part,2), int32)
    ntype_total   = int(size(mass), int32)

    call workspace%ensure_sizes(ntype_tracked, ntype_total, nproc)
    call workspace%clear()

    accepted_total = 0_int32

    ! Build legacy-style global np_mx(ttype).
    ! Charged species: max reduced density.
    ! Neutrals: physical neutral density passed from mod_simulation.
    do ptype = 1, ntype_tracked
      workspace%np_mx(ptype) = maxval(np_red(1:n(1)+1,1:n(2)+1,1:n(3)+1,ptype))
      workspace%np_mx_all(ptype) = workspace%np_mx(ptype)
    end do

    do ttype = ntype_tracked + 1_int32, ntype_total
      ineu = ttype - ntype_tracked
      if (ineu >= 1_int32 .and. ineu <= size(neutral_density)) then
        workspace%np_mx_all(ttype) = neutral_density(ineu)
      else
        workspace%np_mx_all(ttype) = 0.0_real64
      end if
    end do

    if (DEBUG_COLLISIONS .and. mpi_rank == 0) then
      write(*,*) "COLL neutral_density = ", neutral_density
      write(*,*) "COLL np_mx_all       = ", workspace%np_mx_all
      write(*,*) "COLL ntype_tracked   = ", ntype_tracked
      write(*,*) "COLL ntype_total     = ", ntype_total
      write(*,*) "COLL p_ncol          = ", p_ncol(1:min(size(p_ncol),ntype_total))
    end if

    ! ------------------------------------------------------------
    ! Legacy null-collision setup: compute nu_max and Nc.
    ! ------------------------------------------------------------
    do ptype = 1, ntype_tracked
      if (ptype > size(p_ncol)) cycle
      if (p_ncol(ptype) == 0_int32) cycle
      if (mass(ptype) <= 0.0_real64) cycle

      workspace%nu_max(ptype) = 0.0_real64

      do icol = 1, min(p_ncol(ptype), MAX_COL)
        ind_col = sig_list(ptype,icol)
        n_re    = col_info(ind_col,1)
        ttype   = col_info(ind_col, 2 + n_re)

        if (ttype < 1_int32 .or. ttype > ntype_total) cycle

        np_for_nu = workspace%np_mx_all(ttype)
        workspace%nu_max(ptype) = workspace%nu_max(ptype) + &
             np_for_nu * sigv_mx(ptype,ind_col)
      end do

      if (ptype <= size(nu_uplim)) then
        workspace%nu_max(ptype) = min(workspace%nu_max(ptype), nu_uplim(ptype))
      end if

      sum_np_tot_local = 0_int32
      do iproc = 1, nproc
        sum_np_tot_local = sum_np_tot_local + part(ptype,iproc)%n
      end do

      if (nproc_mpi > 1_int32) then
        call MPI_Allreduce(sum_np_tot_local, sum_np_tot_global, 1, &
             MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      else
        sum_np_tot_global = sum_np_tot_local
      end if

      Pmax = workspace%nu_max(ptype) * real(ns_coll, real64) * dt

      if (Pmax > 1.0_real64) then
        if (mpi_rank == 0) then
          write(*,*) 'Collision probability greater than 1.'
          write(*,*) 'ptype, Pmax = ', ptype, Pmax
        end if
        call stop_calculation
      end if

      ! Legacy:
      ! dNc = sum_np_tot * Pmax / (nproc_mpi*nproc)
      dNc = real(sum_np_tot_global, real64) * Pmax / &
            real(max(1_int32,nproc_mpi*nproc), real64)

      Nc_tmp = int(dNc, int32)
      rnd = ran2(iseed(1))
      if (rnd <= (dNc - real(Nc_tmp, real64))) Nc_tmp = Nc_tmp + 1_int32

      workspace%Nc(ptype,:) = Nc_tmp
    end do


    ! ------------------------------------------------------------
    ! Prebuild electron-neutral reaction list for fast path.
    ! ------------------------------------------------------------
    e_neu_count = 0_int32

    if (ntype_tracked >= 1_int32) then
      do icol = 1, min(p_ncol(1), MAX_COL)
        ind_col = sig_list(1,icol)
        n_re    = col_info(ind_col,1)
        ttype   = col_info(ind_col, 2 + n_re)

        if (ttype > ntype_tracked .and. ttype <= ntype_total) then
          e_neu_count = e_neu_count + 1_int32
          e_neu_col(e_neu_count)   = ind_col
          e_neu_ttype(e_neu_count) = ttype
        end if
      end do
    end if


    ! ------------------------------------------------------------
    ! Perform collisions.
    ! ------------------------------------------------------------
    
    !$omp parallel do schedule(static) private( &
    !$omp iproc, ptype, icol, ind_col, ttype, n_re, n_by, ic, ip, ix, iy, iz, ici, &
    !$omp itarget, c_ind, i_by, i_re, btype, ctype, ib, chosen_icol, append_index, &
    !$omp found_target, flag_coll, vx1, vy1, vz1, vx2, vy2, vz2, &
    !$omp mu, vr, Ekr, Eth, Ee, sum_mass_inv, vx_cm, vy_cm, vz_cm, np_t, sig_p, &
    !$omp Ek_L, Ek_R, sig_L, sig_R, nu, sum_nu, sort_arr, indx, saved_target_ip, &
    !$omp target_selected, target_found_by_type, saved_target_vx, saved_target_vy, &
    !$omp saved_target_vz, saved_vr, saved_Ekr, saved_mu, saved_vx_cm, saved_vy_cm, &
    !$omp saved_vz_cm, ipt, ipt_L, ipt_R, ipt_M, npt, th_add, costh, phi, ex, ey, ez, &
    !$omp ex1, ey1, ez1, costh_s, phi_s, vp, st, x0, y0, z0, nu_max_local, v2old, &
    !$omp v2new, flag_add, use_fast_electron_neutral, k )
    do iproc = 1, nproc
      do ptype = 1, ntype_tracked
        if (ptype > size(p_ncol)) cycle
        if (p_ncol(ptype) == 0_int32) cycle
        if (mass(ptype) <= 0.0_real64) cycle
        if (workspace%Nc(ptype,iproc) == 0_int32) cycle
        if (part(ptype,iproc)%n == 0_int32) cycle

        ! Legacy correction:
        ! nu_max_OMP = Nc / np_tot / (ns_coll*dt)
        nu_max_local = real(workspace%Nc(ptype,iproc), real64) / &
             real(max(1_int32, part(ptype,iproc)%n), real64) / &
             (real(ns_coll, real64) * dt)

        if (nu_max_local <= 0.0_real64) cycle

        use_fast_electron_neutral = (ptype == 1_int32)

        do ic = 1, workspace%Nc(ptype,iproc)

          rnd = ran2(iseed(iproc))
          ip  = int(real(part(ptype,iproc)%n, real64) * rnd, int32) + 1_int32
          if (ip > part(ptype,iproc)%n) ip = part(ptype,iproc)%n

          if (allocated(part(ptype,iproc)%flag_dead)) then
            if (part(ptype,iproc)%flag_dead(ip) /= 0) cycle
          end if

          ix = int(part(ptype,iproc)%x(ip) / h(1), int32) + 1_int32
          iy = int(part(ptype,iproc)%y(ip) / h(2), int32) + 1_int32
          iz = int(part(ptype,iproc)%z(ip) / h(3), int32) + 1_int32

          ix = max(1_int32, min(n(1), ix))
          iy = max(1_int32, min(n(2), iy))
          iz = max(1_int32, min(n(3), iz))

          ici = legacy_linear_cell(ix, iy, iz, n)

          vx1 = part(ptype,iproc)%vx(ip)
          vy1 = part(ptype,iproc)%vy(ip)
          vz1 = part(ptype,iproc)%vz(ip)

          ! -----------------------------------------------------
          ! Fast path for electron-neutral collisions.
          ! This avoids the charged-target DSMC search and the full
          ! general target bookkeeping for ptype=1.  It handles only
          ! reactions whose target ttype is neutral.  Electron-charged
          ! reactions are skipped in this fast branch.
          ! -----------------------------------------------------
          if (use_fast_electron_neutral) then
            nu       = 0.0_real64
            sum_nu   = 0.0_real64
            sort_arr = 0.0_real64

            target_selected      = .false.
            target_found_by_type = .false.

            do k = 1, e_neu_count
              ind_col = e_neu_col(k)
              ttype   = e_neu_ttype(k)
              n_re    = col_info(ind_col,1)
              n_by    = col_info(ind_col,2)
              icol    = k

              ! ! Fast branch handles electron-neutral reactions only.
              ! if (ttype <= ntype_tracked) cycle
              ! if (ttype < 1_int32 .or. ttype > ntype_total) cycle

              if (.not. target_selected(ttype)) then
                call neutral_velocity(vx2, vy2, vz2, vt0(ttype), iseed(iproc))

                saved_target_ip(ttype) = 0_int32
                saved_target_vx(ttype) = vx2
                saved_target_vy(ttype) = vy2
                saved_target_vz(ttype) = vz2

                saved_mu(ttype) = abs(mass(ptype))*abs(mass(ttype)) / &
                     (abs(mass(ptype)) + abs(mass(ttype)))

                saved_vr(ttype) = sqrt((vx1-vx2)**2 + (vy1-vy2)**2 + (vz1-vz2)**2)
                saved_Ekr(ttype) = 0.5_real64 * saved_mu(ttype) * saved_vr(ttype)**2 / qe_si

                saved_vx_cm(ttype) = (abs(mass(ptype))*vx1 + abs(mass(ttype))*vx2) / &
                     (abs(mass(ptype)) + abs(mass(ttype)))
                saved_vy_cm(ttype) = (abs(mass(ptype))*vy1 + abs(mass(ttype))*vy2) / &
                     (abs(mass(ptype)) + abs(mass(ttype)))
                saved_vz_cm(ttype) = (abs(mass(ptype))*vz1 + abs(mass(ttype))*vz2) / &
                     (abs(mass(ptype)) + abs(mass(ttype)))

                target_selected(ttype)      = .true.
                target_found_by_type(ttype) = .true.
              end if

              vr  = saved_vr(ttype)
              Ekr = saved_Ekr(ttype)

              npt   = size(sig_Er)
              ipt_L = 1_int32
              ipt_R = npt

              do ipt = 1, 5
                ipt_M = int(real(ipt_R - ipt_L, real64)/2.0_real64, int32) + ipt_L
                Ek_L  = sig_Er(ipt_L)
                Ek_R  = sig_Er(ipt_R)
                if (Ekr >= Ek_L .and. Ekr <= sig_Er(ipt_M)) ipt_R = ipt_M
                if (Ekr >  sig_Er(ipt_M) .and. Ekr <= Ek_R) ipt_L = ipt_M
              end do

              do ipt = ipt_L + 1, ipt_R
                Ek_L = sig_Er(ipt-1)
                Ek_R = sig_Er(ipt)

                if (Ekr > Ek_L .and. Ekr <= Ek_R) then
                  sig_L = sig(ipt-1,ind_col)
                  sig_R = sig(ipt  ,ind_col)

                  sig_p = sig_L + (Ekr - Ek_L) * (sig_R - sig_L) / (Ek_R - Ek_L)

                  np_t = workspace%np_mx_all(ttype)

                  ! Legacy normalized frequency:
                  ! nu = np_t * sigma * vr / nu_max_OMP
                  nu(icol) = np_t * sig_p * vr / nu_max_local
                  sum_nu = sum_nu + nu(icol)
                  sort_arr(icol) = nu(icol)

                  exit
                end if
              end do
            end do

            ! call indexx_small(min(p_ncol(ptype), MAX_COL), sort_arr, indx)
            call indexx_small(e_neu_count, sort_arr, indx)

            rnd       = ran2(iseed(iproc))
            flag_coll = .false.
            c_ind     = 0_int32
            chosen_icol = 0_int32

            if (rnd <= sum_nu) then
              sum_nu = 0.0_real64
              do icol = 1, e_neu_count
                sum_nu = sum_nu + nu(indx(icol))
                if (rnd <= sum_nu) then
                  flag_coll   = .true.
                  chosen_icol = indx(icol)
                  c_ind       = e_neu_col(chosen_icol)
                  exit
                end if
              end do
            end if

            if (.not. flag_coll) cycle

            ! Jump to common product-application block.
            goto 777
          end if

          nu       = 0.0_real64
          sum_nu   = 0.0_real64
          
          !!! MODIFICATION
          sort_arr = 0.0_real64

          target_selected      = .false.
          target_found_by_type = .false.
          saved_target_ip      = 0_int32
          saved_target_vx      = 0.0_real64
          saved_target_vy      = 0.0_real64
          saved_target_vz      = 0.0_real64
          saved_vr             = 0.0_real64
          saved_Ekr            = 0.0_real64
          saved_mu             = 0.0_real64
          saved_vx_cm          = 0.0_real64
          saved_vy_cm          = 0.0_real64
          saved_vz_cm          = 0.0_real64

          ! -----------------------------------------------------
          ! Extract sigma for every reaction of this projectile.
          ! Reuse target by ttype, like legacy flag_ttype(ttype).
          ! -----------------------------------------------------
          do icol = 1, min(p_ncol(ptype), MAX_COL)
            ind_col = sig_list(ptype,icol)
            n_re    = col_info(ind_col,1)
            n_by    = col_info(ind_col,2)
            ttype   = col_info(ind_col, 2 + n_re)

            if (ttype < 1_int32 .or. ttype > ntype_total) cycle

            if (.not. target_selected(ttype)) then

              if (ttype > ntype_tracked) then
                ! Neutral target: no DSMC particle search.
                call neutral_velocity(vx2, vy2, vz2, vt0(ttype), iseed(iproc))
                saved_target_ip(ttype) = 0_int32

              else
                ! Charged target: find a real nearby target macroparticle.
                call pick_target_same_or_neighbor_cell(part(ttype,iproc), n, ici, &
                     iseed(iproc), itarget, found_target)

                if (.not. found_target) then
                  workspace%err_coll(ptype,iproc) = workspace%err_coll(ptype,iproc) + 1_int32
                  target_selected(ttype)      = .true.
                  target_found_by_type(ttype) = .false.
                  cycle
                end if

                saved_target_ip(ttype) = itarget
                vx2 = part(ttype,iproc)%vx(itarget)
                vy2 = part(ttype,iproc)%vy(itarget)
                vz2 = part(ttype,iproc)%vz(itarget)
              end if

              saved_target_vx(ttype) = vx2
              saved_target_vy(ttype) = vy2
              saved_target_vz(ttype) = vz2

              saved_mu(ttype) = abs(mass(ptype))*abs(mass(ttype)) / &
                   (abs(mass(ptype)) + abs(mass(ttype)))

              saved_vr(ttype) = sqrt((vx1-vx2)**2 + (vy1-vy2)**2 + (vz1-vz2)**2)
              saved_Ekr(ttype) = 0.5_real64 * saved_mu(ttype) * saved_vr(ttype)**2 / qe_si

              saved_vx_cm(ttype) = (abs(mass(ptype))*vx1 + abs(mass(ttype))*vx2) / &
                   (abs(mass(ptype)) + abs(mass(ttype)))
              saved_vy_cm(ttype) = (abs(mass(ptype))*vy1 + abs(mass(ttype))*vy2) / &
                   (abs(mass(ptype)) + abs(mass(ttype)))
              saved_vz_cm(ttype) = (abs(mass(ptype))*vz1 + abs(mass(ttype))*vz2) / &
                   (abs(mass(ptype)) + abs(mass(ttype)))

              target_selected(ttype)      = .true.
              target_found_by_type(ttype) = .true.
            end if

            if (.not. target_found_by_type(ttype)) cycle

            vr  = saved_vr(ttype)
            Ekr = saved_Ekr(ttype)

            npt   = size(sig_Er)
            ipt_L = 1_int32
            ipt_R = npt

            do ipt = 1, 5
              ipt_M = int(real(ipt_R - ipt_L, real64)/2.0_real64, int32) + ipt_L
              Ek_L  = sig_Er(ipt_L)
              Ek_R  = sig_Er(ipt_R)
              if (Ekr >= Ek_L .and. Ekr <= sig_Er(ipt_M)) ipt_R = ipt_M
              if (Ekr >  sig_Er(ipt_M) .and. Ekr <= Ek_R) ipt_L = ipt_M
            end do

            do ipt = ipt_L + 1, ipt_R
              Ek_L = sig_Er(ipt-1)
              Ek_R = sig_Er(ipt)

              if (Ekr > Ek_L .and. Ekr <= Ek_R) then
                sig_L = sig(ipt-1,ind_col)
                sig_R = sig(ipt  ,ind_col)

                sig_p = sig_L + (Ekr - Ek_L) * (sig_R - sig_L) / (Ek_R - Ek_L)

                np_t = target_density(ttype, np_red, workspace%np_mx_all, ix, iy, iz, ntype_tracked)

                ! Legacy normalized frequency:
                ! nu = np_t * sigma * vr / nu_max_OMP
                nu(icol) = np_t * sig_p * vr / nu_max_local
                sum_nu = sum_nu + nu(icol)
                
                !!! MODIFICATION
                sort_arr(icol) = nu(icol)

                exit
              end if
            end do
          end do
          

          !!! MODIFICATION
          call indexx_small(min(p_ncol(ptype), MAX_COL), sort_arr, indx)


          !!! MODIFICATION
          rnd       = ran2(iseed(iproc))
          flag_coll = .false.
          c_ind     = 0_int32
          chosen_icol = 0_int32

          if (rnd <= sum_nu) then
            sum_nu = 0.0_real64
            do icol = 1, min(p_ncol(ptype), MAX_COL)
              sum_nu = sum_nu + nu(indx(icol))
              if (rnd <= sum_nu) then
                flag_coll   = .true.
                chosen_icol = indx(icol)
                c_ind       = sig_list(ptype, chosen_icol)
                exit
              end if
            end do
          end if

777       continue
          if (.not. flag_coll) cycle
          
          !$omp atomic
          accepted_total = accepted_total + 1_int32

          n_re  = col_info(c_ind,1)
          n_by  = col_info(c_ind,2)
          ttype = col_info(c_ind, 2 + n_re)

          if (ttype < 1_int32 .or. ttype > ntype_total) cycle
          if (.not. target_found_by_type(ttype)) cycle

          vx2 = saved_target_vx(ttype)
          vy2 = saved_target_vy(ttype)
          vz2 = saved_target_vz(ttype)

          mu    = saved_mu(ttype)
          vr    = saved_vr(ttype)
          Ekr   = saved_Ekr(ttype)
          vx_cm = saved_vx_cm(ttype)
          vy_cm = saved_vy_cm(ttype)
          vz_cm = saved_vz_cm(ttype)

          Eth = sig_Eex(c_ind, ind_Eth)
          Ee  = 0.5_real64 * mu * vr*vr - Eth * qe_si

          if (Ee < 0.0_real64) cycle

          sum_mass_inv = 0.0_real64
          do i_by = 1, n_by
            btype = col_info(c_ind, 2 + n_re + i_by)
            if (btype >= 1_int32 .and. btype <= ntype_total) then
              if (mass(btype) > 0.0_real64) then
                sum_mass_inv = sum_mass_inv + 1.0_real64 / abs(mass(btype))
              end if
            end if
          end do
          if (sum_mass_inv <= 0.0_real64) cycle

          costh   = 1.0_real64 - 2.0_real64 * ran2(iseed(iproc))
          phi     = 2.0_real64 * pi * ran2(iseed(iproc))
          costh_s = 1.0_real64 - 2.0_real64 * ran2(iseed(iproc))
          phi_s   = 2.0_real64 * pi * ran2(iseed(iproc))
          th_add  = 2.0_real64 * pi / real(max(1_int32,n_by), real64)

          x0 = part(ptype,iproc)%x(ip)
          y0 = part(ptype,iproc)%y(ip)
          z0 = part(ptype,iproc)%z(ip)

          flag_add = 0_int32
          flag_add(1:n_re) = -1_int32
          flag_add(n_re+1:n_re+n_by) = 0_int32

          ! -----------------------------------------------------
          ! Create/update byproducts.
          ! -----------------------------------------------------
          do i_by = 1, n_by
            btype = col_info(c_ind, 2 + n_re + i_by)

            if (btype < 1_int32 .or. btype > ntype_tracked) cycle
            if (mass(btype) <= 0.0_real64) cycle

            st = sqrt(max(0.0_real64, 1.0_real64 - costh*costh))
            ex = costh
            ey = st * sin(phi + th_add*real(i_by-1,real64))
            ez = st * cos(phi + th_add*real(i_by-1,real64))

            call scatter_vector(ex1, ey1, ez1, ex, ey, ez, costh_s, phi_s)

            vp = sqrt(2.0_real64 * (Ee / sum_mass_inv) / abs(mass(btype))**2)

            ! Legacy-like reuse rules.
            if (i_by == 1_int32 .and. btype == ptype) then
              ib = ip
              flag_add(1) = 0_int32
            else if (btype == ttype .and. ttype <= ntype_tracked) then
              ib = saved_target_ip(ttype)
              flag_add(2) = 0_int32
            else
              append_index = part(btype,iproc)%n + 1_int32
              ib = append_index
              flag_add(n_re+i_by) = 1_int32
            end if

            if (ib <= 0_int32) cycle

            if (ib <= part(btype,iproc)%n .and. allocated(part(btype,iproc)%vx)) then
              v2old = part(btype,iproc)%vx(ib)**2 + &
                      part(btype,iproc)%vy(ib)**2 + &
                      part(btype,iproc)%vz(ib)**2
            else
              v2old = 0.0_real64
            end if

            call append_or_overwrite_product(part, btype, iproc, ib, x0, y0, z0, &
                 vx_cm + vp*ex1, vy_cm + vp*ey1, vz_cm + vp*ez1)

            v2new = part(btype,iproc)%vx(ib)**2 + &
                    part(btype,iproc)%vy(ib)**2 + &
                    part(btype,iproc)%vz(ib)**2

            ! Legacy-style P_loss(3,btype,bproc), except local Pcoll slice is Pcoll(btype,iproc).
            Pcoll(btype,iproc) = Pcoll(btype,iproc) + &
                 0.5_real64 * Nm(btype) * mass(btype) * (v2new - v2old)
          end do

          ! -----------------------------------------------------
          ! Kill reactants not reused among products.
          ! -----------------------------------------------------
          do i_re = 1, n_re
            ctype = col_info(c_ind, 2 + i_re)

            if (ctype < 1_int32 .or. ctype > ntype_tracked) cycle
            if (mass(ctype) <= 0.0_real64) cycle

            if (flag_add(i_re) == -1_int32) then
              if (ctype == ptype) then
                ib = ip
                if (allocated(part(ctype,iproc)%flag_dead)) &
                  part(ctype,iproc)%flag_dead(ib) = 1
              else if (ctype == ttype .and. ttype <= ntype_tracked) then
                ib = saved_target_ip(ttype)
                if (ib > 0_int32) then
                  if (allocated(part(ctype,iproc)%flag_dead)) &
                    part(ctype,iproc)%flag_dead(ib) = 1
                end if
              end if
            end if
          end do

          ! Threshold energy loss, legacy dEk_lost equivalent.
          !Pcoll(ptype,iproc) = Pcoll(ptype,iproc) - Nm(ptype) * Eth * qe_si

          if (allocated(part(ptype,iproc)%cell_id))    part(ptype,iproc)%cell_id    = 0_int32
          if (allocated(part(ptype,iproc)%cell_count)) part(ptype,iproc)%cell_count = 0_int32
          if (allocated(part(ptype,iproc)%cell_start)) part(ptype,iproc)%cell_start = 0_int32
        end do
      end do
    end do
    !$omp end parallel do


    if (DEBUG_COLLISIONS .and. mpi_rank == 0) then
      write(*,*) "COLL SUMMARY:"
      write(*,*) "  Nc total       = ", sum(workspace%Nc)
      write(*,*) "  accepted total = ", accepted_total
      write(*,*) "  Pcoll raw      = ", sum(Pcoll)
      write(*,*) "  err_coll total = ", sum(workspace%err_coll)
    end if
  end subroutine perform_collisions_step


  subroutine box_muller(v1, v2, vt, r1_in, r2_in)
    real(real64), intent(out) :: v1, v2
    real(real64), intent(in)  :: vt, r1_in, r2_in

    real(real64) :: r1, r2, fac

    r1  = max(r1_in, 1.0e-14_real64)
    r2  = r2_in
    fac = vt * sqrt(-log(1.0_real64 - r1))
    v1  = fac * cos(2.0_real64 * pi * r2)
    v2  = fac * sin(2.0_real64 * pi * r2)
  end subroutine box_muller

  !!! MODIFICATION
  subroutine indexx_small(n, arr, indx)
    integer(int32), intent(in)  :: n
    real(real64),   intent(in)  :: arr(*)
    integer(int32), intent(out) :: indx(*)

    integer(int32) :: i, j, key

    do i = 1, n
      indx(i) = i
    end do

    do i = 2, n
      key = indx(i)
      j = i - 1

      do while (j >= 1 .and. arr(indx(j)) > arr(key))
        indx(j+1) = indx(j)
        j = j - 1
      end do

      indx(j+1) = key
    end do
  end subroutine indexx_small

end module mod_collisions