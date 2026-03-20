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

  type :: CollisionWorkspace
    real(real64), allocatable :: np_mx(:)
    real(real64), allocatable :: nu_max(:)
    integer(int32), allocatable :: Nc(:,:)
    integer(int32), allocatable :: np_add(:,:)
    integer(int32), allocatable :: err_coll(:,:)
  contains
    procedure :: ensure_sizes => collision_workspace_ensure_sizes
    procedure :: clear        => collision_workspace_clear
  end type CollisionWorkspace

contains

  subroutine collision_workspace_ensure_sizes(self, ntype, nproc)
    class(CollisionWorkspace), intent(inout) :: self
    integer(int32),            intent(in)    :: ntype, nproc

    if (.not. allocated(self%np_mx))   allocate(self%np_mx(ntype))
    if (.not. allocated(self%nu_max))  allocate(self%nu_max(ntype))
    if (.not. allocated(self%Nc))      allocate(self%Nc(ntype,nproc))
    if (.not. allocated(self%np_add))  allocate(self%np_add(ntype,nproc))
    if (.not. allocated(self%err_coll)) allocate(self%err_coll(ntype,nproc))
  end subroutine collision_workspace_ensure_sizes


  subroutine collision_workspace_clear(self)
    class(CollisionWorkspace), intent(inout) :: self

    if (allocated(self%np_mx))    self%np_mx    = 0.0_real64
    if (allocated(self%nu_max))   self%nu_max   = 0.0_real64
    if (allocated(self%Nc))       self%Nc       = 0_int32
    if (allocated(self%np_add))   self%np_add   = 0_int32
    if (allocated(self%err_coll)) self%err_coll = 0_int32
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


  subroutine pick_target_same_or_neighbor_cell(part, n, ici, tproc, iseed, itarget, found)
    class(ParticleSet), intent(in)    :: part
    integer(int32),     intent(in)    :: n(3)
    integer(int32),     intent(in)    :: ici, tproc
    integer(int32),     intent(inout) :: iseed
    integer(int32),     intent(out)   :: itarget
    logical,            intent(out)   :: found

    integer(int32) :: cells_to_try(7)
    integer(int32) :: k, icell, i0, i1, count
    real(real64)   :: rnd

    found   = .false.
    itarget = 0_int32

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


  subroutine ensure_particle_capacity(part, needed)
    class(ParticleSet), intent(inout) :: part
    integer(int32),     intent(in)    :: needed

    if (needed <= 0_int32) return
    call part%ensure_capacity(needed)
  end subroutine ensure_particle_capacity


  subroutine perform_collisions_step( &
        part, n, h, np_red, mass, charge, vt0, Nm, &
        p_ncol, sig, sig_Er, sig_list, sig_Eex, col_info, sigv_mx, &
        ns_coll, dt, nu_uplim, iseed, nproc_mpi, mpi_rank, workspace)
    !
    ! Modernized collision driver following legacy logic:
    !   - compute np_mx
    !   - compute nu_max and Nc
    !   - perform null-collision Monte Carlo
    !   - choose targets from same / neighbor cells
    !   - update velocities and create byproducts
    !
    ! Notes:
    !   * This is designed for your current part(ptype,iproc) layout.
    !   * It assumes particles are already sorted by cell.
    !   * It assumes one MPI rank owns a local set of ParticleSet buckets.
    !
    type(ParticleSet),  intent(inout) :: part(:,:)
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: np_red(:,:,:,:)
    real(real64),       intent(in)    :: mass(:), charge(:), vt0(:), Nm(:)
    integer(int32),     intent(in)    :: p_ncol(:)
    real(real64),       intent(in)    :: sig(:,:), sig_Er(:), sig_Eex(:,:)
    integer(int32),     intent(in)    :: sig_list(:,:), col_info(:,:)
    real(real64),       intent(in)    :: sigv_mx(:,:)
    integer(int32),     intent(in)    :: ns_coll
    real(real64),       intent(in)    :: dt, nu_uplim(:)
    integer(int32),     intent(inout) :: iseed(:)
    integer(int32),     intent(in)    :: nproc_mpi, mpi_rank
    type(CollisionWorkspace), intent(inout) :: workspace

    integer(int32) :: ntype, nproc
    integer(int32) :: ptype, iproc, icol, ind_col, ttype
    integer(int32) :: n_re, sum_np_tot_local, sum_np_tot_global, ierr
    integer(int32) :: Nc_tmp
    real(real64)   :: Pmax, dNc, rnd

    ! collision loop locals
    integer(int32) :: ic, ip, ix, iy, iz, ici
    integer(int32) :: itarget, tproc
    logical        :: found_target, flag_coll
    integer(int32) :: c_ind, n_by, i_by, btype, ib
    real(real64)   :: vx1, vy1, vz1, vx2, vy2, vz2
    real(real64)   :: mu, vr, Ekr, Eth, Ee, sum_mass_inv
    real(real64)   :: vx_cm, vy_cm, vz_cm
    real(real64)   :: np_t, sig_p, Ek_L, Ek_R, sig_L, sig_R
    real(real64)   :: nu(128), sum_nu
    real(real64)   :: sort_arr(128)
    integer(int32) :: indx(128)
    integer(int32) :: ipt, ipt_L, ipt_R, ipt_M, npt
    real(real64)   :: th_add, costh, phi, ex, ey, ez, ex1, ey1, ez1
    real(real64)   :: costh_s, phi_s, vp
    integer(int32) :: by_count_needed

    ntype = int(size(part,1), int32)
    nproc = int(size(part,2), int32)

    call workspace%ensure_sizes(ntype, nproc)
    call workspace%clear()

    ! ------------------------------------------------------------
    ! 1) Compute np_mx per species from current reduced density
    ! ------------------------------------------------------------
    do ptype = 1, ntype
      workspace%np_mx(ptype) = maxval(np_red(1:n(1)+1,1:n(2)+1,1:n(3)+1,ptype))
    end do

    ! ------------------------------------------------------------
    ! 2) Legacy-style setup of nu_max and Nc
    ! ------------------------------------------------------------
    do ptype = 1, ntype

      if (p_ncol(ptype) == 0_int32) cycle
      if (mass(ptype) <= 0.0_real64) cycle

      workspace%nu_max(ptype) = 0.0_real64

      do icol = 1, p_ncol(ptype)
        ind_col = sig_list(ptype,icol)
        n_re    = col_info(ind_col,1)
        ttype   = col_info(ind_col, 2 + n_re + col_info(ind_col,2))
        if (ttype < 1 .or. ttype > ntype) cycle
        workspace%nu_max(ptype) = workspace%nu_max(ptype) + workspace%np_mx(ttype) * sigv_mx(ptype,ind_col)
      end do

      workspace%nu_max(ptype) = min(workspace%nu_max(ptype), nu_uplim(ptype))

      sum_np_tot_local = 0_int32
      do iproc = 1, nproc
        sum_np_tot_local = sum_np_tot_local + part(ptype,iproc)%n
      end do

      if (nproc_mpi > 1_int32) then
        call MPI_Allreduce(sum_np_tot_local, sum_np_tot_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      else
        sum_np_tot_global = sum_np_tot_local
      end if

      Pmax = workspace%nu_max(ptype) * real(ns_coll, real64) * dt
      if (Pmax > 1.0_real64) then
        if (mpi_rank == 0) then
          write(*,*) 'Collision probability greater than 1, please correct ...'
          write(*,*) 'Origin: mod_collisions::perform_collisions_step'
        end if
        call stop_calculation
      end if

      dNc   = real(sum_np_tot_global, real64) * Pmax / real(max(1_int32,nproc_mpi*nproc), real64)
      Nc_tmp = int(dNc, int32)

      rnd = ran2(iseed(1))
      if (rnd <= (dNc - real(Nc_tmp, real64))) Nc_tmp = Nc_tmp + 1_int32

      do iproc = 1, nproc
        workspace%Nc(ptype,iproc) = Nc_tmp
      end do
    end do

    ! ------------------------------------------------------------
    ! 3) Collision Monte Carlo, close to legacy collision_OMP
    ! ------------------------------------------------------------
    do iproc = 1, nproc
      do ptype = 1, ntype

        if (p_ncol(ptype) == 0_int32) cycle
        if (mass(ptype) <= 0.0_real64) cycle
        if (workspace%Nc(ptype,iproc) == 0_int32) cycle
        if (part(ptype,iproc)%n == 0_int32) cycle

        do ic = 1, workspace%Nc(ptype,iproc)

          rnd = ran2(iseed(iproc))
          ip  = int(real(part(ptype,iproc)%n, real64) * rnd, int32) + 1_int32
          if (ip > part(ptype,iproc)%n) ip = part(ptype,iproc)%n

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

          nu      = 0.0_real64
          sum_nu  = 0.0_real64
          sort_arr= 0.0_real64

          do icol = 1, p_ncol(ptype)
            ind_col = sig_list(ptype,icol)
            n_re    = col_info(ind_col,1)
            ttype   = col_info(ind_col, 2 + n_re + col_info(ind_col,2))

            if (ttype < 1 .or. ttype > ntype) cycle

            tproc = iproc
            call pick_target_same_or_neighbor_cell(part(ttype,tproc), n, ici, tproc, iseed(iproc), itarget, found_target)

            if (.not. found_target) then
              workspace%err_coll(ptype,iproc) = workspace%err_coll(ptype,iproc) + 1_int32
              cycle
            end if

            if (charge(ttype) == 0.0_real64) then
              call box_muller(vx2, vy2, vt0(ttype), [ran2(iseed(iproc)), ran2(iseed(iproc))])
              call box_muller(vz2, vz2, vt0(ttype), [ran2(iseed(iproc)), ran2(iseed(iproc))])
            else
              vx2 = part(ttype,tproc)%vx(itarget)
              vy2 = part(ttype,tproc)%vy(itarget)
              vz2 = part(ttype,tproc)%vz(itarget)
            end if

            mu = abs(mass(ptype))*abs(mass(ttype)) / (abs(mass(ptype)) + abs(mass(ttype)))
            vr = sqrt((vx1-vx2)**2 + (vy1-vy2)**2 + (vz1-vz2)**2)
            Ekr = 0.5_real64 * mu * vr*vr / 1.60217646e-19_real64

            ! Legacy-style dichotomy on sig_Er
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

                if (workspace%nu_max(ptype) > 0.0_real64) then
                  np_t = density_average_at_particle_cell(np_red, ix, iy, iz, ttype)
                  nu(icol) = np_t * sig_p * vr / workspace%nu_max(ptype)
                  sum_nu   = sum_nu + nu(icol)
                  sort_arr(icol) = nu(icol)
                end if
                exit
              end if
            end do
          end do

          call indexx_small(p_ncol(ptype), sort_arr, indx)

          rnd = ran2(iseed(iproc))
          flag_coll = .false.
          c_ind     = 0_int32

          if (rnd <= sum_nu) then
            sum_nu = 0.0_real64
            do icol = 1, p_ncol(ptype)
              sum_nu = sum_nu + nu(indx(icol))
              if (rnd <= sum_nu) then
                flag_coll = .true.
                c_ind     = sig_list(ptype, indx(icol))
                exit
              end if
            end do
          end if

          if (.not. flag_coll) cycle

          n_re = col_info(c_ind,1)
          n_by = col_info(c_ind,2)
          ttype = col_info(c_ind, 2 + n_re + n_by)
          Eth   = sig_Eex(c_ind, ind_Eth)

          ! Pick target again deterministically from same local structure
          tproc = iproc
          call pick_target_same_or_neighbor_cell(part(ttype,tproc), n, ici, tproc, iseed(iproc), itarget, found_target)
          if (.not. found_target) cycle

          vx2 = part(ttype,tproc)%vx(itarget)
          vy2 = part(ttype,tproc)%vy(itarget)
          vz2 = part(ttype,tproc)%vz(itarget)

          mu = abs(mass(ptype))*abs(mass(ttype)) / (abs(mass(ptype)) + abs(mass(ttype)))
          vr = sqrt((vx1-vx2)**2 + (vy1-vy2)**2 + (vz1-vz2)**2)

          Ee = 0.5_real64 * mu * vr*vr - Eth * 1.60217646e-19_real64
          if (Ee < 0.0_real64) cycle

          vx_cm = (abs(mass(ptype))*vx1 + abs(mass(ttype))*vx2) / (abs(mass(ptype))+abs(mass(ttype)))
          vy_cm = (abs(mass(ptype))*vy1 + abs(mass(ttype))*vy2) / (abs(mass(ptype))+abs(mass(ttype)))
          vz_cm = (abs(mass(ptype))*vz1 + abs(mass(ttype))*vz2) / (abs(mass(ptype))+abs(mass(ttype)))

          sum_mass_inv = 0.0_real64
          do i_by = 1, n_by
            btype = col_info(c_ind, 2 + n_re + i_by)
            if (btype >= 1 .and. btype <= ntype) then
              if (mass(btype) > 0.0_real64) sum_mass_inv = sum_mass_inv + 1.0_real64 / abs(mass(btype))
            end if
          end do
          if (sum_mass_inv <= 0.0_real64) cycle

          th_add = 2.0_real64 * acos(-1.0_real64) / real(max(1_int32,n_by), real64)
          costh  = 1.0_real64 - 2.0_real64 * ran2(iseed(iproc))
          phi    = 2.0_real64 * acos(-1.0_real64) * ran2(iseed(iproc))

          costh_s = 1.0_real64 - 2.0_real64 * ran2(iseed(iproc))
          phi_s   = 2.0_real64 * acos(-1.0_real64) * ran2(iseed(iproc))

          ! remove/overwrite incident particle with first matching byproduct when possible
          by_count_needed = 0_int32

          do i_by = 1, n_by
            btype = col_info(c_ind, 2 + n_re + i_by)
            if (btype < 1 .or. btype > ntype) cycle
            if (mass(btype) <= 0.0_real64) cycle

            ex = cos(acos(costh))
            ey = sin(acos(costh)) * sin(phi)
            ez = sin(acos(costh)) * cos(phi)
            call scatter_vector(ex1, ey1, ez1, ex, ey, ez, costh_s, phi_s)

            vp = sqrt(2.0_real64 * (Ee / sum_mass_inv) / abs(mass(btype))**2)

            if (i_by == 1 .and. btype == ptype) then
              ib = ip
            else if (btype == ttype) then
              ib = itarget
            else
              by_count_needed = by_count_needed + 1_int32
              call ensure_particle_capacity(part(btype,iproc), part(btype,iproc)%n + by_count_needed)
              ib = part(btype,iproc)%n + by_count_needed
            end if

            if (ib > part(btype,iproc)%n) then
              if (ib == part(btype,iproc)%n + 1_int32) then
                part(btype,iproc)%n = ib
              end if
            end if

            part(btype,iproc)%x(ib)  = part(ptype,iproc)%x(ip)
            part(btype,iproc)%y(ib)  = part(ptype,iproc)%y(ip)
            part(btype,iproc)%z(ib)  = part(ptype,iproc)%z(ip)
            part(btype,iproc)%vx(ib) = vx_cm + vp*ex1
            part(btype,iproc)%vy(ib) = vy_cm + vp*ey1
            part(btype,iproc)%vz(ib) = vz_cm + vp*ez1
          end do

          ! Legacy kills extra reactants by marking dead and letting mover remove them.
          ! In your current architecture you do not yet carry a dead-mask in ParticleSet,
          ! so for now we conservatively overwrite incident and target via byproducts only.
          ! When you add a dead mask to ParticleSet, extend this part directly.
        end do
      end do
    end do
  end subroutine perform_collisions_step

    subroutine box_muller(v1, v2, vt, rnd)
    real(real64), intent(out) :: v1, v2
    real(real64), intent(in)  :: vt
    real(real64), intent(in)  :: rnd(2)

    real(real64) :: r1, r2, fac

    r1 = max(rnd(1), 1.0e-14_real64)
    r2 = rnd(2)

    fac = vt * sqrt(-log(1.0_real64 - r1))
    v1  = fac * cos(2.0_real64 * acos(-1.0_real64) * r2)
    v2  = fac * sin(2.0_real64 * acos(-1.0_real64) * r2)
    end subroutine box_muller


  subroutine indexx_small(n, arr, indx)
    integer(int32), intent(in)  :: n
    real(real64),   intent(in)  :: arr(*)
    integer(int32), intent(out) :: indx(*)

    integer(int32) :: i, j, itmp
    real(real64)   :: atmp

    do i = 1, n
      indx(i) = i
    end do

    do i = 1, n-1
      do j = i+1, n
        if (arr(indx(j)) < arr(indx(i))) then
          itmp    = indx(i)
          indx(i) = indx(j)
          indx(j) = itmp
        end if
      end do
    end do
  end subroutine indexx_small

end module mod_collisions