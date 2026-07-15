module mod_simulation
  use iso_fortran_env, only: int32, real64, int8

  use mod_state,                 only: State
  use mod_density,               only: reduce_species_density, build_rho_from_np
  use mod_chargeDeposition,      only: deposit_particle_set_to_np_thread
  use mod_restart,               only: write_restart_modular
  use mod_injection,             only: inject_particles_volume, inject_flux_particles
  use mod_PoissonSolver_legacy,  only: solve_poisson_legacy, print_poisson_breakdown, reset_poisson_breakdown
  use mod_electricField,         only: calc_Efield_modular
  use mod_output_2d,             only: write_density_planes, write_scalar_planes, &
                                       write_vector_component_planes, &
                                       write_plane_xy_scalar_2d, write_plane_xz_scalar_2d, &
                                       write_plane_yz_scalar_2d
  use mod_constants,             only: eps0
  use mod_collisionDiagnostics, only: print_rxn_counts, reset_rxn_counts, &
                                       print_debug_diagnostics, reset_debug_diagnostics
  use mod_collisions, only: perform_collisions_step
  use mpi

  implicit none
  private
  public :: Simulation

  logical, parameter :: DEBUG_SIMULATION = .false.

  type :: Simulation
    type(State) :: state
    !type(CollisionWorkspace) :: coll_ws

    real(real64) :: t_Erho    = 0.0_real64
    real(real64) :: t_poisson = 0.0_real64
    real(real64) :: t_sort    = 0.0_real64
    real(real64) :: t_avg     = 0.0_real64
    real(real64) :: t_MC      = 0.0_real64
    real(real64) :: t_mover   = 0.0_real64
    real(real64) :: t_bck     = 0.0_real64
    real(real64) :: t_total   = 0.0_real64
    integer(int32) :: t_count = 0_int32

    ! Sub-timers breaking down mover and E/rho (deposit) into their actual
    ! component phases, added to localize the remaining mover/E-rho cost
    ! after class(ParticleSet)->type(ParticleSet) had no measurable effect.
    real(real64) :: t_mover_push = 0.0_real64
    real(real64) :: t_mover_bc   = 0.0_real64
    real(real64) :: t_dep_clear  = 0.0_real64
    real(real64) :: t_dep_loop   = 0.0_real64
    real(real64) :: t_dep_reduce = 0.0_real64

    ! Toggles DATA.BAK/ <-> DATA.BAK2/ on each backup write, matching
    ! legacy's flag_wrt (so a crash mid-write never destroys the only
    ! valid backup - the other directory still has the previous one).
    integer(int32) :: flag_wrt = 0_int32

  contains
    procedure :: init
    procedure :: build_initial_fields
    procedure :: write_initial_diagnostics
    procedure :: deposit_all_particles
    procedure :: collisions_step
    procedure :: output_step
    procedure :: reset_2d_averages
    procedure :: advance_one_step
    procedure :: print_debug_state
    procedure :: print_diagnostics
    procedure :: run
    procedure :: finalize
  end type Simulation

contains

  subroutine init(self, comm_in)
    class(Simulation), intent(inout) :: self
    integer, intent(in) :: comm_in

    call self%state%init(comm_in)
    ! if (self%state%mpi_rank == 0) then
    !   write(*,*) "DEBUG max p_ncol =", maxval(self%state%chem%p_ncol)
    !   write(*,*) "DEBUG p_ncol =", self%state%chem%p_ncol(1:self%state%rxn%ntype)
    !   write(*,*) "DEBUG size(sig_list,2) =", size(self%state%rxn%sig_list,2)
    ! end if
  end subroutine init


  subroutine build_initial_fields(self)
    class(Simulation), intent(inout) :: self

    call build_rho_from_np( &
         n        = int(self%state%dom%n, int32), &
         np_red   = self%state%fld%np, &
         charge   = self%state%chem%charge(1:self%state%ntype), &
         ntype    = int(self%state%ntype), &
         rho      = self%state%fld%rho, &
         bcnd     = self%state%dom%bcnd, &
         flag_pbc = int(self%state%dom%flag_pbc, int32))

    call self%state%apply_dielectric_bc_to_phi()

    call solve_poisson_legacy( &
         pdec        = self%state%pdec, &
         phi_global  = self%state%fld%phi, &
         bcnd_global = self%state%dom%bcnd, &
         rhs_global  = self%state%fld%rho(0:self%state%dom%n(1)+1, &
                                          0:self%state%dom%n(2)+1, &
                                          0:self%state%dom%n(3)+1), &
         h           = self%state%dom%h, &
         n_in        = int(self%state%dom%n, int32), &
         ncycl       = 20000, &
         eps         = self%state%cfg%eps, &
         omega       = self%state%cfg%omega, &
         ng          = self%state%cfg%ng, &
         flag_pbc_in = self%state%dom%flag_pbc, &
         flag_nmn_in = self%state%dom%flag_nmn )

    call calc_Efield_modular( &
         n    = int(self%state%dom%n, int32), &
         h    = self%state%dom%h, &
         phi  = self%state%fld%phi, &
         E    = self%state%fld%E, &
         bcnd = self%state%dom%bcnd )
  end subroutine build_initial_fields


  subroutine write_initial_diagnostics(self)
    class(Simulation), intent(in) :: self

    character(len=128) :: prefix
    character(len=10)  :: s
    integer(int32) :: i
    integer(int32) :: iy_plane_phi, iy_plane_density, iy_plane_E

    if (self%state%mpi_rank /= 0) return

    iy_plane_density = int(self%state%dom%n(2)/2 + 1, int32)
    iy_plane_phi     = int(self%state%dom%n(2)/2,     int32)
    iy_plane_E       = int(self%state%dom%n(2)/2 + 1, int32)

    do i = 1, self%state%ntype
      write(s,'(i0)') i
      prefix = './Output/Output_2D/n' // trim(s)

      call write_density_planes( &
        np       = self%state%fld%np, &
        n        = int(self%state%dom%n, int32), &
        ptype    = i, &
        ix_plane = self%state%params%ix_plot_plane, &
        iy_plane = iy_plane_density, &
        iz_plane = self%state%params%iz_plot_plane, &
        every    = 1_int32, &
        prefix   = prefix )
    end do

    call write_scalar_planes( &
      f        = self%state%fld%phi, &
      n        = int(self%state%dom%n, int32), &
      ix_plane = self%state%params%ix_plot_plane, &
      iy_plane = iy_plane_phi, &
      iz_plane = self%state%params%iz_plot_plane, &
      every    = 1_int32, &
      prefix   = './Output/Output_2D/phi1' )

    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      1_int32, self%state%params%ix_plot_plane, &
                                      iy_plane_E, self%state%params%iz_plot_plane, &
                                      1_int32, './Output/Output_2D/Ex')

    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      2_int32, self%state%params%ix_plot_plane, &
                                      iy_plane_E, self%state%params%iz_plot_plane, &
                                      1_int32, './Output/Output_2D/Ey')

    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      3_int32, self%state%params%ix_plot_plane, &
                                      iy_plane_E, self%state%params%iz_plot_plane, &
                                      1_int32, './Output/Output_2D/Ez')
  end subroutine write_initial_diagnostics


  subroutine deposit_all_particles(self)
    class(Simulation), intent(inout) :: self

    integer(int32) :: ptype, iproc
    real(real64)   :: dt0, dt1

    ! No separate whole-array clear pass: each thread now zeros only its
    ! own (:,:,:,ptype,iproc) slice immediately before depositing into it,
    ! inside the parallel loop below - parallel and NUMA-first-touch
    ! friendly, like legacy's per-thread `np(:,:,:,ptype,iproc)=0.d0`
    ! inside its own parallel mover region. Previously this was one
    ! serial np_thread = 0.0_real64 assignment over the WHOLE ~4.86 GB
    ! array (n1+3)*(n2+3)*(n3+3)*ntype*nproc - single-threaded and the
    ! dominant cost of E/rho (~47 ms/step measured). t_dep_clear is kept
    ! at zero now (folded into t_dep_loop) rather than removed, so old
    ! log-parsing scripts expecting that field don't break.
    self%t_dep_clear = self%t_dep_clear

    dt0 = MPI_Wtime()
    !$omp parallel do collapse(2) private(ptype,iproc) schedule(static)
    do ptype = 1, self%state%ntype
      do iproc = 1, self%state%nproc
        self%state%np_thread(:,:,:,ptype,iproc) = 0.0_real64

        if (.not. allocated(self%state%part(ptype,iproc)%x)) cycle
        if (self%state%part(ptype,iproc)%n <= 0_int32) cycle

        call deposit_particle_set_to_np_thread( &
          part       = self%state%part(ptype,iproc), &
          n          = int(self%state%dom%n, int32), &
          h          = self%state%dom%h, &
          kq         = self%state%fld%kq, &
          Nm_species = self%state%params%Nm(ptype), &
          np_local   = self%state%np_thread(:,:,:,ptype,iproc) )

      end do
    end do
    !$omp end parallel do
    dt1 = MPI_Wtime()
    self%t_dep_loop = self%t_dep_loop + (dt1 - dt0)

    dt0 = MPI_Wtime()
    call reduce_species_density( &
      n         = int(self%state%dom%n, int32), &
      bcnd      = self%state%dom%bcnd, &
      np_thread = self%state%np_thread, &
      ntype     = int(self%state%ntype), &
      nproc     = int(self%state%nproc), &
      mpi_comm  = self%state%comm, &
      np_red    = self%state%fld%np )
    dt1 = MPI_Wtime()
    self%t_dep_reduce = self%t_dep_reduce + (dt1 - dt0)
  end subroutine deposit_all_particles


  subroutine collisions_step(self)
    class(Simulation), intent(inout) :: self

  call perform_collisions_step( &
      part          = self%state%part, &
      n             = self%state%dom%n, &
      h             = self%state%dom%h, &
      ntype_tracked = self%state%ntype, &
      ntype_all     = self%state%rxn%ntype, &
      mass          = self%state%chem%mass(1:self%state%rxn%ntype), &
      charge        = self%state%chem%charge(1:self%state%rxn%ntype), &
      Ti            = self%state%chem%Ti(1:self%state%rxn%ntype), &
      Nm            = self%state%params%Nm(1:self%state%rxn%ntype), &
      p_ncol        = self%state%chem%p_ncol(1:self%state%rxn%ntype), &
      sig_list      = self%state%rxn%sig_list, &
      col_info      = self%state%rxn%col_info, &
      sigv_mx       = self%state%rxn%sigv_mx, &
      sig           = self%state%rxn%sig, &
      sig_Er        = self%state%rxn%sig_Er, &
      sig_Eex       = self%state%rxn%sig_Eex, &
      ni0           = self%state%chem%ni0(1:self%state%rxn%ntype), &
      ns_coll       = self%state%params%nb_step_collisions, &
      dt            = self%state%params%dt, &
      nu_uplim      = self%state%params%nu_uplim(1:self%state%rxn%ntype), &
      iseed         = self%state%params%iseed, &
      mpi_rank      = int(self%state%mpi_rank, int32), &
      Pcoll         = self%state%P_loss(3,:,:), &
      dom_volume    = self%state%dom%h(1) * self%state%dom%h(2) * self%state%dom%h(3) * &
                      real(self%state%dom%n(1)*self%state%dom%n(2)*self%state%dom%n(3), real64), &
      np_red        = self%state%fld%np, &
      bcnd          = self%state%dom%bcnd, &
      flag_coulomb  = self%state%cfg%flag_coulomb, &
      n_e           = sum(self%state%fld%np(:,:,:,1)) / &
                      real(product(self%state%dom%n), real64) )

  end subroutine collisions_step


  subroutine reset_2d_averages(self)
    class(Simulation), intent(inout) :: self

    if (allocated(self%state%np_avg_xy))  self%state%np_avg_xy  = 0.0_real64
    if (allocated(self%state%np_avg_xz))  self%state%np_avg_xz  = 0.0_real64
    if (allocated(self%state%np_avg_yz))  self%state%np_avg_yz  = 0.0_real64

    if (allocated(self%state%phi_avg_xy)) self%state%phi_avg_xy = 0.0_real64
    if (allocated(self%state%phi_avg_xz)) self%state%phi_avg_xz = 0.0_real64
    if (allocated(self%state%phi_avg_yz)) self%state%phi_avg_yz = 0.0_real64

    if (allocated(self%state%data_pavg_xy)) self%state%data_pavg_xy = 0.0_real64
    if (allocated(self%state%data_pavg_xz)) self%state%data_pavg_xz = 0.0_real64
    if (allocated(self%state%data_pavg_yz)) self%state%data_pavg_yz = 0.0_real64

    self%state%cnt_avg = 0_int32
  end subroutine reset_2d_averages


  subroutine output_step(self, istep)
    class(Simulation), intent(inout) :: self
    integer(int32), intent(in) :: istep

    integer(int32) :: i
    integer(int32) :: ix_plane, iy_plane_E, iz_plane
    character(len=256) :: prefix
    character(len=16)  :: sstep, sspecies
    real(real64) :: avg_factor
    real(real64), allocatable :: tmp_xy(:,:), tmp_xz(:,:), tmp_yz(:,:)

    if (self%state%mpi_rank /= 0) return

    ix_plane   = self%state%params%ix_plot_plane
    iz_plane   = self%state%params%iz_plot_plane
    iy_plane_E = int(self%state%dom%n(2)/2 + 1, int32)

    avg_factor = real(max(1_int32, self%state%cnt_avg), real64)

    allocate(tmp_xy(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2))
    allocate(tmp_xz(0:self%state%dom%n(1)+2,0:self%state%dom%n(3)+2))
    allocate(tmp_yz(0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))

    write(sstep,'(i0)') istep

    do i = 1, self%state%ntype
      write(sspecies,'(i0)') i

      tmp_xy = self%state%np_avg_xy(:,:,i) / avg_factor
      tmp_xz = self%state%np_avg_xz(:,:,i) / avg_factor
      tmp_yz = self%state%np_avg_yz(:,:,i) / avg_factor

      prefix = './Output/Output_2D/it' // trim(sstep) // '_n' // trim(sspecies)

      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

      tmp_xy = self%state%data_pavg_xy(2,:,:,i) / avg_factor
      tmp_xz = self%state%data_pavg_xz(2,:,:,i) / avg_factor
      tmp_yz = self%state%data_pavg_yz(2,:,:,i) / avg_factor

      prefix = './Output/Output_2D/it' // trim(sstep) // '_T' // trim(sspecies)

      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)
    end do

    tmp_xy = self%state%phi_avg_xy / avg_factor
    tmp_xz = self%state%phi_avg_xz / avg_factor
    tmp_yz = self%state%phi_avg_yz / avg_factor

    prefix = './Output/Output_2D/it' // trim(sstep) // '_phi'

    call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
    call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
    call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

    prefix = './Output/Output_2D/it' // trim(sstep) // '_Ex'
    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      1_int32, ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)

    prefix = './Output/Output_2D/it' // trim(sstep) // '_Ey'
    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      2_int32, ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)

    prefix = './Output/Output_2D/it' // trim(sstep) // '_Ez'
    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      3_int32, ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)

    if (DEBUG_SIMULATION) then
      write(*,*) "MOD output_step: istep, cnt_avg = ", istep, self%state%cnt_avg
    end if

    deallocate(tmp_xy, tmp_xz, tmp_yz)
  end subroutine output_step


  subroutine print_debug_state(self, istep, label)
    class(Simulation), intent(in) :: self
    integer(int32), intent(in) :: istep
    character(len=*), intent(in) :: label

    if (.not. DEBUG_SIMULATION) return
    if (self%state%mpi_rank /= 0) return

    write(*,*) "MOD DEBUG ", trim(label), " istep = ", istep
    write(*,*) "rho min/max/sum = ", minval(self%state%fld%rho), &
                                  maxval(self%state%fld%rho), &
                                  sum(self%state%fld%rho)
    write(*,*) "phi min/max/sum = ", minval(self%state%fld%phi), &
                                  maxval(self%state%fld%phi), &
                                  sum(self%state%fld%phi)
  end subroutine print_debug_state


  subroutine advance_one_step(self, istep)
    class(Simulation), intent(inout) :: self
    integer(int32), intent(in) :: istep

    real(real64)   :: vt_heat
    real(real64)   :: t0, t1, tstep0
    integer(int32) :: iproc_h, i_h, ix_h, ixl_h, ixr_h
    real(real64)   :: xp_h, yp_h, zp_h, Eki_h, yhalf_h, zhalf_h

    tstep0 = MPI_Wtime()

    ! Sorting
    ! t0 = MPI_Wtime()
    ! if (mod(istep, self%state%params%nb_step_sort) == 1_int32) then
    !   call self%state%sort_particles_local()
    ! end if
    ! t1 = MPI_Wtime()
    ! self%t_sort = self%t_sort + (t1 - t0)

    ! E/rho
    t0 = MPI_Wtime()
    call reduce_species_density( &
      n         = int(self%state%dom%n, int32), &
      bcnd      = self%state%dom%bcnd, &
      np_thread = self%state%np_thread, &
      ntype     = int(self%state%ntype), &
      nproc     = int(self%state%nproc), &
      mpi_comm  = self%state%comm, &
      np_red    = self%state%fld%np )

    call build_rho_from_np( &
      n        = int(self%state%dom%n, int32), &
      np_red   = self%state%fld%np, &
      charge   = self%state%chem%charge(1:self%state%ntype), &
      ntype    = int(self%state%ntype), &
      rho      = self%state%fld%rho, &
      bcnd     = self%state%dom%bcnd, &
      flag_pbc = int(self%state%dom%flag_pbc, int32) )
    t1 = MPI_Wtime()
    self%t_Erho = self%t_Erho + (t1 - t0)

    ! Poisson + E field
    t0 = MPI_Wtime()
    call self%state%apply_dielectric_bc_to_phi()

    if (self%state%mpi_rank == 0 .and. istep == 1) then
      write(*,*) "BEFORE POISSON it=", istep
      write(*,*) "rho min/max/sum = ", minval(self%state%fld%rho), maxval(self%state%fld%rho), sum(self%state%fld%rho)
      write(*,*) "phi min/max/sum = ", minval(self%state%fld%phi), maxval(self%state%fld%phi), sum(self%state%fld%phi)
    end if
    call solve_poisson_legacy( &
      pdec        = self%state%pdec, &
      phi_global  = self%state%fld%phi, &
      bcnd_global = self%state%dom%bcnd, &
      rhs_global  = self%state%fld%rho(0:self%state%dom%n(1)+1, &
                                      0:self%state%dom%n(2)+1, &
                                      0:self%state%dom%n(3)+1), &
      h           = self%state%dom%h, &
      n_in        = int(self%state%dom%n, int32), &
      ncycl       = 20000, &
      eps         = self%state%cfg%eps, &
      omega       = self%state%cfg%omega, &
      ng          = self%state%cfg%ng, &
      flag_pbc_in = self%state%dom%flag_pbc, &
      flag_nmn_in = self%state%dom%flag_nmn )

      if (self%state%mpi_rank == 0 .and. istep == 1) then
        write(*,*) "AFTER POISSON it=", istep
        write(*,*) "rho min/max/sum = ", minval(self%state%fld%rho), maxval(self%state%fld%rho), sum(self%state%fld%rho)
        write(*,*) "phi min/max/sum = ", minval(self%state%fld%phi), maxval(self%state%fld%phi), sum(self%state%fld%phi)
      end if

    call calc_Efield_modular( &
      n    = int(self%state%dom%n, int32), &
      h    = self%state%dom%h, &
      phi  = self%state%fld%phi, &
      E    = self%state%fld%E, &
      bcnd = self%state%dom%bcnd )
    t1 = MPI_Wtime()
    self%t_poisson = self%t_poisson + (t1 - t0)

    ! Averaging
    t0 = MPI_Wtime()
    if (self%state%params%nb_step_averaging > 0_int32) then
      if (mod(istep, self%state%params%nb_step_averaging) == 1_int32) then
        call self%state%compute_plane_moments_local()
        call self%state%accumulate_2d_averages()
      end if
    end if
    t1 = MPI_Wtime()
    self%t_avg = self%t_avg + (t1 - t0)

    ! Reset dielectric charge counters
    if (allocated(self%state%sum_q_xz)) self%state%sum_q_xz = 0.0_real64
    if (allocated(self%state%sum_q_yz)) self%state%sum_q_yz = 0.0_real64

    ! Electron heating — fires BEFORE the push, matching legacy's sequence:
    !   [outside OMP] read Nh/sum_dEk from previous deposit → compute vt
    !   [OMP]         reset Nh/sum_dEk; call eheating(vt); push; deposit
    ! In modular the deposit/accumulation stays after BC (below), but the
    ! vt computation and heating application happen here, before the push.
    if (self%state%cfg%Pabs > 0.0_real64 .and. self%state%cfg%flag_heat == 1) then
      if (self%state%params%nb_step_heating > 0_int32) then
        if (mod(istep, self%state%params%nb_step_heating) == 0_int32) then
          call self%state%update_heating_vt(vt_heat)
          if (self%state%cfg%flag_inj == 1) vt_heat = self%state%params%vt0(1)
          if (allocated(self%state%Nh))      self%state%Nh      = 0_int32
          if (allocated(self%state%sum_dEk)) self%state%sum_dEk = 0.0_real64
          call self%state%apply_electron_heating_local(vt_heat)
        end if
      end if
    end if

    ! Mover + particle BC
    t0 = MPI_Wtime()
    call self%state%move_particles_local()
    t1 = MPI_Wtime()
    self%t_mover_push = self%t_mover_push + (t1 - t0)

    t0 = MPI_Wtime()
    call self%state%apply_particle_bc_local()
    t1 = MPI_Wtime()
    self%t_mover_bc = self%t_mover_bc + (t1 - t0)

    ! t_mover kept as the combined total for the existing summary line
    self%t_mover = self%t_mover_push + self%t_mover_bc

    ! Accumulate Nh and sum_dEk from post-push post-BC electrons, matching
    ! legacy's charge_deposition timing. On non-heating steps, also reset
    ! here (legacy resets at OMP block start every step). On heating steps
    ! the reset already happened above before apply_electron_heating_local.
    if (self%state%cfg%Pabs > 0.0_real64 .and. self%state%cfg%flag_heat == 1) then
      if (self%state%params%nb_step_heating <= 0_int32 .or. &
          mod(istep, self%state%params%nb_step_heating) /= 0_int32) then
        if (allocated(self%state%Nh))      self%state%Nh      = 0_int32
        if (allocated(self%state%sum_dEk)) self%state%sum_dEk = 0.0_real64
      end if

      !$omp parallel do private(iproc_h,i_h,ix_h,ixl_h,ixr_h,xp_h,yp_h,zp_h,Eki_h,yhalf_h,zhalf_h) schedule(static)
      do iproc_h = 1, self%state%nproc
        if (.not. allocated(self%state%part(1,iproc_h)%x)) cycle
        if (self%state%part(1,iproc_h)%n <= 0_int32) cycle

        ixl_h  = int(self%state%cfg%xl_pow / self%state%dom%h(1), int32) + 1_int32
        ixr_h  = int(self%state%cfg%xr_pow / self%state%dom%h(1), int32) + 1_int32
        yhalf_h = self%state%dom%ymax / 2.0_real64
        zhalf_h = self%state%dom%zmax / 2.0_real64

        do i_h = 1_int32, self%state%part(1,iproc_h)%n
          if (self%state%part(1,iproc_h)%flag_dead(i_h) /= 0_int8) cycle
          xp_h = self%state%part(1,iproc_h)%x(i_h)
          ix_h = int(xp_h / self%state%dom%h(1), int32) + 1_int32
          if (ix_h < ixl_h .or. ix_h > ixr_h) cycle

          if (self%state%cfg%flag_circxh == 1_int32) then
            yp_h = self%state%part(1,iproc_h)%y(i_h)
            zp_h = self%state%part(1,iproc_h)%z(i_h)
            if (self%state%cfg%flag_ahp == 0_int32) then
              if (((yp_h-yhalf_h)**2 + (zp_h-zhalf_h)**2) > self%state%cfg%R_ahp**2) cycle
            else
              if (((yp_h-yhalf_h)**2 + (zp_h-zhalf_h)**2) < self%state%cfg%R_ahp**2) cycle
            end if
          end if

          Eki_h = 0.5_real64 * self%state%params%Nm(1) * self%state%chem%mass(1) * &
                (self%state%part(1,iproc_h)%vx(i_h)**2 + &
                 self%state%part(1,iproc_h)%vy(i_h)**2 + &
                 self%state%part(1,iproc_h)%vz(i_h)**2)

          self%state%sum_dEk(iproc_h) = self%state%sum_dEk(iproc_h) + Eki_h
          self%state%Nh(iproc_h)      = self%state%Nh(iproc_h) + 1_int32
        end do
      end do
      !$omp end parallel do
    end if

    ! Deposit particles BEFORE collisions (legacy ordering: calc_rho is
    ! called before collisions in Src/main.f90, so the density collisions
    ! read - np_red/np_mx - is always synchronized with the SAME particle
    ! positions collisions are about to act on, never a step stale. This
    ! single deposit still also feeds next iteration's E/rho/Poisson, same
    ! as before - nothing else deposits again before then.)
    t0 = MPI_Wtime()
    call self%deposit_all_particles()
    t1 = MPI_Wtime()
    self%t_Erho = self%t_Erho + (t1 - t0)

    ! Collisions
    t0 = MPI_Wtime()
    if (self%state%params%nb_step_collisions > 0_int32) then
      if (mod(istep, self%state%params%nb_step_collisions) == 1_int32) then

        ! Legacy-like: build Plist/cell lists from current post-mover particles
        call self%state%sort_particles_local()

        t1 = MPI_Wtime()
        self%t_sort = self%t_sort + (t1 - t0)

        t0 = MPI_Wtime()
        call self%collisions_step()

      end if
    end if
    t1 = MPI_Wtime()
    self%t_MC = self%t_MC + (t1 - t0)

    ! Volume particle injection (legacy: part_injection, inside OMP per iproc)
    ! Fires every step when flag_inj==1 (ns_inj=1 hardcoded, same as legacy).
    if (self%state%cfg%flag_inj == 1_int32) then
      !$omp parallel do private(iproc_h) schedule(static)
      do iproc_h = 1, self%state%nproc
        call inject_particles_volume( &
          part         = self%state%part, &
          iproc        = iproc_h, &
          cfg          = self%state%cfg, &
          n            = int(self%state%dom%n, int32), &
          h            = self%state%dom%h, &
          bcnd         = self%state%dom%bcnd, &
          phi          = self%state%fld%phi, &
          vt0          = self%state%params%vt0, &
          Nm           = self%state%params%Nm, &
          mass         = self%state%chem%mass, &
          ni0          = self%state%chem%ni0, &
          charge       = self%state%chem%charge, &
          ntype        = self%state%ntype, &
          nproc        = self%state%nproc, &
          dt           = self%state%params%dt, &
          zg_sec       = self%state%dom%zg_sec, &
          n_cath       = self%state%dom%n_cath, &
          dir_sec_inout = self%state%dom%dir_sec, &
          tag_neg      = int(self%state%chem%tag_neg,  int32), &
          tag_beam     = int(self%state%chem%tag_beam, int32), &
          flag_pbc     = int(self%state%dom%flag_pbc,  int32), &
          flag_pbcz    = int(self%state%dom%flag_pbcz, int32), &
          iseed        = self%state%params%iseed(iproc_h), &
          N_inj        = self%state%N_inj, &
          sour_xy      = self%state%sour_xy, &
          sour_xz      = self%state%sour_xz, &
          sour_yz      = self%state%sour_yz, &
          iz_pl        = self%state%params%iz_plot_plane, &
          ix_pl        = self%state%params%ix_plot_plane, &
          P_loss       = self%state%P_loss )
      end do
      !$omp end parallel do
      self%state%N_inj = 0_int32
    end if

    ! Flux injection of negative ions off extraction electrode
    ! (legacy: part_flux_injection, every step when jne>0)
    if (self%state%cfg%jne > 0.0_real64) then
      call inject_flux_particles( &
        part         = self%state%part, &
        cfg          = self%state%cfg, &
        n            = int(self%state%dom%n, int32), &
        h            = self%state%dom%h, &
        bcnd         = self%state%dom%bcnd, &
        Nm           = self%state%params%Nm, &
        mass         = self%state%chem%mass, &
        ntype        = self%state%ntype, &
        nproc        = self%state%nproc, &
        mpi_rank     = self%state%mpi_rank, &
        mpi_size     = self%state%mpi_size, &
        dt           = self%state%params%dt, &
        istep        = istep, &
        nsav         = self%state%cfg%nsav, &
        xg1          = self%state%dom%xg1, &
        Lgy          = self%state%dom%Lgy, &
        Lgz          = self%state%dom%Lgz, &
        ymax         = self%state%dom%ymax, &
        zmax         = self%state%dom%zmax, &
        Sg           = self%state%dom%Sg, &
        tag_neg      = int(self%state%chem%tag_neg, int32), &
        tag_neu      = int(self%state%chem%tag_neu, int32), &
        iseed        = self%state%params%iseed, &
        sour_xy      = self%state%sour_xy, &
        sour_xz      = self%state%sour_xz, &
        sour_fx_yz   = self%state%sour_fx_yz, &
        iz_pl        = self%state%params%iz_plot_plane, &
        ix_pl        = self%state%params%ix_plot_plane, &
        P_loss       = self%state%P_loss, &
        N_flx        = self%state%N_flx )
      self%state%N_flx = 0_int32
    end if

    ! Output/write
    t0 = MPI_Wtime()
    if (mod(istep, self%state%cfg%nsav) == 1_int32) then
      call self%output_step(istep)
      call self%reset_2d_averages()
    end if
    t1 = MPI_Wtime()
    self%t_avg = self%t_avg + (t1 - t0)

    ! Backup (legacy: it.gt.1 .and. MOD(it,nbak).eq.1)
    t0 = MPI_Wtime()
    if (self%state%cfg%nbak > 0_int32) then
      if (istep > 1_int32 .and. mod(istep, self%state%cfg%nbak) == 1_int32) then
        if (self%state%mpi_rank == 0_int32) write(*,*) 'Backing up simulation data ...'
        call write_restart_modular( &
          mpi_rank = self%state%mpi_rank, &
          ntype    = self%state%ntype, &
          nproc    = self%state%nproc, &
          tag_neg  = int(self%state%chem%tag_neg, int32), &
          part     = self%state%part, &
          time     = real(istep, real64) * self%state%params%dt, &
          flag_wrt = self%flag_wrt )
        if (self%state%mpi_rank == 0_int32) write(*,*) 'Done!'
      end if
    end if
    t1 = MPI_Wtime()
    self%t_bck = self%t_bck + (t1 - t0)

    ! Total wall time
    t1 = MPI_Wtime()
    self%t_total = self%t_total + (t1 - tstep0)
    self%t_count = self%t_count + 1_int32

  end subroutine advance_one_step


  subroutine run(self, nsteps)
    class(Simulation), intent(inout) :: self
    integer, intent(in) :: nsteps

    integer(int32) :: istep, istep_offset, istep_end

    if (self%state%mpi_rank == 0) then
      write(*,*) " "
      write(*,*) " "
      write(*,"(a)") ">>> Entering the PIC loop <<<"
      write(*,*) " "
    end if

    ! On restart, self%state%time is set from the backup file's stored
    ! physical time (= istep_original * dt). Convert back to a step offset
    ! so istep numbering, timing diagnostics, and backup cadence all
    ! continue from where the original run left off instead of resetting
    ! to 1. For a fresh run self%state%time == 0 so istep_offset == 0.
    istep_offset = nint(self%state%time / self%state%params%dt, int32)
    istep_end    = istep_offset + int(nsteps, int32)

    do istep = istep_offset + 1_int32, istep_end

      if (mod(istep, self%state%cfg%nsav) == 1_int32) then
        call self%print_diagnostics(istep)
      end if

      call self%advance_one_step(istep)

    end do
  end subroutine run


  subroutine finalize(self)
    class(Simulation), intent(inout) :: self

    call self%state%finalize()
  end subroutine finalize

  subroutine print_diagnostics(self,istep)
    class(Simulation), intent(inout) :: self
    integer(int32), intent(in)   :: istep
    integer(int32)               :: ptype
    real(real64)                 :: simulation_time
    real(real64) :: Pwall, Pabs, Pcoll, Pinj
    real(real64) :: Iw1, Iw2
    real(real64) :: Ekw1, Ekw2

    call print_rxn_counts()
    call reset_rxn_counts()
    call print_debug_diagnostics()
    call reset_debug_diagnostics()

    Pabs  = sum(self%state%P_loss(2,:,:)) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Pcoll = sum(self%state%P_loss(3,:,:)) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Pinj  = sum(self%state%P_loss(4,:,:))

    Pwall = -sum(self%state%P_loss(1,:,:)) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Iw1 = sum(self%state%p_mac(1,1,:,:)) / &
          (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Iw2 = sum(self%state%p_mac(2:self%state%ntype,1,:,:)) / &
          (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Ekw1 = 0.0_real64
    if (abs(sum(self%state%p_mac(1,1,:,:))) > 0.0_real64) then
      Ekw1 = sum(self%state%p_mac(1,2,:,:)) / &
            abs(sum(self%state%p_mac(1,1,:,:)))
    end if

    Ekw2 = 0.0_real64
    if (sum(abs(self%state%p_mac(2:self%state%ntype,1,:,:))) > 0.0_real64) then
      Ekw2 = sum(self%state%p_mac(2:self%state%ntype,2,:,:)) / &
            sum(abs(self%state%p_mac(2:self%state%ntype,1,:,:)))
    end if

    write(*,'(a)') " "
    write(*,'(a)') " "

    write(*,'(a,i0,a,F8.3,a)') &
      ' TIME STEP ', istep, ' --> ', &
      self%state%params%dt*istep*1.d6, " us"

    write(*,'(a)') " -------------------------"
    write(*,'(a)') " Particles diagnostics: "

    do ptype = 1_int32, self%state%rxn%ntype-self%state%rxn%n_neu
      write(*,'(a,a,a,i0)') &
        " npart ", self%state%chem%pname(ptype) , &
        " = ", sum(self%state%part(ptype,:)%n)
    end do

    write(*,'(a,i12,a)') " ===== MODULAR DIAGNOSTIC it = ", istep, " ====="
    write(*,'(a,3ES24.15)') " phi min/max/sum = ", &
      minval(self%state%fld%phi), &
      maxval(self%state%fld%phi), &
      sum(self%state%fld%phi)
    write(*,'(a)') " ==========================================="

    write(*,'(a)') " -------------------------"

    write(*,"(a)") " Power diagnostics: "

    write(*,'(A,ES10.2)') ' Pwall (W)  = ', Pwall
    write(*,'(A,ES10.2)') ' Pabs  (W)  = ', Pabs
    write(*,'(A,ES10.2)') ' Pcoll (W)  = ', Pcoll
    write(*,'(A,ES10.2)') ' Pinj  (W)  = ', Pinj

    write(*,'(a,ES10.2,ES10.2)') ' I_w  (A)  = ', Iw1, Iw2
    write(*,'(a,ES10.2,ES10.2)') ' Ek_w (eV) = ', Ekw1, Ekw2

    write(*,'(a)') " -------------------------"

    if (self%t_count > 0_int32) then
      write(*,'(a)') " Timing diagnostics: "

      write(*,'(A,F7.1)') ' E/rho      (ms) = ', &
        1000.0_real64*self%t_Erho/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' poisson    (ms) = ', &
        1000.0_real64*self%t_poisson/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' sorting    (ms) = ', &
        1000.0_real64*self%t_sort/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' avg/write  (ms) = ', &
        1000.0_real64*self%t_avg/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' MC         (ms) = ', &
        1000.0_real64*self%t_MC/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' mover      (ms) = ', &
        1000.0_real64*self%t_mover/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' bck        (ms) = ', &
        1000.0_real64*self%t_bck/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' total      (ms) = ', &
        1000.0_real64*self%t_total/real(self%t_count,real64)

      write(*,'(a)') " -------------------------"
    end if

    call print_poisson_breakdown(self%t_count)

    if (self%t_count > 0_int32) then
      write(*,'(a)') " ----- mover / deposit timing breakdown -----"
      write(*,'(a,f8.2,a,f8.2)') &
        "  mover_push(ms)=", 1000.0_real64*self%t_mover_push/real(self%t_count,real64), &
        "  mover_bc(ms)=",   1000.0_real64*self%t_mover_bc/real(self%t_count,real64)
      write(*,'(a,f8.2,a,f8.2,a,f8.2)') &
        "  dep_clear(ms)=",  1000.0_real64*self%t_dep_clear/real(self%t_count,real64), &
        "  dep_loop(ms)=",   1000.0_real64*self%t_dep_loop/real(self%t_count,real64), &
        "  dep_reduce(ms)=", 1000.0_real64*self%t_dep_reduce/real(self%t_count,real64)
      write(*,'(a)') " ---------------------------------------------"
    end if


    

    write(*,'(a)') "  "
    write(*,*) "  "

    ! Reset legacy-style power accumulators after printing
    if (allocated(self%state%P_loss) .or. allocated(self%state%p_mac)) then
      self%state%P_loss = 0.0_real64
      self%state%p_mac  = 0.0_real64
    end if

    ! Reset timers
    self%t_Erho    = 0.0_real64
    self%t_poisson = 0.0_real64
    call reset_poisson_breakdown()
    self%t_sort    = 0.0_real64
    self%t_avg     = 0.0_real64
    self%t_MC      = 0.0_real64
    self%t_mover   = 0.0_real64
    self%t_mover_push = 0.0_real64
    self%t_mover_bc   = 0.0_real64
    self%t_dep_clear  = 0.0_real64
    self%t_dep_loop   = 0.0_real64
    self%t_dep_reduce = 0.0_real64
    self%t_bck     = 0.0_real64
    self%t_total   = 0.0_real64
    self%t_count   = 0_int32

  end subroutine print_diagnostics

end module mod_simulation