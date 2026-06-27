module mod_simulation
  use iso_fortran_env, only: int32, real64

  use mod_state,                 only: State
  use mod_density,               only: reduce_species_density, build_rho_from_np
  use mod_chargeDeposition,      only: clear_np_thread, deposit_particle_set_to_np_thread
  use mod_PoissonSolver_legacy,  only: solve_poisson_legacy
  use mod_electricField,         only: calc_Efield_modular
  use mod_output_2d,             only: write_density_planes, write_scalar_planes, &
                                       write_vector_component_planes, &
                                       write_plane_xy_scalar_2d, write_plane_xz_scalar_2d, &
                                       write_plane_yz_scalar_2d
  use mod_collisions,            only: CollisionWorkspace, perform_collisions_step
  use mod_constants,             only: eps0
  use mpi

  implicit none
  private
  public :: Simulation

  logical, parameter :: DEBUG_SIMULATION = .false.

  type :: Simulation
    type(State) :: state
    type(CollisionWorkspace) :: coll_ws

    real(real64) :: t_Erho    = 0.0_real64
    real(real64) :: t_poisson = 0.0_real64
    real(real64) :: t_sort    = 0.0_real64
    real(real64) :: t_avg     = 0.0_real64
    real(real64) :: t_MC      = 0.0_real64
    real(real64) :: t_mover   = 0.0_real64
    real(real64) :: t_bck     = 0.0_real64
    real(real64) :: t_total   = 0.0_real64
    integer(int32) :: t_count = 0_int32

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
      prefix = '../Output/Output_2D/n' // trim(s)

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
      prefix   = '../Output/Output_2D/phi1' )

    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      1_int32, self%state%params%ix_plot_plane, &
                                      iy_plane_E, self%state%params%iz_plot_plane, &
                                      1_int32, '../Output/Output_2D/Ex')

    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      2_int32, self%state%params%ix_plot_plane, &
                                      iy_plane_E, self%state%params%iz_plot_plane, &
                                      1_int32, '../Output/Output_2D/Ey')

    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      3_int32, self%state%params%ix_plot_plane, &
                                      iy_plane_E, self%state%params%iz_plot_plane, &
                                      1_int32, '../Output/Output_2D/Ez')
  end subroutine write_initial_diagnostics


  subroutine deposit_all_particles(self)
    class(Simulation), intent(inout) :: self

    integer(int32) :: ptype, iproc

    call clear_np_thread( &
      int(self%state%dom%n, int32), &
      self%state%ntype, &
      self%state%nproc, &
      self%state%np_thread )

    !$omp parallel do collapse(2) private(ptype,iproc) schedule(static)
    do ptype = 1, self%state%ntype
      do iproc = 1, self%state%nproc
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

    call reduce_species_density( &
      n         = int(self%state%dom%n, int32), &
      bcnd      = self%state%dom%bcnd, &
      np_thread = self%state%np_thread, &
      ntype     = int(self%state%ntype), &
      nproc     = int(self%state%nproc), &
      mpi_comm  = self%state%comm, &
      np_red    = self%state%fld%np )
  end subroutine deposit_all_particles


  subroutine collisions_step(self)
    class(Simulation), intent(inout) :: self

    call perform_collisions_step( &
      part            = self%state%part, &
      n               = int(self%state%dom%n, int32), &
      h               = self%state%dom%h, &
      np_red          = self%state%fld%np, &
      mass            = self%state%chem%mass(1:self%state%rxn%ntype), &
      charge          = self%state%chem%charge(1:self%state%rxn%ntype), &
      vt0             = self%state%params%vt0(1:self%state%rxn%ntype), &
      Nm              = self%state%params%Nm(1:self%state%rxn%ntype), &
      neutral_density = self%state%chem%ni0(self%state%ntype+1:self%state%rxn%ntype), &
      p_ncol          = self%state%chem%p_ncol(1:self%state%rxn%ntype), &
      sig             = self%state%rxn%sig, &
      sig_Er          = self%state%rxn%sig_Er, &
      sig_list        = self%state%rxn%sig_list, &
      sig_Eex         = self%state%rxn%sig_Eex, &
      col_info        = self%state%rxn%col_info, &
      sigv_mx         = self%state%rxn%sigv_mx, &
      ns_coll         = self%state%params%nb_step_collisions, &
      dt              = self%state%params%dt, &
      nu_uplim        = self%state%params%nu_uplim(1:self%state%rxn%ntype), &
      iseed           = self%state%params%iseed, &
      nproc_mpi       = int(self%state%mpi_size, int32), &
      mpi_rank        = int(self%state%mpi_rank, int32), &
      workspace       = self%coll_ws, & 
      Pcoll           = self%state%P_loss(3,:,:))
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

      prefix = '../Output/Output_2D/it' // trim(sstep) // '_n' // trim(sspecies)

      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

      tmp_xy = self%state%data_pavg_xy(2,:,:,i) / avg_factor
      tmp_xz = self%state%data_pavg_xz(2,:,:,i) / avg_factor
      tmp_yz = self%state%data_pavg_yz(2,:,:,i) / avg_factor

      prefix = '../Output/Output_2D/it' // trim(sstep) // '_T' // trim(sspecies)

      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)
    end do

    tmp_xy = self%state%phi_avg_xy / avg_factor
    tmp_xz = self%state%phi_avg_xz / avg_factor
    tmp_yz = self%state%phi_avg_yz / avg_factor

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_phi'

    call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
    call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
    call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_Ex'
    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      1_int32, ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_Ey'
    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      2_int32, ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_Ez'
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

    real(real64) :: vt_heat
    real(real64) :: t0, t1, tstep0

    tstep0 = MPI_Wtime()

    ! Sorting
    t0 = MPI_Wtime()
    if (mod(istep, self%state%params%nb_step_sort) == 1_int32) then
      call self%state%sort_particles_local()
    end if
    t1 = MPI_Wtime()
    self%t_sort = self%t_sort + (t1 - t0)

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

    ! Mover + particle BC
    t0 = MPI_Wtime()
    call self%state%move_particles_local()
    call self%state%apply_particle_bc_local()
    t1 = MPI_Wtime()
    self%t_mover = self%t_mover + (t1 - t0)

    ! Collisions
    t0 = MPI_Wtime()
    if (self%state%params%nb_step_collisions > 0_int32) then
      if (mod(istep, self%state%params%nb_step_collisions) == 0_int32) then
        call self%collisions_step()
      end if
    end if
    t1 = MPI_Wtime()
    self%t_MC = self%t_MC + (t1 - t0)

    ! Electron heating
    if (self%state%params%nb_step_heating > 0_int32) then
      if (mod(istep, self%state%params%nb_step_heating) == 0_int32) then
        call self%state%compute_heating_region_moments()
        call self%state%update_heating_vt(vt_heat)
        call self%state%apply_electron_heating_local(vt_heat)
      end if
    end if

    ! Deposit particles after all particle operations
    t0 = MPI_Wtime()
    call self%deposit_all_particles()
    t1 = MPI_Wtime()
    self%t_Erho = self%t_Erho + (t1 - t0)

    ! Output/write
    t0 = MPI_Wtime()
    if (mod(istep, self%state%cfg%nsav) == 1_int32) then
      call self%output_step(istep)
      call self%reset_2d_averages()
    end if
    t1 = MPI_Wtime()
    self%t_avg = self%t_avg + (t1 - t0)

    ! Backup not implemented yet
    self%t_bck = self%t_bck + 0.0_real64

    ! Total wall time
    t1 = MPI_Wtime()
    self%t_total = self%t_total + (t1 - tstep0)
    self%t_count = self%t_count + 1_int32

  end subroutine advance_one_step


  subroutine run(self, nsteps)
    class(Simulation), intent(inout) :: self
    integer, intent(in) :: nsteps

    integer(int32) :: istep,ptype

    if (self%state%mpi_rank == 0) then
      write(*,*) " "
      write(*,*) " "
      write(*,"(a)") ">>> Entering the PIC loop <<<"
      write(*,*) " "
    end if

    do istep = 1_int32, int(nsteps, int32)
      if (mod(istep,self%state%cfg%nsav).eq.1_int32) then
        call self%print_diagnostics(istep)
      end if
      call self%advance_one_step(istep)
    end do
  end subroutine run


  subroutine finalize(self)
    class(Simulation), intent(inout) :: self

    call self%state%finalize()
  end subroutine finalize

  ! subroutine print_diagnostics(self,istep)
  !   class(Simulation), intent(inout) :: self
  !   integer(int32), intent(in)   :: istep
  !   integer(int32)               :: ptype
  !   real(real64)                 :: simulation_time
  !   real(real64) :: Pwall, Pabs, Pcoll, Pinj
  !   real(real64) :: Iw1, Iw2
  !   real(real64) :: Ekw1, Ekw2

  !   Pabs  = sum(self%state%P_loss(2,:,:)) / &
  !       (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)
  !   Pcoll = sum(self%state%P_loss(3,:,:)) * self%state%params%Nm(1) / &
  !       (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)
  !   Pinj  = sum(self%state%P_loss(4,:,:))
  !   Pwall = -sum(self%state%P_loss(1,:,:)) / &
  !       (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)



  !   Iw1 = sum(self%state%p_mac(1,1,:,:)) / &
  !         (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

  !   Iw2 = sum(self%state%p_mac(2:self%state%ntype,1,:,:)) / &
  !         (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

  !   Ekw1 = 0.0_real64
  !   if (abs(sum(self%state%p_mac(1,1,:,:))) > 0.0_real64) then
  !     Ekw1 = sum(self%state%p_mac(1,2,:,:)) / &
  !           abs(sum(self%state%p_mac(1,1,:,:)))
  !   end if

  !   Ekw2 = 0.0_real64
  !   if (sum(abs(self%state%p_mac(2:self%state%ntype,1,:,:))) > 0.0_real64) then
  !     Ekw2 = sum(self%state%p_mac(2:self%state%ntype,2,:,:)) / &
  !           sum(abs(self%state%p_mac(2:self%state%ntype,1,:,:)))
  !   end if


  !   write(*,'(a)') " "
  !   write(*,'(a)') " "
  !   write(*,'(a,i0,a,F8.3,a)') ' TIME STEP ', istep, ' --> ', self%state%params%dt*istep*1.d6, " us"

  !   write(*,'(a)') " -------------------------"
  !   write(*,'(a)') " Particles diagnostics: "
  !   do ptype = 1_int32, self%state%rxn%ntype-self%state%rxn%n_neu
  !     write(*,'(a,a,a,i0)') " npart ", self%state%chem%pname(ptype) , " = ", sum(self%state%part(ptype,:)%n)
  !   end do

  !   write(*,'(a)') " -------------------------"
  !   write(*,"(a)") " Power diagnostics: "
  !   write(*,'(A,ES10.2)')  ' Pwall (W)  = ', Pwall
  !   write(*,'(A,ES10.2)')  ' Pabs  (W)  = ', Pabs
  !   write(*,'(A,ES10.2)')  ' Pcoll (W)  = ', Pcoll
  !   write(*,'(A,ES10.2)')  ' Pinj  (W)  = ', Pinj

  !   write(*,'(a,ES10.2,ES10.2)')  ' I_w  (A)  = ', Iw1, Iw2
  !   write(*,'(a,ES10.2,ES10.2)')  ' Ek_w (eV) = ', Ekw1, Ekw2
  !   write(*,'(a)') " -------------------------"


  !   write(*,'(a)') "  "


  !   ! Reset legacy-style power accumulators after printing
  !   if (allocated(self%state%P_loss) .or. allocated(self%state%p_mac)) then
  !     self%state%P_loss = 0.0_real64
  !     self%state%p_mac  = 0.0_real64
  !   end if

  ! end subroutine print_diagnostics

  subroutine print_diagnostics(self,istep)
    class(Simulation), intent(inout) :: self
    integer(int32), intent(in)   :: istep
    integer(int32)               :: ptype
    real(real64)                 :: simulation_time
    real(real64) :: Pwall, Pabs, Pcoll, Pinj
    real(real64) :: Iw1, Iw2
    real(real64) :: Ekw1, Ekw2

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
    self%t_sort    = 0.0_real64
    self%t_avg     = 0.0_real64
    self%t_MC      = 0.0_real64
    self%t_mover   = 0.0_real64
    self%t_bck     = 0.0_real64
    self%t_total   = 0.0_real64
    self%t_count   = 0_int32

  end subroutine print_diagnostics

end module mod_simulation
