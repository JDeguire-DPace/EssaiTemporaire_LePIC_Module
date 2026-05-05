module mod_simulation
  use iso_fortran_env, only: int32, real64, int64

  use mod_state,                 only: State
  use mod_density,               only: reduce_species_density, build_rho_from_np, build_rho_from_np_thread
  use mod_chargeDeposition,      only: clear_np_thread, deposit_particle_set_to_np_thread
  use mod_PoissonSolver_legacy,  only: solve_poisson_legacy
  use mod_electricField,         only: calc_Efield_modular
  use mod_output_2d, only: write_density_planes, write_scalar_planes, &
                           write_plane_xy_scalar_2d, write_plane_xz_scalar_2d, write_plane_yz_scalar_2d
  use mod_collisions, only: CollisionWorkspace, perform_collisions_step
  use mpi
  use omp_lib, only: omp_get_wtime
  use mod_constants, only: eps0
  use mod_particleBC, only: dbg_loss_xright_s1, dbg_loss_xright_s2

  implicit none
  private
  public :: Simulation

  type :: Simulation
    type(State) :: state
    type(CollisionWorkspace) :: coll_ws
  contains
    procedure :: init
    procedure :: build_initial_fields
    procedure :: write_initial_diagnostics
    procedure :: deposit_all_particles
    procedure :: collisions_step
    procedure :: output_step
    procedure :: reset_2d_averages
    procedure :: advance_one_step
    procedure :: run
    procedure :: finalize
  end type Simulation

contains

  subroutine init(self, comm_in)
    class(Simulation), intent(inout) :: self
    integer,           intent(in)    :: comm_in

    call self%state%init(comm_in)
  end subroutine init


  subroutine build_initial_fields(self)
    class(Simulation), intent(inout) :: self

    ! Keep this available for optional initial diagnostics. The main PIC loop
    ! rebuilds rho from np_thread before every Poisson solve.
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

    if (self%state%mpi_rank == 0) then
      print *, "MOD charges = ", self%state%chem%charge(1:self%state%ntype)
    end if
  end subroutine build_initial_fields


  subroutine write_initial_diagnostics(self)
    class(Simulation), intent(in) :: self
    character(len=128) :: prefix
    character(len=10)  :: s
    integer(int32) :: i
    integer(int32) :: iy_plane_phi, iy_plane_density, iy_plane_E
    real(real64), allocatable :: Ex3(:,:,:), Ey3(:,:,:), Ez3(:,:,:)

    iy_plane_density = int(self%state%dom%n(2)/2 + 1, int32)
    iy_plane_phi     = int(self%state%dom%n(2)/2,     int32)
    iy_plane_E       = int(self%state%dom%n(2)/2 + 1, int32)

    if (self%state%mpi_rank == 0) then
      do i = 1, self%state%ntype
        write(s,'(i0)') i
        prefix = '../Output/Output_2D/n'//trim(s)

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

      allocate(Ex3(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))
      allocate(Ey3(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))
      allocate(Ez3(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))

      Ex3 = self%state%fld%E(1,:,:,:)
      Ey3 = self%state%fld%E(2,:,:,:)
      Ez3 = self%state%fld%E(3,:,:,:)

      call write_scalar_planes(Ex3, int(self%state%dom%n, int32), &
                               self%state%params%ix_plot_plane, iy_plane_E, &
                               self%state%params%iz_plot_plane, 1_int32, '../Output/Output_2D/Ex')

      call write_scalar_planes(Ey3, int(self%state%dom%n, int32), &
                               self%state%params%ix_plot_plane, iy_plane_E, &
                               self%state%params%iz_plot_plane, 1_int32, '../Output/Output_2D/Ey')

      call write_scalar_planes(Ez3, int(self%state%dom%n, int32), &
                               self%state%params%ix_plot_plane, iy_plane_E, &
                               self%state%params%iz_plot_plane, 1_int32, '../Output/Output_2D/Ez')

      deallocate(Ex3, Ey3, Ez3)
    end if
  end subroutine write_initial_diagnostics


  subroutine deposit_all_particles(self)
    class(Simulation), intent(inout) :: self
    integer(int32) :: ptype, iproc
    integer :: ix

    call clear_np_thread(int(self%state%dom%n, int32), &
                         self%state%ntype, &
                         self%state%nproc, &
                         self%state%np_thread)

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
        part       = self%state%part, &
        n          = int(self%state%dom%n, int32), &
        h          = self%state%dom%h, &
        np_red     = self%state%fld%np, &
        mass       = self%state%chem%mass(1:self%state%ntype), &
        charge     = self%state%chem%charge(1:self%state%ntype), &
        vt0        = self%state%params%vt0(1:self%state%ntype), &
        Nm         = self%state%params%Nm(1:self%state%ntype), &
        p_ncol     = self%state%chem%p_ncol(1:self%state%ntype), &
        sig        = self%state%rxn%sig, &
        sig_Er     = self%state%rxn%sig_Er, &
        sig_list   = self%state%rxn%sig_list, &
        sig_Eex    = self%state%rxn%sig_Eex, &
        col_info   = self%state%rxn%col_info, &
        sigv_mx    = self%state%rxn%sigv_mx, &
        ns_coll    = self%state%params%nb_step_collisions, &
        dt         = self%state%params%dt, &
        nu_uplim   = self%state%params%nu_uplim(1:self%state%ntype), &
        iseed      = self%state%params%iseed, &
        nproc_mpi  = int(self%state%mpi_size, int32), &
        mpi_rank   = int(self%state%mpi_rank, int32), &
        workspace  = self%coll_ws )
  end subroutine collisions_step


  subroutine reset_2d_averages(self)
    class(Simulation), intent(inout) :: self

    if (allocated(self%state%np_avg_xy))  self%state%np_avg_xy  = 0.0_real64
    if (allocated(self%state%np_avg_xz))  self%state%np_avg_xz  = 0.0_real64
    if (allocated(self%state%np_avg_yz))  self%state%np_avg_yz  = 0.0_real64

    if (allocated(self%state%phi_avg_xy)) self%state%phi_avg_xy = 0.0_real64
    if (allocated(self%state%phi_avg_xz)) self%state%phi_avg_xz = 0.0_real64
    if (allocated(self%state%phi_avg_yz)) self%state%phi_avg_yz = 0.0_real64

    self%state%cnt_avg = 0_int32
  end subroutine reset_2d_averages


  subroutine output_step(self, istep)
    class(Simulation), intent(inout) :: self
    integer(int32),    intent(in)    :: istep

    integer(int32) :: i
    integer(int32) :: ix_plane, iy_plane_density, iy_plane_phi, iy_plane_E, iz_plane
    character(len=256) :: prefix
    character(len=16)  :: sstep, sspecies
    real(real64) :: avg_factor
    real(real64), allocatable :: Ex3(:,:,:), Ey3(:,:,:), Ez3(:,:,:)
    real(real64), allocatable :: tmp_xy(:,:), tmp_xz(:,:), tmp_yz(:,:)

    if (self%state%mpi_rank /= 0) return

    ix_plane         = self%state%params%ix_plot_plane
    iz_plane         = self%state%params%iz_plot_plane
    iy_plane_density = int(self%state%dom%n(2)/2 + 1, int32)
    iy_plane_phi     = int(self%state%dom%n(2)/2,     int32)
    iy_plane_E       = int(self%state%dom%n(2)/2 + 1, int32)

    avg_factor = real(max(1_int32, self%state%cnt_avg), real64)

    allocate(tmp_xy(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2))
    allocate(tmp_xz(0:self%state%dom%n(1)+2,0:self%state%dom%n(3)+2))
    allocate(tmp_yz(0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))

    write(sstep,'(i0)') istep

    do i = 1, self%state%ntype
      write(sspecies,'(i0)') i

      ! Legacy-style averaged density planes.
      tmp_xy = self%state%np_avg_xy(:,:,i) / avg_factor
      tmp_xz = self%state%np_avg_xz(:,:,i) / avg_factor
      tmp_yz = self%state%np_avg_yz(:,:,i) / avg_factor

      prefix = '../Output/Output_2D/it' // trim(sstep) // '_n' // trim(sspecies)
      if (istep == 2001_int32 .and. self%state%mpi_rank == 0) then
        print *, "OUTPUT AVG n1 x tail xy = ", tmp_xy(125:129,9)
        print *, "OUTPUT AVG n1 x tail xz = ", tmp_xz(125:129,9)
      end if
      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

      write(prefix,'(a,"/it",i0,"_T",i0,"_xy.mco")') '../Output/Output_2D', istep, i
      call write_plane_xy_scalar_2d(prefix, self%state%data_pavg_xy(2,:,:,i), int(self%state%dom%n, int32), 1_int32)

      write(prefix,'(a,"/it",i0,"_T",i0,"_xz.mco")') '../Output/Output_2D', istep, i
      call write_plane_xz_scalar_2d(prefix, self%state%data_pavg_xz(2,:,:,i), int(self%state%dom%n, int32), 1_int32)

      write(prefix,'(a,"/it",i0,"_T",i0,"_yz.mco")') '../Output/Output_2D', istep, i
      call write_plane_yz_scalar_2d(prefix, self%state%data_pavg_yz(2,:,:,i), int(self%state%dom%n, int32), 1_int32)
    end do

    ! Legacy-style averaged potential planes.
    tmp_xy = self%state%phi_avg_xy / avg_factor
    tmp_xz = self%state%phi_avg_xz / avg_factor
    tmp_yz = self%state%phi_avg_yz / avg_factor

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_phi'

    call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
    call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
    call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

    allocate(Ex3(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))
    allocate(Ey3(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))
    allocate(Ez3(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))

    Ex3 = self%state%fld%E(1,:,:,:)
    Ey3 = self%state%fld%E(2,:,:,:)
    Ez3 = self%state%fld%E(3,:,:,:)

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_Ex'
    call write_scalar_planes(Ex3, int(self%state%dom%n, int32), &
                            ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_Ey'
    call write_scalar_planes(Ey3, int(self%state%dom%n, int32), &
                            ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_Ez'
    call write_scalar_planes(Ez3, int(self%state%dom%n, int32), &
                            ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)

    write(*,*) "cnt_avg MOD = ", self%state%cnt_avg

    deallocate(Ex3, Ey3, Ez3)
    deallocate(tmp_xy, tmp_xz, tmp_yz)

  end subroutine output_step


  subroutine advance_one_step(self, istep)
    class(Simulation), intent(inout) :: self
    integer(int32),    intent(in)    :: istep
    integer :: imax(3), ix

    real(real64) :: vt_heat

    if (mod(istep, self%state%params%nb_step_sort) == 1_int32) then
      call self%state%sort_particles_local()
    end if

    ! if (istep > 1_int32) then
    !   if (mod(istep, self%state%params%nb_step_collisions) == 1_int32) then
    !     call self%collisions_step()
    !   end if
    ! end if

    ! if (mod(istep, self%state%params%nb_step_heating) == 0_int32) then
    !   call self%state%compute_heating_region_moments()
    !   call self%state%update_heating_vt(vt_heat)
    !   call self%state%apply_electron_heating_local(vt_heat)
    ! end if

    call reduce_species_density( &
        n         = int(self%state%dom%n, int32), &
        bcnd      = self%state%dom%bcnd, &
        np_thread = self%state%np_thread, &
        ntype     = int(self%state%ntype), &
        nproc     = int(self%state%nproc), &
        mpi_comm  = self%state%comm, &
        np_red    = self%state%fld%np )

    call build_rho_from_np_thread( &
        n         = int(self%state%dom%n, int32), &
        np_thread = self%state%np_thread, &
        charge    = self%state%chem%charge(1:self%state%ntype), &
        ntype     = int(self%state%ntype), &
        nproc     = int(self%state%nproc), &
        rho       = self%state%fld%rho, &
        bcnd      = self%state%dom%bcnd, &
        flag_pbc  = int(self%state%dom%flag_pbc, int32) )

    if (istep == 1_int32 .and. self%state%mpi_rank == 0) then
      print *, "MOD rho min/max/sum = ", minval(self%state%fld%rho), &
                                        maxval(self%state%fld%rho), &
                                        sum(self%state%fld%rho)
      print *, "MOD rho center = ", self%state%fld%rho(65,9,9)
    end if

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

    if (istep == 1_int32 .and. self%state%mpi_rank == 0) then
      write(*,*) "phi min/max/sum = ", minval(self%state%fld%phi), &
                                        maxval(self%state%fld%phi), &
                                        sum(self%state%fld%phi)
    end if

    call calc_Efield_modular( &
        n    = int(self%state%dom%n, int32), &
        h    = self%state%dom%h, &
        phi  = self%state%fld%phi, &
        E    = self%state%fld%E, &
        bcnd = self%state%dom%bcnd )

    if (self%state%params%nb_step_averaging > 0) then
      if (mod(istep, self%state%params%nb_step_averaging) == 1_int32) then
        call self%state%accumulate_2d_averages()
      end if
    end if

    if (allocated(self%state%sum_q_xz)) self%state%sum_q_xz = 0.0_real64
    if (allocated(self%state%sum_q_yz)) self%state%sum_q_yz = 0.0_real64

    call self%state%move_particles_local()
    call self%state%apply_particle_bc_local()

    if (mod(istep, 250_int32) == 1_int32 .and. self%state%mpi_rank == 0) then
      print *, 'CHECK NPART after BC:'
      print *, ' electrons = ', sum(self%state%part(1,:)%n)
      print *, ' ions      = ', sum(self%state%part(2,:)%n)
    end if

    call self%deposit_all_particles()

    if (istep == 2001_int32 .and. self%state%mpi_rank == 0) then
      write(*,*) "==== XMAX DENSITY CHECK AFTER DEPOSIT ===="
      print *, "MOD loss xright e/i = ", dbg_loss_xright_s1, dbg_loss_xright_s2

      write(*,*) "nx = ", self%state%dom%n(1)
      write(*,*) "xmax = ", self%state%dom%xmax
      write(*,*) "h(1) = ", self%state%dom%h(1)

      write(*,*) "electron max x = ", maxval(self%state%part(1,1)%x(1:self%state%part(1,1)%n))
      write(*,*) "ion max x      = ", maxval(self%state%part(2,1)%x(1:self%state%part(2,1)%n))

      write(*,*) "electron count x > xmax-h = ", &
          count(self%state%part(1,1)%x(1:self%state%part(1,1)%n) > &
          self%state%dom%xmax - self%state%dom%h(1))

      write(*,*) "ion count x > xmax-h = ", &
          count(self%state%part(2,1)%x(1:self%state%part(2,1)%n) > &
          self%state%dom%xmax - self%state%dom%h(1))

      write(*,*) "ix, sum np1(ix,:,:), sum np2(ix,:,:), kq(ix,9,9), bcnd(ix,9,9)"
      do ix = self%state%dom%n(1)-3, self%state%dom%n(1)+2
        write(*,'(i6,2(1x,es16.8),1x,es16.8,1x,i6)') ix, &
          sum(self%state%fld%np(ix,:,:,1)), &
          sum(self%state%fld%np(ix,:,:,2)), &
          self%state%fld%kq(ix,9,9), &
          self%state%dom%bcnd(ix,9,9)
      end do

      write(*,*) "=========================================="
    end if

    if (mod(istep, self%state%cfg%nsav) == 1_int32) then
      call self%output_step(istep)
      call self%reset_2d_averages()
    end if

  end subroutine advance_one_step


  subroutine run(self, nsteps)
    class(Simulation), intent(inout) :: self
    integer,           intent(in)    :: nsteps
    integer(int32) :: istep,ix

    ! Optional initial fields/planes. For parity runs, you may comment these out
    ! if legacy only writes after the PIC loop.
    ! call self%build_initial_fields()
    ! call self%write_initial_diagnostics()

    if (self%state%mpi_rank == 0) then
      print *, 'MOD eps0 = ', eps0
      print *, 'MOD k_eps0 = ', self%state%cfg%k_eps0
      print *, "MOD ix_plot_plane = ", self%state%params%ix_plot_plane
      print *, "MOD iz_plot_plane = ", self%state%params%iz_plot_plane
      write (*,*) " "
      write (*,'(a)') 'Entering the PIC loop'
    end if
    print *, "DEBUG eps omega ng nsav navg = ", self%state%cfg%eps, self%state%cfg%omega, self%state%cfg%ng, self%state%cfg%nsav, self%state%params%nb_step_averaging
    do istep = 1_int32, int(nsteps, int32)
      call self%advance_one_step(istep)
      if (istep == 2001) then
        write(*,*) "phi xmax line MOD:"
        do ix = self%state%dom%n(1)-3, self%state%dom%n(1)+1
          write(*,'(i6,1x,es16.8)') ix, self%state%fld%phi(ix,9,9)
        end do

        write(*,*) "Ex xmax line MOD:"
        do ix = self%state%dom%n(1)-3, self%state%dom%n(1)+1
          write(*,'(i6,1x,es16.8)') ix, self%state%fld%E(1,ix,9,9)
        end do

        write(*,*) "rho xmax line MOD:"
        do ix = self%state%dom%n(1)-3, self%state%dom%n(1)+1
          write(*,'(i6,1x,es16.8)') ix, self%state%fld%rho(ix,9,9)
        end do

        write(*,*) "np1/np2 center xmax line MOD:"
        do ix = self%state%dom%n(1)-3, self%state%dom%n(1)+1
          write(*,'(i6,2(1x,es16.8))') ix, &
            self%state%fld%np(ix,9,9,1), &
            self%state%fld%np(ix,9,9,2)
        end do

      end if
    end do
  end subroutine run


  subroutine finalize(self)
    class(Simulation), intent(inout) :: self
    call self%state%finalize()
  end subroutine finalize

end module mod_simulation