module mod_simulation
  use iso_fortran_env, only: int32, real64, int64

  use mod_state,                 only: State
  use mod_density,               only: reduce_species_density, build_rho_from_np
  use mod_chargeDeposition,      only: clear_np_thread, deposit_particle_set_to_np_thread
  use mod_PoissonSolver_legacy,  only: solve_poisson_legacy
  use mod_electricField,         only: calc_Efield_modular
  use mod_output_2d,             only: write_density_planes, write_scalar_planes
  use mod_collisions, only: CollisionWorkspace, perform_collisions_step
  use mpi

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

    call build_rho_from_np( &
         n      = int(self%state%dom%n, int32), &
         np_red = self%state%fld%np, &
         charge = self%state%chem%charge(1:self%state%ntype), &
         ntype  = int(self%state%ntype), &
         rho    = self%state%fld%rho )

    call solve_poisson_legacy( &
         pdec        = self%state%pdec, &
         phi_global  = self%state%fld%phi, &
         bcnd_global = self%state%dom%bcnd, &
         rhs_global  = self%state%fld%rho(0:self%state%dom%n(1)+1, &
                                          0:self%state%dom%n(2)+1, &
                                          0:self%state%dom%n(3)+1), &
         h           = self%state%dom%h, &
         n_in        = int(self%state%dom%n, int32), &
         ncycl       = 100, &
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

      allocate(Ex3(0:self%state%dom%n(1)+2, 0:self%state%dom%n(2)+2, 0:self%state%dom%n(3)+2))
      allocate(Ey3(0:self%state%dom%n(1)+2, 0:self%state%dom%n(2)+2, 0:self%state%dom%n(3)+2))
      allocate(Ez3(0:self%state%dom%n(1)+2, 0:self%state%dom%n(2)+2, 0:self%state%dom%n(3)+2))

      Ex3 = self%state%fld%E(1,:,:,:)
      Ey3 = self%state%fld%E(2,:,:,:)
      Ez3 = self%state%fld%E(3,:,:,:)

      call write_scalar_planes( &
        f        = Ex3, &
        n        = int(self%state%dom%n, int32), &
        ix_plane = self%state%params%ix_plot_plane, &
        iy_plane = iy_plane_E, &
        iz_plane = self%state%params%iz_plot_plane, &
        every    = 1_int32, &
        prefix   = '../Output/Output_2D/Ex' )

      call write_scalar_planes( &
        f        = Ey3, &
        n        = int(self%state%dom%n, int32), &
        ix_plane = self%state%params%ix_plot_plane, &
        iy_plane = iy_plane_E, &
        iz_plane = self%state%params%iz_plot_plane, &
        every    = 1_int32, &
        prefix   = '../Output/Output_2D/Ey' )

      call write_scalar_planes( &
        f        = Ez3, &
        n        = int(self%state%dom%n, int32), &
        ix_plane = self%state%params%ix_plot_plane, &
        iy_plane = iy_plane_E, &
        iz_plane = self%state%params%iz_plot_plane, &
        every    = 1_int32, &
        prefix   = '../Output/Output_2D/Ez' )

      deallocate(Ex3, Ey3, Ez3)
    end if
  end subroutine write_initial_diagnostics


  subroutine deposit_all_particles(self)
    class(Simulation), intent(inout) :: self

    integer(int32) :: ptype, iproc
    real(real64), allocatable :: np_thread(:,:,:,:,:)

    allocate(np_thread(0:self%state%dom%n(1)+2, &
                       0:self%state%dom%n(2)+2, &
                       0:self%state%dom%n(3)+2, &
                       self%state%ntype, self%state%nproc))

    call clear_np_thread(int(self%state%dom%n, int32), self%state%ntype, self%state%nproc, np_thread)

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
             np_local   = np_thread(:,:,:,ptype,iproc) )
      end do
    end do

    call reduce_species_density( &
         n         = int(self%state%dom%n, int32), &
         bcnd      = self%state%dom%bcnd, &
         np_thread = np_thread, &
         ntype     = int(self%state%ntype), &
         nproc     = int(self%state%nproc), &
         mpi_comm  = self%state%comm, &
         np_red    = self%state%fld%np )

    deallocate(np_thread)
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

  subroutine output_step(self, istep)
    class(Simulation), intent(inout) :: self
    integer(int32),    intent(in)    :: istep

    integer(int32) :: i
    integer(int32) :: iy_plane_density, iy_plane_phi, iy_plane_E
    character(len=256) :: prefix
    character(len=16)  :: sstep, sspecies
    real(real64), allocatable :: Ex3(:,:,:), Ey3(:,:,:), Ez3(:,:,:)

    if (self%state%mpi_rank /= 0) return

    iy_plane_density = int(self%state%dom%n(2)/2 + 1, int32)
    iy_plane_phi     = int(self%state%dom%n(2)/2,     int32)
    iy_plane_E       = int(self%state%dom%n(2)/2 + 1, int32)

    write(sstep,'(i0)') istep

    ! ------------------------------------------------------------
    ! Species densities
    ! ------------------------------------------------------------
    do i = 1, self%state%ntype
      write(sspecies,'(i0)') i
      prefix = '../Output/Output_2D/it' // trim(sstep) // '_n' // trim(sspecies)

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

    ! ------------------------------------------------------------
    ! Potential
    ! ------------------------------------------------------------
    prefix = '../Output/Output_2D/it' // trim(sstep) // '_phi'

    call write_scalar_planes( &
      f        = self%state%fld%phi, &
      n        = int(self%state%dom%n, int32), &
      ix_plane = self%state%params%ix_plot_plane, &
      iy_plane = iy_plane_phi, &
      iz_plane = self%state%params%iz_plot_plane, &
      every    = 1_int32, &
      prefix   = prefix )

    ! ------------------------------------------------------------
    ! Electric field components
    ! ------------------------------------------------------------
    allocate(Ex3(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))
    allocate(Ey3(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))
    allocate(Ez3(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))

    Ex3 = self%state%fld%E(1,:,:,:)
    Ey3 = self%state%fld%E(2,:,:,:)
    Ez3 = self%state%fld%E(3,:,:,:)

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_Ex'
    call write_scalar_planes( &
      f        = Ex3, &
      n        = int(self%state%dom%n, int32), &
      ix_plane = self%state%params%ix_plot_plane, &
      iy_plane = iy_plane_E, &
      iz_plane = self%state%params%iz_plot_plane, &
      every    = 1_int32, &
      prefix   = prefix )

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_Ey'
    call write_scalar_planes( &
      f        = Ey3, &
      n        = int(self%state%dom%n, int32), &
      ix_plane = self%state%params%ix_plot_plane, &
      iy_plane = iy_plane_E, &
      iz_plane = self%state%params%iz_plot_plane, &
      every    = 1_int32, &
      prefix   = prefix )

    prefix = '../Output/Output_2D/it' // trim(sstep) // '_Ez'
    call write_scalar_planes( &
      f        = Ez3, &
      n        = int(self%state%dom%n, int32), &
      ix_plane = self%state%params%ix_plot_plane, &
      iy_plane = iy_plane_E, &
      iz_plane = self%state%params%iz_plot_plane, &
      every    = 1_int32, &
      prefix   = prefix )

    deallocate(Ex3, Ey3, Ez3)

  end subroutine output_step

  subroutine advance_one_step(self, istep)
    class(Simulation), intent(inout) :: self
    integer(int32),    intent(in)    :: istep

    ! ------------------------------------------------------------
    ! 1) Sort particles on legacy cadence
    ! Legacy: MOD(it,nsort) == 1
    ! ------------------------------------------------------------
    if (mod(istep, self%state%params%nb_step_sort) == 1_int32) then
      call self%state%sort_particles_local()
    end if

    ! ------------------------------------------------------------
    ! 2) Collisions before mover
    ! Legacy: it > 1 and MOD(it,ns_coll) == 1
    ! ------------------------------------------------------------
    if (istep > 1_int32) then

      if (mod(istep, self%state%params%nb_step_collisions) == 1_int32) then
        if(self%state%mpi_rank == 0) write (*,'(a)') "Performing collisions"
        call self%collisions_step()
      end if
    end if

    if (mod(istep, self%state%cfg%nsav) == 1_int32) then
      call self%output_step(istep)
    end if
    
    call self%state%move_particles_local()
    call self%state%apply_particle_bc_local()
    call self%deposit_all_particles()

    call build_rho_from_np( &
        n      = int(self%state%dom%n, int32), &
        np_red = self%state%fld%np, &
        charge = self%state%chem%charge(1:self%state%ntype), &
        ntype  = int(self%state%ntype), &
        rho    = self%state%fld%rho )
        
    call solve_poisson_legacy( &
        pdec        = self%state%pdec, &
        phi_global  = self%state%fld%phi, &
        bcnd_global = self%state%dom%bcnd, &
        rhs_global  = self%state%fld%rho(0:self%state%dom%n(1)+1, &
                                          0:self%state%dom%n(2)+1, &
                                          0:self%state%dom%n(3)+1), &
        h           = self%state%dom%h, &
        n_in        = int(self%state%dom%n, int32), &
        ncycl       = 100, &
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

    if (mod(istep, self%state%cfg%nsav) == 1_int32) then
      call self%output_step(istep)
    end if

  end subroutine advance_one_step


  subroutine run(self, nsteps)
    class(Simulation), intent(inout) :: self
    integer,           intent(in)    :: nsteps
    integer(int32) :: istep

    call self%build_initial_fields()
    call self%write_initial_diagnostics()
    if(self%state%mpi_rank == 0) then
      write (*,*) " "
      write (*,'(a)') 'Entering the PIC loop'
    end if
    do istep = 1_int32, int(nsteps, int32)
      if(self%state%mpi_rank == 0) write(*,'(a,i0)') "Iteration # ", istep
      call self%advance_one_step(istep)
    end do
  end subroutine run


  subroutine finalize(self)
    class(Simulation), intent(inout) :: self
    call self%state%finalize()
  end subroutine finalize

end module mod_simulation