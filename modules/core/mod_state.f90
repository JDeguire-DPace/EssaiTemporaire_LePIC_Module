module mod_state
  use iso_fortran_env, only: int32, real64
  use mpi

  use mod_config,         only: Config
  use mod_readConditions, only: read_input

  use mod_constants,      only: eps0
  use mod_domain,         only: Domain
  use mod_boundary,       only: build_boundary

  use mod_fields,         only: Fields
  use mod_poisson_decomp, only: PoissonDecomp
  use mod_PoissonSolver_legacy, only: solve_poisson_legacy
  use mod_charge_weights, only: build_kq
  use mod_debug_checks,   only: checkpoint_poisson_decomp, checkpoint_kq
  use mod_density,        only: reduce_species_density, build_rho_from_np

  use mod_chemistryState, only: ChemistryState
  use mod_reactionsDB,    only: ReactionsDB
  use mod_part_info,      only: npart
  use mod_particles,      only: ParticleSet
  use mod_particle_loader, only: load_particles_modular
  use mod_magneticField,  only: MagneticField
  use mod_simParams,      only: SimParams
  use mod_output_2d,      only: write_density_planes, write_scalar_planes
  use mod_electricField,  only: calc_Efield_modular

  implicit none
  private
  public :: State

  type :: State
    type(Config)         :: cfg
    type(Domain)         :: dom
    type(Fields)         :: fld
    type(PoissonDecomp)  :: pdec

    type(ChemistryState) :: chem
    type(ReactionsDB)    :: rxn
    type(MagneticField)  :: magField
    type(SimParams)      :: params
    type(ParticleSet), allocatable :: part(:,:)

    integer :: mpi_rank = -1
    integer :: mpi_size = -1
    integer :: comm     = MPI_COMM_NULL
    integer(int32) :: ntype = 0
    integer(int32) :: nproc = 1
  contains
    procedure :: init
    procedure :: build_boundary_only
    procedure :: init_chemistry
    procedure :: init_particles
    procedure :: finalize
  end type State

contains

  subroutine init(self, comm_in)
    class(State), intent(inout) :: self
    integer,      intent(in)    :: comm_in
    integer :: ierr

    self%comm = comm_in
    call MPI_Comm_rank(self%comm, self%mpi_rank, ierr)
    call MPI_Comm_size(self%comm, self%mpi_size, ierr)

    call read_input(self%cfg, self%mpi_rank)

    eps0 = self%cfg%k_eps0 * 8.854187817d-12
    if (self%mpi_rank == 0 .and. self%cfg%k_eps0 /= 1.0_real64) then
      write(*,*) 'eps0 HAS BEEN RE-SCALED: k_eps0 = ', self%cfg%k_eps0
    end if

    call self%dom%init_from_config(self%cfg)
    call self%dom%allocate_masks_domain()

    call self%fld%allocate_from_domain(self%dom)
    call self%fld%zero()

    call self%init_chemistry()
  end subroutine init


  subroutine build_boundary_only(self)
    class(State), intent(inout) :: self
    integer :: ierr
    real(real64), allocatable :: phi_rt(:,:,:)
    real(real64) :: lmax, gmax

    call build_boundary(self%dom, self%cfg, self%fld, self%mpi_rank)
    call self%magField%build_from_cfg(self%cfg, self%dom, self%mpi_rank)
    call self%magField%write_macho_planes('../Output/Output_2D', 1, self%mpi_rank)

    call self%pdec%init(int(self%dom%n(1),int32), int(self%dom%n(2),int32), int(self%dom%n(3),int32), self%comm)
    call self%pdec%scatter_from_global(self%fld%phi, self%dom%bcnd)

    call checkpoint_poisson_decomp(self%mpi_rank, self%comm, &
                                   int(self%dom%n(3),int32), self%pdec%k0, self%pdec%m, &
                                   self%pdec%phi_dom, self%pdec%bcnd_dom)

    call build_kq(self%dom%bcnd, self%fld%kq)
    call checkpoint_kq(self%fld%kq, self%comm, self%mpi_rank, 'after build_kq', self%pdec)

    allocate(phi_rt(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2))
    phi_rt = -999.0_real64

    call self%pdec%gather_phi_to_global(phi_rt)

    lmax = maxval(abs(phi_rt - self%fld%phi))
    call MPI_Allreduce(lmax, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, self%comm, ierr)

    deallocate(phi_rt)

    if (self%mpi_rank == 0) then
      write(*,'(a)') '  '
      write(*,'(a)') 'Building boundary...'
      write(*,'(a,3(i0,1x))')       'n = ', self%dom%n(1), self%dom%n(2), self%dom%n(3)
      write(*,'(a,3(1p,e12.4,1x))') 'h(m) = ', self%dom%h(1), self%dom%h(2), self%dom%h(3)
      write(*,'(a,3(1p,e12.4,1x))') 'box(m) = ', self%dom%xmax, self%dom%ymax, self%dom%zmax
      write(*,'(a,1p,e12.4)')       'Sg(m^2) = ', self%dom%Sg
      write(*,'(a,4(i0,1x))')       'flags(pbc,pbcz,nmn,die) = ', &
           self%dom%flag_pbc, self%dom%flag_pbcz, self%dom%flag_nmn, self%dom%flag_die
      write(*,*) '  '
    end if

    call self%params%build(self%cfg, self%dom, self%chem, self%rxn, self%magField, self%mpi_rank)

    self%nproc = max(1_int32, int(self%cfg%omp_rank_max, int32))
    call self%params%init_seeds(self%nproc, self%mpi_rank)
    call self%params%print_summary(self%mpi_rank, self%cfg%nsav)

    call self%init_particles()
  end subroutine build_boundary_only


  subroutine init_chemistry(self)
    class(State), intent(inout) :: self

    call self%chem%init(npart, self%rxn%ncol_mx, self%rxn%npt_mx)
    self%chem%ngas = self%cfg%ngas
    call self%rxn%load(self%chem, trim(self%cfg%rname), self%mpi_rank)

    if (self%mpi_rank == 0) then
      write(*,'(i0,a,i0,a,i0,a)') self%rxn%ntype, ' species:  (', &
           self%rxn%ntype-self%rxn%n_neu, ' charged and ', self%rxn%n_neu, ' neutral)'
      write(*,'(*(a,1x))') self%chem%pname
      write(*,'(a,i0)') 'Total number of reactions: ', self%chem%ncol
      write(*,*) '  '
    end if
  end subroutine init_chemistry


  subroutine init_particles(self)
    class(State), intent(inout) :: self
    character(len=128) :: prefix
    character(len=10)  :: s
    integer(int32) :: ntype_trk, i
    integer(int32) :: iy_plane_phi, iy_plane_density, iy_plane_E
    real(real64), allocatable :: np_thread(:,:,:,:,:)
    real(real64), allocatable :: Ex3(:,:,:), Ey3(:,:,:), Ez3(:,:,:)

    ntype_trk = int(self%rxn%ntype - self%rxn%n_neu, int32)
    if (ntype_trk < 1_int32) ntype_trk = 1_int32

    self%ntype = ntype_trk
    self%nproc = max(1_int32, int(self%cfg%omp_rank_max, int32))

    if (allocated(self%part)) deallocate(self%part)
    allocate(self%part(self%ntype, self%nproc))

    call self%fld%allocate_species_density(self%dom, self%ntype)

    allocate(np_thread(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2, self%ntype, self%nproc))
    np_thread = 0.0_real64

    call load_particles_modular( &
          cfg       = self%cfg, &
          mpi_rank  = self%mpi_rank, &
          mpi_size  = self%mpi_size, &
          n         = int(self%dom%n, int32), &
          h         = self%dom%h, &
          bcnd      = self%dom%bcnd, &
          kq        = self%fld%kq, &
          vt0       = self%params%vt0, &
          Nm        = self%params%Nm, &
          ni0       = self%chem%ni0(1:self%ntype), &
          iseed     = self%params%iseed, &
          ntype_trk = ntype_trk, &
          part      = self%part, &
          np_thread = np_thread )

    call reduce_species_density( &
         n         = int(self%dom%n, int32), &
         bcnd      = self%dom%bcnd, &
         np_thread = np_thread, &
         ntype     = int(self%ntype), &
         nproc     = int(self%nproc), &
         mpi_comm  = self%comm, &
         np_red    = self%fld%np )

    iy_plane_density = int(self%dom%n(2)/2 + 1, int32)

    if (self%mpi_rank == 0) then
      do i = 1, self%ntype
        write(s,'(i0)') i
        prefix = '../Output/Output_2D/n'//trim(s)

        call write_density_planes( &
          np       = self%fld%np, &
          n        = int(self%dom%n, int32), &
          ptype    = i, &
          ix_plane = self%params%ix_plot_plane, &
          iy_plane = iy_plane_density, &
          iz_plane = self%params%iz_plot_plane, &
          every    = 1_int32, &
          prefix   = prefix )
      end do
    end if

    call build_rho_from_np( &
         n      = int(self%dom%n, int32), &
         np_red = self%fld%np, &
         charge = self%chem%charge(1:self%ntype), &
         ntype  = int(self%ntype), &
         rho    = self%fld%rho )

    call solve_poisson_legacy( &
         pdec        = self%pdec, &
         phi_global  = self%fld%phi, &
         bcnd_global = self%dom%bcnd, &
         rhs_global  = self%fld%rho(0:self%dom%n(1)+1,0:self%dom%n(2)+1,0:self%dom%n(3)+1), &
         h           = self%dom%h, &
         n_in        = int(self%dom%n, int32), &
         ncycl       = 100, &
         eps         = self%cfg%eps, &
         omega       = self%cfg%omega, &
         ng          = self%cfg%ng, &
         flag_pbc_in = self%dom%flag_pbc, &
         flag_nmn_in = self%dom%flag_nmn )

    iy_plane_phi = int(self%dom%n(2)/2, int32)

    if (self%mpi_rank == 0) then
      call write_scalar_planes( &
        f        = self%fld%phi, &
        n        = int(self%dom%n, int32), &
        ix_plane = self%params%ix_plot_plane, &
        iy_plane = iy_plane_phi, &
        iz_plane = self%params%iz_plot_plane, &
        every    = 1_int32, &
        prefix   = '../Output/Output_2D/phi1' )
    end if

    call calc_Efield_modular( &
         n    = int(self%dom%n, int32), &
         h    = self%dom%h, &
         phi  = self%fld%phi, &
         E    = self%fld%E, &
         bcnd = self%dom%bcnd )

    iy_plane_E = int(self%dom%n(2)/2 + 1, int32)

    if (self%mpi_rank == 0) then
      allocate(Ex3(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2))
      allocate(Ey3(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2))
      allocate(Ez3(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2))

      Ex3 = self%fld%E(1,:,:,:)
      Ey3 = self%fld%E(2,:,:,:)
      Ez3 = self%fld%E(3,:,:,:)

      call write_scalar_planes( &
        f        = Ex3, &
        n        = int(self%dom%n, int32), &
        ix_plane = self%params%ix_plot_plane, &
        iy_plane = iy_plane_E, &
        iz_plane = self%params%iz_plot_plane, &
        every    = 1_int32, &
        prefix   = '../Output/Output_2D/Ex' )

      call write_scalar_planes( &
        f        = Ey3, &
        n        = int(self%dom%n, int32), &
        ix_plane = self%params%ix_plot_plane, &
        iy_plane = iy_plane_E, &
        iz_plane = self%params%iz_plot_plane, &
        every    = 1_int32, &
        prefix   = '../Output/Output_2D/Ey' )

      call write_scalar_planes( &
        f        = Ez3, &
        n        = int(self%dom%n, int32), &
        ix_plane = self%params%ix_plot_plane, &
        iy_plane = iy_plane_E, &
        iz_plane = self%params%iz_plot_plane, &
        every    = 1_int32, &
        prefix   = '../Output/Output_2D/Ez' )

      deallocate(Ex3, Ey3, Ez3)
    end if

    deallocate(np_thread)
  end subroutine init_particles


  subroutine finalize(self)
    class(State), intent(inout) :: self
    call self%pdec%destroy()
    call self%fld%destroy()
    call self%rxn%destroy()
    call self%magField%destroy()
    if (allocated(self%part)) deallocate(self%part)
    if (allocated(self%params%iseed)) deallocate(self%params%iseed)
  end subroutine finalize

end module mod_state