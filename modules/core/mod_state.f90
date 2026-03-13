module mod_state
  use iso_fortran_env, only: int32, real64
  use mpi

  use mod_config,         only: Config
  use mod_readConditions, only: read_input

  use mod_domain,         only: Domain
  use mod_boundary,       only: build_boundary

  use mod_fields,         only: Fields
  use mod_poisson_decomp, only: PoissonDecomp
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
  use mod_output_2d, only: write_density_planes

  implicit none
  private
  public :: State

  type :: State
    type(Config)        :: cfg
    type(Domain)        :: dom
    type(Fields)        :: fld
    type(PoissonDecomp) :: pdec

    type(ChemistryState) :: chem
    type(ReactionsDB)    :: rxn
    type(MagneticField)  :: magField
    type(SimParams)      :: params
    type(ParticleSet), allocatable :: part(:,:)   ! per (tracked species, nproc)

    integer :: mpi_rank = -1
    integer :: mpi_size = -1
    integer :: comm     = MPI_COMM_NULL
    integer(int32) :: ntype = 0   ! tracked species count
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

    ! Read config
    call read_input(self%cfg, self%mpi_rank)

    ! Domain setup
    call self%dom%init_from_config(self%cfg)
    call self%dom%allocate_masks_domain()

    ! Fields setup
    call self%fld%allocate_from_domain(self%dom)
    call self%fld%zero()

    ! Chemistry / reactions
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

    ! ------------------------------------------------------------
    ! Poisson Z-slab decomposition (legacy structure)
    ! ------------------------------------------------------------
    call self%pdec%init(int(self%dom%n(1),int32), int(self%dom%n(2),int32), int(self%dom%n(3),int32), self%comm)
    call self%pdec%scatter_from_global(self%fld%phi, self%dom%bcnd)

    call checkpoint_poisson_decomp(self%mpi_rank, self%comm, &
                                   int(self%dom%n(3),int32), self%pdec%k0, self%pdec%m, &
                                   self%pdec%phi_dom, self%pdec%bcnd_dom)

    ! ------------------------------------------------------------
    ! Charge weights kq (global, exact legacy logic)
    ! ------------------------------------------------------------
    call build_kq(self%dom%bcnd, self%fld%kq)
    call checkpoint_kq(self%fld%kq, self%comm, self%mpi_rank, "after build_kq", self%pdec)

    ! ------------------------------------------------------------
    ! Tiny checkpoint: phi round-trip (scatter -> gather)
    ! ------------------------------------------------------------
    allocate(phi_rt(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2))
    phi_rt = -999.0_real64

    call self%pdec%gather_phi_to_global(phi_rt)

    lmax = maxval(abs(phi_rt - self%fld%phi))
    call MPI_Allreduce(lmax, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, self%comm, ierr)

    deallocate(phi_rt)

    ! ------------------------------------------------------------
    ! Print domain summary (rank 0)
    ! ------------------------------------------------------------
    if (self%mpi_rank == 0) then
      write(*,'(a)') "  "
      write(*,'(a)') "Building boundary..."
      write(*,'(a,3(i0,1x))')         "n = ", self%dom%n(1), self%dom%n(2), self%dom%n(3)
      write(*,'(a,3(1p,e12.4,1x))')   "h(m) = ", self%dom%h(1), self%dom%h(2), self%dom%h(3)
      write(*,'(a,3(1p,e12.4,1x))')   "box(m) = ", self%dom%xmax, self%dom%ymax, self%dom%zmax
      write(*,'(a,1p,e12.4)')         "Sg(m^2) = ", self%dom%Sg
      write(*,'(a,4(i0,1x))')         "flags(pbc,pbcz,nmn,die) = ", &
           self%dom%flag_pbc, self%dom%flag_pbcz, self%dom%flag_nmn, self%dom%flag_die
      write(*,*) "  "
    end if

    call self%params%build(self%cfg, self%dom, self%chem, self%rxn, self%magField, self%mpi_rank)

    ! Legacy-equivalent OpenMP lane count from config
    self%nproc = max(1_int32, int(self%cfg%omp_rank_max, int32))

    ! Initialize per-lane seeds in SimParams
    call self%params%init_seeds(self%nproc, self%mpi_rank)

    call self%params%print_summary(self%mpi_rank, self%cfg%nsav)

    call self%init_particles()
  end subroutine build_boundary_only


  subroutine init_chemistry(self)
    class(State), intent(inout) :: self

    ! 1) init chemistry container (species properties)
    call self%chem%init(npart, self%rxn%ncol_mx, self%rxn%npt_mx)

    ! legacy reader expects ngas available via ChemistryState
    self%chem%ngas = self%cfg%ngas

    ! 2) load reaction tables into ReactionsDB
    call self%rxn%load(self%chem, trim(self%cfg%rname), self%mpi_rank)

    if (self%mpi_rank == 0) then
      write(*,'(i0, a,i0,a,i0,a)') self%rxn%ntype, " species:  (", &
           self%rxn%ntype-self%rxn%n_neu, " charged and ", self%rxn%n_neu, " neutral)"
      write(*,'(*(a,1x))') self%chem%pname
      write(*,'(a,i0)') "Total number of reactions: ", self%chem%ncol
      write(*,*) "  "
    end if
  end subroutine init_chemistry


  subroutine init_particles(self)
    class(State), intent(inout) :: self
    character(len=128) :: prefix
    character(len=10)  :: s
    integer(int32) :: ntype_trk,i
    real(real64), allocatable :: np_thread(:,:,:,:,:)

    ! tracked species = charged only
    ntype_trk = int(self%rxn%ntype - self%rxn%n_neu, int32)
    if (ntype_trk < 1_int32) ntype_trk = 1_int32

    ! store tracked count in State
    self%ntype = ntype_trk

    ! nproc already set in build_boundary_only, but keep safe
    self%nproc = max(1_int32, int(self%cfg%omp_rank_max, int32))

    ! allocate tracked particle containers
    if (allocated(self%part)) deallocate(self%part)
    allocate(self%part(self%ntype, self%nproc))

    ! allocate species density field storage
    call self%fld%allocate_species_density(self%dom, self%ntype)

    ! thread-local deposited density
    allocate(np_thread(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2, self%ntype, self%nproc))
    np_thread = 0.0_real64

    ! load particles and build thread-local species density
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

    ! reduce species density over OpenMP threads (+ MPI if needed)
    call reduce_species_density( &
         n         = int(self%dom%n, int32), &
         bcnd      = self%dom%bcnd, &
         np_thread = np_thread, &
         ntype     = int(self%ntype), &
         nproc     = int(self%nproc), &
         mpi_comm  = self%comm, &
         np_red    = self%fld%np )


    if (self%mpi_rank == 0) then

      do i = 1, self%ntype

        write(s,'(i0)') i
        prefix = '../Output/Output_2D/n'//trim(s)

        call write_density_planes( &
          np       = self%fld%np, &
          n        = int(self%dom%n, int32), &
          ptype    = i, &
          ix_plane = self%params%ix_plot_plane, &
          iy_plane = int(self%dom%n(2)/2 + 1, int32), &
          iz_plane = self%params%iz_plot_plane, &
          every    = 1_int32, &
          prefix   = prefix )
      end do

    end if

    ! build charge density rho from species density
    call build_rho_from_np( &
         n      = int(self%dom%n, int32), &
         np_red = self%fld%np, &
         charge = self%chem%charge(1:self%ntype), &
         ntype  = int(self%ntype), &
         rho    = self%fld%rho )
         
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
