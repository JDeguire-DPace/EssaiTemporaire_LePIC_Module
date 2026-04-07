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
  use mod_charge_weights, only: build_kq
  use mod_debug_checks,   only: checkpoint_poisson_decomp, checkpoint_kq
  use mod_density,        only: reduce_species_density

  use mod_chemistryState, only: ChemistryState
  use mod_reactionsDB,    only: ReactionsDB
  use mod_part_info,      only: npart
  use mod_particles,      only: ParticleSet
  use mod_particle_loader,  only: load_particles_modular
  use mod_particle_sorting, only: sort_particles_by_cell, check_particles_are_sorted, check_cell_indexing
  use mod_particleMover,    only: move_particles_electrostatic
  use mod_particleBC,       only: apply_particle_bc_legacy
  use mod_magneticField,    only: MagneticField
  use mod_simParams,        only: SimParams
  use mod_heating,         only: apply_electron_heating

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

    real(real64), allocatable :: np_thread(:,:,:,:,:) 

    integer :: mpi_rank = -1
    integer :: mpi_size = -1
    integer :: comm     = MPI_COMM_NULL

    integer(int32) :: ntype = 0
    integer(int32) :: nproc = 1


    ! P_loss(1,:,:) wall losses
    ! P_loss(2,:,:) heating power
    ! P_loss(3,:,:) collision power exchange
    ! P_loss(4,:,:) injected/secondary power
    real(real64) :: p_loss_heating_local
    real(real64), allocatable :: P_loss(:,:,:)
  contains
    procedure :: init
    procedure :: init_chemistry
    procedure :: init_domain_and_boundary
    procedure :: init_particles

    procedure :: sort_particles_local
    procedure :: move_particles_local
    procedure :: apply_particle_bc_local
    procedure :: apply_electron_heating_local

    procedure :: finalize_particles_only
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
    call self%init_domain_and_boundary()

    call self%params%build(self%cfg, self%dom, self%chem, self%rxn, self%magField, self%mpi_rank)

    self%nproc = max(1_int32, int(self%cfg%omp_rank_max, int32))
    call self%params%init_seeds(self%nproc, self%mpi_rank)
    call self%params%print_summary(self%mpi_rank, self%cfg%nsav)

    call self%init_particles()
    allocate(self%np_thread(0:self%dom%n(1)+2, &
                        0:self%dom%n(2)+2, &
                        0:self%dom%n(3)+2, &
                        self%ntype, self%nproc))
    self%np_thread = 0.0_real64
    allocate(self%P_loss(4, self%ntype, self%nproc))
    self%P_loss = 0.0_real64
  end subroutine init


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


  subroutine init_domain_and_boundary(self)
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
  end subroutine init_domain_and_boundary


  subroutine init_particles(self)
    class(State), intent(inout) :: self
    integer(int32) :: ntype_trk
    real(real64), allocatable :: np_thread(:,:,:,:,:)

    ntype_trk = int(self%rxn%ntype - self%rxn%n_neu, int32)
    if (ntype_trk < 1_int32) ntype_trk = 1_int32

    self%ntype = ntype_trk
    self%nproc = max(1_int32, int(self%cfg%omp_rank_max, int32))

    if (allocated(self%part)) then
      call self%finalize_particles_only()
    end if
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

    call self%sort_particles_local()

    call reduce_species_density( &
         n         = int(self%dom%n, int32), &
         bcnd      = self%dom%bcnd, &
         np_thread = np_thread, &
         ntype     = int(self%ntype), &
         nproc     = int(self%nproc), &
         mpi_comm  = self%comm, &
         np_red    = self%fld%np )

    deallocate(np_thread)
  end subroutine init_particles


  subroutine sort_particles_local(self)
    class(State), intent(inout) :: self
    integer(int32) :: ptype, iproc
    logical :: ok_sorted, ok_cells

    if (.not. allocated(self%part)) return

    do ptype = 1, self%ntype
      do iproc = 1, self%nproc
        if (.not. allocated(self%part(ptype,iproc)%x)) cycle
        if (self%part(ptype,iproc)%n <= 1_int32) cycle

        call sort_particles_by_cell( &
             part = self%part(ptype,iproc), &
             n    = int(self%dom%n, int32), &
             h    = self%dom%h )

        ok_sorted = check_particles_are_sorted(self%part(ptype,iproc))
        if (.not. ok_sorted) then
          write(*,'(a,3(i0,1x))') 'Sorting failed on rank, ptype, iproc = ', &
               self%mpi_rank, ptype, iproc
          error stop 'mod_state%sort_particles_local: particle sorting failed'
        end if

        ok_cells = check_cell_indexing(self%part(ptype,iproc), int(self%dom%n, int32))
        if (.not. ok_cells) then
          write(*,'(a,3(i0,1x))') 'Cell indexing failed on rank, ptype, iproc = ', &
               self%mpi_rank, ptype, iproc
          error stop 'mod_state%sort_particles_local: cell indexing failed'
        end if
      end do
    end do
  end subroutine sort_particles_local

  subroutine apply_electron_heating_local(self, vt)
    class(State), intent(inout) :: self
    real(real64), intent(in)    :: vt

    integer(int32) :: iproc

    if (.not. allocated(self%part)) return
    if (self%ntype < 1) return
    if (.not. allocated(self%P_loss)) return
    if(self%cfg%flag_circxh.eq.1) self%cfg%R_ahp= self%cfg%yr_pow-self%dom%ymax/2.d0
    !$omp parallel do private(iproc) schedule(static)
    do iproc = 1, self%nproc
      if (.not. allocated(self%part(1,iproc)%x)) cycle
      if (self%part(1,iproc)%n <= 0_int32) cycle

      call apply_electron_heating( &
          part           = self%part(1,iproc), &
          h              = self%dom%h, &
          vt             = vt, &
          iseed          = self%params%iseed(iproc), &
          nudt           = self%params%nudt, &
          xl_pow        = self%cfg%xl_pow, &
          xr_pow        = self%cfg%xr_pow, &
          flag_circxh    = self%cfg%flag_circxh, &
          flag_ahp       = self%cfg%flag_ahp, &
          R_ahp          = self%cfg%R_ahp, &
          ymax           = self%dom%ymax, &
          zmax           = self%dom%zmax, &
          Nm_e           = self%params%Nm(1), &
          mass_e         = self%chem%mass(1), &
          p_loss_heating = self%P_loss(2,1,iproc) )
    end do
    !$omp end parallel do

  end subroutine apply_electron_heating_local


  subroutine move_particles_local(self)
    use omp_lib, only: omp_get_max_threads
    class(State), intent(inout) :: self
    integer(int32) :: ptype, iproc
    real(real64)   :: q_species, m_species, dt_local

    if (.not. allocated(self%part)) return

    dt_local = self%params%dt

    !$omp parallel do collapse(2) private(ptype,iproc,q_species,m_species) schedule(static)
    do ptype = 1, self%ntype
      do iproc = 1, self%nproc
        q_species = self%chem%charge(ptype)
        m_species = self%chem%mass(ptype)

        if (m_species <= 0.0_real64) cycle
        if (.not. allocated(self%part(ptype,iproc)%x)) cycle
        if (self%part(ptype,iproc)%n <= 0_int32) cycle

        call move_particles_electrostatic( &
            part = self%part(ptype,iproc), &
            n    = int(self%dom%n, int32), &
            h    = self%dom%h, &
            E    = self%fld%E, &
            q    = q_species, &
            m    = m_species, &
            dt   = dt_local )
      end do
    end do
    !$omp end parallel do
  end subroutine move_particles_local


  subroutine apply_particle_bc_local(self)
    class(State), intent(inout) :: self
    integer(int32) :: ptype, iproc
    integer(int32) :: tag_neg_local
    integer(int32) :: ispec

    if (.not. allocated(self%part)) return

    tag_neg_local = -1_int32
    do ispec = 1_int32, self%ntype
      if (self%chem%charge(ispec) < 0.0_real64 .and. ispec /= 1_int32) then
        tag_neg_local = ispec
        exit
      end if
    end do

    !$omp parallel do collapse(2) private(ptype,iproc) schedule(static)
    do ptype = 1, self%ntype
      do iproc = 1, self%nproc
        if (.not. allocated(self%part(ptype,iproc)%x)) cycle
        if (self%part(ptype,iproc)%n <= 0_int32) cycle

        call apply_particle_bc_legacy( &
            part     = self%part(ptype,iproc), &
            n        = int(self%dom%n, int32), &
            h        = self%dom%h, &
            bcnd     = self%dom%bcnd, &
            xmax     = self%dom%xmax, &
            ymax     = self%dom%ymax, &
            zmax     = self%dom%zmax, &
            flag_pbc = int(self%dom%flag_pbc, int32), &
            flag_nmn = int(self%dom%flag_nmn, int32), &
            ptype    = ptype, &
            tag_neg  = tag_neg_local )
      end do
    end do
    !$omp end parallel do
  end subroutine apply_particle_bc_local


  subroutine finalize_particles_only(self)
    class(State), intent(inout) :: self
    integer(int32) :: ptype, iproc

    if (.not. allocated(self%part)) return

    do ptype = 1, size(self%part,1)
      do iproc = 1, size(self%part,2)
        call self%part(ptype,iproc)%destroy()
      end do
    end do

    deallocate(self%part)
  end subroutine finalize_particles_only


  subroutine finalize(self)
    class(State), intent(inout) :: self

    call self%pdec%destroy()
    call self%fld%destroy()
    call self%rxn%destroy()
    call self%magField%destroy()
    call self%finalize_particles_only()

    if (allocated(self%params%iseed)) deallocate(self%params%iseed)
  end subroutine finalize

end module mod_state