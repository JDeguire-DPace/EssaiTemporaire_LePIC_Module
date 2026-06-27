module mod_state
  use iso_fortran_env, only: int32, real64
  use omp_lib, only: omp_get_max_threads, omp_get_num_procs
  use mpi

  use mod_config,         only: Config
  use mod_readConditions, only: read_input

  use mod_constants,      only: eps0, qe
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
  use mod_heating,          only: apply_electron_heating
  use mod_planeMoments,     only: compute_particle_plane_moments_species

  implicit none
  private
  public :: State

  ! Production defaults. Keep expensive consistency checks available in the
  ! source, but do not run them in normal legacy-comparison/production runs.
  logical, parameter :: DEBUG_INIT_CHECKS = .false.
  logical, parameter :: VALIDATE_SORTING  = .false.

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
    real(real64), allocatable :: sum_q_xz(:,:,:,:)
    real(real64), allocatable :: sum_q_yz(:,:,:,:,:)

    integer :: mpi_rank = -1
    integer :: mpi_size = -1
    integer :: comm     = MPI_COMM_NULL

    integer(int32) :: ntype = 0
    integer(int32) :: nproc = 1
    integer(int32), allocatable :: Nh(:)
    real(real64), allocatable :: sum_dEk(:)

    real(real64) :: p_loss_heating_local
    real(real64), allocatable :: P_loss(:,:,:)
    real(real64), allocatable :: p_mac(:,:,:,:)

    real(real64), allocatable :: data_pavg_xy(:,:,:,:)
    real(real64), allocatable :: data_pavg_xz(:,:,:,:)
    real(real64), allocatable :: data_pavg_yz(:,:,:,:)

    real(real64), allocatable :: np_avg_xy(:,:,:)
    real(real64), allocatable :: np_avg_xz(:,:,:)
    real(real64), allocatable :: np_avg_yz(:,:,:)

    real(real64), allocatable :: phi_avg_xy(:,:)
    real(real64), allocatable :: phi_avg_xz(:,:)
    real(real64), allocatable :: phi_avg_yz(:,:)

    integer(int32) :: cnt_avg = 0_int32

  contains
    procedure :: init
    procedure :: init_chemistry
    procedure :: init_domain_and_boundary
    procedure :: init_particles

    procedure :: sort_particles_local
    procedure :: move_particles_local
    procedure :: apply_particle_bc_local
    procedure :: apply_electron_heating_local
    procedure :: compute_heating_region_moments
    procedure :: update_heating_vt
    procedure :: compute_plane_moments_local
    procedure :: apply_dielectric_bc_to_phi
    procedure :: accumulate_2d_averages

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

    self%ntype = int(self%rxn%ntype - self%rxn%n_neu, int32)
    if (self%ntype < 1_int32) self%ntype = 1_int32

    self%nproc = max(1_int32, int(self%cfg%omp_rank_max, int32))
    call self%params%init_seeds(self%nproc, self%mpi_rank)
    call self%params%print_summary(self%mpi_rank, self%cfg%nsav)

    allocate(self%np_thread(0:self%dom%n(1)+2, &
                            0:self%dom%n(2)+2, &
                            0:self%dom%n(3)+2, &
                            self%ntype, self%nproc))
    self%np_thread = 0.0_real64

    call self%init_particles()

    allocate(self%sum_q_xz(0:self%dom%n(1)+2, &
                           0:self%dom%n(3)+2, &
                           self%ntype, self%nproc))
    self%sum_q_xz = 0.0_real64

    allocate(self%sum_q_yz(2, &
                           0:self%dom%n(2)+2, &
                           0:self%dom%n(3)+2, &
                           self%ntype, self%nproc))
    self%sum_q_yz = 0.0_real64

    allocate(self%P_loss(4, self%ntype, self%nproc))
    self%P_loss = 0.0_real64

    allocate(self%p_mac(self%ntype, 2, 0:self%cfg%ngrid, self%nproc))
    self%p_mac = 0.0_real64

    allocate(self%Nh(self%nproc))
    allocate(self%sum_dEk(self%nproc))
    self%Nh      = 0_int32
    self%sum_dEk = 0.0_real64

    allocate(self%np_avg_xy(0:self%dom%n(1)+2,0:self%dom%n(2)+2,self%ntype))
    allocate(self%np_avg_xz(0:self%dom%n(1)+2,0:self%dom%n(3)+2,self%ntype))
    allocate(self%np_avg_yz(0:self%dom%n(2)+2,0:self%dom%n(3)+2,self%ntype))

    allocate(self%phi_avg_xy(0:self%dom%n(1)+2,0:self%dom%n(2)+2))
    allocate(self%phi_avg_xz(0:self%dom%n(1)+2,0:self%dom%n(3)+2))
    allocate(self%phi_avg_yz(0:self%dom%n(2)+2,0:self%dom%n(3)+2))

    self%np_avg_xy  = 0.0_real64
    self%np_avg_xz  = 0.0_real64
    self%np_avg_yz  = 0.0_real64
    self%phi_avg_xy = 0.0_real64
    self%phi_avg_xz = 0.0_real64
    self%phi_avg_yz = 0.0_real64
    self%cnt_avg    = 0_int32

    allocate(self%data_pavg_xy(5, 0:self%dom%n(1)+2, 0:self%dom%n(2)+2, self%ntype))
    allocate(self%data_pavg_xz(5, 0:self%dom%n(1)+2, 0:self%dom%n(3)+2, self%ntype))
    allocate(self%data_pavg_yz(5, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2, self%ntype))

    self%data_pavg_xy = 0.0_real64
    self%data_pavg_xz = 0.0_real64
    self%data_pavg_yz = 0.0_real64

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

    if (DEBUG_INIT_CHECKS) then
      call checkpoint_poisson_decomp(self%mpi_rank, self%comm, &
                                     int(self%dom%n(3),int32), self%pdec%k0, self%pdec%m, &
                                     self%pdec%phi_dom, self%pdec%bcnd_dom)
    end if

    call build_kq(self%dom%bcnd, self%fld%kq)

    if (DEBUG_INIT_CHECKS) then
      call checkpoint_kq(self%fld%kq, self%comm, self%mpi_rank, 'after build_kq', self%pdec)

      allocate(phi_rt(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2))
      phi_rt = -999.0_real64

      call self%pdec%gather_phi_to_global(phi_rt)

      lmax = maxval(abs(phi_rt - self%fld%phi))
      call MPI_Allreduce(lmax, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, self%comm, ierr)

      deallocate(phi_rt)
    end if

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

    ntype_trk = self%ntype

    if (allocated(self%part)) then
      call self%finalize_particles_only()
    end if
    allocate(self%part(self%ntype, self%nproc))

    call self%fld%allocate_species_density(self%dom, self%ntype)

    self%np_thread = 0.0_real64

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
          np_thread = self%np_thread )
    call reduce_species_density( &
         n         = int(self%dom%n, int32), &
         bcnd      = self%dom%bcnd, &
         np_thread = self%np_thread, &
         ntype     = int(self%ntype), &
         nproc     = int(self%nproc), &
         mpi_comm  = self%comm, &
         np_red    = self%fld%np )

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

        if (VALIDATE_SORTING) then
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
    if (self%cfg%flag_circxh == 1) self%cfg%R_ahp = self%cfg%yr_pow - self%dom%ymax/2.0_real64

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
          xl_pow         = self%cfg%xl_pow, &
          xr_pow         = self%cfg%xr_pow, &
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


  subroutine apply_dielectric_bc_to_phi(self)
    use iso_fortran_env, only: int32, real64
    class(State), intent(inout) :: self

    integer(int32) :: ix, iy, iz, ptype, iproc
    integer(int32) :: igrid
    real(real64)   :: Ca
    real(real64)   :: qsum_xz, qsum_yz_1, qsum_yz_2

    if (self%dom%flag_die /= 1) return
    if (.not. allocated(self%sum_q_xz)) return
    if (.not. allocated(self%sum_q_yz)) return

    Ca = self%params%Ca_cell

    do iz = 0, self%dom%n(3)+2
      do iy = 0, self%dom%n(2)+2
        do ix = 0, self%dom%n(1)+2

          igrid = self%dom%bcnd(ix,iy,iz)
          if (igrid <= 0) cycle

          do ptype = 1, self%ntype

            if (self%dom%dtype(igrid) == 2) then
              qsum_xz = 0.0_real64
              do iproc = 1, self%nproc
                qsum_xz = qsum_xz + self%sum_q_xz(ix,iz,ptype,iproc)
              end do
              self%fld%phi(ix,iy,iz) = self%fld%phi(ix,iy,iz) + qsum_xz / (2.0_real64 * Ca)
            end if

            if (self%dom%dtype(igrid) == 3) then
              qsum_yz_1 = 0.0_real64
              do iproc = 1, self%nproc
                qsum_yz_1 = qsum_yz_1 + self%sum_q_yz(1,iy,iz,ptype,iproc)
              end do
              self%fld%phi(ix,iy,iz) = self%fld%phi(ix,iy,iz) + qsum_yz_1 / Ca
            end if

            if (self%dom%dtype(igrid) == 4) then
              qsum_yz_2 = 0.0_real64
              do iproc = 1, self%nproc
                qsum_yz_2 = qsum_yz_2 + self%sum_q_yz(2,iy,iz,ptype,iproc)
              end do
              self%fld%phi(ix,iy,iz) = self%fld%phi(ix,iy,iz) + qsum_yz_2 / Ca
            end if

          end do
        end do
      end do
    end do
  end subroutine apply_dielectric_bc_to_phi


  subroutine move_particles_local(self)
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
    use iso_fortran_env, only: int32, real64
    class(State), intent(inout) :: self

    integer(int32) :: ptype, iproc
    integer(int32) :: tag_neg_local
    integer(int32) :: ispec
    real(real64)   :: qmacro

    if (.not. allocated(self%part)) return

    tag_neg_local = -1_int32
    do ispec = 1_int32, self%ntype
      if (self%chem%charge(ispec) < 0.0_real64 .and. ispec /= 1_int32) then
        tag_neg_local = ispec
        exit
      end if
    end do

    !$omp parallel do collapse(2) private(ptype,iproc,qmacro) schedule(static)
    do ptype = 1, self%ntype
      do iproc = 1, self%nproc

        if (.not. allocated(self%part(ptype,iproc)%x)) cycle
        if (self%part(ptype,iproc)%n <= 0_int32) cycle

        ! The caller zeros sum_q_xz/sum_q_yz before applying particle BCs.
        ! Each OpenMP iteration writes to a unique (ptype,iproc) slice, so no
        ! temporary allocation or post-copy is needed here.
        qmacro = self%params%Nm(ptype) * self%chem%charge(ptype)

        call apply_particle_bc_legacy( &
            part           = self%part(ptype,iproc), &
            n              = int(self%dom%n, int32), &
            h              = self%dom%h, &
            bcnd           = self%dom%bcnd, &
            xmax           = self%dom%xmax, &
            ymax           = self%dom%ymax, &
            zmax           = self%dom%zmax, &
            flag_pbc       = int(self%dom%flag_pbc, int32), &
            flag_nmn       = int(self%dom%flag_nmn, int32), &
            ptype          = ptype, &
            tag_neg        = tag_neg_local, &
            flag_die       = int(self%dom%flag_die, int32), &
            dtype          = self%dom%dtype, &
            qmacro         = qmacro, &
            sum_q_xz_local = self%sum_q_xz(:,:,ptype,iproc), &
            sum_q_yz_local = self%sum_q_yz(:,:,:,ptype,iproc), &
            p_mac_boundary = self%p_mac(ptype,:,:,iproc), &
            mass_species   = self%chem%mass(ptype), &
            P_loss_wall    = self%P_loss(1,ptype,iproc), &
            Nm_species     = self%params%Nm(ptype))
      end do
    end do
    !$omp end parallel do
  end subroutine apply_particle_bc_local


  subroutine compute_heating_region_moments(self)
    class(State), intent(inout) :: self
    integer(int32) :: iproc, i, ix
    real(real64)   :: x, y, z, v2
    real(real64)   :: ymax_half, zmax_half

    if (.not. allocated(self%part)) return
    if (.not. allocated(self%Nh)) return
    if (.not. allocated(self%sum_dEk)) return

    self%Nh      = 0_int32
    self%sum_dEk = 0.0_real64

    ymax_half = self%dom%ymax / 2.0_real64
    zmax_half = self%dom%zmax / 2.0_real64

    !$omp parallel do private(iproc,i,ix,x,y,z,v2) schedule(static)
    do iproc = 1, self%nproc
      if (.not. allocated(self%part(1,iproc)%x)) cycle
      if (self%part(1,iproc)%n <= 0_int32) cycle

      do i = 1, self%part(1,iproc)%n
        x  = self%part(1,iproc)%x(i)
        ix = int(x / self%dom%h(1), int32) + 1_int32

        if (ix < self%cfg%xl_pow/self%dom%h(1) .or. ix > self%cfg%xr_pow/self%dom%h(1)) cycle

        if (self%cfg%flag_circxh == 1) then
          y = self%part(1,iproc)%y(i)
          z = self%part(1,iproc)%z(i)

          if (self%cfg%flag_ahp == 0) then
            if (((y-ymax_half)**2 + (z-zmax_half)**2) > self%cfg%R_ahp**2) cycle
          else
            if (((y-ymax_half)**2 + (z-zmax_half)**2) < self%cfg%R_ahp**2) cycle
          end if
        end if

        v2 = self%part(1,iproc)%vx(i)**2 + &
             self%part(1,iproc)%vy(i)**2 + &
             self%part(1,iproc)%vz(i)**2

        self%sum_dEk(iproc) = self%sum_dEk(iproc) + &
            0.5_real64 * self%params%Nm(1) * self%chem%mass(1) * v2

        self%Nh(iproc) = self%Nh(iproc) + 1_int32
      end do
    end do
    !$omp end parallel do
  end subroutine compute_heating_region_moments


  subroutine update_heating_vt(self, vt_heat)
    use mpi
    class(State), intent(inout) :: self
    real(real64), intent(out)   :: vt_heat

    integer :: ierr
    integer(int32) :: sum_Nh, sum_Nh_global
    real(real64)   :: sum_dEk_tot, sum_dEk_global
    real(real64)   :: Te

    vt_heat = self%params%vt0(1)

    if (.not. allocated(self%Nh)) return
    if (.not. allocated(self%sum_dEk)) return

    sum_Nh       = sum(self%Nh)
    sum_dEk_tot  = sum(self%sum_dEk)

    sum_Nh_global  = sum_Nh
    sum_dEk_global = sum_dEk_tot

    if (self%mpi_size > 1) then
      call MPI_Allreduce(sum_Nh, sum_Nh_global, 1, MPI_INTEGER, MPI_SUM, self%comm, ierr)
      call MPI_Allreduce(sum_dEk_tot, sum_dEk_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
    end if

    if (sum_Nh_global <= 0) return

    Te = (2.0_real64/3.0_real64) * &
         (sum_dEk_global + self%cfg%Pabs / self%cfg%nu_h) / &
         (self%params%Nm(1) * qe * real(sum_Nh_global, real64))

    if (Te > 0.0_real64) then
      vt_heat = sqrt(2.0_real64 * qe * Te / abs(self%chem%mass(1)))
    end if
    
  end subroutine update_heating_vt


  subroutine compute_plane_moments_local(self)
    class(State), intent(inout) :: self
    integer(int32) :: ptype

    if (.not. allocated(self%part)) return

    do ptype = 1, self%ntype
      call compute_particle_plane_moments_species( &
          part         = self%part(ptype,:), &
          nproc        = self%nproc, &
          n            = int(self%dom%n, int32), &
          h            = self%dom%h, &
          ptype        = ptype, &
          mass_species = self%chem%mass(ptype), &
          ix_plane     = self%params%ix_plot_plane , &
          iy_plane     = int(self%dom%n(2)/2 +1, int32), &
          iz_plane     = self%params%iz_plot_plane, &
          data_xy      = self%data_pavg_xy(:,:,:,ptype), &
          data_xz      = self%data_pavg_xz(:,:,:,ptype), &
          data_yz      = self%data_pavg_yz(:,:,:,ptype), &
          params       = self%params)
    end do
  end subroutine compute_plane_moments_local

  subroutine accumulate_2d_averages(self)
    class(State), intent(inout) :: self

    integer(int32) :: ptype
    integer(int32) :: ix_plane, iy_density, iy_phi, iz_plane
    integer(int32) :: ny

    ny = self%dom%n(2)

    ix_plane   = self%params%ix_plot_plane
    iy_density = int(ny/2 + 1, int32)
    iy_phi     = int(ny/2,     int32)
    iz_plane   = self%params%iz_plot_plane

    do ptype = 1, self%ntype
      self%np_avg_xy(:,:,ptype) = self%np_avg_xy(:,:,ptype) + &
          self%fld%np(:,:,iz_plane,ptype)

      self%np_avg_xz(:,:,ptype) = self%np_avg_xz(:,:,ptype) + &
          self%fld%np(:,iy_density,:,ptype)

      self%np_avg_yz(:,:,ptype) = self%np_avg_yz(:,:,ptype) + &
          self%fld%np(ix_plane,:,:,ptype)
    end do

    self%phi_avg_xy(:,:) = self%phi_avg_xy(:,:) + self%fld%phi(:,:,iz_plane)
    self%phi_avg_xz(:,:) = self%phi_avg_xz(:,:) + self%fld%phi(:,iy_phi,:)
    self%phi_avg_yz(:,:) = self%phi_avg_yz(:,:) + self%fld%phi(ix_plane,:,:)

    self%cnt_avg = self%cnt_avg + 1_int32

  end subroutine accumulate_2d_averages


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

    if (allocated(self%np_thread))    deallocate(self%np_thread)
    if (allocated(self%sum_q_xz))     deallocate(self%sum_q_xz)
    if (allocated(self%sum_q_yz))     deallocate(self%sum_q_yz)
    if (allocated(self%P_loss))       deallocate(self%P_loss)
    if (allocated(self%p_mac))        deallocate(self%p_mac)
    if (allocated(self%Nh))           deallocate(self%Nh)
    if (allocated(self%sum_dEk))      deallocate(self%sum_dEk)
    if (allocated(self%np_avg_xy))    deallocate(self%np_avg_xy)
    if (allocated(self%np_avg_xz))    deallocate(self%np_avg_xz)
    if (allocated(self%np_avg_yz))    deallocate(self%np_avg_yz)
    if (allocated(self%phi_avg_xy))   deallocate(self%phi_avg_xy)
    if (allocated(self%phi_avg_xz))   deallocate(self%phi_avg_xz)
    if (allocated(self%phi_avg_yz))   deallocate(self%phi_avg_yz)
    if (allocated(self%data_pavg_xy)) deallocate(self%data_pavg_xy)
    if (allocated(self%data_pavg_xz)) deallocate(self%data_pavg_xz)
    if (allocated(self%data_pavg_yz)) deallocate(self%data_pavg_yz)

    if (allocated(self%params%iseed)) deallocate(self%params%iseed)
  end subroutine finalize

end module mod_state