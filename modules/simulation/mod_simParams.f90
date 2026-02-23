module mod_simParams
  !!
  !! Derived “simulation parameters” module (legacy Step 7)
  !!
  !! Computes:
  !!  - Debye length lbd_d
  !!  - plasma frequency wp
  !!  - thermal velocities vt0(ptype)
  !!  - macroparticle weights Nm(ptype)
  !!  - dt, nb_step_sort, nb_step_collisions, nb_step_heating, nb_step_averaging, nseq, tseq corrected, nudt, nu_uplim
  !!  - plot plane indices ix_plot_plane, iz_plot_plane (and optional thruster ix_thr)
  !!  - plt_src
  !!  - random seeds iseed(iproc)
  !!
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use mod_config,         only: Config, npart
  use mod_domain,         only: Domain
  use mod_utils,          only: stop_calculation
  use mod_chemistryState, only: ChemistryState
  use mod_reactionsDB,    only: ReactionsDB
  use mod_magneticField, only: MagneticField
  use mod_constants,      only: qe, eps0
  implicit none
  private
  public :: SimParams

  type :: SimParams
    ! --- physical scales ---
    real(real64) :: lbd_d = 0.0_real64
    real(real64) :: wp    = 0.0_real64

    ! --- species derived arrays ---
    real(real64) :: vt0(npart) = 0.0_real64   ! thermal velocity
    real(real64) :: Nm (npart) = 0.0_real64   ! macro weight

    ! --- timestep / frequencies ---
    real(real64)   :: dt    = 0.0_real64
    integer(int32) :: nb_step_sort         = 10
    integer(int32) :: nb_step_collisions   = 10
    integer(int32) :: nb_step_heating      = 4
    integer(int32) :: nb_step_averaging    = 0
    integer(int32) :: nb_step_Neg_PE       = 1
    integer(int32) :: nb_step_Part_Injection = 1

    integer(int32) :: nseq    = 0
    real(real64)   :: tseq    = 0.0_real64

    real(real64) :: nudt = 0.0_real64
    real(real64) :: nu_uplim(npart) = 0.0_real64

    ! --- plotting / diagnostics ---
    integer(int32) :: ix_plot_plane = 0
    integer(int32) :: iz_plot_plane = 0
    integer(int32) :: plt_src = 1

    ! --- dielectric convenience (don’t mutate cfg) ---
    real(real64) :: Ca_cell = 0.0_real64

    ! --- thruster convenience (optional) ---
    real(real64)    :: x_thr  = 0.0_real64
    integer(int32)  :: ix_thr = 0
  contains
    procedure, public :: build         => build_sim_params
    procedure, public :: init_seeds
    procedure, public :: print_summary
  end type SimParams

contains

  subroutine build_sim_params(self, cfg, dom, chem, rxn, magF, mpi_rank)
    class(SimParams),    intent(inout) :: self
    type(Config),        intent(in)    :: cfg
    type(Domain),        intent(in)    :: dom
    type(ChemistryState),intent(in)    :: chem
    type(ReactionsDB),   intent(in)    :: rxn
    type(MagneticField), intent(in)    :: magF
    integer,             intent(in)    :: mpi_rank

    integer :: ptype
    integer :: ntype, n_neu
    real(real64) :: hmin
    real(real64) :: qabs, mabs
    logical :: bak_mode

    ntype = rxn%ntype
    n_neu = rxn%n_neu

    ! legacy: bak mode depends on cfg%nbak
    bak_mode = (cfg%nbak /= 0)

    ! -------------------------------
    ! Debye length / plasma frequency
    ! -------------------------------
    if (cfg%n0 <= 0.0_real64) then
      if (mpi_rank == 0) write(*,*) 'ERROR: n0 must be > 0 to compute lbd_d/wp'
      call stop_calculation
    end if

    if (abs(chem%charge(1)) <= 0.0_real64 .or. chem%mass(1) <= 0.0_real64) then
      if (mpi_rank == 0) write(*,*) 'ERROR: invalid chem%charge(1) or chem%mass(1)'
      call stop_calculation
    end if

    self%lbd_d = sqrt( (eps0 * cfg%Ti(1)) / (cfg%n0 * abs(chem%charge(1))) )
    self%wp    = sqrt( cfg%n0 * chem%charge(1)**2 / (eps0 * chem%mass(1)) )

    ! -------------------------------
    ! vt0(ptype) and Nm(ptype)
    ! -------------------------------
    self%vt0 = 0.0_real64
    self%Nm  = 0.0_real64

    if (cfg%np_cell == 0) then
      if (mpi_rank == 0) write(*,*) 'ERROR: np_cell must be nonzero for Nm()'
      call stop_calculation
    end if

    do ptype = 1, min(ntype + n_neu, npart)

      qabs = abs(chem%charge(ptype))
      mabs = abs(chem%mass(ptype))

      if (mabs > 0.0_real64) then
        ! Thermal velocities [sqrt(2*kT/m)]
        self%vt0(ptype) = sqrt( 2.0_real64 * qe * cfg%Ti(ptype) / mabs )
      else
        self%vt0(ptype) = 0.0_real64
      end if

      ! Macroparticle weight
      if (qabs > 0.0_real64) then
        self%Nm(ptype) = cfg%n0   * dom%h(1) * dom%h(2) * dom%h(3) / abs(real(cfg%np_cell, real64))
      else
        self%Nm(ptype) = cfg%ngas * dom%h(1) * dom%h(2) * dom%h(3) / abs(real(cfg%np_cell, real64))
      end if
    end do

    ! -------------------------------
    ! dt = kt*min(h)/vt0(1)
    ! -------------------------------
    hmin = min(dom%h(1), min(dom%h(2), dom%h(3)))
    if (self%vt0(1) <= 0.0_real64) then
      if (mpi_rank == 0) write(*,*) 'ERROR: vt0(1) <= 0, cannot compute dt'
      call stop_calculation
    end if
    self%dt = cfg%kt * hmin / self%vt0(1)

    ! -------------------------------
    ! call frequencies / averaging
    ! -------------------------------
    self%nb_step_sort       = 10
    self%nb_step_collisions = 1 * self%nb_step_sort

    self%nb_step_heating = 4
    if (cfg%kt >= 0.05_real64 .and. cfg%kt < 0.1_real64) self%nb_step_heating = 20
    if (cfg%kt <  0.05_real64)                          self%nb_step_heating = 40

    if (.not. bak_mode) then
      self%nb_step_averaging = 5 * self%nb_step_sort
    else
      self%nb_step_averaging = cfg%nsav
    end if

    ! sequence saving frequency (tseq corrected)
    self%tseq = cfg%tseq
    self%nseq = 0
    if (self%tseq > 0.0_real64) then
      self%nseq = max(nint((self%tseq / self%dt) / real(cfg%nsav, real64)), 1) * cfg%nsav
      self%tseq = real(self%nseq, real64) * self%dt
    end if

    ! -------------------------------
    ! nudt and nu_uplim
    ! -------------------------------
    self%nudt = cfg%nu_h * ( real(self%nb_step_heating, real64) * self%dt )

    self%nu_uplim = 0.0_real64
    self%nu_uplim(1) = 5.0e8_real64
    if (ntype >= 2) self%nu_uplim(2:ntype) = 1.0e7_real64

    self%nb_step_Neg_PE         = 1
    self%nb_step_Part_Injection = 1

    ! -------------------------------
    ! plot planes
    ! -------------------------------
    self%ix_plot_plane = int(dom%n(1)/2 + 1, int32)
    self%iz_plot_plane = int(dom%n(3)/2 + 1, int32)

    if (cfg%flag_grd == 1) then
      ! legacy: ix_pl = ixg (grid location)
      ! you currently store some grid index in cfg%ind_g; swap if you store ixg elsewhere
      self%ix_plot_plane = int(cfg%ind_g, int32)
    end if

    ! dielectric: Ca per cell (don’t change cfg%Ca)
    self%Ca_cell = 0.0_real64
    if (dom%flag_die == 1) self%Ca_cell = cfg%Ca * dom%h(1) * dom%h(3)

    ! thruster convenience
    self%x_thr  = 0.0_real64
    self%ix_thr = 0
    if (magF%flag_thr == 1) then
      self%x_thr  = dom%xmax - 1.0e-3_real64
      self%ix_thr = int(self%x_thr / dom%h(1), int32) + 1_int32
      self%ix_plot_plane = int( 0.5_real64*(cfg%xr_pow + cfg%xl_pow) / dom%h(1), int32 ) + 1_int32
    end if

    ! source plotting flag (legacy)
    self%plt_src = 1
    if ( (cfg%Pabs < 0.0_real64 .and. abs(cfg%opt_inj) < 3) .or. (cfg%flag_avg3D == 0) ) self%plt_src = 0

    ! warning like legacy
    if (mpi_rank == 0) then
      if (self%nudt > 1.0_real64) then
        write(*,*) 'Warning, nu*dt > 1, please correct ...'
        write(*,*) 'nu*dt = ', self%nudt
        call stop_calculation
      end if
    end if

  end subroutine build_sim_params


  subroutine init_seeds(self, iseed, nproc, mpi_rank)
    class(SimParams), intent(in)    :: self
    integer,          intent(inout) :: iseed(:)
    integer,          intent(in)    :: nproc, mpi_rank
    integer :: iproc

    if (size(iseed) < nproc) then
      if (mpi_rank == 0) write(*,*) 'ERROR: iseed(:) too small for nproc'
      call stop_calculation
    end if

    do iproc = 1, nproc
      iseed(iproc) = 123456 * iproc * (10*mpi_rank + 1)
    end do
  end subroutine init_seeds


  subroutine print_summary(self, mpi_rank, nsav)
    class(SimParams), intent(in) :: self
    integer,          intent(in) :: mpi_rank
    integer,          intent(in) :: nsav

    if (mpi_rank /= 0) return

    write(*,'(1x,"Frequency of calls to collision subroutine= ",i0)') self%nb_step_collisions
    write(*,'(1x,"Frequency of calls to electron (Maxwellian) heating subroutine= ",i0)') self%nb_step_heating
    write(*,'(1x,"Frequency of calls to sort subroutine= ",i0)') self%nb_step_sort

    if (self%nb_step_averaging == nsav) then
      write(*,'(1x,"Data will not be averaged!")')
    else
      write(*,'(1x,"Frequency for averaging= ",i0)') self%nb_step_averaging
    end if
  end subroutine print_summary

end module mod_simParams