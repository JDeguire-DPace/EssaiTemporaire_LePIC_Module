module mod_config
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: Config, NBMAX, nm_rg, npart

  integer, parameter :: NBMAX = 5
  integer, parameter :: nm_rg = 5
  integer, parameter :: npart = 32

  type :: Config
    integer           :: n(3)      = 0
    real(real64)      :: tmax      = 0.0_real64
    real(real64)      :: xl_pow    = 0.0_real64, xr_pow = 0.0_real64
    real(real64)      :: yl_pow    = 0.0_real64, yr_pow = 0.0_real64
    integer           :: flag_circxh = 0
    real(real64)      :: R_ahp    = 0.0_real64
    real(real64)      :: zl_pow    = 0.0_real64, zr_pow = 0.0_real64
    integer           :: nsav      = 0
    real(real64)      :: eps       = 0.0_real64
    real(real64)      :: omega     = 0.0_real64
    real(real64)      :: kt        = 0.0_real64
    character(len=20) :: rname     = ""
    integer           :: ngrid     = 0
    integer           :: ng        = 0

    integer           :: n_B(3)    = 1
    integer           :: flag_avg3D = 1
    real(real64)      :: np_dup    = 1.0_real64

    real(real64)      :: xl_rg(nm_rg) = 0.0_real64
    real(real64)      :: xr_rg(nm_rg) = 0.0_real64
    real(real64)      :: yl_rg     = 0.0_real64, yr_rg = 0.0_real64
    real(real64)      :: zl_rg     = 0.0_real64, zr_rg = 0.0_real64

    real(real64)      :: I_inj     = 0.0_real64
    real(real64)      :: Ca        = 0.0_real64

    real(real64)      :: phi0_RF   = 0.0_real64, f0_RF = 0.0_real64
    real(real64)      :: phi1_RF   = 0.0_real64, f1_RF = 0.0_real64

    real(real64)      :: Ti(npart)     = 0.0_real64
    real(real64)      :: x_load    = 0.0_real64
    integer           :: flag_ahp  = 0

    integer           :: flag_B_pos = 0
    integer           :: flag_B     = 0
    integer           :: nB         = 0

    character(len=1)  :: B_file(NBMAX)  = " "
    character(len=20) :: B_name(NBMAX)  = ""
    character(len=2)  :: B_info(NBMAX)  = ""
    real(real64)      :: B_scale(NBMAX) = 0.0_real64
    real(real64)      :: B0(NBMAX)      = 0.0_real64
    real(real64)      :: dL(NBMAX)      = 0.0_real64
    real(real64)      :: x0(NBMAX)      = 0.0_real64
    real(real64)      :: y0(NBMAX)      = 0.0_real64
    real(real64)      :: z0(NBMAX)      = 0.0_real64

    real(real64)      :: Pabs      = 0.0_real64
    real(real64)      :: nu_h      = 0.0_real64
    integer           :: opt_inj   = 0
    integer           :: flag_heat = 0
    integer           :: flag_inj  = 0

    real(real64)      :: n0        = 0.0_real64
    real(real64)      :: ngas      = 0.0_real64
    integer           :: np_cell   = 0

    integer           :: nbak      = 0
    integer           :: cnt_plt   = 0
    real(real64)      :: tseq_init = 0.0_real64
    real(real64)      :: tseq_final= 0.0_real64
    real(real64)      :: tseq      = 0.0_real64

    real(real64)      :: k_eps0    = 0.0_real64
    integer           :: flag_convP= 0

    integer           :: flag_pdf  = 0
    integer           :: ptype_pdf = 0
    integer           :: n_rg      = 0

    integer           :: omp_rank_max = 1
    integer           :: mpi_rank_max = 1
    integer           :: flag_restart = 0

    real(real64)      :: Lhy = 0.0_real64
    real(real64)      :: Lhz = 0.0_real64
    integer           :: nhy = 0
    integer           :: nhz = 0
    integer           :: ind_g = 0

    real(real64)      :: jne = 0.0_real64
    real(real64)      :: THm = 0.0_real64
    integer           :: num_grd = 0

    real(real64)      :: gam_sec  = 0.0_real64
    integer           :: igrid_sec = 0

    integer           :: flag_RFpot = 0
    integer           :: flag_grd   = 0
  end type Config

end module mod_config
