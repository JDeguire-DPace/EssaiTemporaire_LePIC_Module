module mod_part_info
  ! Exact replacement for Src/particle_info.h (legacy-global style)
  implicit none
  private

  ! ---- Public API: export everything (legacy style) ----
  public :: npart, np_cell, flag_pdf, flag_restart, n_cell, nbak, n_rg, &
            flag_pbc, cnt_seg, tag_neg, num_grd, tag_neu, flag_nmn, flag_heat, flag_inj, &
            opt_inj, flag_B_pos, omp_rank_max, mpi_rank_max, flag_die, flag_thr, ix_pl, iz_pl, ix_thr, flag_convP, &
            ig_die, flag_pbcz, tag_beam, every, flag_bak, flag_circxh, flag_ahp, flag_gridB, flag_RFpot, &
            ncol, p_ncol, sig_npt_mx, nscol, p_nscol, np_mx, nu_max, nu_uplim, &
            np_loss, P_w, ind_nre, ind_nby, ind_Eth, ind_dE, cnt_neg, charge, mass, vt0, Ti, &
            nB, B_scale, B0, dL, x0, y0, z0, B_file, B_name, B_info, &
            ixl_pow, ixr_pow, ns_heat, ns_coll, nm_rg, n_holes, flag_grd, ind_g, ptype_pdf, ns_flx, &
            flag_B, nhy, nhz, ns_inj, ixg, ix_PE, igrid_sec, dir_sec, n_cath, plt_src, &
            ixl_rg, ixr_rg, iyl_rg, iyr_rg, izl_rg, izr_rg, cnt_plt, &
            Pabs, n0, ngas, Nm, dt, xmax, ymax, zmax, x_load, nu_h, Lgy, Lgz, Sg, Lhy, Lhz, xg1, Bmax, nudt, R_ahp, &
            k_eps0, jne, RNeta, THm, x_thr, y_thr, tseq, tseq_init, tseq_final, gam_sec, zg_sec, pname, &
            ! --- Added for modular read_input/generate_boundary pipeline ---
            n, ngrid, ng, n_B, mpi_rank, flag_avg3D, np_dup, Ca, rname, &
            xl_pow, xr_pow, yl_pow, yr_pow, zl_pow, zr_pow, &
            xl_rg, xr_rg, yl_rg, yr_rg, zl_rg, zr_rg, &
            phi0_RF, f0_RF, phi1_RF, f1_RF, tmax, eps, omega, kt, I_inj, nsav
  


  ! -------------------------
  ! Counters and seed number
  ! -------------------------
  integer, parameter :: npart = 10

  integer :: np_cell, flag_pdf, &
             flag_restart, n_cell, nbak, n_rg, &
             flag_pbc, cnt_seg, tag_neg, num_grd, &
             tag_neu, flag_nmn, flag_heat, flag_inj, &
             opt_inj, flag_B_pos, omp_rank_max, mpi_rank_max, &
             flag_die, flag_thr, ix_pl, iz_pl, ix_thr, flag_convP, &
             ig_die(5), flag_pbcz, tag_beam, every, flag_bak, &
             flag_circxh, flag_ahp, flag_gridB, flag_RFpot

  ! ----------- Added: MPI / grid sizing / MG -----------
  integer :: mpi_rank
  integer :: n(3)        ! legacy exponent/size array (as in read_input signature)
  integer :: ngrid       ! number of wall labels
  integer :: ng          ! number of MG levels (as in read_input)
  integer :: flag_avg3D  ! saved from conditions.inp (-ngrid trick)
  integer :: n_B(3)      ! magnetic-field map grid
  real(kind=8) :: np_dup ! duplication factor for restart logic

  ! ----------- Added: file-driven config values -----------
  character(len=20) :: rname
  real(kind=8) :: Ca

  ! power deposition region (meters in legacy after conversion)
  real(kind=8) :: xl_pow, xr_pow, yl_pow, yr_pow, zl_pow, zr_pow

  ! PDF regions (meters after conversion)
  integer, parameter :: nm_rg = 5
  real(kind=8) :: xl_rg(nm_rg), xr_rg(nm_rg)
  real(kind=8) :: yl_rg, yr_rg, zl_rg, zr_rg

  ! RF params (as per read_input signature)
  real(kind=8) :: phi0_RF, f0_RF, phi1_RF, f1_RF

  ! ----------- Existing: Collisions -----------
  integer :: ncol, p_ncol(npart), sig_npt_mx, nscol, p_nscol(npart)
  real(kind=8) :: np_mx(npart), nu_max(npart), nu_uplim(npart)

  ! ------------------------
  ! Particle arrays & flags
  ! ------------------------
  integer, parameter :: np_loss=1, P_w=2, ind_nre=1, ind_nby=2, ind_Eth=1, ind_dE=2
  integer :: cnt_neg
  real(kind=8) :: charge(npart), mass(npart), vt0(npart), Ti(npart)

  ! ---------------
  ! Magnetic field
  ! ---------------
  integer :: nB
  real(kind=8) :: B_scale(5), B0(5), dL(5), x0(5), y0(5), z0(5)
  character(len=1)  :: B_file(5)
  character(len=20) :: B_name(5)
  character(len=2)  :: B_info(5)

  ! ----------------------------
  ! Particle and simulation info
  ! ----------------------------
  integer :: ixl_pow, ixr_pow, ns_heat, ns_coll, n_holes, flag_grd, ind_g, &
             ptype_pdf, ns_flx, flag_B, nhy, nhz, ns_inj, ixg, ix_PE, igrid_sec, dir_sec, &
             n_cath, plt_src

  integer :: ixl_rg(nm_rg), ixr_rg(nm_rg), iyl_rg, iyr_rg, izl_rg, izr_rg, cnt_plt

  real(kind=8) :: Pabs, n0, ngas, Nm(npart), dt, xmax, ymax, zmax, x_load, &
                  nu_h, Lgy, Lgz, Sg, Lhy, Lhz, xg1, Bmax, nudt, R_ahp

  real(kind=8) :: k_eps0, jne, RNeta, THm, x_thr, y_thr, tseq, tseq_init, &
                  tseq_final, gam_sec, zg_sec(2)

  character(len=6) :: pname(npart)
  integer :: nsav
  real(kind=8) :: tmax, eps, omega, kt, I_inj
end module mod_part_info
