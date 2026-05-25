module mod_particleBC
  use iso_fortran_env, only: int32, int64, int8, real64
  use mod_particles,   only: ParticleSet
  use mod_constants,   only: qe

  implicit none
  private

  public :: apply_particle_bc_legacy
  public :: particle_is_lost_legacy

  public :: dbg_loss_xright_s1, dbg_loss_xright_s2
  public :: dbg_loss_zlow_s1,   dbg_loss_zlow_s2
  public :: dbg_loss_zhigh_s1,  dbg_loss_zhigh_s2

  integer(int64), save :: dbg_loss_xright_s1 = 0_int64
  integer(int64), save :: dbg_loss_xright_s2 = 0_int64
  integer(int64), save :: dbg_loss_zlow_s1   = 0_int64
  integer(int64), save :: dbg_loss_zlow_s2   = 0_int64
  integer(int64), save :: dbg_loss_zhigh_s1  = 0_int64
  integer(int64), save :: dbg_loss_zhigh_s2  = 0_int64

contains

  logical function particle_is_lost_legacy(bcnd, ix, iy, iz, n) result(is_lost)
    integer(int32), intent(in) :: ix, iy, iz
    integer(int32), intent(in) :: n(3)
    integer(int32), intent(in) :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    is_lost = &
         bcnd(ix  ,iy  ,iz  ) >= 1_int32 .and. &
         bcnd(ix+1,iy  ,iz  ) >= 1_int32 .and. &
         bcnd(ix+1,iy+1,iz  ) >= 1_int32 .and. &
         bcnd(ix  ,iy+1,iz  ) >= 1_int32 .and. &
         bcnd(ix  ,iy  ,iz+1) >= 1_int32 .and. &
         bcnd(ix+1,iy  ,iz+1) >= 1_int32 .and. &
         bcnd(ix+1,iy+1,iz+1) >= 1_int32 .and. &
         bcnd(ix  ,iy+1,iz+1) >= 1_int32
  end function particle_is_lost_legacy


  subroutine apply_particle_bc_legacy( part, n, h, bcnd, xmax, ymax, zmax, &
                                       flag_pbc, flag_nmn, ptype, tag_neg, &
                                       flag_die, dtype, qmacro, &
                                       sum_q_xz_local, sum_q_yz_local, &
                                       p_mac_boundary, mass_species, P_loss_wall, Nm_species )

    class(ParticleSet), intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    integer(int32),     intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(in)    :: xmax, ymax, zmax
    integer(int32),     intent(in)    :: flag_pbc, flag_nmn, ptype, tag_neg
    integer(int32),     intent(in)    :: flag_die
    integer(int32),     intent(in)    :: dtype(:)
    real(real64),       intent(in)    :: qmacro
    real(real64),       intent(inout) :: sum_q_xz_local(0:n(1)+2,0:n(3)+2)
    real(real64),       intent(inout) :: sum_q_yz_local(2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(inout) :: p_mac_boundary(:,:)
    real(real64),       intent(in)    :: mass_species
    real(real64),       intent(inout) :: P_loss_wall
    real(real64),       intent(in)    :: Nm_species

    integer(int32) :: i, i_shift
    integer(int32) :: ix, iy, iz
    integer(int32) :: flag_lost, np_lost
    integer(int32) :: igrid, d_ind

    real(real64) :: xp_new, yp_new, zp_new
    real(real64) :: vpx_new, vpy_new, vpz_new
    real(real64) :: px, py, pz
    real(real64) :: ki4(4)
    real(real64) :: Ek_eV, Ek_J

    if (.not. allocated(part%x)) return
    if (part%n <= 0_int32) return

    np_lost = 0_int32

    do i = 1, part%n

      xp_new  = part%x(i)
      yp_new  = part%y(i)
      zp_new  = part%z(i)
      vpx_new = part%vx(i)
      vpy_new = part%vy(i)
      vpz_new = part%vz(i)

      if (allocated(part%flag_dead)) then
        if (part%flag_dead(i) == 1_int8) then
          np_lost = np_lost + 1_int32
          cycle
        end if
      end if

      ! Legacy: compute post-move cell BEFORE periodic wrapping.
      ix = floor(xp_new / h(1)) + 1_int32
      iy = floor(yp_new / h(2)) + 1_int32
      iz = floor(zp_new / h(3)) + 1_int32

      if (ix < 0_int32) ix = 0_int32
      if (ix > n(1)+1_int32) ix = n(1)+1_int32
      if (iy < 0_int32) iy = 0_int32
      if (iy > n(2)+1_int32) iy = n(2)+1_int32
      if (iz < 0_int32) iz = 0_int32
      if (iz > n(3)+1_int32) iz = n(3)+1_int32

      flag_lost = 0_int32

      ! Legacy loss test: all 8 surrounding nodes are solid/wall.
      if (particle_is_lost_legacy(bcnd, ix, iy, iz, n)) flag_lost = 1_int32

      ! Negative ions with Neumann BC: no specular reflection.
      if (ptype == tag_neg) then
        if (xp_new < 0.0_real64 .and. flag_nmn == 1_int32) flag_lost = 2_int32
      end if

      if (flag_lost >= 1_int32) then

        if (ptype == 1_int32) then
          if (xp_new > xmax - h(1)) dbg_loss_xright_s1 = dbg_loss_xright_s1 + 1_int64
          if (zp_new < h(3))        dbg_loss_zlow_s1   = dbg_loss_zlow_s1   + 1_int64
          if (zp_new > zmax-h(3))   dbg_loss_zhigh_s1  = dbg_loss_zhigh_s1  + 1_int64
        else if (ptype == 2_int32) then
          if (xp_new > xmax - h(1)) dbg_loss_xright_s2 = dbg_loss_xright_s2 + 1_int64
          if (zp_new < h(3))        dbg_loss_zlow_s2   = dbg_loss_zlow_s2   + 1_int64
          if (zp_new > zmax-h(3))   dbg_loss_zhigh_s2  = dbg_loss_zhigh_s2  + 1_int64
        end if

        igrid = bcnd(ix,iy,iz)

        if (flag_lost == 2_int32 .and. ptype == tag_neg) igrid = 0_int32

        if (flag_die == 1_int32 .and. igrid > 0_int32) then
          if (dtype(igrid) > 1_int32) then

            ! Legacy only corrects z for dielectric + periodic case here.
            if (flag_pbc == 1_int32) then
              if (zp_new >= zmax) then
                zp_new = zp_new - zmax
                iz = floor(zp_new / h(3), int32) + 1_int32
              end if
              if (zp_new <= 0.0_real64) then
                zp_new = zmax + zp_new
                iz = floor(zp_new / h(3), int32) + 1_int32
              end if
            end if

            if (ix < 0_int32)      ix = 0_int32
            if (ix > n(1)+1_int32) ix = n(1)+1_int32
            if (iy < 0_int32)      iy = 0_int32
            if (iy > n(2)+1_int32) iy = n(2)+1_int32
            if (iz < 0_int32)      iz = 0_int32
            if (iz > n(3)+1_int32) iz = n(3)+1_int32

            pz = (real(iz,real64)*h(3) - zp_new) / h(3)

            if (dtype(igrid) == 2_int32) then
              px = (real(ix,real64)*h(1) - xp_new) / h(1)

              ki4(1) = qmacro * px               * pz
              ki4(2) = qmacro * (1.0_real64-px) * pz
              ki4(3) = qmacro * (1.0_real64-px) * (1.0_real64-pz)
              ki4(4) = qmacro * px               * (1.0_real64-pz)

              sum_q_xz_local(ix  ,iz  ) = sum_q_xz_local(ix  ,iz  ) + ki4(1)
              sum_q_xz_local(ix+1,iz  ) = sum_q_xz_local(ix+1,iz  ) + ki4(2)
              sum_q_xz_local(ix+1,iz+1) = sum_q_xz_local(ix+1,iz+1) + ki4(3)
              sum_q_xz_local(ix  ,iz+1) = sum_q_xz_local(ix  ,iz+1) + ki4(4)
            end if

            if (dtype(igrid) == 3_int32 .or. dtype(igrid) == 4_int32) then
              py = (real(iy,real64)*h(2) - yp_new) / h(2)

              ki4(1) = qmacro * py               * pz
              ki4(2) = qmacro * (1.0_real64-py) * pz
              ki4(3) = qmacro * (1.0_real64-py) * (1.0_real64-pz)
              ki4(4) = qmacro * py               * (1.0_real64-pz)

              d_ind = dtype(igrid)

              sum_q_yz_local(d_ind-2,iy  ,iz  ) = sum_q_yz_local(d_ind-2,iy  ,iz  ) + ki4(1)
              sum_q_yz_local(d_ind-2,iy+1,iz  ) = sum_q_yz_local(d_ind-2,iy+1,iz  ) + ki4(2)
              sum_q_yz_local(d_ind-2,iy+1,iz+1) = sum_q_yz_local(d_ind-2,iy+1,iz+1) + ki4(3)
              sum_q_yz_local(d_ind-2,iy  ,iz+1) = sum_q_yz_local(d_ind-2,iy  ,iz+1) + ki4(4)
            end if

          end if
        end if

        ! ------------------------------------------
        ! Legacy wall diagnostics accumulation
        ! ------------------------------------------

        if (igrid < 0_int32) igrid = 0_int32

        Ek_eV = 0.5_real64 * mass_species * &
          (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new) / qe

        p_mac_boundary(1,igrid) = p_mac_boundary(1,igrid) + qmacro
        p_mac_boundary(2,igrid) = p_mac_boundary(2,igrid) + abs(qmacro) * Ek_eV
        
        Ek_J = 0.5_real64 * mass_species * &
            (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new)

        P_loss_wall = P_loss_wall + Nm_species * Ek_J


        np_lost = np_lost + 1_int32
        cycle
      end if

      ! Survivors: Neumann reflection only on LHS.
      if (flag_nmn == 1_int32) then
        if (xp_new <= 0.0_real64) then
          xp_new  = -xp_new
          vpx_new = -vpx_new
        end if
      end if

      ! Legacy periodic wrap happens only for surviving particles.
      if (flag_pbc == 1_int32) then
        if (yp_new >= ymax) yp_new = yp_new - ymax
        if (yp_new <= 0.0_real64) yp_new = ymax + yp_new
        if (zp_new >= zmax) zp_new = zp_new - zmax
        if (zp_new <= 0.0_real64) zp_new = zmax + zp_new
      end if

      i_shift = i - np_lost

      part%x(i_shift)  = xp_new
      part%y(i_shift)  = yp_new
      part%z(i_shift)  = zp_new
      part%vx(i_shift) = vpx_new
      part%vy(i_shift) = vpy_new
      part%vz(i_shift) = vpz_new

      if (allocated(part%w))         part%w(i_shift)         = part%w(i)
      if (allocated(part%sp))        part%sp(i_shift)        = part%sp(i)
      if (allocated(part%flag_dead)) part%flag_dead(i_shift) = 0_int8
      if (allocated(part%flag_cex))  part%flag_cex(i_shift)  = part%flag_cex(i)

    end do

    part%n = part%n - np_lost
    if (part%n < 0_int32) part%n = 0_int32

    if (allocated(part%flag_dead)) then
      if (part%n < part%nmax) part%flag_dead(part%n+1:part%nmax) = 0_int8
    end if
    if (allocated(part%flag_cex)) then
      if (part%n < part%nmax) part%flag_cex(part%n+1:part%nmax) = 0_int32
    end if

    if (allocated(part%cell_id))    part%cell_id    = 0_int32
    if (allocated(part%cell_count)) part%cell_count = 0_int32
    if (allocated(part%cell_start)) part%cell_start = 0_int32

  end subroutine apply_particle_bc_legacy

    pure real(real64) function self_charge_current(qmacro) result(val)
    real(real64), intent(in) :: qmacro
    val = qmacro
  end function self_charge_current


  ! pure real(real64) function particle_energy_eV(vx,vy,vz) result(Ek)
  !   use mod_constants, only: qe
  !   real(real64), intent(in) :: vx,vy,vz

  !   real(real64) :: v2

  !   v2 = vx*vx + vy*vy + vz*vz

  !   ! electron mass convention handled outside
  !   Ek = 0.5_real64 * mass_species * v2 / qe
  ! end function particle_energy_eV

end module mod_particleBC