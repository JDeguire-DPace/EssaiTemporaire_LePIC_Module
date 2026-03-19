module mod_particleBC
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  implicit none
  private

  public :: apply_particle_bc_legacy

contains

  pure integer(int32) function clamp_index(i, ilo, ihi) result(ic)
    integer(int32), intent(in) :: i, ilo, ihi
    ic = max(ilo, min(ihi, i))
  end function clamp_index


  subroutine apply_particle_bc_legacy(part, n, h, bcnd, xmax, ymax, zmax, flag_pbc, flag_nmn, ptype, tag_neg)
    !===========================================================
    ! First boundary-condition layer adapted from legacy
    ! part_expmover.f90.
    !
    ! Implemented:
    !   - lost particle detection using 8-corner bcnd test
    !   - special negative-ion loss rule at x < 0 if flag_nmn=1
    !   - Neumann reflection at x=0 for particles that remain alive
    !   - periodic wrap in y and z if flag_pbc=1
    !   - compaction of surviving particles
    !
    ! Not yet implemented:
    !   - dielectric charge deposition
    !   - secondary emission
    !   - wall power accounting
    !   - cex flags / beam diagnostics
    !   - MPI exchange
    !===========================================================
    class(ParticleSet), intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    integer(int32),     intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(in)    :: xmax, ymax, zmax
    integer(int32),     intent(in)    :: flag_pbc
    integer(int32),     intent(in)    :: flag_nmn
    integer(int32),     intent(in)    :: ptype
    integer(int32),     intent(in)    :: tag_neg

    integer(int32) :: i, i_shift
    integer(int32) :: ix, iy, iz
    integer(int32) :: flag_lost
    integer(int32) :: np_lost

    real(real64) :: xp_new, yp_new, zp_new
    real(real64) :: vpx_new, vpy_new, vpz_new

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

      !-------------------------------------------
      ! Get particle grid location from new position
      !-------------------------------------------
      ix = floor(xp_new / h(1)) + 1_int32
      iy = floor(yp_new / h(2)) + 1_int32
      iz = floor(zp_new / h(3)) + 1_int32

      if (ix < 0_int32)        ix = 0_int32
      if (ix > n(1)+1_int32)   ix = n(1)+1_int32
      if (iy < 0_int32)        iy = 0_int32
      if (iy > n(2)+1_int32)   iy = n(2)+1_int32
      if (iz < 0_int32)        iz = 0_int32
      if (iz > n(3)+1_int32)   iz = n(3)+1_int32

      !-------------------------------------------
      ! Check if particle is inside a solid region
      !-------------------------------------------
      flag_lost = 0_int32
      if ( bcnd(ix  ,iy  ,iz  ) >= 1_int32 .and. &
           bcnd(ix+1,iy  ,iz  ) >= 1_int32 .and. &
           bcnd(ix+1,iy+1,iz  ) >= 1_int32 .and. &
           bcnd(ix  ,iy+1,iz  ) >= 1_int32 .and. &
           bcnd(ix  ,iy  ,iz+1) >= 1_int32 .and. &
           bcnd(ix+1,iy  ,iz+1) >= 1_int32 .and. &
           bcnd(ix+1,iy+1,iz+1) >= 1_int32 .and. &
           bcnd(ix  ,iy+1,iz+1) >= 1_int32 ) then
        flag_lost = 1_int32
      end if

      !-------------------------------------------
      ! Legacy special case:
      ! negative ions + Neumann BC at LHS => lost
      !-------------------------------------------
      if (tag_neg > 0_int32) then
        if (ptype == tag_neg) then
          if (xp_new < 0.0_real64 .and. flag_nmn == 1_int32) flag_lost = 2_int32
        end if
      end if

      !-------------------------------------------
      ! Surviving particles: apply BCs and compact
      !-------------------------------------------
      if (flag_lost == 0_int32) then

        ! Neumann BC on LHS only: specular reflection
        if (flag_nmn == 1_int32) then
          if (xp_new <= 0.0_real64) then
            xp_new  = -xp_new
            vpx_new = -vpx_new
          end if
        end if

        ! Periodic BCs in y and z
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

        if (allocated(part%w))  part%w(i_shift)  = part%w(i)
        if (allocated(part%sp)) part%sp(i_shift) = part%sp(i)

      else
        np_lost = np_lost + 1_int32
      end if

    end do

    part%n = part%n - np_lost
    if (part%n < 0_int32) part%n = 0_int32

    ! Invalidate sorting metadata after compaction / BC handling
    if (allocated(part%cell_id))    part%cell_id    = 0_int32
    if (allocated(part%cell_count)) part%cell_count = 0_int32
    if (allocated(part%cell_start)) part%cell_start = 0_int32

  end subroutine apply_particle_bc_legacy

end module mod_particleBC