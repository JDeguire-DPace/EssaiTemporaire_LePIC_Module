module mod_particleBC
  use iso_fortran_env, only: int32, int8, real64
  use mod_particles,   only: ParticleSet
  implicit none
  private

  public :: apply_particle_bc_legacy

contains

  subroutine apply_particle_bc_legacy(part, n, h, bcnd, xmax, ymax, zmax, flag_pbc, flag_nmn, ptype, tag_neg)
    !===========================================================
    ! Legacy-style particle boundary handling:
    !
    ! Implemented:
    !   - explicit box-crossing loss detection
    !   - solid-region loss detection using 8-corner bcnd test
    !   - special negative-ion loss rule at x < 0 if flag_nmn=1
    !   - Neumann reflection at x=0 for surviving particles
    !   - periodic wrap in y and z if flag_pbc=1
    !   - compaction of surviving particles
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

      flag_lost = 0_int32

      ! --------------------------------------------------------
      ! 1) Explicit box-crossing tests
      ! --------------------------------------------------------

      ! X boundaries
      if (xp_new < 0.0_real64) then
        if (flag_nmn == 1_int32) then
          ! For negative ions, legacy can force loss at x<0
          if (tag_neg > 0_int32 .and. ptype == tag_neg) then
            flag_lost = 2_int32
          end if
        else
          flag_lost = 1_int32
        end if
      end if

      if (xp_new > xmax) then
        flag_lost = 1_int32
      end if

      ! Y/Z boundaries
      if (flag_pbc /= 1_int32) then
        if (yp_new < 0.0_real64 .or. yp_new > ymax) flag_lost = 1_int32
        if (zp_new < 0.0_real64 .or. zp_new > zmax) flag_lost = 1_int32
      end if

      ! --------------------------------------------------------
      ! 2) Legacy solid-region test with clamped indices
      ! Only do this if particle was not already lost by box exit
      ! --------------------------------------------------------
      if (flag_lost == 0_int32) then

        ix = floor(xp_new / h(1)) + 1_int32
        iy = floor(yp_new / h(2)) + 1_int32
        iz = floor(zp_new / h(3)) + 1_int32

        if (ix < 0_int32)      ix = 0_int32
        if (ix > n(1)+1_int32) ix = n(1)+1_int32
        if (iy < 0_int32)      iy = 0_int32
        if (iy > n(2)+1_int32) iy = n(2)+1_int32
        if (iz < 0_int32)      iz = 0_int32
        if (iz > n(3)+1_int32) iz = n(3)+1_int32

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

      end if

      ! --------------------------------------------------------
      ! 3) Surviving particles: apply BCs and compact
      ! --------------------------------------------------------
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
          if (yp_new <  0.0_real64) yp_new = ymax + yp_new
          if (zp_new >= zmax) zp_new = zp_new - zmax
          if (zp_new <  0.0_real64) zp_new = zmax + zp_new
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
        if (allocated(part%flag_dead)) part%flag_dead(i_shift) = part%flag_dead(i)
        if (allocated(part%flag_cex))  part%flag_cex(i_shift)  = part%flag_cex(i)

      else
        np_lost = np_lost + 1_int32
      end if

    end do

    part%n = part%n - np_lost
    if (part%n < 0_int32) part%n = 0_int32

    ! Optional cleanup of inactive tail
    if (allocated(part%flag_dead)) then
      if (part%n < part%nmax) part%flag_dead(part%n+1:part%nmax) = 0_int8
    end if
    if (allocated(part%flag_cex)) then
      if (part%n < part%nmax) part%flag_cex(part%n+1:part%nmax) = 0_int32
    end if

    ! Invalidate sorting metadata after compaction / BC handling
    if (allocated(part%cell_id))    part%cell_id    = 0_int32
    if (allocated(part%cell_count)) part%cell_count = 0_int32
    if (allocated(part%cell_start)) part%cell_start = 0_int32

  end subroutine apply_particle_bc_legacy

end module mod_particleBC