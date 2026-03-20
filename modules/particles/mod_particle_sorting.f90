module mod_particle_sorting
  use iso_fortran_env, only: int32, real64, int8
  use mod_particles,   only: ParticleSet
  implicit none
  private

  public :: sort_particles_by_cell
  public :: compute_particle_cell_ids
  public :: cell_index_from_position
  public :: check_particles_are_sorted
  public :: check_cell_indexing
  public :: get_cell_particle_range

contains

  pure integer(int32) function cell_index_from_position(x, y, z, h, n) result(ic)
    !=============================================================
    ! Return the 1D physical cell index corresponding to particle
    ! position, using the same convention as the legacy code:
    !
    !   ic = (ix-1) + n(1)*((iy-1) + n(2)*(iz-1)) + 1
    !
    ! where:
    !   ix = INT(x/hx) + 1
    !   iy = INT(y/hy) + 1
    !   iz = INT(z/hz) + 1
    !
    ! Valid cell ids are therefore:
    !   1 ... n(1)*n(2)*n(3)
    !
    ! Particle positions are clamped into the physical cell range.
    !=============================================================
    real(real64),   intent(in) :: x, y, z
    real(real64),   intent(in) :: h(3)
    integer(int32), intent(in) :: n(3)

    integer(int32) :: ix, iy, iz

    ix = int(x / h(1), int32) + 1_int32
    iy = int(y / h(2), int32) + 1_int32
    iz = int(z / h(3), int32) + 1_int32

    ix = max(1_int32, min(n(1), ix))
    iy = max(1_int32, min(n(2), iy))
    iz = max(1_int32, min(n(3), iz))

    ic = (ix - 1_int32) + n(1) * ((iy - 1_int32) + n(2) * (iz - 1_int32)) + 1_int32
  end function cell_index_from_position


  subroutine compute_particle_cell_ids(part, n, h)
    !=============================================================
    ! Compute cell_id(i) for each active particle in ParticleSet.
    ! Also fills:
    !   - cell_count(ic): number of particles in cell ic
    !   - cell_start(ic): first particle index of cell ic
    !
    ! This routine does NOT reorder the particle arrays.
    !=============================================================
    class(ParticleSet), intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)

    integer(int32) :: i, ic, ncells

    ncells = n(1) * n(2) * n(3)

    call part%ensure_cell_storage(ncells)

    part%cell_count = 0_int32
    part%cell_start = 0_int32

    if (part%n <= 0_int32) return

    do i = 1, part%n
      ic = cell_index_from_position(part%x(i), part%y(i), part%z(i), h, n)
      part%cell_id(i)     = ic
      part%cell_count(ic) = part%cell_count(ic) + 1_int32
    end do

    part%cell_start(1) = 1_int32
    do ic = 2, ncells
      part%cell_start(ic) = part%cell_start(ic-1) + part%cell_count(ic-1)
    end do
  end subroutine compute_particle_cell_ids


  subroutine sort_particles_by_cell(part, n, h)
    !=============================================================
    ! Reorder particles so that particles belonging to the same cell
    ! are contiguous in memory.
    !
    ! On exit:
    !   - x/y/z, vx/vy/vz, w, sp are sorted by cell
    !   - cell_id(i) is sorted and nondecreasing
    !   - cell_count(ic) is valid
    !   - cell_start(ic) is valid
    !
    ! This is the modern replacement of the legacy particle sorting
    ! + Plist construction, at the level of one ParticleSet.
    !=============================================================
    class(ParticleSet), intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)

    integer(int32) :: np, ncells
    integer(int32) :: i, ic, pos
    integer(int32), allocatable :: next_slot(:)
    integer(int32), allocatable :: cell_id_new(:)
    integer(int32), allocatable :: sp_new(:)

    real(real64), allocatable :: x_new(:), y_new(:), z_new(:)
    real(real64), allocatable :: vx_new(:), vy_new(:), vz_new(:)
    real(real64), allocatable :: w_new(:)

    integer(int8),  allocatable :: flag_dead_new(:)
    integer(int32), allocatable :: flag_cex_new(:)

    if (.not. allocated(part%x)) return

    np     = part%n
    ncells = n(1) * n(2) * n(3)

    call part%ensure_cell_storage(ncells)

    if (np <= 1_int32) then
      call compute_particle_cell_ids(part, n, h)
      return
    end if

    call compute_particle_cell_ids(part, n, h)

    allocate(next_slot(ncells))
    next_slot = part%cell_start

    allocate(x_new(np), y_new(np), z_new(np))
    allocate(vx_new(np), vy_new(np), vz_new(np))
    allocate(w_new(np))
    allocate(sp_new(np))
    allocate(cell_id_new(np))

    allocate(flag_dead_new(np))
    allocate(flag_cex_new(np))

    do i = 1, np
      ic  = part%cell_id(i)
      pos = next_slot(ic)

      x_new(pos)       = part%x(i)
      y_new(pos)       = part%y(i)
      z_new(pos)       = part%z(i)
      vx_new(pos)      = part%vx(i)
      vy_new(pos)      = part%vy(i)
      vz_new(pos)      = part%vz(i)
      w_new(pos)       = part%w(i)
      sp_new(pos)      = part%sp(i)
      cell_id_new(pos) = ic
      flag_dead_new(pos) = part%flag_dead(i)
      flag_cex_new(pos)  = part%flag_cex(i)

      next_slot(ic) = next_slot(ic) + 1_int32
    end do

    part%x(1:np)       = x_new
    part%y(1:np)       = y_new
    part%z(1:np)       = z_new
    part%vx(1:np)      = vx_new
    part%vy(1:np)      = vy_new
    part%vz(1:np)      = vz_new
    part%w(1:np)       = w_new
    part%sp(1:np)      = sp_new
    part%cell_id(1:np) = cell_id_new
    part%flag_dead(1:np) = flag_dead_new
    part%flag_cex(1:np)  = flag_cex_new

    deallocate(next_slot)
    deallocate(x_new, y_new, z_new)
    deallocate(vx_new, vy_new, vz_new)
    deallocate(w_new, sp_new, cell_id_new)
    deallocate(flag_dead_new, flag_cex_new)
  end subroutine sort_particles_by_cell


  logical function check_particles_are_sorted(part) result(ok)
    !=============================================================
    ! Return .true. if cell_id(:) is nondecreasing over active
    ! particles.
    !=============================================================
    class(ParticleSet), intent(in) :: part
    integer(int32) :: i

    ok = .true.

    if (.not. allocated(part%x)) return
    if (part%n <= 1_int32) return

    if (.not. allocated(part%cell_id)) then
      ok = .false.
      return
    end if

    do i = 2, part%n
      if (part%cell_id(i) < part%cell_id(i-1)) then
        ok = .false.
        return
      end if
    end do
  end function check_particles_are_sorted


  logical function check_cell_indexing(part, n) result(ok)
    !=============================================================
    ! Strong consistency check for the Plist-equivalent metadata:
    !
    !   cell_id(i)
    !   cell_count(ic)
    !   cell_start(ic)
    !
    ! Checks:
    !   1) all cell ids lie in [1,ncells]
    !   2) cell_count matches a direct scan over cell_id
    !   3) cell_start/cell_count define valid contiguous ranges
    !   4) each particle in a cell range really belongs to that cell
    !=============================================================
    class(ParticleSet), intent(in) :: part
    integer(int32),     intent(in) :: n(3)

    integer(int32) :: ncells
    integer(int32) :: icell, i, i0, i1
    integer(int32), allocatable :: counts_scan(:)

    ok = .true.

    if (.not. allocated(part%x)) return
    if (part%n <= 0_int32) return

    if (.not. allocated(part%cell_id)) then
      ok = .false.
      return
    end if
    if (.not. allocated(part%cell_start)) then
      ok = .false.
      return
    end if
    if (.not. allocated(part%cell_count)) then
      ok = .false.
      return
    end if

    ncells = n(1) * n(2) * n(3)

    if (size(part%cell_count) /= ncells) then
      ok = .false.
      return
    end if

    if (size(part%cell_start) /= ncells) then
      ok = .false.
      return
    end if

    allocate(counts_scan(ncells))
    counts_scan = 0_int32

    do i = 1, part%n
      if (part%cell_id(i) < 1_int32 .or. part%cell_id(i) > ncells) then
        ok = .false.
        deallocate(counts_scan)
        return
      end if
      counts_scan(part%cell_id(i)) = counts_scan(part%cell_id(i)) + 1_int32
    end do

    do icell = 1, ncells
      if (part%cell_count(icell) /= counts_scan(icell)) then
        ok = .false.
        deallocate(counts_scan)
        return
      end if
    end do

    do icell = 1, ncells
      if (part%cell_count(icell) > 0_int32) then
        i0 = part%cell_start(icell)
        i1 = i0 + part%cell_count(icell) - 1_int32

        if (i0 < 1_int32 .or. i1 > part%n) then
          ok = .false.
          deallocate(counts_scan)
          return
        end if

        do i = i0, i1
          if (part%cell_id(i) /= icell) then
            ok = .false.
            deallocate(counts_scan)
            return
          end if
        end do
      end if
    end do

    deallocate(counts_scan)
  end function check_cell_indexing


  subroutine get_cell_particle_range(part, icell, i0, i1, count)
    !=============================================================
    ! Helper for collisions:
    ! return the contiguous particle range for one cell.
    !
    ! If the cell is empty:
    !   count = 0
    !   i0 = 0
    !   i1 = -1
    !=============================================================
    class(ParticleSet), intent(in)  :: part
    integer(int32),     intent(in)  :: icell
    integer(int32),     intent(out) :: i0, i1, count

    count = 0_int32
    i0    = 0_int32
    i1    = -1_int32

    if (.not. allocated(part%cell_count)) return
    if (.not. allocated(part%cell_start)) return

    if (icell < 1_int32 .or. icell > size(part%cell_count)) return

    count = part%cell_count(icell)
    if (count <= 0_int32) return

    i0 = part%cell_start(icell)
    i1 = i0 + count - 1_int32
  end subroutine get_cell_particle_range

end module mod_particle_sorting