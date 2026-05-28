module mod_collision_flat_backend
  use iso_fortran_env, only: int8, int32, real64
  use mod_particles,   only: ParticleSet

  implicit none
  private

  public :: FlatCollisionBuffers
  public :: init_flat_collision_buffers
  public :: destroy_flat_collision_buffers
  public :: pack_particles_to_flat_sorted
  public :: unpack_flat_to_particles
  public :: flat_capacity_report

  type :: FlatCollisionBuffers
    integer(int32) :: ntype  = 0_int32
    integer(int32) :: nproc  = 0_int32
    integer(int32) :: nmax   = 0_int32
    integer(int32) :: pl_max = 0_int32

    ! Legacy-style flat particle arrays:
    !   vxp(1,:,:,:) = x
    !   vxp(2,:,:,:) = y
    !   vxp(3,:,:,:) = z
    !   vxp(4,:,:,:) = vx
    !   vxp(5,:,:,:) = vy
    !   vxp(6,:,:,:) = vz
    real(real64), allocatable :: vxp(:,:,:,:)

    ! Number of active particles per species/thread.
    integer(int32), allocatable :: np_tot(:,:)

    ! Legacy cumulative particle list by cell:
    !   Plist(icell,ptype,iproc) = number of particles with cell_id <= icell
    !   particles in cell icell are from Plist(icell-1)+1 to Plist(icell)
    integer(int32), allocatable :: Plist(:,:,:)

    ! Legacy dead flags, one per flat particle.
    integer(int8), allocatable :: flag_dead(:,:,:)

    ! Legacy CEX flag. In the old code this is indexed as flag_cex(ib,bproc).
    ! It is mostly relevant for negative ions; kept here for compatibility.
    integer(int32), allocatable :: flag_cex(:,:)

  end type FlatCollisionBuffers

contains

  subroutine init_flat_collision_buffers(buf, part, n, extra_fraction, extra_min)
    type(FlatCollisionBuffers), intent(inout) :: buf
    type(ParticleSet), intent(in)             :: part(:,:)
    integer(int32), intent(in)                :: n(3)
    real(real64), intent(in), optional        :: extra_fraction
    integer(int32), intent(in), optional      :: extra_min

    integer(int32) :: ptype, iproc
    integer(int32) :: ntype, nproc, nmax_needed, extra_i
    real(real64)   :: frac

    call destroy_flat_collision_buffers(buf)

    ntype = int(size(part,1), int32)
    nproc = int(size(part,2), int32)

    frac = 0.50_real64
    if (present(extra_fraction)) frac = extra_fraction

    extra_i = 100000_int32
    if (present(extra_min)) extra_i = extra_min

    nmax_needed = 0_int32
    do iproc = 1_int32, nproc
      do ptype = 1_int32, ntype
        nmax_needed = max(nmax_needed, part(ptype,iproc)%n)
      end do
    end do

    ! Same idea as legacy nmax: one common capacity for all species/procs.
    buf%nmax = nmax_needed + max(extra_i, int(frac * real(max(1_int32,nmax_needed), real64), int32))
    buf%ntype = ntype
    buf%nproc = nproc
    buf%pl_max = n(1)*n(2)*n(3)

    allocate(buf%vxp(6, buf%nmax, ntype, nproc))
    allocate(buf%np_tot(ntype, nproc))
    allocate(buf%Plist(0:buf%pl_max, ntype, nproc))
    allocate(buf%flag_dead(buf%nmax, ntype, nproc))
    allocate(buf%flag_cex(buf%nmax, nproc))

    buf%vxp       = 0.0_real64
    buf%np_tot    = 0_int32
    buf%Plist     = 0_int32
    buf%flag_dead = 0_int8
    buf%flag_cex  = 0_int32
  end subroutine init_flat_collision_buffers


  subroutine destroy_flat_collision_buffers(buf)
    type(FlatCollisionBuffers), intent(inout) :: buf

    if (allocated(buf%vxp))       deallocate(buf%vxp)
    if (allocated(buf%np_tot))    deallocate(buf%np_tot)
    if (allocated(buf%Plist))     deallocate(buf%Plist)
    if (allocated(buf%flag_dead)) deallocate(buf%flag_dead)
    if (allocated(buf%flag_cex))  deallocate(buf%flag_cex)

    buf%ntype  = 0_int32
    buf%nproc  = 0_int32
    buf%nmax   = 0_int32
    buf%pl_max = 0_int32
  end subroutine destroy_flat_collision_buffers


  subroutine pack_particles_to_flat_sorted(buf, part, n, h)
    type(FlatCollisionBuffers), intent(inout) :: buf
    type(ParticleSet), intent(in)             :: part(:,:)
    integer(int32), intent(in)                :: n(3)
    real(real64), intent(in)                  :: h(3)

    integer(int32) :: ptype, iproc, i, k
    integer(int32) :: ix, iy, iz, icell
    integer(int32), allocatable :: cell_count(:)
    integer(int32), allocatable :: cursor(:)

    if (.not. allocated(buf%vxp)) error stop "pack_particles_to_flat_sorted: buffers not allocated"
    if (int(size(part,1),int32) /= buf%ntype) error stop "pack_particles_to_flat_sorted: wrong ntype"
    if (int(size(part,2),int32) /= buf%nproc) error stop "pack_particles_to_flat_sorted: wrong nproc"

    buf%vxp       = 0.0_real64
    buf%np_tot    = 0_int32
    buf%Plist     = 0_int32
    buf%flag_dead = 0_int8
    buf%flag_cex  = 0_int32

    allocate(cell_count(buf%pl_max))
    allocate(cursor(buf%pl_max))

    do iproc = 1_int32, buf%nproc
      do ptype = 1_int32, buf%ntype

        if (.not. allocated(part(ptype,iproc)%x)) cycle

        if (part(ptype,iproc)%n > buf%nmax) then
          write(*,*) "pack_particles_to_flat_sorted: n > nmax"
          write(*,*) "ptype, iproc, n, nmax = ", ptype, iproc, part(ptype,iproc)%n, buf%nmax
          error stop
        end if

        buf%np_tot(ptype,iproc) = part(ptype,iproc)%n

        cell_count = 0_int32

        ! Count particles per cell.
        do i = 1_int32, part(ptype,iproc)%n
          ! if (allocated(part(ptype,iproc)%flag_dead)) then
          !   if (part(ptype,iproc)%flag_dead(i) /= 0) cycle
          ! end if

          icell = particle_cell_id(part(ptype,iproc)%x(i), &
                                   part(ptype,iproc)%y(i), &
                                   part(ptype,iproc)%z(i), n, h)
          cell_count(icell) = cell_count(icell) + 1_int32
        end do

        ! Build cumulative Plist.
        buf%Plist(0,ptype,iproc) = 0_int32
        do icell = 1_int32, buf%pl_max
          buf%Plist(icell,ptype,iproc) = buf%Plist(icell-1,ptype,iproc) + cell_count(icell)
        end do

        if (ptype == 1_int32) then
          write(*,*) "DEBUG pack iproc=", iproc, &
                    " n=", part(ptype,iproc)%n, &
                    " plist=", buf%Plist(buf%pl_max,ptype,iproc), &
                    " diff=", part(ptype,iproc)%n - buf%Plist(buf%pl_max,ptype,iproc)
        end if

        cursor = 0_int32

        ! Fill flat arrays sorted by cell.
        do i = 1_int32, part(ptype,iproc)%n
          ! if (allocated(part(ptype,iproc)%flag_dead)) then
          !   if (part(ptype,iproc)%flag_dead(i) /= 0) cycle
          ! end if

          icell = particle_cell_id(part(ptype,iproc)%x(i), &
                                   part(ptype,iproc)%y(i), &
                                   part(ptype,iproc)%z(i), n, h)

          cursor(icell) = cursor(icell) + 1_int32
          k = buf%Plist(icell-1,ptype,iproc) + cursor(icell)

          if (k < 1_int32 .or. k > buf%nmax) then
            write(*,*) "pack_particles_to_flat_sorted: bad packed index"
            write(*,*) "ptype, iproc, i, k, nmax = ", ptype, iproc, i, k, buf%nmax
            error stop
          end if

          buf%vxp(1,k,ptype,iproc) = part(ptype,iproc)%x(i)
          buf%vxp(2,k,ptype,iproc) = part(ptype,iproc)%y(i)
          buf%vxp(3,k,ptype,iproc) = part(ptype,iproc)%z(i)
          buf%vxp(4,k,ptype,iproc) = part(ptype,iproc)%vx(i)
          buf%vxp(5,k,ptype,iproc) = part(ptype,iproc)%vy(i)
          buf%vxp(6,k,ptype,iproc) = part(ptype,iproc)%vz(i)

          if (allocated(part(ptype,iproc)%flag_dead)) then
            buf%flag_dead(k,ptype,iproc) = int(part(ptype,iproc)%flag_dead(i), int8)
          else
            buf%flag_dead(k,ptype,iproc) = 0_int8
          end if

          if (allocated(part(ptype,iproc)%flag_cex)) then
            if (ptype == 1_int32) then
              ! Keep this conservative. Real tag_neg-specific handling can be added later.
              buf%flag_cex(k,iproc) = part(ptype,iproc)%flag_cex(i)
            end if
          end if
        end do

        ! Important: after removing dead particles during packing, update active count.
        ! buf%np_tot(ptype,iproc) = buf%Plist(buf%pl_max,ptype,iproc)

      end do
    end do

    deallocate(cell_count)
    deallocate(cursor)
  end subroutine pack_particles_to_flat_sorted


  subroutine unpack_flat_to_particles(buf, part)
    type(FlatCollisionBuffers), intent(in) :: buf
    type(ParticleSet), intent(inout)       :: part(:,:)

    integer(int32) :: ptype, iproc, i, nnew

    if (.not. allocated(buf%vxp)) error stop "unpack_flat_to_particles: buffers not allocated"

    do iproc = 1_int32, buf%nproc
      do ptype = 1_int32, buf%ntype

        nnew = buf%np_tot(ptype,iproc)
        call part(ptype,iproc)%ensure_capacity(nnew)

        part(ptype,iproc)%n = nnew

        if (nnew <= 0_int32) cycle

        part(ptype,iproc)%x (1:nnew) = buf%vxp(1,1:nnew,ptype,iproc)
        part(ptype,iproc)%y (1:nnew) = buf%vxp(2,1:nnew,ptype,iproc)
        part(ptype,iproc)%z (1:nnew) = buf%vxp(3,1:nnew,ptype,iproc)
        part(ptype,iproc)%vx(1:nnew) = buf%vxp(4,1:nnew,ptype,iproc)
        part(ptype,iproc)%vy(1:nnew) = buf%vxp(5,1:nnew,ptype,iproc)
        part(ptype,iproc)%vz(1:nnew) = buf%vxp(6,1:nnew,ptype,iproc)

        if (allocated(part(ptype,iproc)%flag_dead)) then
          do i = 1_int32, nnew
            part(ptype,iproc)%flag_dead(i) = int(buf%flag_dead(i,ptype,iproc), kind(part(ptype,iproc)%flag_dead(i)))
          end do
        end if

        if (allocated(part(ptype,iproc)%flag_cex)) then
          if (ptype == 1_int32) then
            part(ptype,iproc)%flag_cex(1:nnew) = buf%flag_cex(1:nnew,iproc)
          else
            part(ptype,iproc)%flag_cex(1:nnew) = 0
          end if
        end if

        if (allocated(part(ptype,iproc)%sp)) part(ptype,iproc)%sp(1:nnew) = ptype
        if (allocated(part(ptype,iproc)%w))  part(ptype,iproc)%w(1:nnew)  = 1.0_real64

        ! Sorting/cell metadata will be rebuilt by your normal sorter.
        if (allocated(part(ptype,iproc)%cell_id))    part(ptype,iproc)%cell_id    = 0_int32
        if (allocated(part(ptype,iproc)%cell_count)) part(ptype,iproc)%cell_count = 0_int32
        if (allocated(part(ptype,iproc)%cell_start)) part(ptype,iproc)%cell_start = 0_int32

      end do
    end do
  end subroutine unpack_flat_to_particles


  subroutine flat_capacity_report(buf)
    type(FlatCollisionBuffers), intent(in) :: buf

    if (.not. allocated(buf%np_tot)) then
      write(*,*) "FlatCollisionBuffers not allocated."
      return
    end if

    write(*,*) "FlatCollisionBuffers:"
    write(*,*) "  ntype  = ", buf%ntype
    write(*,*) "  nproc  = ", buf%nproc
    write(*,*) "  nmax   = ", buf%nmax
    write(*,*) "  pl_max = ", buf%pl_max
    write(*,*) "  max(np_tot) = ", maxval(buf%np_tot)
  end subroutine flat_capacity_report


  pure integer(int32) function particle_cell_id(x, y, z, n, h) result(icell)
    real(real64), intent(in)    :: x, y, z
    integer(int32), intent(in)  :: n(3)
    real(real64), intent(in)    :: h(3)

    integer(int32) :: ix, iy, iz

    ix = int(x / h(1), int32) + 1_int32
    iy = int(y / h(2), int32) + 1_int32
    iz = int(z / h(3), int32) + 1_int32

    ix = max(1_int32, min(n(1), ix))
    iy = max(1_int32, min(n(2), iy))
    iz = max(1_int32, min(n(3), iz))

    icell = (ix - 1_int32) + n(1) * ((iy - 1_int32) + n(2) * (iz - 1_int32)) + 1_int32
  end function particle_cell_id

end module mod_collision_flat_backend