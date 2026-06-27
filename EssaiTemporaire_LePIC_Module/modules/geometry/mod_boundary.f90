module mod_boundary
  use iso_fortran_env, only: real64
  use mod_domain,      only: Domain
  use mod_config,      only: Config
  use mod_fields,      only: Fields
  use mod_generateBoundary
  implicit none
  private
  public :: build_boundary, GeometryData, read_geometry_file, read_boundary_potentials

  type :: GeometryData
    real(real64) :: xmax=0.0_real64, ymax=0.0_real64, zmax=0.0_real64
    integer      :: nseg = 0
    real(real64), allocatable :: xl(:), yl(:), zl(:), xr(:), yr(:), zr(:)
    integer,      allocatable :: ind(:)
  end type GeometryData
contains

  subroutine read_geometry_file(fname, geo, dom)
    character(len=*), intent(in)  :: fname
    type(GeometrySegments), intent(out) :: geo
    type(Domain), intent(inout) :: dom
    integer :: i

    open(10,file=fname)
    read(10,*) geo%xmax, geo%ymax, geo%zmax
    dom%xmax = geo%xmax/100.0_real64
    dom%ymax = geo%ymax/100.0_real64
    dom%zmax = geo%zmax/100.0_real64
    ! allocate max 100 like legacy
    geo%nseg = 0
    allocate(geo%xl(100),geo%yl(100),geo%zl(100),geo%xr(100),geo%yr(100),geo%zr(100),geo%ind(100))

    do i=1,100
      read(10,*,end=20) geo%xl(i),geo%yl(i),geo%zl(i),geo%xr(i),geo%yr(i),geo%zr(i),geo%ind(i)
      geo%nseg = geo%nseg + 1
    end do
   20 continue
    close(10)
  end subroutine read_geometry_file


  subroutine read_boundary_potentials(fname, V, ngrid)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: ngrid
    real(real64), intent(out) :: V(ngrid)
    integer :: ig
    V = 0.0_real64
    open(10,file=fname)
    do ig=1,ngrid
      read(10,*,end=10) V(ig)
    end do
   10 continue
    close(10)
  end subroutine read_boundary_potentials




  subroutine build_boundary(dom, cfg, fld, mpi_rank)
    implicit none
    type(Domain), intent(inout) :: dom
    type(Config), intent(in)    :: cfg
    type(Fields), intent(inout) :: fld  
    integer,      intent(in)    :: mpi_rank
    type(GeometrySegments) :: geo
    type(BoundaryInputs)  :: inp
    type(BoundaryOutputs) :: outp
    real(real64), allocatable :: V(:)
    integer,      allocatable :: dtype(:)

    ! --- allocate temps ---
    allocate(V(cfg%ngrid))
    allocate(dtype(0:cfg%ngrid))
    V     = 0.0_real64
    dtype = 0


    call read_boundary_potentials('../input_dir/boundary.inp', V, cfg%ngrid)
    call read_geometry_file('../input_dir/geometry.inp', geo, dom)

    ! --- ensure dom%dtype bounds ---
    if (allocated(dom%dtype)) then
      if (lbound(dom%dtype,1) /= 0 .or. ubound(dom%dtype,1) /= cfg%ngrid) then
        deallocate(dom%dtype)
      end if
    end if
    if (.not. allocated(dom%dtype)) allocate(dom%dtype(0:cfg%ngrid))

    inp%flag_restart = cfg%flag_restart
    inp%flag_convP   = cfg%flag_convP
    inp%flag_grd     = cfg%flag_grd

    inp%ind_g = cfg%ind_g
    inp%nhy   = cfg%nhy
    inp%nhz   = cfg%nhz
    inp%Lhy   = cfg%Lhy
    inp%Lhz   = cfg%Lhz

    inp%gam_sec   = cfg%gam_sec
    inp%opt_inj   = cfg%opt_inj
    inp%igrid_sec = cfg%igrid_sec

    ! ------------------------------------------------------------
    ! call kernel wrapper
    ! ------------------------------------------------------------
    call generate_boundary(fld%phi , dom%n, dom%h, dom%bcnd, V, cfg%ngrid, dtype, &
                           cfg%xl_pow, cfg%xr_pow, mpi_rank, inp, outp, geo)

    ! ------------------------------------------------------------
    ! outp -> dom (explicit outputs)
    ! ------------------------------------------------------------
    dom%xmax = outp%xmax
    dom%ymax = outp%ymax
    dom%zmax = outp%zmax

    dom%Lgy  = outp%Lgy
    dom%Lgz  = outp%Lgz
    dom%Sg   = outp%Sg
    dom%xg1  = outp%xg1
    dom%ixg  = outp%ixg

    dom%flag_pbc  = outp%flag_pbc
    dom%flag_pbcz = outp%flag_pbcz
    dom%flag_nmn  = outp%flag_nmn
    dom%flag_die  = outp%flag_die

    ! dtype -> dom
    dom%dtype = dtype

    deallocate(V, dtype)
  end subroutine build_boundary

end module mod_boundary
