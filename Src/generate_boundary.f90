subroutine generate_boundary(u,n,h,bcnd,V,ngrid,dtype,xl_pow,xr_pow,mpi_rank)
!     ==============================================================
!     VERSION:         0.6
!     LAST MOD:       Dec/23
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Generates boundary condition matrix
!     NOTE:          dtype is the wall structure type
!                    dtype=1 is for metals
!                         =2 dielectric surface in XZ plane
!                         =3-4 dielectric surfaces in YZ plane
!     --------------------------------------------------------------
  implicit none
  integer:: ig,ix,iy,iz,ind,ngrid,ixl,ixr,iyl,iyr,izl,izr,iyh,izh, &
       igl,igr,cnt_hy,cnt_hz,mpi_rank,dtype(0:ngrid),cnt_yz_planes,&
       ind_val,flag_circz,flag_circx,flag_circx_tmp
  integer:: n(3),bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
  integer :: udon
  real(kind=8):: h(3),u(0:n(1)+2,0:n(2)+2,0:n(3)+2),V(ngrid), &
       xl,xr,yl,yr,zl,zr,x,y,z,yd,zd,R,Sh,xl_pow,xr_pow
  include 'particle_info.h'
  include 'constants.h'

  ! Initialization
  iyh=0
  izh=0
  Lgy=0.d0
  Lgz=0.d0
  bcnd=0
  u=0.d0
  flag_pbc=0
  flag_pbcz=0
  flag_nmn=0
  flag_die=0
  cnt_yz_planes=0
  ig_die=0
  zg_sec= 0.d0
  n_cath=0
  R=0.d0
  every=1
  if(MAXVAL(n).gt.512) every=2
  if(MAXVAL(n).gt.1024) every=4
  if(MAXVAL(n).gt.2048) every=8
  ! Read potential values
  if( flag_restart.eq.1 .and. flag_convP.eq.1 ) then
     ! When Vgrd is modified iteratively
     open(41,file='DATA.BAK/Vgrd.bak',form='UNFORMATTED')
     read(41) V(1:ngrid)
     close(41)
  else
     open(10,file='input_dir/boundary.inp')
     do ig=1,ngrid
        read(10,*,end=90) V(ig)
     enddo
  endif
  
90 continue
  close(10)

  ! Open input files
  open(10,file='input_dir/geometry.inp')

  ! Define box size
  read(10,*) xmax,ymax,zmax

  ! Neumann BC on the LHS
  if(xmax.lt.0) flag_nmn=1
  xmax= ABS(xmax)

  ! Option to draw cylindrical segments along (Ox)
  flag_circx=0
  flag_circx_tmp=0
  if(ymax.lt.0) then
     flag_circx_tmp=1
  endif
  ymax= ABS(ymax)

  ! Draw cylindrical segments along (Oz)
  flag_circz=0
  if(zmax.lt.0) then
     flag_circz=1
  endif
  zmax= ABS(zmax)
  
  ! Define dx,dy,dz
  h(1)=xmax/n(1)
  h(2)=ymax/n(2)   
  h(3)=zmax/n(3)   

  !
  ! Draw walls and grid segments
  !
  do ig=1,100

     read(10,*,end=50) xl,yl,zl,xr,yr,zr,ind
     dtype(ABS(ind))= 1 ! Metallic wall structure type

     if(flag_circx_tmp.eq.1) then
        flag_circx=1
     else
        flag_circx=0
     endif
     if(yl.lt.0) then
        ! Reset flag
        if(flag_circx_tmp.eq.1) flag_circx=0
        if(flag_circx_tmp.eq.0) flag_circx=1
     endif
     yl= ABS(yl)

     if(flag_circx.eq.1) then
        ! Radius of segment
        R= MIN((yr-yl)/2.d0,(zr-zl)/2.d0)
        ! Cylindrical heating power deposition profile is set
        if( xl.le.xl_pow .and. xr.ge.xr_pow ) flag_circxh=1
     endif
     
     if(flag_circz.eq.1) R= MIN((xr-xl)/2.d0,(yr-yl)/2.d0)
     
     if(ind.eq.0) flag_pbc=1
     if(xr.lt.0) then ! Dielectric along X=cste planes
        xr= ABS(xr)
        if(cnt_yz_planes.eq.0) then 
           dtype(ABS(ind))= 3 ! plane #1
        else
           dtype(ABS(ind))= 4 ! plane #2
        endif
        cnt_yz_planes= cnt_yz_planes + 1
        flag_die=1
        ig_die(dtype(ABS(ind)))= INT( xl/h(1) ) + 1
     endif
     if(yr.lt.0) then ! Dielectric along Y=cste planes
        yr= ABS(yr)
        dtype(ABS(ind))=2
        flag_die=1
        ig_die(2)= INT( yl/h(2) ) + 1
     endif
     if(zr.lt.0) then  ! Draw one surface as periodic
        zr= ABS(zr)
        dtype(ABS(ind))= -dtype(ABS(ind))
        flag_pbc=1
        flag_pbcz=1
     endif

     ! Warning
     if( ABS(ind).gt.ngrid ) then
        if(mpi_rank.eq.0) then
           print*, 'Insufficient number of wall labels found in file boundary.inp '
           print*, 'please correct ...'
        endif
        call stop_calculation
     endif

     ixl=INT( xl/h(1) ) + 1
     ! LHS simulation box wall
     if(ixl.le.1) ixl=0

     ixr=INT( xr/h(1) ) + 1
     ! RHS simulation box wall
     if(ixr.ge.n(1)) ixr= n(1)+2
 
     iyl=INT( yl/h(2) ) + 1
     ! South side
     if(iyl.le.1) iyl=0

     iyr=INT( yr/h(2) ) + 1
     ! North side
     if(iyr.ge.n(2)) iyr= n(2)+2

     izl=INT( zl/h(3) ) + 1
     ! Bottom
     if(izl.le.1) izl=0

     izr=INT( zr/h(3) ) + 1
     ! Top
     if(izr.ge.n(3)) izr= n(3)+2

     ! Secondary particle emission along (Oz)
     dir_sec=1
     if( (gam_sec.gt.0.d0 .or. ABS(opt_inj).eq.4) .and. igrid_sec.eq.ABS(ind) ) then
        ! Warning
        if(zr.ne.zl) then
           print*, 'Warning: secondary particle emission model only along (Oz), please correct...'
           call stop_calculation
        endif
           
        if(zl.eq.0.d0) then ! @ Bottom of the simulation domain
           zg_sec(1)= zr
           n_cath= n_cath + 1
        endif
        
        if(zr.eq.zmax) then ! Top 
           zg_sec(2)= zl
           dir_sec=-1
           n_cath= n_cath + 1
        endif
     endif
     
     ! Grid index where holes will be drawn
     if(ABS(ind).eq.ind_g) then
        igl= ixl
        igr= ixr
        Lgy= MIN((iyr-iyl)*h(2),ymax)
        Lgz= MIN((izr-izl)*h(3),zmax)
        ixg= NINT(real(igl+igr)/2.d0)
        xg1= (ixl-1)*h(1)
     endif

     if(ind.eq.0) then ! Periodic BC's
        ixl=MAX(ixl,2)
        ixr=MIN(ixr,n(1))
     endif

     ! Draw walls
     if(ind.ge.0) then
        !$OMP PARALLEL
        !$OMP DO
        do iz=0,n(3)+2
           do iy=0,n(2)+2
              do ix=ixl,ixr
                 bcnd(ix,iy,iz)=ind
                 if(ind.gt.0) u(ix,iy,iz)=V(ind)
              enddo
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL   
     endif
     
     ! Interior domain
     if(ind.ge.0) then
        ixl=MAX(ixl,2)
        ixr=MIN(ixr,n(1))
        iyl=MAX(iyl,2)
        iyr=MIN(iyr,n(2))
        izl=MAX(izl,2)
        izr=MIN(izr,n(3))
     endif
     
     !$OMP PARALLEL
     !$OMP DO
     do iz=izl,izr
        z= (iz-1)*h(3)
        do iy=iyl,iyr
           y= (iy-1)*h(2)
           do ix=ixl,ixr
              x= (ix-1)*h(1)
              if(ind.ge.0) then
                 ! Draw cylinder with axis along (Ox) or (Oz)
                 if(flag_circz.eq.1 .and. ind.gt.0 .and.  &
                      ((x-xmax/2.d0)**2+(y-ymax/2.d0)**2).gt.R**2 ) goto 40
                 
                 if(flag_circx.eq.1 .and. ind.gt.0 .and. &
                      ((y-ymax/2.d0)**2+(z-zmax/2.d0)**2).gt.R**2 ) goto 40
                 
                 bcnd(ix,iy,iz)=-1
                 u(ix,iy,iz)= 0.d0
40               continue
              else
                 if( flag_circz.eq.1 ) then
                    ! Draw disk in Z=cste plane
                    if( ((x-xmax/2.d0)**2+(y-ymax/2.d0)**2).gt.R**2 ) goto 45
                 endif

                 if( flag_circx.eq.1 ) then
                    ! Draw disk in X=cste plane
                    if( ((y-ymax/2.d0)**2+(z-zmax/2.d0)**2).gt.R**2 ) goto 45
                 endif

                 bcnd(ix,iy,iz)=ABS(ind)
                 if(ABS(ind).gt.0) u(ix,iy,iz)=V(ABS(ind))
45               continue
              endif
           enddo
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL   

  enddo

  ! Neumann BC's (LHS only)
50  if(flag_nmn.eq.1) then
     bcnd(0:1,:,:)=-2
     u(0:1,:,:)=0
  endif

  ! Combine a dielectric or metal surface with periodic BC's in the Z=1 & Z=n+1 planes
  if(flag_pbc.eq.1) then
     !$OMP PARALLEL
     !$OMP DO
     do iy=0,n(2)+2
        do ix=0,n(1)+2
           ig= bcnd(ix,iy,1)
           if(dtype(ig).lt.0) then
              if( bcnd(ix,iy,2).eq.-1 ) then
                 bcnd(ix,iy,0:1)= 0
                 bcnd(ix,iy,n(3)+1:n(3)+2)= 0
              endif
           endif
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL   
  endif

  ! Reset dtype
  dtype= ABS(dtype)

  !
  ! Draw Holes
  !
  if(flag_grd.eq.0) goto 100

  iyh= NINT( Lhy/2.d0/h(2) ) + 1 
  izh= NINT( Lhz/2.d0/h(3) ) + 1 
     
  if(Lhy.gt.ymax) Lhy=ymax
  if(Lhz.gt.zmax) Lhz=zmax

  if(Lhy.eq.Lhz) then ! Disk
     R= Lhy/2.d0
     Sh= pi*R**2
  else ! Slit
     Sh= Lhy*Lhz
  endif

  if(mpi_rank.eq.0) then
     print'(1x,"Corrected dimensions of holes, Ly(cm)= ",f6.2,", Lz(cm)= ",f6.2)', &
          (2*iyh+1)*h(2),(2*izh+1)*h(3) ! Size is +1 larger
  endif
     
  do cnt_hz=1,nhz

     ! z-location of aperture
     zd= -(nhz-1)*3.d0*Lhz/4.d0 + (cnt_hz-1)*3.d0*Lhz/2.d0 + zmax/2.d0
     izl= NINT( zd/h(3) ) + 1 - izh
     izr= izl + 2*izh

     if( (izl.lt.1 .or. izr.gt.n(3)+1) .and. nhz.gt.1 ) then
           if(mpi_rank.eq.0) then
              print*, 'Warning: too many apertures along Oz'
              print*, 'Please correct ...'
           endif
           call stop_calculation
     endif

     if(izl.lt.1) izl=1
     if(izr.gt.n(3)+1) izr=n(3)+1

     do cnt_hy=1,nhy

        ! y-location of aperture
        yd= -(nhy-1)*3.d0*Lhy/4.d0 + (cnt_hy-1)*3.d0*Lhy/2.d0 + ymax/2.d0
        iyl= NINT( yd/h(2) ) + 1 - iyh
        iyr= iyl + 2*iyh

        if( (iyl.lt.1 .or. iyr.gt.n(2)+1) .and. nhy.gt.1 ) then 
           if(mpi_rank.eq.0) then
              print*, 'Warning: too many apertures along Oy'
              print*, 'Please correct ...'
           endif
           call stop_calculation
        endif

        if(iyl.lt.1) iyl=1
        if(iyr.gt.n(2)+1) iyr=n(2)+1

        do iz=izl,izr
           do iy=iyl,iyr
              y= (iy-1)*h(2) 
              z= (iz-1)*h(3)
              ind_val= -1
              if( iy.eq.1 .or. iy.eq.n(2)+1 .or. iz.eq.1 .or. iz.eq.n(3)+1 ) &
                   ind_val= 0
              do ix=igl,igr
                 if( Lhy.eq.Lhz ) then ! Draw a disk
                    if( ((y-yd)**2+(z-zd)**2).le.R**2 ) then                       
                       bcnd(ix,iy,iz)= ind_val 
                       u(ix,iy,iz)= 0.d0
                    endif
                 else ! Draw a rectangle
                    bcnd(ix,iy,iz)= ind_val 
                    u(ix,iy,iz)= 0.d0
                 endif
              enddo
           enddo
        enddo

     enddo
  enddo

  ! Convert units in meters
100 h= h*1.d-2
  xmax= xmax*1.d-2
  ymax= ymax*1.d-2
  zmax= zmax*1.d-2

  xg1= xg1*1.d-2
  Lgy= Lgy*1.d-2
  Lgz= Lgz*1.d-2

  Sg= Lgy*Lgz - (cnt_hy-1)*(cnt_hz-1)*Sh*1.d-4

  zg_sec= zg_sec*1.d-2

  if(mpi_rank.eq.0) then 
     print'(1x,"Type of wall surfaces (1=metal, 2 & 3= dielectric):",10(1x,i2))',dtype(1:ngrid)
     if(flag_grd.eq.1) then
        print'(1x,"Dimension of grid #",i2," : Lgy(cm)= ",f6.2,", Lgz(cm)= ",f6.2)', &
             ind_g,Lgy*1.d2,Lgz*1.d2
        print'(1x,"Grid surface (cm2): ",f6.2,", surface occupied by apertures (cm2): ",f6.2)',&
             Sg*1.d4,(cnt_hy-1)*(cnt_hz-1)*Sh
     endif
  endif

  close(10)

  !  
  ! Write boundaries in files
  !
  if(mpi_rank.eq.0) then

     ! 2D XY (Z=0) plane
     izl= n(3)/2+1
     open(12,file='DATA/DATA_2D/bcnd_xy.mco')
     write(12,*) n(1)/every,n(2)/every
     do iy=n(2)+1,1,-1*every
        write(12,101) ( bcnd(ix,iy,izl), ix=1,n(1)+1,every )
    101     format(800(i3,1x)) 
     enddo
     close(12)

     ! 2D XZ (Y=0) plane
     open(12,file='DATA/DATA_2D/bcnd_xz.mco')
     write(12,*) n(1)/every,n(3)/every
     do iz=n(3)+1,1,-1*every
        write(12,101) ( bcnd(ix,n(2)/2+1,iz), ix=1,n(1)+1,every )
     enddo
     close(12)

     ! 2D YZ plane
     ix=n(1)/2+1
     if(flag_grd.eq.1) ix= ixg
     open(12,file='DATA/DATA_2D/bcnd_yz.mco')
     write(12,*) n(2)/every,n(3)/every
     do iz=n(3)+1,1,-1*every 
        write(12,101) ( bcnd(ix,iy,iz), iy=1,n(2)+1,every ) 
     enddo
     close(12)
        
  endif
  
  
   ! Full 3D boundary condition file
   if(mpi_rank.eq.0) then
      write(*,*) "nx, ny, nz = ", n(1), n(2), n(3)
      open(newunit=udon,file="DATA/bcnd_type_full.dat",status="replace",action="write")
      write(udon,*) n(1), n(2), n(3)
      do iz=0,n(3)+2
      do iy=0,n(2)+2
         write(udon,*) (bcnd, ix=0,n(1)+2)
      end do
      end do
      close(udon)
   endif
  return
end subroutine generate_boundary

