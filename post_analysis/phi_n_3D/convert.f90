program convert
  implicit none
  integer:: ix,iy,iz,ix0,iy0,iz0,nx,ny,nz,it,imin,imax,i_pl,flag_plt,flag_sav
  real(kind=8), allocatable:: phi_n(:,:,:,:),Ei(:,:,:,:)
  real(kind=8):: dx,dy,dz
  character:: plnum*3,corrnum*3

  open(12,file='input.dat')
  read(12,*) ix0,iy0,iz0
  read(12,*) imin,imax
  read(12,*) dx,dy,dz
  close(12)
  
  flag_plt=0
  flag_sav=0
  if(dz.lt.0) then
     flag_sav=1
     dz= ABS(dz)
  endif

  ! Save time vs Ez(ix0,iy0,iz)
  open(13,file='Ez_vs_time.mco')
  open(14,file='n2_vs_time.mco')

  if(imin.eq.0 .or. imax.eq.0) then
     imin=0
     imax=0
  endif
  
  do it=imin,imax

     !
     ! Read file in binary format
     !
     if( imin.eq.0 .or. imax.eq.0 ) then
        open(12,file='../../DATA/phi_n_3D.dat',form='UNFORMATTED')
     else
        flag_plt=1

        if(it.lt.10) then 
           write (plnum,'(i1)'),it
           i_pl=1
        endif
        if( it.ge.10 .and. it.lt.100 ) then 
           write (plnum,'(i2)'),it
           i_pl=2
        endif
        if( it.ge.100 .and. it.lt.1000 ) then 
           write (plnum,'(i3)'),it
           i_pl=3
        endif
                 
        if(i_pl.eq.1) corrnum= '_00'
        if(i_pl.eq.2) corrnum= '_0'
        if(i_pl.eq.3) corrnum= '_'
        
        open(12,file='../../DATA/phi_n_3D'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.dat',form='UNFORMATTED')        
     endif
     read(12) nx,ny,nz
     print*, 'it=',it,', nx=',nx,', ny=',ny,', nz=',nz
     allocate( phi_n(2,nx+1,ny+1,nz+1), Ei(3,nx+1,ny+1,nz+1) )
     read(12) phi_n(:,1:nx+1,1:ny+1,1:nz+1)           
     close(12)

     if(it.eq.imin) write(13,*) nz,imax-imin+1
     
     !
     ! Calculate E-field
     !
     Ei=0.d0
     
     do iz= 2,nz ! Interior points
        do iy= 2,ny
           do ix= 2,nx
              if(phi_n(2,ix,iy,iz).gt.0.d0) then
                 Ei(1,ix,iy,iz)= -( phi_n(1,ix+1,iy,iz) - phi_n(1,ix-1,iy,iz) )/(2.d0*dx)
                 Ei(2,ix,iy,iz)= -( phi_n(1,ix,iy+1,iz) - phi_n(1,ix,iy-1,iz) )/(2.d0*dy)
                 Ei(3,ix,iy,iz)= -( phi_n(1,ix,iy,iz+1) - phi_n(1,ix,iy,iz-1) )/(2.d0*dz)
              endif
           enddo
        enddo
     enddo
     
     !
     ! Print in Macho format
     !

     if(flag_sav.eq.1) goto 1000
     
     ! XY plane
     if(flag_plt.eq.0) then
        open(12,file='phi_xy.mco')
     else
        open(12,file='phi_xy'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) nx,ny
     do iy=ny+1,1,-1
        write(12,100) ( phi_n(1,ix,iy,iz0), ix=1,nx+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='n2_xy.mco')
     else
        open(12,file='n2_xy'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) nx,ny
     do iy=ny+1,1,-1
        write(12,100) ( phi_n(2,ix,iy,iz0), ix=1,nx+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='Ex_xy.mco')
     else
        open(12,file='Ex_xy'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) nx,ny
     do iy=ny+1,1,-1
        write(12,100) ( Ei(1,ix,iy,iz0), ix=1,nx+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='Ey_xy.mco')
     else
        open(12,file='Ey_xy'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) nx,ny
     do iy=ny+1,1,-1
        write(12,100) ( Ei(2,ix,iy,iz0), ix=1,nx+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='Ez_xy.mco')
     else
        open(12,file='Ez_xy'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) nx,ny
     do iy=ny+1,1,-1
        write(12,100) ( Ei(3,ix,iy,iz0), ix=1,nx+1 ) 
     enddo
     close(12)
     
     ! XZ plane
     if(flag_plt.eq.0) then
        open(12,file='phi_xz.mco')
     else
        open(12,file='phi_xz'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) nx,nz
     do iz=nz+1,1,-1
        write(12,100) ( phi_n(1,ix,iy0,iz), ix=1,nx+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='n2_xz.mco')
     else
        open(12,file='n2_xz'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) nx,nz
     do iz=nz+1,1,-1
        write(12,100) ( phi_n(2,ix,iy0,iz), ix=1,nx+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='Ex_xz.mco')
     else
        open(12,file='Ex_xz'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) nx,nz
     do iz=nz+1,1,-1
        write(12,100) ( Ei(1,ix,iy0,iz), ix=1,nx+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='Ey_xz.mco')
     else
        open(12,file='Ey_xz'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) nx,nz
     do iz=nz+1,1,-1
        write(12,100) ( Ei(2,ix,iy0,iz), ix=1,nx+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='Ez_xz.mco')
     else
        open(12,file='Ez_xz'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) nx,nz
     do iz=nz+1,1,-1
        write(12,100) ( Ei(3,ix,iy0,iz), ix=1,nx+1 ) 
     enddo
     close(12)
     
     ! YZ plane
     if(flag_plt.eq.0) then
        open(12,file='phi_yz.mco')
     else
        open(12,file='phi_yz'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) ny,nz
     do iz=nz+1,1,-1
        write(12,100) ( phi_n(1,ix0,iy,iz), iy=1,ny+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='n2_yz.mco')
     else
        open(12,file='n2_yz'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) ny,nz
     do iz=nz+1,1,-1
        write(12,100) ( phi_n(2,ix0,iy,iz), iy=1,ny+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='Ex_yz.mco')
     else
        open(12,file='Ex_yz'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) ny,nz
     do iz=nz+1,1,-1
        write(12,100) ( Ei(1,ix0,iy,iz), iy=1,ny+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='Ey_yz.mco')
     else
        open(12,file='Ey_yz'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) ny,nz
     do iz=nz+1,1,-1
        write(12,100) ( Ei(2,ix0,iy,iz), iy=1,ny+1 ) 
     enddo
     close(12)

     if(flag_plt.eq.0) then
        open(12,file='Ez_yz.mco')
     else
        open(12,file='Ez_yz'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco')
     endif
     write(12,*) ny,nz
     do iz=nz+1,1,-1
        write(12,100) ( Ei(3,ix0,iy,iz), iy=1,ny+1 ) 
     enddo
     close(12)

     ! Save time vs Ez(ix0,iy0,iz)
1000 write(13,100) ( Ei(3,ix0,iy0,iz), iz=1,nz+1 )
     write(14,100) ( phi_n(2,ix0,iy0,iz), iz=1,nz+1 ) 
     
     deallocate(phi_n,Ei)

  enddo

  close(13)
  close(14)
  
  100 format(1200(e18.6,1x))

end program convert
