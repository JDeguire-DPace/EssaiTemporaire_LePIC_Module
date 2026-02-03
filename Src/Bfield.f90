subroutine read_Bfield_map(Bi,n,name,namlen,B_dir,scaling,mpi_rank)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      Oct/11
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  integer:: i,ix,iy,iz,nx,ny,nz,n(3),B_dir,namlen,nmax,mpi_rank
  real(kind=8):: Bi(4,0:n(1)+2,0:n(2)+2,0:n(3)+2),B,scaling
  character:: name*20

  nmax= (n(1)+3)*(n(2)+3)*(n(3)+3)
  
  if(mpi_rank.eq.0) print*, 'Opening: ',name(1:namlen)
  open(30,file='mag_field_files/'//name(1:namlen),status='OLD')

  nx= 0
  ny= 0
  nz= 0

  ! Generate B-field map 
  do i=1,nmax
     read(30,*,end=999) ix,iy,iz,B
     B= B*scaling
     Bi(B_dir,ix,iy,iz)= Bi(B_dir,ix,iy,iz) + B   
     nx= MAX(nx,ix)
     ny= MAX(ny,iy)
     nz= MAX(nz,iz)
  enddo
  
999 continue

  if(mpi_rank.eq.0) &
       print'(1x,"nx= ",i4,", ny= ",i4,", nz= ",i4,", scaling= ",f5.2)', nx,ny,nz,scaling
  
  close(30)

  return

end subroutine read_Bfield_map

subroutine gaussian_Bfield(Bi,n,h,B0,x0,Lx,B_dir)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      Oct/11
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  implicit none
  integer:: ix,iy,iz,n(3),B_dir
  real(kind=8):: h(3),Bi(4,0:n(1)+2,0:n(2)+2,0:n(3)+2),x,B0,x0,Lx

  ! Generate B-field map 
  do ix=0,n(1)+2
     x= (ix-1)*h(1)
     do iy=0,n(2)+2
        do iz=0,n(3)+2
           Bi(B_dir,ix,iy,iz)= Bi(B_dir,ix,iy,iz) + & ! update B-field array
                B0*dexp( -(x-x0)**2/(2*Lx**2) )
        enddo
     enddo
  enddo

  return

end subroutine gaussian_Bfield

subroutine PG_current(Bi,n,h,B0,B_dir)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Jun/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:    1D fit for the magnetic filter field profile 
!                  in ELISE device, IPP Garching.
!                  Franzen et al., Plasma Phys. Control. Fusion 56 
!                  (2014) p. 025007, Fig 5.
!     --------------------------------------------------------------
  implicit none
  integer:: ix,iy,iz,n(3),B_dir
  real(kind=8):: h(3),Bi(4,0:n(1)+2,0:n(2)+2,0:n(3)+2), &
       a(5),b(2),x,B0,B_tmp
     
  a(1) = 0.0175596       
  a(2) = 0.0227362       
  a(3) = 0.00193113      
  a(4) = -8.51235e-05    
  a(5) = 9.58673e-07     

  b(1) = 15.7416
  b(2) = -0.377132

  ! Generate B-field map 
  do ix=0,n(1)+2

     x= (ix-1)*h(1)*1.d2 ! in cm
     if(x.le.39.d0) then 
        B_tmp=  a(1) + a(2)*x + a(3)*x**2 + a(4)*x**3 + &
             a(5)*x**4
     else
        B_tmp=  b(1) + b(2)*x
     endif

     do iy=0,n(2)+2
        do iz=0,n(3)+2
           Bi(B_dir,ix,iy,iz)= Bi(B_dir,ix,iy,iz) + B0*B_tmp ! update B-field array
        enddo
     enddo
  enddo

  return

end subroutine PG_current

subroutine EE_magnets(Bi,n,h,B0,x0,z0,d)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Jan/16
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      d is the distance between the extraction magnets
!                    x0 is the position of the magnets
!                    z0 is the location in the middle of the magnets
!     --------------------------------------------------------------
  implicit none
  integer:: ix,iy,iz,n(3)
  real(kind=8):: h(3),Bi(4,0:n(1)+2,0:n(2)+2,0:n(3)+2),x,z,B0,x0,z0,d
  include 'constants.h'

  ! Generate B-field map 
  do ix=0,n(1)+2
     x= (ix-1)*h(1)
     do iy=0,n(2)+2
        do iz=0,n(3)+2
           z= (iz-1)*h(3)
           Bi(1,ix,iy,iz)= Bi(1,ix,iy,iz) + & ! Bx
                B0*dsin(pi*(z-z0)/d)*dexp(-pi*(x0-x)/d)
           Bi(3,ix,iy,iz)= Bi(3,ix,iy,iz) + & ! Bz
                B0*dcos(pi*(z-z0)/d)*dexp(-pi*(x0-x)/d)
        enddo
     enddo
  enddo

  return

end subroutine EE_magnets

subroutine Bfield_Hall_thruster(Bi,n,h,B0,x0,a)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      July/18
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:   B0= 100 G (max at x=x0)
!                 x0= 0.75 cm
!                 y0= 1.25 cm (channel axis)
!                 d= 1 cm (channel length)
!                 a= 1 cm (radial distance for the magnets)
!     --------------------------------------------------------------
  implicit none
  integer:: ix,iy,iz,n(3)
  real(kind=8):: h(3),Bi(4,0:n(1)+2,0:n(2)+2,0:n(3)+2),x,y,x0,a,d,&
       B0,B1,B2,alpha,y0
  real(kind=8):: Bx_magnet,By_magnet

  d= 1.5d-2
  y0= 2.d-2
  alpha=0.75d0 ! B2/B1
  B1= ABS(B0)/(1-alpha*a**2/(4.d0*d**2+a**2)) 
  B2=alpha*B1

  ! Pseudo 1D configuration
  if(B0.lt.0) y0= 0.d0
  
  ! Generate B-field map 
  do iz=0,n(3)+2
     do iy=0,n(2)+2
        do ix=0,n(1)+2
           x= (ix-1)*h(1)
           y= (iy-1)*h(2)
           ! Bx
           Bi(1,ix,iy,iz)= Bi(1,ix,iy,iz) + & 
                B1*Bx_magnet(x-x0,y-y0,a) - B2*Bx_magnet(x-x0+2.d0*d,y-y0,a)
           ! By
           Bi(2,ix,iy,iz)= Bi(2,ix,iy,iz) + & 
                B1*By_magnet(x-x0,y-y0,a) - B2*By_magnet(x-x0+2.d0*d,y-y0,a)
        enddo
     enddo
  enddo

  if(B0.lt.0) then
     Bi(1,:,:,:)=0.d0 ! Bx=0
     do ix=0,n(1)+2
        Bi(2,ix,:,:)= Bi(2,ix,1,1) ! By uniform along (0y) 
     enddo
  endif
  
  return

end subroutine Bfield_Hall_thruster

function Bx_magnet(x,y,a)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      July/18
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: x component of the analytically expression for the Bfield 
!               profile of a Hall thruster
!     --------------------------------------------------------------
  implicit none
  real(kind=8):: x,y,a,Bx_magnet

  Bx_magnet= 0.5d0*( a*x/(x**2+(y+a)**2) - a*x/(x**2+(y-a)**2) )

  return
end function Bx_magnet

function By_magnet(x,y,a)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      July/18
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: y component of the analytically expression for the Bfield 
!               profile of a Hall thruster
!     --------------------------------------------------------------
  implicit none
  real(kind=8):: x,y,a,By_magnet

  By_magnet= 0.5d0*( a*(y+a)/(x**2+(y+a)**2) - a*(y-a)/(x**2+(y-a)**2) )

  return
end function By_magnet

subroutine find_B_dir(B_info,B_dir)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      Oct/11
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  implicit none
  integer:: B_dir
  character:: B_info*2

  ! Find direction of magnetic field
  if( B_info.eq.'Bx' .or. B_info.eq.'BX' .or. &
       B_info.eq.'bx' .or. B_info.eq.'bX' ) B_dir=1
  if( B_info.eq.'By' .or. B_info.eq.'BY' .or. &
       B_info.eq.'by' .or. B_info.eq.'bY' ) B_dir=2
  if( B_info.eq.'Bz' .or. B_info.eq.'BZ' .or. &
       B_info.eq.'bz' .or. B_info.eq.'bZ' ) B_dir=3
  
  return

end subroutine find_B_dir
