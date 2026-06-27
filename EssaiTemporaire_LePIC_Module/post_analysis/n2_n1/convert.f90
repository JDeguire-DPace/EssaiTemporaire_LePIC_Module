program convert
  implicit none
  integer:: ix,iy,nx,ny,n
  parameter (n=500)
  real(kind=8):: n2_n1(n+1,n+1),np(2,n+1,n+1)

  !
  ! Read Macho format
  !
  open(12,file='../../DATA/n1_xy.mco')
  read(12,*) nx,ny
  print*, 'nx=',nx,' ny=',ny
  do iy=ny+1,1,-1
     read(12,*) ( np(1,ix,iy), ix=1,nx+1 ) 
  enddo
  close(12)

  open(12,file='../../DATA/n2_xy.mco')
  read(12,*) nx,ny
  do iy=ny+1,1,-1
     read(12,*) ( np(2,ix,iy), ix=1,nx+1 ) 
  enddo
  close(12)

  do iy=2,ny
     do ix=2,nx
        n2_n1(ix,iy)=0.d0
        if(np(1,ix,iy).gt.0.d0) n2_n1(ix,iy)= np(2,ix,iy)/np(1,ix,iy)-1.d0        
     enddo
  enddo

  !
  ! Print in Macho format
  !
  open(12,file='n2_n1.mco')
  write(12,*) nx,ny
  do iy=ny+1,1,-1
     write(12,100) ( n2_n1(ix,iy), ix=1,nx+1 ) 
  enddo
  close(12)

  100 format(800(e18.6,1x))

end program convert
