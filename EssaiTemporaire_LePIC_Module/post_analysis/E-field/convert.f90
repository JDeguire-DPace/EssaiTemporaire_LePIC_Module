program convert
  implicit none
  integer:: ix,iy,nx,ny,n
  parameter (n=512)
  real(kind=8):: phi(n+1,n+1),E(2,n+1,n+1),dx,dy

  !
  ! Read Macho format
  !
  open(12,file='../../DATA/phi_xz.mco')
  read(12,*) nx,ny
  print*, 'nx=',nx,' ny=',ny
  do iy=ny+1,1,-1
     read(12,*) ( phi(ix,iy), ix=1,nx+1 ) 
  enddo
  close(12)

  open(13,file='input.dat')
  read(13,*) dx,dy
  close(13)
  dx= dx*1.d-3
  dy= dy*1.d-3

  do iy=2,ny
     do ix=2,nx
        E(1,ix,iy)= -( phi(ix+1,iy) - phi(ix-1,iy) )/(2.d0*dx)
        E(2,ix,iy)= -( phi(ix,iy+1) - phi(ix,iy-1) )/(2.d0*dy)
     enddo
  enddo

  !
  ! Print in Macho format
  !
  open(12,file='Ex.mco')
  write(12,*) nx,ny
  do iy=ny+1,1,-1
     write(12,100) ( E(1,ix,iy), ix=1,nx+1 ) 
  enddo
  close(12)

  open(12,file='Ey.mco')
  write(12,*) nx,ny
  do iy=ny+1,1,-1
     write(12,100) ( E(2,ix,iy), ix=1,nx+1 ) 
  enddo
  close(12)

  100 format(800(e18.6,1x))

end program convert
