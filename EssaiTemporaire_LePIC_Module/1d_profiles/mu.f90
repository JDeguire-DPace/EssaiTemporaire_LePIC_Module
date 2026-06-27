  program mu

  implicit none
  integer:: i,nmax,cnt,flux,Ex,ne,Te
  parameter (nmax=1000,flux=1,Ex=2,ne=3,Te=4)
  real(kind=8):: data(4,nmax),gradP_n(1000),dx,qe

  cnt=0
  dx= 78.125d-6
  qe=1.6d-19
  
  open(10,file='sum_dens.dat',status='OLD')

  do i=1,nmax
     read(10,*,end=999) data(flux,i),data(Ex,i),data(ne,i),data(Te,i)
     cnt= cnt + 1
  enddo

999  close(10)

  print*, 'cnt=',cnt

  do i=2,cnt-1
     gradP_n(i)= (data(ne,i+1)*data(Te,i+1)-data(ne,i-1)*data(Te,i-1))/(2.d0*dx*data(ne,i))
  enddo

  open(10,file='mu.dat')
  do i=2,cnt-1
     write(10,100) (data(flux,i)/data(ne,i))/(-(data(Ex,i)+gradP_n(i))),data(flux,i)/data(ne,i),data(Ex,i),gradP_n(i)
  enddo
  close(10)

100  format(800(e18.6,1x))
  
end program mu
