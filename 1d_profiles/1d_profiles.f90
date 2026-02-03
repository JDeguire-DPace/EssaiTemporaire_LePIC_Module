program prof
  implicit none
  integer:: ic,opt,ix,iy,ix0,iy0,n,nx,ny,lgh,cnt
  parameter (n=1000)
  real(kind=8):: ni(n,n),ni_avg,ni2_avg
  character:: name*20,pnum*3

  open(9,file='namelist.inp')
  read(9,*,end=1000) name
  print*, 'Reading file: ',name

  do ic=1,20
     if( name(ic:ic).eq.'.' .or. name(ic:ic).eq.' ' ) exit
  enddo

  ! Open file
  open(10,file='../DATA/DATA_2D/'//name(1:ic-1)//'.mco',status='OLD')

  ! Read file
  read(10,*,end=1000) nx,ny
  do iy=ny+1,1,-1
     read(10,100) ( ni(ix,iy), ix=1,nx+1 )
100  format(800(e18.6,1x))
  enddo


  ! Plot options
  print*, 'along Ox=1, Oy=2?'
  read(*,*) opt

  if(opt.eq.1) then
     print*, 'ix=?'
     read(*,*) ix0
     write (pnum,'(i3)'),ix0
     if(ix0.le.9) lgh=3
     if(ix0.ge.10 .and. ix0.le.99 ) lgh=2
     if(ix0.ge.100 .and. ix0.le.999 ) lgh=1        
     open(11,file=name(1:ic-1)//'_ix_'//pnum(lgh:3)//'.dat')
  else
     print*, 'iy=?'
     read(*,*) iy0
     write (pnum,'(i3)'),iy0
     if(iy0.le.9) lgh=3
     if(iy0.ge.10 .and. iy0.le.99 ) lgh=2
     if(iy0.ge.100 .and. iy0.le.999 ) lgh=1        
     open(11,file=name(1:ic-1)//'_iy_'//pnum(lgh:3)//'.dat')
  endif
  
  ! Draw plot
  ni_avg= 0.d0
  ni2_avg= 0.d0
  cnt= 0
  
  do iy=1,ny+1
     do ix=1,nx+1
        if( opt.eq.1 .and. ix.eq.ix0 ) then
           if( ix0-1.ge.1 .or. ix0+1.le.nx+1 ) then
              write(11,*) SUM(ni(ix0-1:ix0+1,iy))/3.d0
           else
              write(11,*) ni(ix0,iy)
           endif

           ni_avg= ni_avg + ni(ix0,iy)
           ni2_avg= ni2_avg + ni(ix0,iy)**2
           if(ABS(ni(ix0,iy)).ne.0) cnt= cnt + 1
        endif
        if( opt.eq.2 .and. iy.eq.iy0 ) then
           if( iy0-1.ge.1 .or. iy0+1.le.ny+1 ) then
              write(11,*) SUM(ni(ix,iy0-1:iy0+1))/3.d0
           else
              write(11,*) ni(ix,iy0)
           endif

           ni_avg= ni_avg + ni(ix,iy0)
           ni2_avg= ni2_avg + ni(ix,iy0)**2
           if(ABS(ni(ix,iy0)).ne.0) cnt= cnt + 1
        endif
     enddo
  enddo

  ni_avg= ni_avg/cnt
  ni2_avg= ni2_avg/cnt
  
  print*, 'cnt=',cnt,', <a>=',ni_avg,', sqrt(<a²>)=',dsqrt(ni2_avg),', a_rms=',dsqrt(ni2_avg-ni_avg**2)
  
  goto 1001

  ! Messages
1000 continue
  print*, 'Error opening file'

1001 continue
  print*, 'End of file reached'

end program prof
