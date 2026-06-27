  implicit none
  integer:: i,i_rg,nmax1,nmax2,cnt,cnt_err,ptype
  parameter (nmax1=10**7,nmax2=250)
  integer:: nvx(-nmax2:nmax2),nvy(-nmax2:nmax2),nvz(-nmax2:nmax2),nEk(0:2*nmax2)
  real(kind=8):: vx(nmax1),vy(nmax1),vz(nmax1),dv,int
  real(kind=8):: vavg(3),v2avg(3),v2rms(3),vt,Tp,Ek_tmp,Tp_tmp
  real(kind=8):: Ek(nmax1),me,qe,dE,Ek_max,int_v,int_E,fE,fE_min,Ek_all,amu,mass
  character:: pnum*1

  !
  ! Simulation parameters
  !
  qe=1.60217646e-19 ! Coulombs
  me=9.10938188d-31
  amu=1.66053886d-27
  
  print*, '1=electrons, 2=ions'
  read(*,*) ptype
  
  mass=0.d0
  if(ptype.eq.1) then 
     mass=me
  else
     print*, 'mass/amu'
     read(*,*) mass
     mass= mass*amu
  endif

  !
  ! Plot parameters
  !
  dv=0.1d0
  cnt=0
  cnt_err=0
  dE=0.1d0
  Ek_all=2*nmax2*dE

  print*, 'Which region index would you like to plot?'
  read(*,*) i_rg
  write (pnum,'(i1)'),i_rg
  open(10,file='../DATA/xv'//pnum//'.dat')
  open(11,file='pdf.dat')
  open(12,file='pedf.dat')

  !
  ! Initialize variables & arrays
  !
  Ek_max=0.d0

  do i=1,3
     vavg(i)=0.d0
     v2avg(i)=0.d0
  enddo

  do i=-nmax2,nmax2
     nvx(i)=0
     nvy(i)=0
     nvz(i)=0
     nEk(i+nmax2)=0
  enddo

  !
  ! Read phase space cooridnates and calculate macroscopic quantities
  !
  do i=1,nmax1

     read(10,*,end=999) vx(i),vy(i),vz(i)
     
     vavg(1)= vavg(1) + vx(i)
     vavg(2)= vavg(2) + vy(i)
     vavg(3)= vavg(3) + vz(i)

     v2avg(1)= v2avg(1) + vx(i)**2.
     v2avg(2)= v2avg(2) + vy(i)**2.
     v2avg(3)= v2avg(3) + vz(i)**2.

     cnt= cnt + 1

  enddo

999 print*, '# of particles=',cnt

  if(cnt.ge.nmax1) then ! Warning
     print*, 'Warning: nmax1 is too small, please correct ...'
     STOP
  endif

  vavg= vavg/cnt
  v2avg= v2avg/cnt

  v2rms(1)= v2avg(1) - vavg(1)**2.
  v2rms(2)= v2avg(2) - vavg(2)**2.
  v2rms(3)= v2avg(3) - vavg(3)**2.
  Tp_tmp= v2rms(1) + v2rms(2) + v2rms(3)
  Tp= (1.d0/3.d0)*mass*Tp_tmp/qe
  vt= dsqrt(2.d0*qe*Tp/mass)

  !
  ! Normalize variables and calculate histograms
  !
  do i=1,cnt

     Ek_tmp= vx(i)**2. + vy(i)**2. + vz(i)**2.
     Ek(i)= 0.5d0*mass*Ek_tmp/qe
     if(Ek(i).gt.Ek_all) then 
        cnt_err= cnt_err + 1
        goto 10
     endif
     Ek_max= MAX(Ek_max,Ek(i))

     nEk(NINT(Ek(i)/dE))=  nEk(NINT(Ek(i)/dE)) +1  

     vx(i)= vx(i)/vt
     nvx(NINT(vx(i)/dv))=  nvx(NINT(vx(i)/dv)) +1

     vy(i)= vy(i)/vt
     nvy(NINT(vy(i)/dv))=  nvy(NINT(vy(i)/dv)) +1

     vz(i)= vz(i)/vt
     nvz(NINT(vz(i)/dv))=  nvz(NINT(vz(i)/dv)) +1

10 continue
  enddo

  print*, 'particules with Ek > Ek_allowed:',cnt_err
  print*, 'vt(m/s)=',vt,', T(eV)=',Tp,', Ek_max(eV)=',Ek_max
  print*, '<vx> (m/s)=',vavg(1),', <vy>=',vavg(2),', <vz>=',vavg(3)
  print*, 'Tx (eV)=',mass*v2rms(1)/qe,', Ty=',mass*v2rms(2)/qe,', Tz=',mass*v2rms(3)/qe


  !
  ! Write to a file
  !

  int_v= cnt*dv
  int_E= (cnt-cnt_err)*dE

  do i=-nmax2,nmax2

     write(11,*) real(i*dv),nvx(i)/int_v,nvy(i)/int_v,nvz(i)/int_v

     fE= real(nEk(i+nmax2))/int_E/real( (i+nmax2+0.5)*dE )**0.5
     if( (i+nmax2).eq. NINT(Ek_max/dE) ) fE_min= fE
     write(12,*) real((i+nmax2)*dE),fE

  enddo


  call generate_gnuplot(Ek_max,nEk(1)/int_E/(3.*dE/2.)**0.5,Tp,fE_min,pnum)

  close(10)
  close(11)
  close(12)

end program

subroutine generate_gnuplot(xmax,x0,Tp,ymin,pnum)
  implicit none
  real(kind=8):: xmax,x0,Tp,ymin,ymax
  character:: pnum*1

  open(10,file='plot.gnu')

  write(10,*) 'set style line 1 lt 1 lw 3 # red'
  write(10,*) 'set style line 2 lt 2 lw 3 # green'
  write(10,*) 'set style line 3 lt 3 lw 3 # blue'
  write(10,*) 'set style line 4 lt 4 lw 3 # purple'
  write(10,*) 'set style line 5 lt 5 lw 3 # light blue'
  write(10,*) 'set style line 6 lt 6 lw 3 # yellow'
  write(10,*) 'set style line 7 lt 7 lw 3 # black'
  write(10,*) 'set style line 8 lt 7 lw 5 # thick black'
  write(10,*) ' '

  write(10,*) 'set size square'
  write(10,*) 'set terminal postscript eps solid color enhanced "Helvetica" 24'
  write(10,*) 'set output "fig'//pnum//'.eps"'
  write(10,*) ' '

  write(10,*) 'set xrange [-4:4]'
  write(10,*) 'set xlabel "v_i/v_t"'
  write(10,*) 'plot "pdf.dat" u 1:2 w histeps ls 1 noti,', & 
       '"pdf.dat" u 1:3 w histeps ls 2 noti, ', &
       '"pdf.dat" u 1:4 w histeps ls 3 noti,exp(-x**2.)/sqrt(pi) ls 7 noti'
  write(10,*) ' '

  write(10,*) 'x0=',x0
  write(10,*) 'xmax=',xmax
  write(10,*) 'ymin=',ymin
!  write(10,*) 'ymax=',1.1*x0
  write(10,*) 'Tp=',Tp
  write(10,*) 'set logscale y'
  write(10,*) 'set format y "%2.0t{/Symbol \327}10^{%L}"'
  write(10,*) 'set xrange[0:xmax]'
  write(10,*) 'set yrange[ymin:*]'
  write(10,*) 'set xlabel "E_k(eV)"'
  write(10,*) 'set ylabel "f (eV^{-3/2})"'
  write(10,*) 'plot "pedf.dat" u 1:2 w l ls 1 noti, x0*exp(-x/Tp) ls 7 noti'

  return
end subroutine generate_gnuplot
