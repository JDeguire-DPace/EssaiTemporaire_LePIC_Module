!     ================================================================
!     VERSION:         3.3.6.2
!     LAST MOD:       Sep/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      A 3D explicit parallel hybrid OpenMP/MPI PIC code
!                   
!     NOTE:         Domain extend from ix(iy)=1 to ix(iy)=nx(ny)+1
!                   ix(iy)=1 & ix(iy)=nx(ny)+1 are boundary conditions, 
!                   unknown interior points are located between ix(iy)=2 
!                   and ix(iy)=nx(ny)
!                   Example: nx=5 & along (Ox)
!                   BC                  BC
!                    |---|---|---|---|---|
!                    1   2   3   4   5   6
!                   (0)                (xmax)  dx=xmax/nx
!
!
! Copyright or © or Copr. Gwenael Fubiani (2022/11/22)
!
! gwenael.fubiani@cnrs.fr
!
! This software is a computer program whose purpose is to model low
! temperature plasmas with a Particle-In-Cell algorihm in 3D-3V
! dimensions.
!
! This software is governed by the CeCILL-C license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL-C
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 

! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
!
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-C license and that you accept its terms.
!     ----------------------------------------------------------------
program main
  use omp_lib
  use mpi
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  ! MPI
  integer ierr,mpi_rank,nproc_mpi
  ! Loop indexes
  integer:: it,ix,iy,iz,i,i_rg,ptype,niter,kl,kr,m
  ! Simulation parameters (integers)
  integer:: n(3),n_B(3),ncycle,nmax,nmax_tmp,ngrid,ntype,ntmax,nsav, &
       ng,iproc,nproc,ctime(0:10),MSTIMER,B_dir,namlen,&
       sum_Nh_tmp,sum_Nh,n_neu,flag_sav,flag_wrt,igrid,flag_avg3D,&
       flag_nopart,nseq,nEf(3),i_pl,flag_diag,flag_updatephi
  real(kind=8):: h(3),h_B(3)
  parameter (nmax_tmp=1*10**8)
  ! Physical scales
  real(kind=8):: lbd_d,wp,kt,xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow,&
       wc,ni0(npart),Te,vt,sum_dEk_tmp,sum_dEk_tot,I_inj
  ! Arrays
  integer, allocatable:: bcnd(:,:,:),iseed(:), &
       bcnd_dom(:,:,:),shift(:),length(:),N_inj(:,:),N_flx(:,:)
  integer(kind=1), allocatable:: flag_dead(:,:,:)
  real(kind=8), allocatable:: phi(:,:,:),sum_dEk(:), &
       vxp(:,:,:,:),np(:,:,:,:,:),Ei(:,:,:,:),p_mac(:,:,:,:),&
       P_loss(:,:,:),kq(:,:,:),Bi(:,:,:,:),&
       np_red(:,:,:,:),phi_dom(:,:,:),rhs_dom(:,:,:),cnt_dead(:),beam_div(:,:),&
       avg3D(:,:,:,:)
  ! Sorting
  integer:: pl_max,nsort,sum_cex(3)
  integer, allocatable:: Plist(:,:,:),flag_cex(:,:),cnt_cex(:,:)
  real(kind=8), allocatable:: sorting(:,:)
  ! Macroscopic parameters
  integer, allocatable:: sour_xz(:,:,:,:),sour_xy(:,:,:,:),sour_fx_yz(:,:,:,:),&
       sour_yz(:,:,:,:),np_tot(:,:),Nh(:),sum_np_tot_tmp(:),sum_np_tot(:),dtype(:)
  real(kind=8), allocatable:: cnt_col(:,:,:),p_mts_xy(:,:,:,:,:),&
       p_mts_xz(:,:,:,:,:),p_mts_yz(:,:,:,:,:),phi_avg_xy(:,:),phi_avg_xz(:,:),&
       phi_avg_yz(:,:),data_pavg_xy(:,:,:,:),data_pavg_xz(:,:,:,:),&
       data_pavg_yz(:,:,:,:),Vgrd(:),j_die_xz(:,:,:),j_die_yz(:,:,:),&
       sum_q_xz(:,:,:,:),sum_q_yz(:,:,:,:,:),sum_q_red_xz(:,:,:),&
       sum_q_red_yz(:,:,:,:),ss2D_xz(:,:,:,:,:),ss2D_xy(:,:,:,:,:),&
       ss2D_yz(:,:,:,:,:)
  ! Collisions
  integer:: ncol_mx,npt_mx
  parameter (ncol_mx=65, npt_mx=1300)
  integer:: sig_list(npart,ncol_mx),col_info(ncol_mx,10), &
       sig_type(ncol_mx),scol_rank(npart,ncol_mx)
  real(kind=8):: sig(npt_mx,ncol_mx),sig_Er(npt_mx),sig_Eex(ncol_mx,2), &
       scol_info(ncol_mx,4),sigv_mx(npart,ncol_mx)
  ! Simulation parameters (reals)
  integer:: cnt_i,cnt_f,cnt_rate,navg,n_mts,cnt_avg(3),lgh,iseed_OMP,sav_np(2),&
       np_pos,np_pos0,ns_convP,iyl_pow,iyr_pow,izl_pow,izr_pow,xBm(3)
  real(kind=8):: res,eps,ksor,time,tmax,sum_time,a1,a2,a3,xl_rg(nm_rg),xr_rg(nm_rg),&
       yl_rg,yr_rg,zl_rg,zr_rg,dtime,omega,cnt_dead_tmp,sum_beam_div(4),Ca,x,np_dup,&
       phi0_RF,f0_RF,phi1_RF,f1_RF
  character:: rname*20,name*20,pnum*1,pnum_bck*3,plnum*3,corrnum*3

  CALL MPI_Init(ierr)                             ! starts MPI
  CALL MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)  ! get current process id
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc_mpi, ierr) ! get # of procs

  call system_clock(cnt_i,cnt_rate)

  !
  ! Write code info on screen
  !
  if(mpi_rank.eq.0) call introduction

  !
  ! Constants
  !
  qe=1.60217646e-19 ! Coulombs
  c=2.99792458d8 ! m/s
  eps0=8.854187817d-12 ! S.I.
  pi=4.d0*datan(1.d0)
  amu=1.66053886d-27


  !
  ! Get number of processors
  !
  nproc= omp_get_max_threads()


  if(mpi_rank.eq.0) then
     print'(1x,"Number of OpenMP threads= ",i3)',nproc
     print'(1x,"Number of MPI threads= ",i3)',nproc_mpi
  endif

  !
  ! Read input parameters
  !
  call read_input(n,tmax,xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow,&
       nsav,eps,omega,kt,rname,ngrid,ng, xl_rg,xr_rg,yl_rg,yr_rg,&
       zl_rg,zr_rg,I_inj,Ca,mpi_rank,flag_avg3D,np_dup,n_B,phi0_RF,f0_RF,&
       phi1_RF,f1_RF)

  flag_sav= 0
  if(nsav.lt.0) then
     flag_sav= 1
     nsav= ABS(nsav)
  endif

  flag_bak= 0
  if(nbak.lt.0) then
     flag_bak= 1
     nbak= ABS(nbak)
  endif

  ! Save every n-points
  every=1
  if(MAXVAL(n).gt.512) every=2
  if(MAXVAL(n).gt.1024) every=4
  if(MAXVAL(n).gt.2048) every=8
  
  ! Correct for nmax
  nmax= NINT(nmax_tmp/real(nproc)/real(nproc_mpi))

  if( mpi_rank.eq.0 .and. k_eps0.ne.1 ) &
       print*, 'eps0 HAS BEEN RE-SCALED: k_eps0=',k_eps0
  eps0=k_eps0*eps0
  
  ! Warning
  if(k_eps0.eq.0) then
     print*, 'Please correct k_eps0'
     call stop_calculation
  endif

  !
  ! Read cross sections from files
  !
  call read_reactions(sig,sig_Er,sig_list,sig_Eex,ncol_mx,sig_type, &
       npt_mx,rname,col_info,scol_rank,scol_info,sigv_mx,ni0,ntype,n_neu,mpi_rank)

  ntype= ntype - n_neu
  ! Neutral density
  np_mx(ntype+1:ntype+n_neu)= ni0(ntype+1:ntype+n_neu)

  ! Electron fraction
  if(tag_neg.gt.0) ni0(1)= ni0(1) - ni0(tag_neg)

  !
  ! Allocate arrays
  !
  allocate (  phi(0:n(1)+2,0:n(2)+2,0:n(3)+2),&
       vxp(6,nmax,ntype,nproc),np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc), &
       Ei(3,0:n(1)+2,0:n(2)+2,0:n(3)+2),p_mac(ntype,2,0:ngrid,nproc),P_loss(4,ntype,nproc), &
       Bi(4,0:n_B(1)+2,0:n_B(2)+2,0:n_B(3)+2), &
       np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype), &
       kq(0:n(1)+2,0:n(2)+2,0:n(3)+2), &
       sum_dEk(nproc) )
  allocate (  rhs_dom(0:n(1)+1,0:n(2)+1,0:n(3)/nproc_mpi+1), &
       phi_dom(0:n(1)+2,0:n(2)+2,-1:n(3)/nproc_mpi+2), &
       bcnd_dom(0:n(1)+2,0:n(2)+2,0:n(3)/nproc_mpi+2), &
       shift(0:nproc_mpi-1),length(0:nproc_mpi-1), &
       N_inj(ntype,nproc),N_flx(ntype,nproc) )
  allocate ( np_tot(ntype,nproc),sour_xy(0:n(1)+2,0:n(2)+2,ntype,nproc), &
       sour_xz(0:n(1)+2,0:n(3)+2,ntype,nproc),sour_fx_yz(0:n(2)+2,0:n(3)+2,ntype,nproc),&
       sour_yz(0:n(2)+2,0:n(3)+2,ntype,nproc),Nh(nproc),sum_np_tot_tmp(ntype),sum_np_tot(ntype) )
  allocate ( cnt_col(ncol_mx,nproc,nm_rg),iseed(nproc) )
  allocate ( bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2) )
  allocate ( flag_dead(nmax,ntype,nproc), flag_cex(nmax,nproc), cnt_cex(3,nproc), &
       cnt_dead(nproc),beam_div(4,nproc),j_die_xz(2,0:n(1)+2,0:n(3)+2),j_die_yz(2,0:n(2)+2,0:n(3)+2),&
       sum_q_xz(0:n(1)+2,0:n(3)+2,ntype,nproc),sum_q_yz(2,0:n(2)+2,0:n(3)+2,ntype,nproc),&
       sum_q_red_xz(ntype,0:n(1)+2,0:n(3)+2),sum_q_red_yz(2,ntype,0:n(2)+2,0:n(3)+2))
  allocate ( phi_avg_xy(0:n(1)+2,0:n(2)+2), &
       phi_avg_xz(0:n(1)+2,0:n(3)+2), &
       phi_avg_yz(0:n(2)+2,0:n(3)+2), &
       data_pavg_xy(8,0:n(1)+2,0:n(2)+2,ntype),&
       data_pavg_xz(8,0:n(1)+2,0:n(3)+2,ntype),&
       data_pavg_yz(8,0:n(2)+2,0:n(3)+2,ntype),&
       Vgrd(ngrid),dtype(0:ngrid),&
       ss2D_xy(2,0:n(1)+2,0:n(2)+2,ntype,nproc),&
       ss2D_xz(2,0:n(1)+2,0:n(3)+2,ntype,nproc),&
       ss2D_yz(2,0:n(2)+2,0:n(3)+2,ntype,nproc) )

  pl_max= n(1) + n(1)*(n(2) + n(2)*n(3)) + 1
  allocate ( Plist(0:pl_max,ntype,nproc) )

  nEf=0
  if(flag_avg3D.eq.1) nEf= n
  allocate( avg3D(5,1:nEf(1)+1,1:nEf(2)+1,1:nEf(3)+1) )
  
  !
  ! Internal code parameters
  !
  ntmax= 5000001
  ncycle= 20000 ! PDE max cycle

  ! Calculate maximum number multigrid levels
  a1= alog10(real(n(1)))/alog10(2.)
  a2= alog10(real(n(2)))/alog10(2.)
  a3= alog10(real(n(3)))/alog10(2.)
  if( ABS(NINT(a1)-a1).le.1d-3 .and. ABS(NINT(a2)-a2).le.1d-3 .and. &
       ABS(NINT(a3)-a3).le.1d-3 ) then
     ng= NINT(MIN(a1,a2,a3) )
  else
     a1= real(n(1))/2**(ng-1)
     a2= real(n(2))/2**(ng-1)
     a3= real(n(3))/2**(ng-1)
     ! Warning
     if( ABS(NINT(a1)-a1).gt.0 .or. ABS(NINT(a2)-a2).gt.0 .or. &
          ABS(NINT(a3)-a3).gt.0 ) then
        print*, 'Incorrect number of grid levels in MG'
        print*, 'Please correct ...'
        call stop_calculation
     endif
  endif

  !
  ! Initialize variables
  !
  cnt_col=0.d0
  time=0.d0
  Ei=0.d0
  phi=0.d0
  sour_xy=0
  sour_xz=0
  sour_fx_yz=0
  sour_yz=0
  ctime=0
  flag_dead=0
  Bmax=0.d0
  cnt_avg=0
  p_mac=0.d0
  P_loss=0.d0
  data_pavg_xy=0.d0
  data_pavg_xz=0.d0
  data_pavg_yz=0.d0
  phi_avg_xy=0.d0
  phi_avg_xz=0.d0
  phi_avg_yz=0.d0
  flag_cex=0
  cnt_cex=0
  cnt_dead=0.d0
  beam_div=0.d0
  N_inj=0
  flag_wrt=0
  sum_q_xz= 0.d0
  sum_q_yz= 0.d0
  N_flx=0
  j_die_xz= 0.d0
  j_die_yz= 0.d0
  dtype= 0
  cnt_plt= 0
  ss2D_xy=0.d0
  ss2D_xz=0.d0
  ss2D_yz=0.d0
  Bi=0.d0
  
  !
  ! Generate boundary counditions
  !
  call generate_boundary(phi,n,h,bcnd,Vgrd,ngrid,dtype,xl_pow,xr_pow,mpi_rank)

  ! Domain decomposition for Poisson solver
  m= n(3)/nproc_mpi
  phi_dom(:,:,0:m+2)= phi(:,:,mpi_rank*m:(mpi_rank+1)*m+2) 
  bcnd_dom(:,:,0:m+2)= bcnd(:,:,mpi_rank*m:(mpi_rank+1)*m+2) 

  if(flag_circxh.eq.1) R_ahp= yr_pow-ymax/2.d0 ! Radius of Antenna heating profile
  
  !
  ! Coordinates of sub-regions
  !
  do i=1,n_rg
     xl_rg(i)= MAX(0.d0,xl_rg(i))
     xr_rg(i)= MIN(xmax,xr_rg(i))
     ixl_rg(i)= NINT(xl_rg(i)/h(1)) + 1
     ixr_rg(i)= NINT(xr_rg(i)/h(1)) + 1
  enddo
  yl_rg= MAX(0.d0,yl_rg)
  yr_rg= MIN(ymax,yr_rg)
  zl_rg= MAX(0.d0,zl_rg)
  zr_rg= MIN(zmax,zr_rg)
  iyl_rg= NINT(yl_rg/h(2)) + 1
  iyr_rg= NINT(yr_rg/h(2)) + 1  
  izl_rg= NINT(zl_rg/h(3)) + 1
  izr_rg= NINT(zr_rg/h(3)) + 1  

  !
  ! Define charge per cell coefficient
  !
  kq= 1.d0
  !$OMP PARALLEL
  !$OMP DO
  do ix=0,n(1)+2
     do iy=0,n(2)+2
        do iz=0,n(3)+2
           if(bcnd(ix,iy,iz).ne.-1) kq(ix,iy,iz)= 2.d0
        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL  

  !
  ! Open (some) output files
  !
  open(12,file='DATA/time.dat')
  if(tag_neg.gt.0) then 
     open(13,file='DATA/neg_ion_impacts.dat')
     write(13,'(A120)') '# time(s), neg. ions generated on the PE surface (%), CEX(%), volume production(%), &
          YP and ZP RMS beam divergence (rad)'
  endif
  open(40,file='DATA/residual.dat')

  !
  ! Magnetic field profile
  !
  if(flag_B.eq.1) then
     h_B(1)= xmax/n_B(1)
     h_B(2)= ymax/n_B(2)
     h_B(3)= zmax/n_B(3)
  
     flag_gridB=0
     if( n_B(1).ne.n(1) .or. n_B(2).ne.n(2) .or. &
          n_B(3).ne.n(3)) flag_gridB=1
  
     do i=1,nB

        ! Skip
        if(B_file(i).eq.'s'.or.B_file(i).eq.'S') goto 10
              
        ! Extract file name
        name= B_name(i)
        do namlen=20,1,-1
           if( name(namlen:namlen).ne.' ' ) goto 20
        enddo
20      continue

        ! B-field map from an external file
        if( B_file(i).eq.'y'.or. B_file(i).eq.'Y' ) then      
           do B_dir=1,3
              if(B_dir.eq.1) name= name(1:namlen)//'Bx.dat'
              if(B_dir.eq.2) name= name(1:namlen)//'By.dat'
              if(B_dir.eq.3) name= name(1:namlen)//'Bz.dat'
              call read_Bfield_map(Bi,n_B,name,namlen+6,B_dir,B_scale(i),mpi_rank)
           enddo
        endif

        ! B-field map from gaussian profile/PG current
        if( B_file(i).eq.'n'.or. B_file(i).eq.'N' ) then 
           call find_B_dir(B_info,B_dir)
           if( name(1:1).eq. 'g' .or. name(1:1).eq. 'G' ) then
              if(mpi_rank.eq.0) print*, 'Gaussian magnetic filter field profile' 
              call gaussian_Bfield(Bi,n_B,h_B,B0(i),x0(i),dL(i),B_dir)
           endif
           if( name(1:1).eq. 'P' .or. name(1:1).eq. 'p' ) then
              if(mpi_rank.eq.0) print*, 'Magnetic filter field generated by a PG current' 
              call PG_current(Bi,n_B,h_B,B0(i),B_dir)
           endif
           if( name(1:1).eq. 'e' .or. name(1:1).eq. 'E' ) then
              if(mpi_rank.eq.0) print*, 'Cusp magnetic field profile is generated' 
              call EE_magnets(Bi,n_B,h_B,B0(i),x0(i),z0(i),dL(i))
           endif
           flag_thr=0
           if( name(1:1).eq. 't' .or. name(1:1).eq. 'T' ) then
              if(mpi_rank.eq.0) print*, 'Magnetic field profile of the Hall thruster' 
              call Bfield_Hall_thruster(Bi,n_B,h_B,B0(i),x0(i),dL(i))
              flag_thr=1
              y_thr= MIN(ymax, 8.d-3) ! Inject cathode e-flux +/- 4 mm radially
              if(B0(i).lt.0.d0) y_thr= ymax
           endif
        endif
        
10      continue
     enddo
 
     ! Calculate ||B||
     !$OMP PARALLEL REDUCTION(MAX:Bmax)
     !$OMP DO
     do iz=0,n_B(3)+2
        do iy=0,n_B(2)+2
           do ix=0,n_B(1)+2
              Bi(4,ix,iy,iz)= dsqrt( Bi(1,ix,iy,iz)*Bi(1,ix,iy,iz) + &
                   Bi(2,ix,iy,iz)*Bi(2,ix,iy,iz) + &
                   Bi(3,ix,iy,iz)*Bi(3,ix,iy,iz) )
              Bmax= MAX(Bmax,Bi(4,ix,iy,iz))                               
           enddo
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL

     xBm= MAXLOC(Bi(4,:,:,:))
     
     ! Write on files
     if( n_B(1).gt.1 .and. mpi_rank.eq.0 ) then     

        print*, 'Bmax @ location:',xBm(1:3)-1
        
        ! XY plane
        open(30,file='DATA/DATA_2D/B_xy.mco')
        write(30,*) n_B(1)/every,n_B(2)/every
        do iy=n_B(2)+1,1,-1*every ! in Gauss
           write(30,99) ( Bi(4,ix,iy,n_B(3)/2+1)*1.d4, ix=1,n_B(1)+1,every )
99         format(800(e18.6,1x)) 
        enddo
        write(30,*) 'vector'
        do iy=n_B(2)+1,1,-1*every ! Theta= arctg(By/Bx)
           write(30,99) ( datan2( Bi(2,ix,iy,n_B(3)/2+1),&
                Bi(1,ix,iy,n_B(3)/2+1) ), ix=1,n_B(1)+1,every )  
        enddo
        close(30)

        open(30,file='DATA/DATA_2D/By_xy.mco')
        write(30,*) n_B(1)/every,n_B(2)/every
        do iy=n_B(2)+1,1,-1*every ! in Gauss
           write(30,99) ( Bi(2,ix,iy,n_B(3)/2+1)*1.d4, ix=1,n_B(1)+1,every )
        enddo
        close(30)

        open(30,file='DATA/DATA_2D/Bx_xy.mco')
        write(30,*) n_B(1)/every,n_B(2)/every
        do iy=n_B(2)+1,1,-1*every ! in Gauss
           write(30,99) ( Bi(1,ix,iy,n_B(3)/2+1)*1.d4, ix=1,n_B(1)+1,every )
        enddo
        close(30)
        
        ! XZ plane
        open(30,file='DATA/DATA_2D/B_xz.mco')
        write(30,*) n_B(1)/every,n_B(3)/every
        do iz=n_B(3)+1,1,-1*every ! in Gauss
           write(30,99) ( Bi(4,ix,n_B(2)/2+1,iz)*1.d4, ix=1,n_B(1)+1,every )
        enddo
        write(30,*) 'vector'
        do iz=n_B(3)+1,1,-1*every ! Theta= arctg(Bz/Bx)
           write(30,99) ( datan2( Bi(3,ix,n_B(2)/2+1,iz),&
                Bi(1,ix,n_B(2)/2+1,iz) ), ix=1,n_B(1)+1,every )  
        enddo
        close(30)

        open(30,file='DATA/DATA_2D/Bx_xz.mco')
        write(30,*) n_B(1)/every,n_B(3)/every
        do iz=n_B(3)+1,1,-1*every ! in Gauss
           write(30,99) ( Bi(1,ix,n_B(2)/2+1,iz)*1.d4, ix=1,n_B(1)+1,every )
        enddo
        close(30)

        open(30,file='DATA/DATA_2D/Bz_xz.mco')
        write(30,*) n_B(1)/every,n_B(3)/every
        do iz=n_B(3)+1,1,-1*every ! in Gauss
           write(30,99) ( Bi(3,ix,n_B(2)/2+1,iz)*1.d4, ix=1,n_B(1)+1,every )
        enddo
        close(30)

        ! YZ plane
        ix_pl= n_B(1)/2+1
        open(30,file='DATA/DATA_2D/B_yz.mco')
        write(30,*) n_B(2)/every,n_B(3)/every
        do iz=n_B(3)+1,1,-1*every ! in Gauss
           write(30,99) ( Bi(4,ix_pl,iy,iz)*1.d4, iy=1,n_B(2)+1,every )
        enddo
        write(30,*) 'vector'
        do iz=n_B(3)+1,1,-1*every ! Theta= arctg(Bz/By)
           write(30,99) ( datan2( Bi(3,ix_pl,iy,iz),&
                Bi(2,ix_pl,iy,iz) ), iy=1,n_B(2)+1,every )  
        enddo
        close(30)

        open(30,file='DATA/DATA_2D/By_yz.mco')
        write(30,*) n_B(2)/every,n_B(3)/every
        do iz=n_B(3)+1,1,-1*every ! in Gauss
           write(30,99) ( Bi(2,ix_pl,iy,iz)*1.d4, iy=1,n_B(2)+1,every )
        enddo
        close(30)

        open(30,file='DATA/DATA_2D/Bz_yz.mco')
        write(30,*) n_B(2)/every,n_B(3)/every
        do iz=n_B(3)+1,1,-1*every ! in Gauss
           write(30,99) ( Bi(3,ix_pl,iy,iz)*1.d4, iy=1,n_B(2)+1,every )
        enddo
        close(30)

     endif
         
  endif

  !
  ! Simulation parameters
  !
  lbd_d=dsqrt( (eps0*Ti(1))/(n0*ABS(charge(1))) ) ! Debye length
  wp=dsqrt( n0*charge(1)**2./(eps0*mass(1)) ) ! Plasma frequency

  do ptype=1,ntype+n_neu
     ! Thermal velocities [sqrt(2*kT/m)]
     vt0(ptype)=dsqrt(2.d0*qe*Ti(ptype)/ABS(mass(ptype)))
     ! Charged macroparticle weight
     if(ABS(charge(ptype)).gt.0) then
        Nm(ptype)=n0*h(1)*h(2)*h(3)/ABS(np_cell)
     else ! Gas macroparticle weight
        Nm(ptype)=ngas*h(1)*h(2)*h(3)/ABS(np_cell)
     endif
  enddo
  dt= kt*MIN(h(1),h(2),h(3))/vt0(1) ! Time step k*dr/vte, k is arbitrary

  nsort=10 ! Frequency of calls to sorting subroutine
  ns_coll= 1*nsort ! Frequency of calls to collision subroutine
  ns_heat=4 ! Frequency of calls to heating subroutine
  if( kt.ge.0.05 .and. kt.lt.0.1 ) ns_heat=20
  if(kt.lt.0.05) ns_heat=40
  if(flag_bak.eq.0) then ! Frequency for averaging
     navg=5*nsort
  else
     navg= nsav
  endif
  if(tseq.gt.0) then ! Save sequence of profiles every tseq seconds
     nseq= MAX(NINT((tseq/dt)/nsav),1)*nsav 
     tseq= nseq*dt ! correct tseq
  endif
  nudt=nu_h*(ns_heat*dt)  ! % of electrons which should be heated
  nu_uplim(1)=5.d8 ! Maximum collision frequency (e-)
  nu_uplim(2:ntype)=1.d7 ! (ions)
  ns_flx=1 ! Negative ion injection frequency on the PE
  ns_inj=1 ! Frequency for particle injection
  ix_pl= n(1)/2+1 ! 2D plots in YZ plane
  iz_pl= n(3)/2+1 ! 2D plots in XY plane
  if(flag_grd.eq.1) ix_pl= ixg ! Grid location
  if(flag_die.eq.1) Ca= Ca*h(1)*h(3) ! Farads per grid cell
  if(flag_thr.eq.1) then 
     x_thr= xmax - 1.d-3 ! Inject cathode e-flux 1 mm from RHS boundary
     ix_thr= INT( x_thr/h(1) ) + 1
     ix_pl= INT( 0.5d0*(xr_pow+xl_pow)/h(1) ) + 1 ! Middle of the ion injection area
  endif
  plt_src=1 ! Flag for plotting source/sink terms (0: count macroparticle production, 1: n*nu is evaluated)
  if((Pabs.lt.0 .and. ABS(opt_inj).lt.3) .or. flag_avg3D.eq.0) plt_src=0

  if(mpi_rank.eq.0) then
     print'(1x,"Frequency of calls to collision subroutine= ",i3)',ns_coll
     print'(1x,"Frequency of calls to electron (Maxwellian) heating subroutine= ",i3)',ns_heat
     print'(1x,"Frequency of calls to sort subroutine= ",i3)',nsort
     if(navg.eq.nsav) then
        print'(1x,"Data will not be averaged!")'
     else
        print'(1x,"Frequency for averaging= ",i5)',navg
     endif

     ! Warning
     if(nudt.gt.1.d0) then
        print*, 'Warning, nu*dt>1, please correct ...'
        print*, 'nu*dt=',nudt
        call stop_calculation
     endif
  endif

  !
  ! Set random seed number for each processor
  !
  do iproc=1,nproc
     iseed(iproc)= 123456*iproc*(10*mpi_rank+1)
  enddo
  
  !
  ! Load particles
  !

  ! Correct for heating boundaries out of range
  if(xr_pow.gt.xmax) xr_pow=xmax
  if(xl_pow.lt.0.d0) xl_pow=0.d0
  ixr_pow= INT( xr_pow/h(1) ) + 1
  ixl_pow= INT( xl_pow/h(1) ) + 1
  if(yr_pow.gt.ymax) yr_pow=ymax
  if(yl_pow.lt.0.d0) yl_pow=0.d0
  iyr_pow= INT( yr_pow/h(2) ) + 1
  iyl_pow= INT( yl_pow/h(2) ) + 1
  if(zr_pow.gt.zmax) zr_pow=zmax
  if(zl_pow.lt.0.d0) zl_pow=0.d0
  izr_pow= INT( zr_pow/h(3) ) + 1
  izl_pow= INT( zl_pow/h(3) ) + 1

  ! Location of the plasma electrode
  if(x_load.gt.xmax) x_load= xmax
  ix_PE= INT( x_load/h(1) ) + 1

  if(mpi_rank.eq.0) then
     print'(1x,"Heating region: ix<=",i3,", ix>=",i5)',ixl_pow,ixr_pow
     print'(1x,"Injection region: ix<=",i3,", ix>=",i5,", iy<=",i3,", iy>=",i5,"&
          , iz<=",i3,", iz>=",i5)',ixl_pow,ixr_pow,iyl_pow,iyr_pow,izl_pow,izr_pow
     print'(1x,"Particle loading region: nx<= ",i3)',ix_PE
  endif

  ! Estimate total number of cells
  n_cell= 0
  do iz=1,n(3)
     do iy=1,n(2)
        do ix=1,ix_PE
           ! calculate # of cells inside simulation domain
           if( bcnd(ix,iy,iz).eq.-1 .or. bcnd(ix+1,iy,iz).eq.-1 .or. &
                bcnd(ix+1,iy+1,iz).eq.-1 .or. bcnd(ix,iy+1,iz).eq.-1 .or. &
                bcnd(ix,iy,iz+1).eq.-1 .or. bcnd(ix+1,iy,iz+1).eq.-1 .or. &
                bcnd(ix+1,iy+1,iz+1).eq.-1 .or. bcnd(ix,iy+1,iz+1).eq.-1 ) &
                n_cell= n_cell + 1           
        enddo
     enddo
  enddo
  
  ! Restart using .bak files
  if(flag_restart.eq.0) then
     if(np_cell.gt.0) then
        call load_part(n,h,bcnd,np,vxp,ntype,nmax,kq,ni0,&
             np_tot,nproc,iseed,sum_dEk,Nh,mpi_rank,nproc_mpi)
     else
        np= 0.d0
        np_tot(1:ntype,:)= 0
     endif
  else
     call restart(n,h,np,vxp,nmax,ntype,kq,time,nproc,iseed,&
          np_tot,sum_dEk,Nh,flag_cex,mpi_rank,nproc_mpi,np_dup)
     ! Load potential on dielectric surfaces
     if( flag_die.eq.1 .and. Ca.gt.0 ) then
        if(flag_restart.eq.2) then
           open(41,file='SAV_DATA/phi.bak',form='UNFORMATTED')
        else
           open(41,file='DATA.BAK/phi.bak',form='UNFORMATTED')
        endif
        read(41) phi(0:n(1)+2,0:n(2)+2,0:n(3)+2)
        close(41)
        phi_dom(:,:,0:m+2)= phi(:,:,mpi_rank*m:(mpi_rank+1)*m+2) 
     endif
     Ca= ABS(Ca)
     if( tseq.gt.0.d0 .and. time.ge.tseq_init .and. time.le.tseq_final ) &
          cnt_plt= NINT((time-tseq_init)/tseq)
  endif

  !
  ! Print info. on screen
  ! 
  if(mpi_rank.eq.0) then
     print*, 'lde (mm)=',lbd_d*1.d3
     print*, 'lpe (mm)=',2*pi*c/wp*1.d3
     print*, 'vt (m/s)=',vt0(1:ntype) ! Electron thermal velocity
     print*, '<|v|> (m/s)=',2.d0/dsqrt(pi)*vt0(1) ! Electron electron speed
     print*, 'dx (mm)=',h(1)*1.d3,'dy (mm)=',h(2)*1.d3,'dz (mm)=',h(3)*1.d3
     print*, 'dx/lde=',h(1)/lbd_d,'dy/lde=',h(2)/lbd_d,'dz/lde=',h(3)/lbd_d
     print*, 'Nm=',Nm(1)
     print*, 'dt (s)=',dt,', wpe*dt=',wp*dt
     print*, 'vte*dt/dx=',vt0(1)*dt/MIN(h(1),h(2),h(3))
     
     if(flag_nmn.eq.1) print'(1x,"Left-hand-side boundary condition set to Neumann type")'
     if(flag_thr.eq.1) then
        print'(1x,"Thruster simulation: the anode corresponds to grid #",i2)',ind_g        
        print'(1x,"Electron flux from the cathode injected at ix=",i5)',ix_thr
     endif

     print'(1x,"Figures in the YZ plane are drawn at the location ix=",i5)',ix_pl
     print'(1x,"Figures in the XY plane are drawn at the location iz=",i5)',iz_pl

     if(Bmax.gt.0.d0) then
        if(n_B(1).eq.1) then
           print'(1x,"B-field is constant")'
           print'(1x,"Bx(T)= ",es10.2,", By(T)= ",es10.2,", Bz(T)= ",es10.2)', Bi(1,1,1,1), Bi(2,1,1,1), Bi(3,1,1,1)
        endif
     
        wc=qe*Bmax/mass(1)
        a1= 2.*datan(wc*dt/2.)/(wc*dt)
        a2= dsqrt(1+(wc*dt/2.)**2)
        print'(1x,"Bm(T)= ",es10.2,", wc*dt= ",f6.2,", wc (/s)= ",es12.4,", & 
             re (mm)= ",f5.2)', Bmax,wc*dt,wc,vt0(1)/wc*1.d3
        print'(1x,"err: w*/wc= ",f5.2,", r*/re= ",f5.2)', a1,a2
     endif
     
     print'(1x,"Multigrid PDE Solver is used, eps= ",es10.2)',eps
     print'(1x,"ng=",i2,", omega=",f6.1)',ng,omega
     print'(1x,"nx=",i4,", ny=",i4,", nz=",i4,", # of cells: ",i8)',n(1),n(2),n(3),n_cell
     if(flag_avg3D.eq.1) print*, 'Save 3D density, potential and E-field maps option on'

     if(opt_inj.lt.0 .and. I_inj.gt.0) print'(1x,"Particle injection profile is a Cosine")'
     if(ABS(opt_inj).eq.3 .and. I_inj.gt.0) print'(1x,"Electron beam injection along +Z direction")'
     if( (ABS(opt_inj).eq.4 .and. I_inj.gt.0 ) .or. gam_sec.gt.0.d0 ) then
        igrid= 1 ! Z=cste plane and @ bottom 
        if(dir_sec.eq.-1) igrid=2 ! Top        
        if(n_cath.eq.1) then
           print'(1x,"Emissive cathode along (Oz) at zg(cm)=",f6.2)',zg_sec(igrid)*1.d2
        else
           print'(1x,"Emissive cathode along (Oz) at zg(cm)=",2(f6.2,1x))',zg_sec(1:n_cath)*1.d2
        endif
     endif
     
     if(flag_sav.eq.0 .and. tseq.gt.0 .and. time.le.tseq_final) &
          print'(1x,"Frequency for saving sequence of profiles= ",i6,", dt(us)=",f5.2,", cnt_plt=",i4)', &
          nseq,tseq*1.d6,cnt_plt
     if(plt_src.eq.1) print'(1x,"Source/sink plotting option on")'
     if(jne.gt.0.d0) &
          print'(1x,"Negative ion flux emmited on grid located at xg=",f6.2,1x,"cm, ig=",i4)',xg1*1.d2,ind_g       
  endif

  ! Iterative process to find Pabs or I_inj : initialization
  if(flag_convP.eq.1) then
     ns_convP=10000
     ! Save total number of targeted positive ions in calculation
     np_pos0= NINT(n0*h(1)*h(2)*h(3)*real(n_cell)/Nm(1)) ! np_cell*n_cell
     if(mpi_rank.eq.0) then
        print*, 'Pabs, I_inj or Vgrd will be calculated through an iterative process'
        print*, 'Total number of targeted positive ions in calculation: ',np_pos0
     endif
  endif

  ! Reset timer        
  ctime(0)= MSTIMER()
  ctime=0

  call system_clock(cnt_f)
  dtime=REAL(cnt_f-cnt_i)/REAL(cnt_rate)
  if(mpi_rank.eq.0) print'(" *** Startup time (s): ",f5.1," ***")', dtime

  if(mpi_rank.eq.0) then
     print*, ' '
     print*, 'Running simulation ...'
  endif

  !
  ! Start iteration 
  !
  do it=1,ntmax

     time= time + dt
     if(time.ge.tmax) exit ! Stop calculation

     ctime(0)= MSTIMER()

     ! Stat without any background plasma
     flag_nopart= 0
     if(SUM(np_tot(1:ntype,1:nproc)).eq.0) flag_nopart= 1

     flag_updatephi= 0
     ! Iterate to find Pabs, I_inj or V
     if( flag_convP.eq.1 .and. MOD(it,ns_convP).eq.0 ) then
        np_pos= SUM(np_tot(2:ntype,1:nproc))
        if(tag_neg.gt.0) np_pos= np_pos - SUM(np_tot(tag_neg,1:nproc))
        if(nproc_mpi.gt.1) then
           call MPI_ALLREDUCE(MPI_IN_PLACE, np_pos,1, &
                MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
        endif
        
        if( ABS(opt_inj).ne.2 .and. I_inj.gt.0.d0 ) then
           ! Inject particle beam
           I_inj=I_inj*real(np_pos0)/real(np_pos)
           if(nproc_mpi.gt.1) call MPI_Bcast(I_inj,1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        endif
        
        ! Inject external power
        if(Pabs.gt.0.d0) then
           Pabs=Pabs*real(np_pos0)/real(np_pos)
           if(nproc_mpi.gt.1) call MPI_Bcast(Pabs,1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)     
        endif
        
        ! Secondary electron emission
        if( Pabs.le.0.d0 .and. I_inj.le.0.d0 .and. gam_sec.gt.0.d0 ) then
           ! Update flag
           flag_updatephi= 1
           ! Iterative Vgrd
           Vgrd(igrid_sec)=Vgrd(igrid_sec)*real(np_pos0)/real(np_pos)       
        endif
     endif

     if(flag_RFpot.eq.1) then
        ! Warning
        if(flag_updatephi.eq.1) then
           print*, 'Cathode potential update used twice : for iterative proceduce and RF. Please correct...'
           call stop_calculation
        endif
        ! Update flag
        flag_updatephi=1
        ! Calculate RF potential
        Vgrd(igrid_sec)= phi0_RF*dsin(2.d0*pi*f0_RF*time) + &
             phi1_RF*dsin(2.d0*pi*f1_RF*time)
     endif

     if(flag_updatephi.eq.1) then
        ! Update cathode potential
        if(nproc_mpi.gt.1) call MPI_Bcast(Vgrd(igrid_sec),1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        !$OMP PARALLEL PRIVATE(igrid)
        !$OMP DO
        do iz=0,m+2
           do iy=0,n(2)+2
              do ix=0,n(1)+2
                 igrid= bcnd_dom(ix,iy,iz)
                 if(igrid.eq.igrid_sec) &
                      phi_dom(ix,iy,iz)= Vgrd(igrid_sec)                
              enddo
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
     endif
        
     !
     ! Calculate rho
     !
     call calc_rho(n,np,rhs_dom,bcnd,ntype,nproc,nproc_mpi,mpi_rank)
     ctime(1)= ctime(1) + MSTIMER()

     !
     ! Dielectric boundary conditions
     !
     if(flag_die.eq.1) then
     
        ! Reduction
        sum_q_red_xz= 0.d0
        sum_q_red_yz= 0.d0
        !$OMP PARALLEL
        !$OMP DO
        do iz=0,n(3)+2 ! XZ planes
           do ix=0,n(1)+2
              do ptype=1,ntype
                 do iproc=1,nproc
                    sum_q_red_xz(ptype,ix,iz)= sum_q_red_xz(ptype,ix,iz) + sum_q_xz(ix,iz,ptype,iproc)
                 enddo
              enddo
           enddo
        enddo
        !$OMP END DO NOWAIT

        !$OMP DO
        do iz=0,n(3)+2 ! YZ planes
           do iy=0,n(2)+2
              do ptype=1,ntype
                 do iproc=1,nproc
                    sum_q_red_yz(:,ptype,iy,iz)= sum_q_red_yz(:,ptype,iy,iz) + sum_q_yz(:,iy,iz,ptype,iproc)
                 enddo
              enddo
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL  
       
        call MPI_ALLREDUCE(MPI_IN_PLACE, sum_q_red_xz, (ntype*(n(1)+3)*(n(3)+3)), &
             MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        call MPI_ALLREDUCE(MPI_IN_PLACE, sum_q_red_yz, (2*ntype*(n(2)+3)*(n(3)+3)), &
             MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! Periodic BC's along (OZ) 
        if(flag_pbc.eq.1) then
           do ptype= 1,ntype
              ! The cell collects half the charge at the boundaries
              sum_q_red_xz(ptype,:,1)= 2.d0*sum_q_red_xz(ptype,:,1)
              sum_q_red_xz(ptype,:,n(3)+1)= 2.d0*sum_q_red_xz(ptype,:,n(3)+1)
              sum_q_red_xz(ptype,:,0)= 2.d0*sum_q_red_xz(ptype,:,n(3))
              sum_q_red_xz(ptype,:,n(3)+2)= 2.d0*sum_q_red_xz(ptype,:,2)

              sum_q_red_yz(:,ptype,:,1)= 2.d0*sum_q_red_yz(:,ptype,:,1)
              sum_q_red_yz(:,ptype,:,n(3)+1)= 2.d0*sum_q_red_yz(:,ptype,:,n(3)+1)
              sum_q_red_yz(:,ptype,:,0)= 2.d0*sum_q_red_yz(:,ptype,:,n(3))
              sum_q_red_yz(:,ptype,:,n(3)+2)= 2.d0*sum_q_red_yz(:,ptype,:,2)
           enddo
        endif

        ! Re-initialize surface charge counters on dielectrics
        sum_q_xz=0.d0
        sum_q_yz=0.d0

        ! Update potential at boundary
        !$OMP PARALLEL PRIVATE(igrid)
        !$OMP DO
        do iz=0,m+2
           do iy=0,n(2)+2
              do ix=0,n(1)+2
                 do ptype= 1,ntype
                    ! dV/dt= dQ/dt/C : V^k+1= V^k + (Q^k+1 - Q^k)/C (1)
                    ! V^1= V_0 + Q^1/C
                    ! V^2 = V^1 + (Q^2-Q^1)/C
                    !    = V_0 + Q^2/C
                    ! ...
                    ! V^n = V_0 + Q^n/C (2) : Q^n is the accumulated charge on the surface at time k= n
                    ! We chose to implement (1) with no loss of generality.
                    igrid= bcnd_dom(ix,iy,iz)
                    if(igrid.gt.0) then
                       if(dtype(igrid).eq.2) & ! XZ planes
                            phi_dom(ix,iy,iz)= phi_dom(ix,iy,iz) +  & ! x2 because of posited symmetry
                            sum_q_red_xz(ptype,ix,iz+mpi_rank*m)/(2.d0*Ca)
                       if(dtype(igrid).eq.3) & ! YZ plane #1
                            phi_dom(ix,iy,iz)= phi_dom(ix,iy,iz) +  &
                            sum_q_red_yz(1,ptype,iy,iz+mpi_rank*m)/Ca
                       if(dtype(igrid).eq.4) & ! YZ plane #2
                            phi_dom(ix,iy,iz)= phi_dom(ix,iy,iz) +  &
                            sum_q_red_yz(2,ptype,iy,iz+mpi_rank*m)/Ca
                    endif
                 enddo
              enddo
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL  
        
     endif
     ctime(1)= ctime(1) + MSTIMER()

     !
     ! Calculate Potential
     !   
     call pdesolver(phi_dom,rhs_dom,bcnd_dom,h,n,ncycle,eps,omega,niter,&
          ksor,res,ng,mpi_rank,nproc_mpi)

     ! Concatenate *_dom()
     if(it.eq.1) then
        shift(0)=0
        do i=0,nproc_mpi-1
           length(i)=m
           if(i.ge.1) shift(i)=i*m+1
        enddo
        length(0)=length(0)+1
        length(nproc_mpi-1)=length(nproc_mpi-1)+2
        
        length= length*(n(1)+3)*(n(2)+3)
        shift= shift*(n(1)+3)*(n(2)+3)
        
        kl= 1
        if(mpi_rank.eq.0) kl=kl-1
        kr= m
        if(mpi_rank.eq.nproc_mpi-1) kr=kr+1
     endif

     call MPI_ALLGATHERV(phi_dom(0:n(1)+2,0:n(2)+2,kl:kr),length(mpi_rank),MPI_REAL8,phi, &
          length,shift,MPI_REAL8,MPI_COMM_WORLD,ierr)

     ! Plot dielectric potential & currents
     if(flag_die.eq.1) then  

        ! Electron flux profiles on dielectric surfaces
        if(ig_die(2).gt.0) &
             j_die_xz(1,:,:)= j_die_xz(1,:,:) + sum_q_red_xz(1,:,:)/(dt*2.d0*h(1)*h(3)*charge(1)) ! XZ planes
        if(ig_die(3).gt.0) &
             j_die_yz(1,:,:)= j_die_yz(1,:,:) + sum_q_red_yz(1,1,:,:)/(dt*h(2)*h(3)*charge(1)) ! YZ plane #1
        ! Ion current profiles
        do ptype=2,ntype
           if(ig_die(2).gt.0) &
                j_die_xz(2,:,:)= j_die_xz(2,:,:) + sum_q_red_xz(ptype,:,:)/(dt*2.d0*h(1)*h(3)*charge(ptype))
           if(ig_die(3).gt.0) &
                j_die_yz(2,:,:)= j_die_yz(2,:,:) + sum_q_red_yz(1,ptype,:,:)/(dt*h(2)*h(3)*charge(ptype))
        enddo

        ! Save in Macho format
        if( mpi_rank.eq.0 .and. MOD(it,nsav).eq.1) then 
           
           
           if(ig_die(2).gt.0) then
              open(14,file='DATA/DATA_2D/phi_die_xz.mco')
              write(14,*) n(1)/every,n(3)/every
              do iz=n(3)+1,1,-1*every 
                 write(14,'(800(e18.6,1x))') ( phi(ix,ig_die(2),iz), ix=1,n(1)+1,every )
              enddo
              close(14)
           endif
           
           if(ig_die(3).gt.0) then
              open(14,file='DATA/DATA_2D/phi_die_yz1.mco')
              write(14,*) n(2)/every,n(3)/every
              do iz=n(3)+1,1,-1*every 
                 write(14,'(800(e18.6,1x))') ( phi(ig_die(3),iy,iz), iy=1,n(2)+1,every )
              enddo
              close(14)
           endif

           if(ig_die(4).gt.0) then
              open(14,file='DATA/DATA_2D/phi_die_yz2.mco')
              write(14,*) n(2)/every,n(3)/every
              do iz=n(3)+1,1,-1*every 
                 write(14,'(800(e18.6,1x))') ( phi(ig_die(4),iy,iz), iy=1,n(2)+1,every )
              enddo
              close(14)
           endif
           
           if(ig_die(2).gt.0) then
              open(14,file='DATA/DATA_2D/je_die_xz.mco')
              write(14,*) n(1)/every,n(3)/every
              do iz=n(3)+1,1,-1*every 
                 write(14,'(800(e18.6,1x))') ( j_die_xz(1,ix,iz)/real(it), ix=1,n(1)+1,every )
              enddo
              close(14)
              
              open(14,file='DATA/DATA_2D/ji_die_xz.mco')
              write(14,*) n(1)/every,n(3)/every
              do iz=n(3)+1,1,-1*every 
                 write(14,'(800(e18.6,1x))') ( j_die_xz(2,ix,iz)/real(it), ix=1,n(1)+1,every )
              enddo
              close(14)
           endif

           if(ig_die(3).gt.0) then
              open(14,file='DATA/DATA_2D/je_die_yz1.mco')
              write(14,*) n(2)/every,n(3)/every
              do iz=n(3)+1,1,-1*every 
                 write(14,'(800(e18.6,1x))') ( j_die_yz(1,iy,iz)/real(it), iy=1,n(2)+1,every )
              enddo
              close(14)
              
              open(14,file='DATA/DATA_2D/ji_die_yz1.mco')
              write(14,*) n(2)/every,n(3)/every
              do iz=n(3)+1,1,-1*every 
                 write(14,'(800(e18.6,1x))') ( j_die_yz(2,iy,iz)/real(it), iy=1,n(2)+1,every )
              enddo
              close(14)
           endif

        endif
     endif

     ctime(2)= ctime(2) + MSTIMER()

     !
     ! Calculate electric field
     ! 
     call calc_Efield(n,h,phi,Ei,bcnd)
     ctime(1)= ctime(1) + MSTIMER()
     
     !
     ! Sort particles inside vxp()
     !
     if( flag_nopart.eq.0 .and. MOD(it,nsort).eq.1 ) then

        ! Allocate array
        allocate ( sorting(6,nmax*nproc) )

        do ptype=1,ntype
           if( SUM(np_tot(ptype,1:nproc)).gt.0 ) then
              call part_sorting(vxp,sorting,nmax,ntype,n,h,Plist,pl_max,nproc,np_tot,ptype)
           endif
        enddo
        
        ! deallocate
        deallocate ( sorting )
        
        ctime(3)= ctime(3) + MSTIMER()
     endif

     !
     ! Calculate averages
     !
     if( flag_nopart.eq.0 .and. MOD(it,navg).eq.1 ) then 

        n_mts=5
        allocate ( p_mts_xy(n_mts,0:n(1)+2,0:n(2)+2,ntype,nproc), &
             p_mts_xz(n_mts,0:n(1)+2,0:n(3)+2,ntype,nproc), &
             p_mts_yz(n_mts,0:n(2)+2,0:n(3)+2,ntype,nproc) )
        
        !$OMP PARALLEL PRIVATE(iproc,ptype)
        iproc= omp_get_thread_num() + 1
        
        ! Particle moments
        do ptype=1,ntype
           call part_moments(n,h,vxp,nmax,ntype,kq,nproc,np_tot, &
                iproc,ptype,p_mts_xy,p_mts_xz,p_mts_yz,n_mts)           
        enddo ! enddo over ptype particle
        !$OMP END PARALLEL

        if(plt_src.eq.1) then
           flag_diag=1
           ss2D_xz= 0.d0
           ss2D_xy= 0.d0
           ss2D_yz= 0.d0

           call collisions(it,vxp,n,h,ntype,nmax,sig,sig_Er,sig_list,sig_Eex,&
                ncol_mx,npt_mx,cnt_col,P_loss,sour_xy,sour_xz,Plist,pl_max,sigv_mx,&
                col_info,np_red,flag_dead,flag_cex,nproc,np_tot,iseed,nproc_mpi,mpi_rank,&
                ss2D_xy,ss2D_xz,ss2D_yz,flag_diag)
        endif
        
        call calc_avg(n,h,np,p_mts_xy,p_mts_xz,p_mts_yz,data_pavg_xy,&
             data_pavg_xz,data_pavg_yz,sour_xy,sour_xz,sour_yz,&
             ntype,n_mts,ss2D_xy,ss2D_xz,ss2D_yz,nproc,nproc_mpi)

        sour_xy= 0
        sour_xz= 0
        sour_yz= 0
        deallocate( p_mts_xy, p_mts_xz, p_mts_yz )

        ! Calculate average value of the potential
        !$OMP PARALLEL
        !$OMP DO
        do iy=0,n(2)+2
           do ix=0,n(1)+2
              phi_avg_xy(ix,iy)= phi_avg_xy(ix,iy) + phi(ix,iy,iz_pl)
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL  

        ! XZ plane
        !$OMP PARALLEL
        !$OMP DO
        do iz=0,n(3)+2
           do ix=0,n(1)+2
              phi_avg_xz(ix,iz)= phi_avg_xz(ix,iz) +  phi(ix,n(2)/2,iz)
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL  

        ! YZ plane
        if(flag_thr.eq.1) x= (ix_pl-1)*h(1)
        !$OMP PARALLEL
        !$OMP DO
        do iz=0,n(3)+2
           do iy=0,n(2)+2
              phi_avg_yz(iy,iz)= phi_avg_yz(iy,iz) +  phi(ix_pl,iy,iz)
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL  

        if(flag_avg3D.eq.1) then
           call dens_red(n,np,np_red,bcnd,ntype,nproc,nproc_mpi)
           
           !$OMP PARALLEL
           !$OMP DO
           do iz=1,n(3)+1
              do iy=1,n(2)+1
                 do ix=1,n(1)+1
                    avg3D(1:3,ix,iy,iz)= avg3D(1:3,ix,iy,iz) + Ei(1:3,ix,iy,iz)
                    avg3D(4,ix,iy,iz)= avg3D(4,ix,iy,iz) + phi(ix,iy,iz)
                    avg3D(5,ix,iy,iz)= avg3D(5,ix,iy,iz) + np_red(ix,iy,iz,2)
                 enddo
              enddo
           enddo
           !$OMP END DO NOWAIT
           !$OMP END PARALLEL 
        endif
        
        ! Update counter
        cnt_avg(1)= cnt_avg(1) + 1 ! Will be reset every nsav
        cnt_avg(3)= cnt_avg(3) + 1 ! Count from the beginning

        ctime(4)= ctime(4) + MSTIMER()
     endif

     !
     ! Collisions
     !
     if( it.gt.1 .and. ncol.gt.0 .and. MOD(it,ns_coll).eq.1 ) then
        call dens_red(n,np,np_red,bcnd,ntype,nproc,nproc_mpi)
        flag_diag=0
        call collisions(it,vxp,n,h,ntype,nmax,sig,sig_Er,sig_list,sig_Eex,&
             ncol_mx,npt_mx,cnt_col,P_loss,sour_xy,sour_xz,Plist,pl_max,sigv_mx,&
             col_info,np_red,flag_dead,flag_cex,nproc,np_tot,iseed,nproc_mpi,mpi_rank,&
             ss2D_xy,ss2D_xz,ss2D_yz,flag_diag)
        ctime(5)= ctime(5) + MSTIMER() 
     endif

     !
     ! Move particles
     !
     if( MOD(it,ns_heat).eq.0 .and. Pabs.gt.0.d0 ) then
        ! Reduction
        sum_dEk_tot= SUM(sum_dEk) ! sum over OMP proc
        if(nproc_mpi.gt.1) then ! sum over MPI proc
           sum_dEk_tmp=0.d0
           call MPI_ALLREDUCE(sum_dEk_tot, sum_dEk_tmp, 1, MPI_REAL8, MPI_SUM, &
                MPI_COMM_WORLD, ierr)
           sum_dEk_tot= sum_dEk_tmp
        endif

        sum_Nh= SUM(Nh)
        if(nproc_mpi.gt.1) then
           sum_Nh_tmp=0
           call MPI_ALLREDUCE(sum_Nh, sum_Nh_tmp, 1, MPI_INTEGER, MPI_SUM, &
                MPI_COMM_WORLD, ierr)
           sum_Nh= sum_Nh_tmp
        endif

        Te= (2.d0/3.d0)*( sum_dEk_tot + Pabs/nu_h )/(Nm(1)*qe*real(sum_Nh))
        ! Calculated thermal velocity
        vt=dsqrt(2.d0*qe*Te/ABS(mass(1)))
     endif
     
     !$OMP PARALLEL PRIVATE(iproc,ptype,iseed_OMP,sav_np)
     iproc= omp_get_thread_num() + 1
     iseed_OMP= iseed(iproc)

     ! Initialize counters and arrays
     Nh(iproc)=0
     sum_dEk(iproc)=0.d0
        
     ! Electron heating
     if( MOD(it,ns_heat).eq.0 .and. flag_heat.eq.1 ) then 
        if(flag_inj.eq.1) vt= vt0(1) ! Constant electron temperature 
        call eheating(h,vxp,nmax,ntype,nproc,iseed_OMP,np_tot,vt,iproc,P_loss)
     endif

     ! Push particles
     if(flag_nopart.eq.0) then
        do ptype=1,ntype
           ! Do not consider neutrals in this subroutine
           if( ABS(charge(ptype)).gt.0.d0 ) then
              if(ABS(opt_inj).eq.2) sav_np(1)= SUM(p_mac(ptype,np_loss,0:ngrid,iproc))
              if(flag_thr.eq.1) sav_np(2)= p_mac(ptype,np_loss,ind_g,iproc)
              call part_mover(n,h,Ei,Bi,p_mac,P_loss,vxp,bcnd,nmax,ntype,ngrid,flag_dead,&
                   nproc,np_tot,iproc,ptype,flag_cex,cnt_cex,cnt_dead,beam_div,sum_q_xz,&
                   sum_q_yz,dtype,iseed_OMP,n_B,h_B)
              if(ABS(opt_inj).eq.2) N_inj(ptype,iproc)= N_inj(ptype,iproc) + &
                   ( SUM(p_mac(ptype,np_loss,0:ngrid,iproc)) - sav_np(1) )
              if(flag_thr.eq.1) N_flx(ptype,iproc)= N_flx(ptype,iproc) + &
                   ( p_mac(ptype,np_loss,ind_g,iproc) - sav_np(2) )
           endif
        enddo
     endif
     
     ! Particle injection
     if( MOD(it,ns_inj).eq.0 .and. flag_inj.eq.1 ) then
        call part_injection(n,h,bcnd,vxp,sour_xy,sour_xz,sour_yz,ntype,nmax,&
             I_inj,xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow,np_tot,nproc,&
             iseed_OMP,nproc_mpi,N_inj,iproc,P_loss,phi,ni0)
        ! Reset counter 
        N_inj(:,iproc)= 0
     endif

     ! Calculate particle density
     do ptype=1,ntype
        ! Initialize particle density array
        np(:,:,:,ptype,iproc)=0.d0
        if(np_tot(ptype,iproc).gt.0) &
             call charge_deposition(n,h,vxp,nmax,ntype,kq,np,nproc,np_tot,sum_dEk,&
             Nh,iproc,ptype)
     enddo
     
     iseed(iproc)= iseed_OMP
     !$OMP END PARALLEL
     ctime(6)= ctime(6) + MSTIMER()

     ! Inject a flux of particles at a specific location
     if( (jne.gt.0.d0 .or. flag_thr.eq.1) .and. MOD(it,ns_flx).eq.0 ) then
        call part_flux_injection(it,vxp,n,h,ntype,nmax,bcnd,sour_xy,sour_xz,sour_fx_yz,nproc,&
             iseed,flag_cex,np_tot,N_flx,nproc_mpi,mpi_rank,nsav,P_loss)
        N_flx= 0  ! Reset counter
        ctime(7)= ctime(7) + MSTIMER()
     endif

     !
     ! Backup simulation data
     !
     if( it.gt.1 .and. MOD(it,nbak).eq.1 ) then 

        if(mpi_rank.eq.0) print*, 'Backing up simulation data ...'
           
        write (pnum_bck,'(i3)'),mpi_rank
        if(mpi_rank.le.9) lgh=3
        if(mpi_rank.ge.10 .and. mpi_rank.le.99 ) lgh=2
        if(mpi_rank.ge.100 .and. mpi_rank.le.999 ) lgh=1

        ! Create files
        if(flag_wrt.eq.0) then
           open(41+mpi_rank,file='DATA.BAK/particles'//pnum_bck(lgh:3)//'.bak',form='UNFORMATTED')
           if( flag_die.eq.1 .and. mpi_rank.eq.0 ) &
                open(43,file='DATA.BAK/phi.bak',form='UNFORMATTED')
           if( flag_convP.eq.1 .and. mpi_rank.eq.0 ) &
                open(44,file='DATA.BAK/Vgrd.bak',form='UNFORMATTED')
           flag_wrt= 1
        else
           open(41+mpi_rank,file='DATA.BAK2/particles'//pnum_bck(lgh:3)//'.bak',form='UNFORMATTED')
           if(flag_die.eq.1 .and. mpi_rank.eq.0 ) &
                open(43,file='DATA.BAK2/phi.bak',form='UNFORMATTED')
           if( flag_convP.eq.1 .and. mpi_rank.eq.0 ) &
                open(44,file='DATA.BAK2/Vgrd.bak',form='UNFORMATTED')
           flag_wrt= 0
        endif

        ! Save particles
        write(41+mpi_rank) time
        write(41+mpi_rank) np_tot(1:ntype,1:nproc)
        do iproc=1,nproc
           do ptype=1,ntype
              write(41+mpi_rank) vxp(:,1:np_tot(ptype,iproc),ptype,iproc)
           enddo
           if(tag_neg.gt.0) write(41+mpi_rank) flag_cex(1:np_tot(tag_neg,iproc),iproc)
        enddo
        close(41+mpi_rank)

        ! Save dielectric potential
        if( flag_die.eq.1 .and. mpi_rank.eq.0 ) then
           write(43) phi(0:n(1)+2,0:n(2)+2,0:n(3)+2)
           close(43)       
        endif

        ! Save Vgrd
        if( flag_convP.eq.1 .and. mpi_rank.eq.0 ) then
           write(44) Vgrd(1:ngrid)
           close(44)           
        endif

        if(mpi_rank.eq.0) print*, 'Done!'
        
        !
        ! Save particle distribution function
        !
        if(flag_pdf.eq.1) then
           if(mpi_rank.eq.0) then 

              print*, 'Saving particle distribution function ...'
           
              do i_rg=1,n_rg
                 write (pnum,'(i1)'),i_rg
                 open(41+i_rg,file='DATA/xv'//pnum//'.dat')
              enddo
              
              do iproc=1,nproc
                 do i=1,np_tot(ptype_pdf,iproc)
                    do i_rg=1,n_rg
                       
                       ! Save velocity for each predefined sub-regions
                       if ( vxp(1,i,ptype_pdf,iproc).ge.xl_rg(i_rg) .and. &
                            vxp(1,i,ptype_pdf,iproc).le.xr_rg(i_rg) .and. &
                            vxp(2,i,ptype_pdf,iproc).ge.yl_rg .and. &
                            vxp(2,i,ptype_pdf,iproc).le.yr_rg .and. &
                            vxp(3,i,ptype_pdf,iproc).ge.zl_rg .and. &
                            vxp(3,i,ptype_pdf,iproc).le.zr_rg ) then  
                          write(41+i_rg,'(3(es14.6,1x))') vxp(4,i,ptype_pdf,iproc),&
                               vxp(5,i,ptype_pdf,iproc),vxp(6,i,ptype_pdf,iproc)
                          goto 30
                       endif
                       
                    enddo
30                  continue
                 enddo
              enddo
           
              do i_rg=1,n_rg
                 write (pnum,'(i1)'),i_rg
                 close(41+i_rg)
              enddo
              
              print*, 'Done!'
           endif
        endif
        
        ctime(8)= ctime(8) + MSTIMER()
        
     endif

     !
     ! Save data
     !
     if( MOD(it,nsav).eq.1 ) then

        ! MPI SUM
        sum_np_tot= SUM(np_tot,DIM=2) ! sum over OMP proc
        if(nproc_mpi.gt.1) then ! sum over MPI proc
           sum_np_tot_tmp=0
           call MPI_ALLREDUCE(sum_np_tot, sum_np_tot_tmp, ntype, MPI_INTEGER, MPI_SUM, &
                MPI_COMM_WORLD, ierr)
           sum_np_tot= sum_np_tot_tmp
        endif
            
        ! Negative ion impacts
        if(tag_neg.gt.0) then 
           sum_cex=0
           call MPI_REDUCE(SUM(cnt_cex,DIM=2), sum_cex, 3, MPI_INTEGER, MPI_SUM, &
                0,MPI_COMM_WORLD, ierr)
           
           sum_beam_div=0.d0
           call MPI_REDUCE(SUM(beam_div,DIM=2), sum_beam_div, 4, MPI_REAL8, MPI_SUM, &
                0,MPI_COMM_WORLD, ierr)
           if(sum_cex(1).gt.0) sum_beam_div= sum_beam_div/sum_cex(1)
        endif

        cnt_avg(2)= nsav
        if( it.eq.1 .or. flag_sav.eq.1 ) cnt_avg(2)= it
        flag_bak=0
        if( flag_sav.eq.0 .and. tseq.gt.0.d0 .and. &
             time.ge.tseq_init .and. time.le.tseq_final ) then
           ! Save sequence of profiles in Macho format.
           if(MOD(it,nseq).eq.1 ) flag_bak=2 
        endif
        call write_data(it,time,n,h,p_mac,P_loss,phi_avg_xy,phi_avg_xz,&
             phi_avg_yz,data_pavg_xy,data_pavg_xz,data_pavg_yz,ntype,ngrid,&
             cnt_col,ncol_mx,sig_list,nproc,cnt_avg,mpi_rank,nproc_mpi,Vgrd,&
             avg3D,flag_avg3D,I_inj,nEf)        
        
        if(mpi_rank.eq.0) then
           
           ! Save energy conservation, particle errors & calculation time
           sum_time= real(SUM(ctime(1:10)))
           if(it.eq.1) write(12,'(a105)') &
                '# time(s), nsteps, ctime{E/rho, poisson, sorting, avg/write, MC, mover, &
                neg. ions, bck}(%), total(ms)'
           write(12,100) time,it,( real(ctime(i))/sum_time*100., i=1,8 ), &
                real(SUM(ctime(1:10)))/real(it)
100        format(es15.8,2x,i8,8(2x,es10.2),(2x,f7.1))
           
           ! Print on screen
           write(*,101) it,time*1.d6, (sum_np_tot(ptype), ptype=1,ntype)
101        format(' it= ',i8,', t (us)= ',f6.1,', np= ',10(i10,2x))
           
           ! Save PDE solver residual & # of iterations
           write(40,*) res,niter,nint(ksor)

           ! Save negative ion impact characteristics on num_grd
           if( tag_neg.gt.0 .and. sum_cex(1).gt.0 ) &
                write(13,102) time,real(sum_cex(2))/real(sum_cex(1))*100.d0,&
                real(sum_cex(3))/real(sum_cex(1))*100.d0, &
                real(sum_cex(1)-sum_cex(2))/real(sum_cex(1))*100.d0,&
                dsqrt( sum_beam_div(2)-sum_beam_div(1)**2  ),&
                dsqrt( sum_beam_div(4)-sum_beam_div(3)**2  )
102        format(6(es12.2,1x))
           
           ! Print on screen
           if(it.gt.1) &
                write(*,103) ( real(ctime(i))/real(cnt_avg(2)), i=1,6 ),&
                real(ctime(8))/real(cnt_avg(2)), &
                real(SUM(ctime(1:10)))/real(cnt_avg(2))
103        format(' <t> (ms): E/rho= ',f6.1,', poisson= ',f6.1,', sorting= ',f6.1, &
                ', avg/write= ',f6.1,', MC= ',f6.1,', mover= ',f7.1,', bck= ',f6.1,', total= ',f7.1)

           ! Save 3D Efield map
           if(flag_avg3D.eq.1) then
              print*, 'Saving 3D potentiel and density maps ...'
              if(flag_bak.eq.2) then
                 ! Save 3D profiles for a given time window
                 cnt_plt= cnt_plt+1
                 if(cnt_plt.lt.10) then 
                    write (plnum,'(i1)'),cnt_plt
                    i_pl=1
                 endif
                 if( cnt_plt.ge.10 .and. cnt_plt.lt.100 ) then 
                    write (plnum,'(i2)'),cnt_plt
                    i_pl=2
                 endif
                 if( cnt_plt.ge.100 .and. cnt_plt.lt.1000 ) then 
                    write (plnum,'(i3)'),cnt_plt
                    i_pl=3
                 endif
                 
                 if(i_pl.eq.1) corrnum= '_00'
                 if(i_pl.eq.2) corrnum= '_0'
                 if(i_pl.eq.3) corrnum= '_'
                 
                 open(41,file='DATA/phi_n_3D'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.dat',form='UNFORMATTED')
              else
                 open(41,file='DATA/phi_n_3D.dat',form='UNFORMATTED')
              endif
              write(41) n(1),n(2),n(3)
              write(41) avg3D(4:5,1:n(1)+1,1:n(2)+1,1:n(3)+1)/real(cnt_avg(1))           
              close(41)
           endif
           
        endif
        
        ! Negative ion current lost through collisions
        if( it.gt.1 .and. tag_neg.gt.0 ) then 
              cnt_dead_tmp=0.d0
              call MPI_REDUCE(SUM(cnt_dead,DIM=1), cnt_dead_tmp, 1, MPI_REAL8, MPI_SUM, &
                   0,MPI_COMM_WORLD, ierr)
              if(mpi_rank.eq.0 .and. cnt_dead_tmp.gt.0) print'(1x,"Neg. ion current lost through collisions (A)=",es10.2)', &
                   cnt_dead_tmp/(real(cnt_avg(2))*dt)
        endif

        ! Reset arrays
        if(flag_sav.eq.0) then
           phi_avg_xy=0.d0
           phi_avg_xz=0.d0
           phi_avg_yz=0.d0
           if(plt_src.eq.1) then
              data_pavg_xy(1:8,:,:,:)=0.d0
              data_pavg_xz(1:8,:,:,:)=0.d0
              data_pavg_yz(1:8,:,:,:)=0.d0
           else
              data_pavg_xy(1:6,:,:,:)=0.d0 ! All except sour_avg
              data_pavg_xz(1:6,:,:,:)=0.d0
              data_pavg_yz(1:6,:,:,:)=0.d0
           endif
           cnt_avg(1)=0
           p_mac=0.d0
           P_loss=0.d0
           cnt_cex=0
           cnt_dead=0.d0
           beam_div=0.d0
           ctime=0
           if(flag_avg3D.eq.1) avg3D= 0.d0
        endif

        ! Time lag for writing data and calculating averages
        ctime(4)= ctime(4) + MSTIMER()
        
     endif

  enddo ! enddo over ntmax iterations

  ! 
  ! Close external files
  !
  close(12)
  close(13)
  close(40)
  close(42)

  call MPI_Finalize(ierr)

end program main

subroutine introduction
!     ==============================================================
!     VERSION:         3.3.6.2
!     LAST MOD:       Sep/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Display code info
!     NOTE:              /
!     --------------------------------------------------------------
  implicit none
  
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) '+                                                                          +'
  write(*,*) '+   3D explicit parallel particle in cell code                             +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ LePIC 3D algorithm version: 3.3.6.2                                      +'
  write(*,*) '+ Last modifications: Sep/2024                                             +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ Copyright or © or Copr. Gwenael Fubiani (2022/11/22)                     +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ gwenael.fubiani@cnrs.fr                                                  +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ This software is a computer program whose purpose is to model low        +'
  write(*,*) '+ temperature plasmas with a Particle-In-Cell algorihm in 3D-3V            +'
  write(*,*) '+ dimensions.                                                              +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ This software is governed by the CeCILL-C license under French law and   +'
  write(*,*) '+ abiding by the rules of distribution of free software.  You can  use,    +'
  write(*,*) '+ modify and/ or redistribute the software under the terms of the CeCILL-C +'
  write(*,*) '+ license as circulated by CEA, CNRS and INRIA at the following URL        +'
  write(*,*) '+ "http://www.cecill.info".                                                +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ As a counterpart to the access to the source code and  rights to copy,   +'
  write(*,*) '+ modify and redistribute granted by the license, users are provided only  +'
  write(*,*) "+ with a limited warranty and the software's author, the holder of the     +"
  write(*,*) '+ economic rights,  and the successive licensors  have only  limited       +'
  write(*,*) '+ liability.                                                               +'
  write(*,*) '+                                                                          +'
  write(*,*) "+ In this respect, the user's attention is drawn to the risks associated   +"
  write(*,*) '+ with loading,  using,  modifying and/or developing or reproducing the    +'
  write(*,*) '+ software by the user in light of its specific status of free software,   +'
  write(*,*) '+ that may mean  that it is complicated to manipulate,  and  that  also    +'
  write(*,*) '+ therefore means  that it is reserved for developers  and  experienced    +'
  write(*,*) '+ professionals having in-depth computer knowledge. Users are therefore    +'
  write(*,*) "+ encouraged to load and test the software's suitability as regards their  +"
  write(*,*) '+ requirements in conditions enabling the security of their systems and/or +'
  write(*,*) '+ data to be ensured and,  more generally, to use and operate it in the    +'
  write(*,*) '+ same conditions as regards security.                                     +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ The fact that you are presently reading this means that you have had     +'
  write(*,*) '+ knowledge of the CeCILL-C license and that you accept its terms.         +'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) '                                                  '
  
  return
end subroutine introduction
