subroutine read_input(n,tmax,xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow,&
     nsav,eps,omega,kt,rname,ngrid,ng,xl_rg,xr_rg,yl_rg,yr_rg,zl_rg,zr_rg,&
     I_inj,Ca,mpi_rank,flag_avg3D,np_dup,n_B,phi0_RF,f0_RF,&
     phi1_RF,f1_RF)
!     ==============================================================
!     VERSION:         0.3
!     LAST MOD:      Dec/23
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  integer:: iB,n(3),nsav,flag_read,ngrid,i_rg,ng,mpi_rank,flag_avg3D,n_B(3)
  real(kind=8):: eps,tmax,kt,xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow,&
       xl_rg(nm_rg),xr_rg(nm_rg),yl_rg,yr_rg,zl_rg,zr_rg,omega,I_inj,Ca,np_dup,&
       phi0_RF,f0_RF,phi1_RF,f1_RF
  character:: end_file*3,ans*1,rname*20

  ! Initialize variables & arrays
  xl_rg=0.d0
  xr_rg=0.d0
  flag_grd=0

  ! Open input file
  open(10,file='input_dir/conditions.inp')

  read(10,*,end=999) rname ! name of file storing reactions
  read(10,*,end=999) Ti(1) ! read electron temperature (eV)

  !
  ! Simulation parameters
  !
  read(10,*,err=999) ng  ! # of levels in MG
  read(10,*,err=999) eps,omega ! Min. precision in residual, w for SOR
  read(10,*,err=999) n(1),n(2),n(3) ! n: nx=2**n, ny=2**n., nz=2**n., Number of grid points
  
  read(10,*,err=999) xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow ! particle heating region [cm]
  xl_pow= xl_pow*1.d-2 ! convert to meters
  xr_pow= xr_pow*1.d-2
  yl_pow= yl_pow*1.d-2
  yr_pow= yr_pow*1.d-2
  zl_pow= zl_pow*1.d-2
  zr_pow= zr_pow*1.d-2
  read(10,*,err=999) x_load ! particle loading region [cm]
  x_load= x_load*1.d-2

  flag_ahp=0
  if(yr_pow.lt.0.d0) then
     flag_ahp=1 ! Option to deposit power @ r > Rahp
     yr_pow= ABS(yr_pow)
  endif
  
  read(10,*,err=999) nB ! # of magnetic field files (<0 pos. ions not magnetized)   
  flag_B_pos=0
  if(nB.lt.0) flag_B_pos=1
  nB= ABS(nB)
  read(10,*,err=999) n_B(1),n_B(2),n_B(3) ! Grid for the magnetic field map
  flag_B=0
  do iB=1,nB
     ! Ext. file (y/n/s), name, scaling, B(G), direction, L/x0(m)
     read(10,*,err=999) B_file(iB),B_name(iB),B_scale(iB),B0(iB),B_info(iB),dL(iB),x0(iB),y0(iB),z0(iB) 
     B0(iB)= B0(iB)*1.d-4 ! convert into Teslas
     if(B_file(iB).ne.'s' .or. B_file(iB).ne.'S') flag_B=1
  enddo
  if( n_B(1).le.1 .or.  n_B(2).le.1 .or.  n_B(3).le.1 .or. flag_B.eq.0 .or. nB.eq.0 ) n_B= 1
  if(flag_B.eq.0) nB=0
  dL= dL*1.d-2
  x0= x0*1.d-2
  y0= y0*1.d-2
  z0= z0*1.d-2

  read(10,*,err=999) Pabs,nu_h ! Pabs(W) (<0 to turn off), Heating freqency (<0 == off)
  flag_heat=1
  if(nu_h.lt.0.d0) flag_heat=0
  read(10,*,err=999) I_inj, opt_inj ! Inject current I(A/m) (<0 == off), options (1=fixed, 2= equal to ion losses)
  flag_inj=1
  if(I_inj.lt.0.d0) flag_inj=0
 
  read(10,*,err=999) tmax ! s
  read(10,*,err=999) kt ! time coeff. kt ; dt= kt*dx/vt
  read(10,*,err=999) n0 ! m-3
  read(10,*,err=999) ngas ! m-3
  read(10,*,err=999) np_cell ! number of particles per cell
  read(10,*,err=999) ngrid ! number of wall labels

  flag_avg3D= 1
  if(ngrid.lt.0) then
     flag_avg3D= 0 ! Save 3D Efield map : 0 no, 1 yes.
     ngrid= ABS(ngrid)
  endif
  
  read(10,*,err=999) nsav ! frequency for data saving
  read(10,*,err=999) nbak,cnt_plt,tseq_init,tseq_final ! Frequency for simulation backup
  tseq= ABS(tseq_final-tseq_init)/real(cnt_plt)
  
  read(10,*,err=999) k_eps0 
  flag_convP=0
  if(k_eps0.lt.0) then
     flag_convP=1
     k_eps0= ABS(k_eps0)
  endif

  read(10,*,err=999) ans,ptype_pdf,yl_rg,yr_rg,zl_rg,zr_rg,n_rg,(xl_rg(i_rg),xr_rg(i_rg),i_rg=1,n_rg) ! Save pdf (y/n), ptype, yl & yr, # of regions n_rg, xl & xr for each n_rg
  zl_rg= zl_rg*1.d-2
  zr_rg= zr_rg*1.d-2
  yl_rg= yl_rg*1.d-2
  yr_rg= yr_rg*1.d-2
  xl_rg= xl_rg*1.d-2
  xr_rg= xr_rg*1.d-2

  if( mpi_rank.eq.0 .and. n_rg.gt.nm_rg ) then ! Warning
     print*, 'Warning: nm_rg is too small, please correct ...'
     call stop_calculation
  endif

  flag_pdf=0
  if(ans.eq.'y'.or.ans.eq.'Y') flag_pdf=1
   
  read(10,*,err=999) ans,omp_rank_max,mpi_rank_max,np_dup ! restart (yes/no/external/gauss),  # of OMP threads, # of MPI, # for particle duplication

  flag_restart=0
  if(ans.eq.'y'.or.ans.eq.'Y') flag_restart=1
  if(ans.eq.'e'.or.ans.eq.'E') flag_restart=2
  if(flag_restart.eq.0) np_dup=1
   
  read(10,*,err=999) ans,Lhy,Lhz,nhy,nhz,ind_g
  if(ans.eq.'y'.or.ans.eq.'Y') flag_grd=1

  read(10,*,err=999) jne,THm,num_grd,Ca ! jH-(mA/cm2) from H on PE surface, TH-(eV), characteristics of neg. ions on grid # [index], dielectric capacitance Ca/S (Fd/m2)

  read(10,*,err=999) gam_sec,igrid_sec,ans,phi0_RF,f0_RF,phi1_RF,f1_RF ! secondary emission coeff, grid index (<0, fixed Is instead), phi0(V), f0(Hz), phi1(V), f1(Hz)
  
  flag_RFpot=0
  if(ans.eq.'y'.or.ans.eq.'Y') flag_RFpot=1
  
  if( I_inj.gt.0.d0 .and. gam_sec.gt.0.d0 ) then
     print*, 'Warning: gam_sec>0 is incompatible with I_inj>0, please correct ...'
     call stop_calculation
  endif
  
  flag_read=0
  do while( flag_read.eq.0 )
     read(10,*,err=999) end_file
     if( end_file.eq.'END' .or. end_file.eq.'end'  ) flag_read=1
  enddo
 
  if(flag_read.eq.1) then
     if(mpi_rank.eq.0) print*, 'Input file was read correctly ...'
     close(10)
     return
  endif

999 if(mpi_rank.eq.0) print*, 'Input file was not read correctly!'  
  call stop_calculation

end subroutine read_input
