subroutine write_data(it,time,n,h,p_mac,P_loss,phi_avg_xy,&
     phi_avg_xz,phi_avg_yz,data_pavg_xy,data_pavg_xz,data_pavg_yz,&
     ntype,ngrid,cnt_col,ncol_mx,sig_list,nproc,cnt_avg,mpi_rank,&
     nproc_mpi,Vgrd,avg3D,flag_avg3D,I_inj,nEf)
!     ==============================================================
!     VERSION:         0.5
!     LAST MOD:      Sep/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:
!     NOTE: 
!     --------------------------------------------------------------
  use omp_lib
  use mpi
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer ierr
  integer:: ix,iy,iz,ptype,n(3),nproc,ntype,i_rg,cnt_avg(3)
  integer:: np_avg,Tp_avg,sour_avg,iproc,dummy, &
       u1_avg,u2_avg,u3_avg,sink_avg,ngrid,igrid,icol,ncol_mx,it,&
       sig_list(npart,ncol_mx),mpi_rank,nproc_mpi,kn,flag_avg3D,nEf(3)
  parameter ( np_avg=1, Tp_avg=2, u1_avg=3, u2_avg=4, u3_avg=5, dummy=6, sour_avg=7,&
       sink_avg=8 ) 
  real(kind=8):: h(3),p_mac(ntype,2,0:ngrid,nproc), &
       phi_avg_xy(0:n(1)+2,0:n(2)+2), &
       phi_avg_xz(0:n(1)+2,0:n(3)+2), &
       phi_avg_yz(0:n(2)+2,0:n(3)+2), &
       data_pavg_xy(8,0:n(1)+2,0:n(2)+2,ntype),&
       data_pavg_xz(8,0:n(1)+2,0:n(3)+2,ntype),&
       data_pavg_yz(8,0:n(2)+2,0:n(3)+2,ntype),&
       avg3D(5,1:nEf(1)+1,1:nEf(2)+1,1:nEf(3)+1),&
       P_loss(4,ntype,nproc),time,cnt_col(ncol_mx,nproc,nm_rg), &
       sum_nu(ntype,nm_rg),P_loss_tmp(4,ntype),cnt_col_tmp(ncol_mx,nm_rg),&
       p_mac_tmp(ntype,2,0:ngrid),dr,Ip,Im,Pwall,Pext,Pinj,Pcoll,I_mw,&
       Ek_ew,Ek_iw,Vgrd(ngrid),Ptmp(ntype),I_inj,Isec(ngrid),Ip_cath,np_tmp
  real(kind=8),allocatable:: ld_xy(:,:),ld_xz(:,:),ld_yz(:,:),cnt_col_red(:,:),&
       p_mac_red(:,:,:),P_loss_red(:,:)       
  character, save:: name(192)*15,pnum*1,name0D(16)*15,strg*3

  allocate( ld_xy(0:n(1)+2,0:n(2)+2), &
       ld_xz(0:n(1)+2,0:n(3)+2), &
       ld_yz(0:n(2)+2,0:n(3)+2) )
  
  ! Initialize
  ld_xy=0.d0
  ld_xz=0.d0
  ld_yz=0.d0

  !
  ! Calculate average collision frequency  
  !
  cnt_col_tmp= 0

  ! Reduction 
  do i_rg=1,n_rg
     do iproc=1,nproc
        do icol=1,ncol
           cnt_col_tmp(icol,i_rg)= cnt_col_tmp(icol,i_rg) + &
                cnt_col(icol,iproc,i_rg)/(real(it)*dt)
        enddo
     enddo
  enddo

  if(nproc_mpi.gt.1) then
     allocate( cnt_col_red(ncol_mx,nm_rg) )
     cnt_col_red=0.d0
     call MPI_REDUCE(cnt_col_tmp(:,:), cnt_col_red, ncol_mx*nm_rg, &
          MPI_REAL8, MPI_SUM, 0,MPI_COMM_WORLD, ierr)
     cnt_col_tmp= cnt_col_red
     deallocate( cnt_col_red )
  endif

  !
  ! Calculate average power & current loss 
  !

  P_loss_tmp= 0.d0
  p_mac_tmp= 0.d0

  ! Reduction
  do iproc=1,nproc
     do ptype=1,ntype
        P_loss_tmp(:,ptype)= P_loss_tmp(:,ptype) + &
             P_loss(:,ptype,iproc)/real(cnt_avg(2))/dt
        p_mac_tmp(ptype,np_loss,:)= p_mac_tmp(ptype,np_loss,:) + &
             p_mac(ptype,np_loss,:,iproc)*charge(ptype)*Nm(ptype)/real(cnt_avg(2))/dt
        p_mac_tmp(ptype,P_w,:)= p_mac_tmp(ptype,P_w,:) + &
             p_mac(ptype,P_w,:,iproc)/real(cnt_avg(2))/dt
     enddo
  enddo

  if(nproc_mpi.gt.1) then
     allocate( p_mac_red(ntype,2,0:ngrid) )
     p_mac_red=0.d0
     call MPI_REDUCE(p_mac_tmp(:,:,:), p_mac_red, ntype*2*(ngrid+1), &
          MPI_REAL8, MPI_SUM, 0,MPI_COMM_WORLD, ierr)
     p_mac_tmp= p_mac_red
     deallocate( p_mac_red )

     allocate( P_loss_red(4,ntype) )
     P_loss_red=0.d0
     call MPI_REDUCE(P_loss_tmp(:,:), P_loss_red, 4*ntype, &
          MPI_REAL8, MPI_SUM, 0,MPI_COMM_WORLD, ierr)
     P_loss_tmp= P_loss_red
     deallocate( P_loss_red )
  endif

  Ip= SUM(p_mac_tmp(2:ntype,np_loss,0:ngrid))
  Ip_cath= SUM(p_mac_tmp(2:ntype,np_loss,igrid_sec))
  if(tag_neg.gt.0) then
     Ip= Ip - SUM(p_mac_tmp(tag_neg,np_loss,0:ngrid))
     Ip_cath= Ip_cath - p_mac_tmp(tag_neg,np_loss,igrid_sec)
  endif
  if(nproc_mpi.gt.1) then
     call MPI_Bcast(Ip, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(Ip_cath, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  endif

  ! Isec >0 as it corresponds to electrons drawn from the power supply
  Isec=0.d0
  ! Secondary electron current emmited on cathode.
  if(gam_sec.gt.0) Isec(igrid_sec)= gam_sec*Ip_cath
  if(I_inj.gt.0 .and. ABS(opt_inj).eq.4) Isec(igrid_sec)= ABS(I_inj)
  ! Note that we consider here an incoming beam the same way as electrons drawn from a biased surface.
  p_mac_tmp(1,np_loss,igrid_sec)= p_mac_tmp(1,np_loss,igrid_sec) + Isec(igrid_sec)
  
  !
  ! Only on rank 0  
  !
  if(mpi_rank.eq.0) then

     dr= (h(1)*h(2)*h(3))**(1./3.)

     ! Calculate Debye length
     do ptype=1,ntype

        ! XY plane
        !$OMP PARALLEL
        !$OMP DO
        do iy=0,n(2)+2
           do ix=0,n(1)+2
              np_tmp= data_pavg_xy(np_avg,ix,iy,ptype)  
              if(np_tmp.gt.0.d0) then
                 ! Electron Debye length
                 if(ptype.eq.1) then 
                    ld_xy(ix,iy)= dsqrt( (eps0*data_pavg_xy(Tp_avg,ix,iy,1))/ &
                         (data_pavg_xy(np_avg,ix,iy,1)*qe) )
                    ! Calculate ratio with dr
                    if(ld_xy(ix,iy).gt.0) ld_xy(ix,iy)=dr/ld_xy(ix,iy)
                 endif
                 ! Calculate average velocity u                   
                 data_pavg_xy(u1_avg,ix,iy,ptype)= data_pavg_xy(u1_avg,ix,iy,ptype)/np_tmp
                 data_pavg_xy(u2_avg,ix,iy,ptype)= data_pavg_xy(u2_avg,ix,iy,ptype)/np_tmp
                 data_pavg_xy(u3_avg,ix,iy,ptype)= data_pavg_xy(u3_avg,ix,iy,ptype)/np_tmp                  
              endif
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL  
        
        ! XZ plane
        !$OMP PARALLEL
        !$OMP DO
        do iz=0,n(3)+2
           do ix=0,n(1)+2
              np_tmp= data_pavg_xz(np_avg,ix,iz,ptype)
              if(np_tmp.gt.0.d0) then
                 ! Electron Debye length
                 if(ptype.eq.1) then 
                    ld_xz(ix,iz)= dsqrt( (eps0*data_pavg_xz(Tp_avg,ix,iz,1))/ &
                         (data_pavg_xz(np_avg,ix,iz,1)*qe) )
                    ! Calculate ratio with dr
                    if(ld_xz(ix,iz).gt.0) ld_xz(ix,iz)=dr/ld_xz(ix,iz)
                 endif
                 ! Calculate average velocity u                   
                 data_pavg_xz(u1_avg,ix,iz,ptype)= data_pavg_xz(u1_avg,ix,iz,ptype)/np_tmp
                 data_pavg_xz(u2_avg,ix,iz,ptype)= data_pavg_xz(u2_avg,ix,iz,ptype)/np_tmp
                 data_pavg_xz(u3_avg,ix,iz,ptype)= data_pavg_xz(u3_avg,ix,iz,ptype)/np_tmp  
              endif
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL  

        ! YZ plane
        !$OMP PARALLEL
        !$OMP DO
        do iz=0,n(3)+2
           do iy=0,n(2)+2
              np_tmp= data_pavg_yz(np_avg,iy,iz,ptype)
              if(np_tmp.gt.0.d0) then
                 ! Electron Debye length
                 if(ptype.eq.1) then 
                    ld_yz(iy,iz)= dsqrt( (eps0*data_pavg_yz(Tp_avg,iy,iz,1))/ &
                         (data_pavg_yz(np_avg,iy,iz,1)*qe) )
                    ! Calculate ratio with dr
                    if(ld_yz(iy,iz).gt.0) ld_yz(iy,iz)=dr/ld_yz(iy,iz)
                 endif
                 ! Calculate average velocity u                   
                 data_pavg_yz(u1_avg,iy,iz,ptype)= data_pavg_yz(u1_avg,iy,iz,ptype)/np_tmp
                 data_pavg_yz(u2_avg,iy,iz,ptype)= data_pavg_yz(u2_avg,iy,iz,ptype)/np_tmp
                 data_pavg_yz(u3_avg,iy,iz,ptype)= data_pavg_yz(u3_avg,iy,iz,ptype)/np_tmp  
              endif
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL  

     enddo

     !
     ! Write data
     !
     Im= 0.d0
     if(tag_neg.gt.0) Im= SUM(p_mac_tmp(tag_neg,np_loss,0:ngrid))
     do ptype=1,ntype ! Include contribution of the bias voltage on the total power 
        Ptmp(ptype)= SUM(p_mac_tmp(ptype,np_loss,1:ngrid)*Vgrd(1:ngrid))
     enddo
     Pwall= -SUM(P_loss_tmp(1,1:ntype))
     Pext= SUM(P_loss_tmp(2,1:ntype)) - SUM(Ptmp(1:ntype)) ! -sum(I)*Vbias
     Pcoll= SUM(P_loss_tmp(3,1:ntype))
     Pinj= SUM(P_loss_tmp(4,1:ntype))
     I_mw= SUM(p_mac_tmp(1,np_loss,0:ngrid))+Im
     Ek_ew=0.d0
     if(ABS(SUM(p_mac_tmp(1,np_loss,0:ngrid))).gt.0.d0) &
          Ek_ew= SUM(p_mac_tmp(1,P_w,0:ngrid))/ABS(SUM(p_mac_tmp(1,np_loss,0:ngrid)))
     Ek_iw=0.d0
     if(SUM(ABS(p_mac_tmp(2:ntype,np_loss,0:ngrid))).gt.0.d0) &
          Ek_iw= SUM(p_mac_tmp(2:ntype,P_w,0:ngrid))/SUM(ABS(p_mac_tmp(2:ntype,np_loss,0:ngrid))) 
     write(*,102) Pwall,Pext,Pcoll,Pinj,I_mw,Ip,Ek_ew,Ek_iw  ! print on screen
     
     open(22,file='DATA/Ptot.dat',access='APPEND')
     if( it.eq.1 .and. flag_restart.eq.0 ) &
          write(22,'("# Time (s), Pwall (W), Pabs, Pcoll, Pinj, I_mw (A), I_pw, Ek_ew (eV), Ek_iw, I_inj(A), Vgrd")')
     write(22,101) time,Pwall,Pext,Pcoll,Pinj,I_mw,Ip,Ek_ew,Ek_iw,I_inj,Vgrd(igrid_sec)  
     close(22)
     
     !
     ! Get data file names
     !
     do ptype=1,ntype
        write (pnum,'(i1)'),ptype

        name(24*ptype-23)='n'//pnum//'_xy.mco'
        name(24*ptype-22)='n'//pnum//'_xz.mco'
        name(24*ptype-21)='n'//pnum//'_yz.mco'
        name(24*ptype-20)='T'//pnum//'_xy.mco'        
        name(24*ptype-19)='T'//pnum//'_xz.mco'
        name(24*ptype-18)='T'//pnum//'_yz.mco'        
        name(24*ptype-17)='sour'//pnum//'_xy.mco'
        name(24*ptype-16)='sour'//pnum//'_xz.mco'
        name(24*ptype-15)='sour'//pnum//'_yz.mco'
        name(24*ptype-14)='j'//pnum//'_xy.mco'
        name(24*ptype-13)='j'//pnum//'_xz.mco'
        name(24*ptype-12)='j'//pnum//'_yz.mco'
        name(24*ptype-11)='sink'//pnum//'_xy.mco'
        name(24*ptype-10)='sink'//pnum//'_xz.mco'
        name(24*ptype-9)='sink'//pnum//'_yz.mco'
        name(24*ptype-8)='u'//pnum//'x_xy.mco'
        name(24*ptype-7)='u'//pnum//'x_xz.mco'
        name(24*ptype-6)='u'//pnum//'x_yz.mco'
        name(24*ptype-5)='u'//pnum//'y_xy.mco'
        name(24*ptype-4)='u'//pnum//'y_xz.mco'
        name(24*ptype-3)='u'//pnum//'y_yz.mco'
        name(24*ptype-2)='u'//pnum//'z_xy.mco'
        name(24*ptype-1)='u'//pnum//'z_xz.mco'
        name(24*ptype)='u'//pnum//'z_yz.mco'
        
        name0D(2*ptype-1)='Iw'//pnum//'.dat' 
        name0D(2*ptype)='Pw'//pnum//'.dat' 
                
     enddo
 
     ! Collision frequencies per specie per sub-regions
     do i_rg=1,n_rg

        write (pnum,'(i1)'),i_rg
        
        do ptype=1,ntype
           sum_nu(ptype,i_rg)= SUM(cnt_col_tmp(sig_list(ptype,1:p_ncol(ptype)),i_rg))

           ! Warning
           if(sum_nu(ptype,i_rg).gt.0.25d0*nu_uplim(ptype)) then
              print*, 'Warning: sum_nu/nu_uplim=',sum_nu(ptype,i_rg)/nu_uplim(ptype),'>25% for ptype=',ptype
              print*, 'nu_max=',nu_max(ptype),', nu_uplim=',nu_uplim(ptype)
              call stop_calculation
           endif
           
           ! Warning
           if( (sum_nu(ptype,i_rg)*real(ns_coll)*dt).gt.0.2d0 ) then
              print*, 'Warning: sum_nu*dt > 20% for ptype=',ptype
           endif
           
        enddo
        
        open(14,file='DATA/frequency'//pnum//'.out')
        
        write(14,*) '# null collision frequency for each species'
        write(14,'(10(1x,es15.8))') ( nu_max(ptype),ptype=1,ntype )
        write(14,*) '# sum of collision frequencies for each species'
        write(14,'(10(1x,es15.8))') ( sum_nu(ptype,i_rg),ptype=1,ntype )
        write(14,*) '# reaction index/frequency for each species'
        
        do icol=1,ncol
           do ptype=1,ntype
              strg='no'
              if(ptype.eq.ntype) strg='yes'
              if(icol.le.p_ncol(ptype)) then
                 write(14,'((i3,1x,es15.8))',advance=strg) sig_list(ptype,icol), &
                      cnt_col_tmp(sig_list(ptype,icol),i_rg)
              else
                 write(14,'((i3,1x,es15.8))',advance=strg) 0,0.
              endif
           enddo
        enddo
        
        close(14)

     enddo
  endif

  ! Save potential and Debye length
  if(mpi_rank.eq.0) then

     open(14,file='DATA/DATA_2D/dr_xy.mco')
     open(15,file='DATA/DATA_2D/phi_xy.mco')
     open(16,file='DATA/DATA_2D/dr_xz.mco')
     open(17,file='DATA/DATA_2D/phi_xz.mco')
     open(18,file='DATA/DATA_2D/dr_yz.mco')
     open(19,file='DATA/DATA_2D/phi_yz.mco')
          
     ! XY plane
     write(14,*) n(1)/every,n(2)/every
     write(15,*) n(1)/every,n(2)/every
     do iy=n(2)+1,1,-1*every
        write(14,103) ( ld_xy(ix,iy), ix=1,n(1)+1,every )
        write(15,103) ( phi_avg_xy(ix,iy)/real(cnt_avg(1)), ix=1,n(1)+1,every )
     enddo
     
     ! XZ plane
     write(16,*) n(1)/every,n(3)/every
     write(17,*) n(1)/every,n(3)/every
     do iz=n(3)+1,1,-1*every
        write(16,103) ( ld_xz(ix,iz), ix=1,n(1)+1,every )
        write(17,103) ( phi_avg_xz(ix,iz)/real(cnt_avg(1)), ix=1,n(1)+1,every )
     enddo

     ! YZ plane
     write(18,*) n(2)/every,n(3)/every
     write(19,*) n(2)/every,n(3)/every
     do iz=n(3)+1,1,-1*every
        write(18,103) ( ld_yz(iy,iz), iy=1,n(2)+1,every )
        write(19,103) ( phi_avg_yz(iy,iz)/real(cnt_avg(1)), iy=1,n(2)+1,every )
     enddo
     
     close(14)
     close(15)
     close(16)
     close(17)
     close(18)
     close(19)


     ! Save E-field
     if(flag_avg3D.eq.1) then

        ! XY plane
        open(14,file='DATA/DATA_2D/Ex_xy.mco')
        write(14,*) n(1)/every,n(2)/every
        do iy=n(2)+1,1,-1*every
           write(14,103) ( avg3D(1,ix,iy,iz_pl)/real(cnt_avg(1)), ix=1,n(1)+1,every ) 
        enddo
        close(14)
        
        open(14,file='DATA/DATA_2D/Ey_xy.mco')
        write(14,*) n(1)/every,n(2)/every
        do iy=n(2)+1,1,-1*every
           write(14,103) ( avg3D(2,ix,iy,iz_pl)/real(cnt_avg(1)), ix=1,n(1)+1,every ) 
        enddo
        close(14)
        
        open(14,file='DATA/DATA_2D/Ez_xy.mco')
        write(14,*) n(1)/every,n(2)/every
        do iy=n(2)+1,1,-1*every
           write(14,103) ( avg3D(3,ix,iy,iz_pl)/real(cnt_avg(1)), ix=1,n(1)+1,every ) 
        enddo
        close(14)
        
        ! XZ plane
        open(14,file='DATA/DATA_2D/Ex_xz.mco')
        write(14,*) n(1)/every,n(3)/every
        do iz=n(3)+1,1,-1*every
           write(14,103) ( avg3D(1,ix,n(2)/2+1,iz)/real(cnt_avg(1)), ix=1,n(1)+1,every ) 
        enddo
        close(14)

        open(14,file='DATA/DATA_2D/Ey_xz.mco')
        write(14,*) n(1)/every,n(3)/every
        do iz=n(3)+1,1,-1*every
           write(14,103) ( avg3D(2,ix,n(2)/2+1,iz)/real(cnt_avg(1)), ix=1,n(1)+1,every ) 
        enddo
        close(14)

        open(14,file='DATA/DATA_2D/Ez_xz.mco')
        write(14,*) n(1)/every,n(3)/every
        do iz=n(3)+1,1,-1*every
           write(14,103) ( avg3D(3,ix,n(2)/2+1,iz)/real(cnt_avg(1)), ix=1,n(1)+1,every ) 
        enddo
        close(14)
        
        ! YZ plane
        open(14,file='DATA/DATA_2D/Ex_yz.mco')
        write(14,*) n(2)/every,n(3)/every
        do iz=n(3)+1,1,-1*every
           write(14,103) ( avg3D(1,ix_pl,iy,iz)/real(cnt_avg(1)), iy=1,n(2)+1,every ) 
        enddo
        close(14)

        open(14,file='DATA/DATA_2D/Ey_yz.mco')
        write(14,*) n(2)/every,n(3)/every
        do iz=n(3)+1,1,-1*every
           write(14,103) ( avg3D(2,ix_pl,iy,iz)/real(cnt_avg(1)), iy=1,n(2)+1,every ) 
        enddo
        close(14)

        open(14,file='DATA/DATA_2D/Ez_yz.mco')
        write(14,*) n(2)/every,n(3)/every
        do iz=n(3)+1,1,-1*every
           write(14,103) ( avg3D(3,ix_pl,iy,iz)/real(cnt_avg(1)), iy=1,n(2)+1,every ) 
        enddo
        close(14)
        
     endif
     
     ! Write density, temperature, source term and currents for each species
     do ptype=1,ntype

        ! XY plane
        open(15,file='DATA/DATA_2D/'//name(24*ptype-23))
        open(16,file='DATA/DATA_2D/'//name(24*ptype-20))
        if(plt_src.eq.0 .or. (plt_src.eq.1 .and. ptype.gt.1)) &
             open(17,file='DATA/DATA_2D/'//name(24*ptype-17))
        open(18,file='DATA/DATA_2D/'//name(24*ptype-14))
        if(plt_src.eq.1 .and. ptype.gt.1) &
             open(19,file='DATA/DATA_2D/'//name(24*ptype-11))
        open(180,file='DATA/DATA_2D/'//name(24*ptype-8))
        open(181,file='DATA/DATA_2D/'//name(24*ptype-5))
        open(182,file='DATA/DATA_2D/'//name(24*ptype-2))
        
        write(15,*) n(1)/every,n(2)/every
        write(16,*) n(1)/every,n(2)/every
        if(plt_src.eq.0 .or. (plt_src.eq.1 .and. ptype.gt.1)) &
             write(17,*) n(1)/every,n(2)/every
        write(18,*) n(1)/every,n(2)/every
        if(plt_src.eq.1 .and. ptype.gt.1) &
             write(19,*) n(1)/every,n(2)/every
        write(180,*) n(1)/every,n(2)/every
        write(181,*) n(1)/every,n(2)/every
        write(182,*) n(1)/every,n(2)/every
        do iy=n(2)+1,1,-1*every
           write(15,103) ( data_pavg_xy(np_avg,ix,iy,ptype)/real(cnt_avg(1)), ix=1,n(1)+1,every )
           write(16,103) ( data_pavg_xy(Tp_avg,ix,iy,ptype)/real(cnt_avg(1)), ix=1,n(1)+1,every )
           if(plt_src.eq.0) then
              write(17,103) ( data_pavg_xy(sour_avg,ix,iy,ptype)/real(cnt_avg(3)*dt), ix=1,n(1)+1,every )
           else
              if(ptype.gt.1) &
                   write(17,103) ( data_pavg_xy(sour_avg,ix,iy,ptype)/real(cnt_avg(1)), ix=1,n(1)+1,every )
           endif
           write(18,103) ( data_pavg_xy(np_avg,ix,iy,ptype)*dsqrt( data_pavg_xy(u1_avg,ix,iy,ptype)**2 + &
                data_pavg_xy(u2_avg,ix,iy,ptype)**2 + &
                data_pavg_xy(u3_avg,ix,iy,ptype)**2 )/real(cnt_avg(1)), ix=1,n(1)+1,every )
           if(plt_src.eq.1 .and. ptype.gt.1) &
                write(19,103) ( data_pavg_xy(sink_avg,ix,iy,ptype)/real(cnt_avg(1)), ix=1,n(1)+1,every )
           write(180,103) ( data_pavg_xy(u1_avg,ix,iy,ptype), ix=1,n(1)+1,every )
           write(181,103) ( data_pavg_xy(u2_avg,ix,iy,ptype), ix=1,n(1)+1,every )
           write(182,103) ( data_pavg_xy(u3_avg,ix,iy,ptype), ix=1,n(1)+1,every )
        enddo
        
        write(18,*) 'vector'
        do iy=n(2)+1,1,-1*every ! Theta= arctg(jy/jx)
           write(18,103) ( datan2(data_pavg_xy(u2_avg,ix,iy,ptype),&
                data_pavg_xy(u1_avg,ix,iy,ptype)), ix=1,n(1)+1,every )  
        enddo
        
        close(15)
        close(16)
        if(plt_src.eq.0 .or. (plt_src.eq.1 .and. ptype.gt.1)) close(17)
        close(18)
        if(plt_src.eq.1 .and. ptype.gt.1) close(19)
        close(180)
        close(181)
        close(182)
        
        ! XZ plane
        open(15,file='DATA/DATA_2D/'//name(24*ptype-22))
        open(16,file='DATA/DATA_2D/'//name(24*ptype-19))
        if(plt_src.eq.0 .or. (plt_src.eq.1 .and. ptype.gt.1)) &
             open(17,file='DATA/DATA_2D/'//name(24*ptype-16))
        open(18,file='DATA/DATA_2D/'//name(24*ptype-13))
        if(plt_src.eq.1 .and. ptype.gt.1) &
             open(19,file='DATA/DATA_2D/'//name(24*ptype-10))
        open(180,file='DATA/DATA_2D/'//name(24*ptype-7))
        open(181,file='DATA/DATA_2D/'//name(24*ptype-4))
        open(182,file='DATA/DATA_2D/'//name(24*ptype-1))        
        
        write(15,*) n(1)/every,n(3)/every
        write(16,*) n(1)/every,n(3)/every
        if(plt_src.eq.0 .or. (plt_src.eq.1 .and. ptype.gt.1)) &
             write(17,*) n(1)/every,n(3)/every
        write(18,*) n(1)/every,n(3)/every
        if(plt_src.eq.1 .and. ptype.gt.1) &
             write(19,*) n(1)/every,n(3)/every
        write(180,*) n(1)/every,n(3)/every
        write(181,*) n(1)/every,n(3)/every
        write(182,*) n(1)/every,n(3)/every
        do iz=n(3)+1,1,-1*every
           write(15,103) ( data_pavg_xz(np_avg,ix,iz,ptype)/real(cnt_avg(1)), ix=1,n(1)+1,every )
           write(16,103) ( data_pavg_xz(Tp_avg,ix,iz,ptype)/real(cnt_avg(1)), ix=1,n(1)+1,every )
           if(plt_src.eq.0) then
              write(17,103) ( data_pavg_xz(sour_avg,ix,iz,ptype)/real(cnt_avg(3)*dt), ix=1,n(1)+1,every )
           else
              if(ptype.gt.1) &
                   write(17,103) ( data_pavg_xz(sour_avg,ix,iz,ptype)/real(cnt_avg(1)), ix=1,n(1)+1,every )
           endif
           write(18,103) ( data_pavg_xz(np_avg,ix,iz,ptype)*dsqrt( data_pavg_xz(u1_avg,ix,iz,ptype)**2 + &
                data_pavg_xz(u2_avg,ix,iz,ptype)**2 + &
                data_pavg_xz(u3_avg,ix,iz,ptype)**2 )/real(cnt_avg(1)), ix=1,n(1)+1,every )
           if(plt_src.eq.1 .and. ptype.gt.1) &
                write(19,103) ( data_pavg_xz(sink_avg,ix,iz,ptype)/real(cnt_avg(1)), ix=1,n(1)+1,every )
           write(180,103) ( data_pavg_xz(u1_avg,ix,iz,ptype), ix=1,n(1)+1,every )
           write(181,103) ( data_pavg_xz(u2_avg,ix,iz,ptype), ix=1,n(1)+1,every )
           write(182,103) ( data_pavg_xz(u3_avg,ix,iz,ptype), ix=1,n(1)+1,every )
        enddo
        
        write(18,*) 'vector'
        do iz=n(3)+1,1,-1*every ! Theta= arctg(jz/jx)
           write(18,103) ( datan2(data_pavg_xz(u3_avg,ix,iz,ptype),&
                data_pavg_xz(u1_avg,ix,iz,ptype)), ix=1,n(1)+1,every )  
        enddo
        
        close(15)
        close(16)
        if(plt_src.eq.0 .or. (plt_src.eq.1 .and. ptype.gt.1)) close(17)
        close(18)
        if(plt_src.eq.1 .and. ptype.gt.1) close(19)
        close(180)
        close(181)
        close(182)

        ! YZ plane
        open(15,file='DATA/DATA_2D/'//name(24*ptype-21))
        open(16,file='DATA/DATA_2D/'//name(24*ptype-18))
        if(plt_src.eq.0 .or. (plt_src.eq.1 .and. ptype.gt.1)) &
             open(17,file='DATA/DATA_2D/'//name(24*ptype-15))
        open(18,file='DATA/DATA_2D/'//name(24*ptype-12))
        if(plt_src.eq.1 .and. ptype.gt.1) &
             open(19,file='DATA/DATA_2D/'//name(24*ptype-9))
        open(180,file='DATA/DATA_2D/'//name(24*ptype-6))
        open(181,file='DATA/DATA_2D/'//name(24*ptype-3))
        open(182,file='DATA/DATA_2D/'//name(24*ptype))
        
        write(15,*) n(2)/every,n(3)/every
        write(16,*) n(2)/every,n(3)/every
        if(plt_src.eq.0 .or. (plt_src.eq.1 .and. ptype.gt.1)) &
             write(17,*) n(2)/every,n(3)/every
        write(18,*) n(2)/every,n(3)/every
        if(plt_src.eq.1 .and. ptype.gt.1) &
             write(19,*) n(2)/every,n(3)/every
        write(180,*) n(2)/every,n(3)/every
        write(181,*) n(2)/every,n(3)/every
        write(182,*) n(2)/every,n(3)/every
        do iz=n(3)+1,1,-1*every
           write(15,103) ( data_pavg_yz(np_avg,iy,iz,ptype)/real(cnt_avg(1)), iy=1,n(2)+1,every )
           write(16,103) ( data_pavg_yz(Tp_avg,iy,iz,ptype)/real(cnt_avg(1)), iy=1,n(2)+1,every )
           if(plt_src.eq.0) then
              write(17,103) ( data_pavg_yz(sour_avg,iy,iz,ptype)/real(cnt_avg(3)*dt), iy=1,n(2)+1,every )
           else
              if(ptype.gt.1) &
                   write(17,103) ( data_pavg_yz(sour_avg,iy,iz,ptype)/real(cnt_avg(1)), iy=1,n(2)+1,every )
           endif
           write(18,103) ( data_pavg_yz(np_avg,iy,iz,ptype)*dsqrt( data_pavg_yz(u1_avg,iy,iz,ptype)**2 + &
                data_pavg_yz(u2_avg,iy,iz,ptype)**2 + &
                data_pavg_yz(u3_avg,iy,iz,ptype)**2 )/real(cnt_avg(1)), iy=1,n(2)+1,every )
           if(plt_src.eq.1 .and. ptype.gt.1) &
                write(19,103) ( data_pavg_yz(sink_avg,iy,iz,ptype)/real(cnt_avg(1)), iy=1,n(2)+1,every )
           write(180,103) ( data_pavg_yz(u1_avg,iy,iz,ptype), iy=1,n(2)+1,every )
           write(181,103) ( data_pavg_yz(u2_avg,iy,iz,ptype), iy=1,n(2)+1,every )
           write(182,103) ( data_pavg_yz(u3_avg,iy,iz,ptype), iy=1,n(2)+1,every )
        enddo
        
        write(18,*) 'vector'
        do iz=n(3)+1,1,-1*every ! Theta= arctg(jz/jy)
           write(18,103) ( datan2(data_pavg_yz(u3_avg,iy,iz,ptype),&
                data_pavg_yz(u2_avg,iy,iz,ptype)), iy=1,n(2)+1,every )  
        enddo
        
        close(15)
        close(16)
        if(plt_src.eq.0 .or. (plt_src.eq.1 .and. ptype.gt.1)) close(17)
        close(18)
        if(plt_src.eq.1 .and. ptype.gt.1) close(19)
        close(180)
        close(181)
        close(182)
        
     enddo ! end loop over ptype
     
     !
     ! Write space integrated data
     !
     kn=1
     if(flag_nmn.eq.1) kn=0
     do ptype=1,ntype
        
        open(22,file='DATA/'//name0D(2*ptype-1),access='APPEND')
        open(23,file='DATA/'//name0D(2*ptype),access='APPEND')
        
        write(22,101) time,( p_mac_tmp(ptype,np_loss,igrid), igrid=kn,ngrid )
        write(23,101) time,( p_mac_tmp(ptype,P_w,igrid), igrid=kn,ngrid )

        close(22)
        close(23)
             
     enddo

     ! Save some local values of phi, ne and Te
     open(24,file='DATA/phi_Te_ne_cntr.dat',access='APPEND')
     if( it.eq.1 .and. flag_restart.eq.0 ) write(24,'("# Time(s), phi(xm/2,ym/2,z_pl), Te(xm/2,ym/2,z_pl), ne at ym/2, z_pl and &
          xm/2, xm/10, 2xm/10, 3xm/10, 4xm/10, 6xm/10, 7xm/10, 8xm/10, 9xm/10 ")')
     write(24,101) time, phi_avg_xy(n(1)/2+1,n(2)/2+1)/real(cnt_avg(1)),&
          data_pavg_xy(Tp_avg,n(1)/2+1,n(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg_xy(np_avg,n(1)/2+1,n(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg_xy(np_avg,NINT(n(1)/10.)+1,n(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg_xy(np_avg,NINT(2*n(1)/10.)+1,n(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg_xy(np_avg,NINT(3*n(1)/10.)+1,n(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg_xy(np_avg,NINT(4*n(1)/10.)+1,n(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg_xy(np_avg,NINT(6*n(1)/10.)+1,n(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg_xy(np_avg,NINT(7*n(1)/10.)+1,n(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg_xy(np_avg,NINT(8*n(1)/10.)+1,n(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg_xy(np_avg,NINT(9*n(1)/10.)+1,n(2)/2+1,1)/real(cnt_avg(1))
     close(24)
          
  endif
  
  deallocate( ld_xy, ld_xz, ld_yz )

101 format(20(1x,es16.8))
102 format(' Pwall (W)= ',es10.2,', Pabs (W)= ',es10.2,', Pcoll (W)= ',es10.2,&
         ', Pinj (W)= ',es10.2,', I_w (A)= ',2(es10.2,2x),', Ek_w (eV)= ',2(es10.2,2x))
103 format(800(e18.6,1x)) 
     
  return

end subroutine write_data
