subroutine part_moments(n,h,vxp,nmax,ntype,kq,nproc,np_tot, &
     iproc,ptype,p_mts_xy,p_mts_xz,p_mts_yz,n_mts)
!     ==============================================================
!     VERSION:         0.4
!     LAST MOD:      Nov/22
!     MOD AUTHOR
!     NOTE:        Indexes in vxp(): 1== x
!                                    2== y    
!                                    3== z
!                                    4== vx
!                                    5== vy
!                                    6== vz                
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  integer:: ix,iy,iz,i,ptype,n(3),nmax,ntype,iproc,nproc,n_mts, &
       Ekg,j1,j2,j3
  parameter ( Ekg=1, j1=2, j2=3, j3=4 ) 
  real(kind=8):: h(3),ki(4),k,k_tmp
  ! Particle arrays
  real(kind=8):: vxp(6,nmax,ntype,nproc),xp,yp,zp,vx,&
       vy,vz,Eki,px,py,pz,p,kq(0:n(1)+2,0:n(2)+2,0:n(3)+2),&
       p_mts_xy(n_mts,0:n(1)+2,0:n(2)+2,ntype,nproc), &
       p_mts_xz(n_mts,0:n(1)+2,0:n(3)+2,ntype,nproc), &
       p_mts_yz(n_mts,0:n(2)+2,0:n(3)+2,ntype,nproc)
  ! Macroscopic parameters
  integer:: np_tot(ntype,nproc)

  ! Initialization
  p_mts_xy(:,:,:,ptype,iproc)=0.d0 
  p_mts_xz(:,:,:,ptype,iproc)=0.d0 
  p_mts_yz(:,:,:,ptype,iproc)=0.d0 
  k=Nm(ptype)/(h(1)*h(2)*h(3))

  ! Loop over ptype particles
  do i=1,np_tot(ptype,iproc)

     vx= vxp(4,i,ptype,iproc)
     vy= vxp(5,i,ptype,iproc)
     vz= vxp(6,i,ptype,iproc)        
     xp= vxp(1,i,ptype,iproc) - vx*dt/2.d0
     yp= vxp(2,i,ptype,iproc) - vy*dt/2.d0
     zp= vxp(3,i,ptype,iproc) - vz*dt/2.d0

     ! Neumann BCs (LHS only)
     if(flag_nmn.eq.1) then
        if( xp.le.0.d0 ) then
           ! Specular reflection
           xp= -xp
           vx= -vx
        endif
     endif

     ! Periodic BCs
     if(flag_pbc.eq.1) then
        if( yp.ge.ymax ) then 
           yp= yp - ymax           
        endif
        if( yp.le.0.d0 ) then 
           yp= ymax + yp           
        endif
        if( zp.ge.zmax ) then 
           zp= zp - zmax           
        endif
        if( zp.le.0.d0 ) then 
           zp= zmax + zp           
        endif
     endif

     if( xp.gt.xmax .or. xp.lt.0.d0 .or. &
          yp.gt.ymax .or. yp.lt.0.d0 .or. &
          zp.gt.zmax .or. zp.lt.0.d0 ) goto 100

     ix= INT( xp/h(1) ) + 1
     iy= INT( yp/h(2) ) + 1
     iz= INT( zp/h(3) ) + 1

     ! Calculate moments only for particles within mid-planes
     if( ix.eq.ix_pl .or. ix.eq.(ix_pl-1) .or.  &
          iy.eq.(n(2)/2+1) .or. iy.eq.(n(2)/2) .or. &
          iz.eq.iz_pl .or. iz.eq.(iz_pl-1) ) then
        
        Eki= vx*vx + vy*vy + vz*vz

        ! Particle fluxes (trilinear interpolation)
        px=( ix*h(1) - xp )/h(1)
        py=( iy*h(2) - yp )/h(2)
        pz=( iz*h(3) - zp )/h(3)
           
        ! XZ plane
        if( iy.eq.(n(2)/2+1) .or. iy.eq.(n(2)/2) ) then 

           if(iy.eq.(n(2)/2+1)) then
              p=py
           else
              p=1-py
           endif
           ki(1)= k*px*pz*p
           ki(2)= k*(1.d0-px)*pz*p
           ki(3)= k*px*(1.d0-pz)*p
           ki(4)= k*(1.d0-px)*(1.d0-pz)*p

           ! Grid node 000
           k_tmp= kq(ix,iy,iz)*ki(1)
           p_mts_xz(j1,ix,iz,ptype,iproc)= p_mts_xz(j1,ix,iz,ptype,iproc) + k_tmp*vx
           p_mts_xz(j2,ix,iz,ptype,iproc)= p_mts_xz(j2,ix,iz,ptype,iproc) + k_tmp*vy
           p_mts_xz(j3,ix,iz,ptype,iproc)= p_mts_xz(j3,ix,iz,ptype,iproc) + k_tmp*vz
           p_mts_xz(Ekg,ix,iz,ptype,iproc)= p_mts_xz(Ekg,ix,iz,ptype,iproc) + k_tmp*Eki
           
           ! Grid node +00
           k_tmp= kq(ix+1,iy,iz)*ki(2)
           p_mts_xz(j1,ix+1,iz,ptype,iproc)= p_mts_xz(j1,ix+1,iz,ptype,iproc) + k_tmp*vx
           p_mts_xz(j2,ix+1,iz,ptype,iproc)= p_mts_xz(j2,ix+1,iz,ptype,iproc) + k_tmp*vy
           p_mts_xz(j3,ix+1,iz,ptype,iproc)= p_mts_xz(j3,ix+1,iz,ptype,iproc) + k_tmp*vz
           p_mts_xz(Ekg,ix+1,iz,ptype,iproc)= p_mts_xz(Ekg,ix+1,iz,ptype,iproc) + k_tmp*Eki
           
           ! Grid node 00+
           k_tmp= kq(ix,iy,iz+1)*ki(3)
           p_mts_xz(j1,ix,iz+1,ptype,iproc)= p_mts_xz(j1,ix,iz+1,ptype,iproc) + k_tmp*vx
           p_mts_xz(j2,ix,iz+1,ptype,iproc)= p_mts_xz(j2,ix,iz+1,ptype,iproc) + k_tmp*vy
           p_mts_xz(j3,ix,iz+1,ptype,iproc)= p_mts_xz(j3,ix,iz+1,ptype,iproc) + k_tmp*vz
           p_mts_xz(Ekg,ix,iz+1,ptype,iproc)= p_mts_xz(Ekg,ix,iz+1,ptype,iproc) + k_tmp*Eki
           
           ! Grid node +0+
           k_tmp= kq(ix+1,iy,iz+1)*ki(4)
           p_mts_xz(j1,ix+1,iz+1,ptype,iproc)= p_mts_xz(j1,ix+1,iz+1,ptype,iproc) + k_tmp*vx
           p_mts_xz(j2,ix+1,iz+1,ptype,iproc)= p_mts_xz(j2,ix+1,iz+1,ptype,iproc) + k_tmp*vy
           p_mts_xz(j3,ix+1,iz+1,ptype,iproc)= p_mts_xz(j3,ix+1,iz+1,ptype,iproc) + k_tmp*vz
           p_mts_xz(Ekg,ix+1,iz+1,ptype,iproc)= p_mts_xz(Ekg,ix+1,iz+1,ptype,iproc) + k_tmp*Eki
           
        endif

        ! XY plane
        if( iz.eq.iz_pl .or. iz.eq.(iz_pl-1) ) then

           if(iz.eq.iz_pl) then
              p=pz
           else
              p=1-pz
           endif
           ki(1)= k*px*py*p
           ki(2)= k*(1.d0-px)*py*p
           ki(3)= k*(1.d0-px)*(1.d0-py)*p
           ki(4)= k*px*(1.d0-py)*p

           ! Grid node 000
           k_tmp= kq(ix,iy,iz)*ki(1)
           p_mts_xy(j1,ix,iy,ptype,iproc)= p_mts_xy(j1,ix,iy,ptype,iproc) + k_tmp*vx
           p_mts_xy(j2,ix,iy,ptype,iproc)= p_mts_xy(j2,ix,iy,ptype,iproc) + k_tmp*vy
           p_mts_xy(j3,ix,iy,ptype,iproc)= p_mts_xy(j3,ix,iy,ptype,iproc) + k_tmp*vz
           p_mts_xy(Ekg,ix,iy,ptype,iproc)= p_mts_xy(Ekg,ix,iy,ptype,iproc) + k_tmp*Eki

           ! Grid node +00
           k_tmp= kq(ix+1,iy,iz)*ki(2)
           p_mts_xy(j1,ix+1,iy,ptype,iproc)= p_mts_xy(j1,ix+1,iy,ptype,iproc) + k_tmp*vx
           p_mts_xy(j2,ix+1,iy,ptype,iproc)= p_mts_xy(j2,ix+1,iy,ptype,iproc) + k_tmp*vy
           p_mts_xy(j3,ix+1,iy,ptype,iproc)= p_mts_xy(j3,ix+1,iy,ptype,iproc) + k_tmp*vz
           p_mts_xy(Ekg,ix+1,iy,ptype,iproc)= p_mts_xy(Ekg,ix+1,iy,ptype,iproc) + k_tmp*Eki

           ! Grid node ++0
           k_tmp= kq(ix+1,iy+1,iz)*ki(3)
           p_mts_xy(j1,ix+1,iy+1,ptype,iproc)= p_mts_xy(j1,ix+1,iy+1,ptype,iproc) + k_tmp*vx
           p_mts_xy(j2,ix+1,iy+1,ptype,iproc)= p_mts_xy(j2,ix+1,iy+1,ptype,iproc) + k_tmp*vy
           p_mts_xy(j3,ix+1,iy+1,ptype,iproc)= p_mts_xy(j3,ix+1,iy+1,ptype,iproc) + k_tmp*vz
           p_mts_xy(Ekg,ix+1,iy+1,ptype,iproc)= p_mts_xy(Ekg,ix+1,iy+1,ptype,iproc) + k_tmp*Eki

           ! Grid node 0+0
           k_tmp= kq(ix,iy+1,iz)*ki(4)
           p_mts_xy(j1,ix,iy+1,ptype,iproc)= p_mts_xy(j1,ix,iy+1,ptype,iproc) + k_tmp*vx
           p_mts_xy(j2,ix,iy+1,ptype,iproc)= p_mts_xy(j2,ix,iy+1,ptype,iproc) + k_tmp*vy
           p_mts_xy(j3,ix,iy+1,ptype,iproc)= p_mts_xy(j3,ix,iy+1,ptype,iproc) + k_tmp*vz
           p_mts_xy(Ekg,ix,iy+1,ptype,iproc)= p_mts_xy(Ekg,ix,iy+1,ptype,iproc) + k_tmp*Eki

        endif

        ! YZ plane
        if( ix.eq.ix_pl .or. ix.eq.(ix_pl-1) ) then
           
           if(ix.eq.ix_pl) then
              p=px
           else
              p=1-px
           endif
           ki(1)= k*py*pz*p
           ki(2)= k*(1.d0-py)*pz*p
           ki(3)= k*py*(1.d0-pz)*p
           ki(4)= k*(1.d0-py)*(1.d0-pz)*p

           ! Grid node 000
           k_tmp= kq(ix,iy,iz)*ki(1)
           p_mts_yz(j1,iy,iz,ptype,iproc)= p_mts_yz(j1,iy,iz,ptype,iproc) + k_tmp*vx
           p_mts_yz(j2,iy,iz,ptype,iproc)= p_mts_yz(j2,iy,iz,ptype,iproc) + k_tmp*vy
           p_mts_yz(j3,iy,iz,ptype,iproc)= p_mts_yz(j3,iy,iz,ptype,iproc) + k_tmp*vz
           p_mts_yz(Ekg,iy,iz,ptype,iproc)= p_mts_yz(Ekg,iy,iz,ptype,iproc) + k_tmp*Eki
           
           ! Grid node +00
           k_tmp= kq(ix,iy+1,iz)*ki(2)
           p_mts_yz(j1,iy+1,iz,ptype,iproc)= p_mts_yz(j1,iy+1,iz,ptype,iproc) + k_tmp*vx
           p_mts_yz(j2,iy+1,iz,ptype,iproc)= p_mts_yz(j2,iy+1,iz,ptype,iproc) + k_tmp*vy
           p_mts_yz(j3,iy+1,iz,ptype,iproc)= p_mts_yz(j3,iy+1,iz,ptype,iproc) + k_tmp*vz
           p_mts_yz(Ekg,iy+1,iz,ptype,iproc)= p_mts_yz(Ekg,iy+1,iz,ptype,iproc) + k_tmp*Eki
           
           ! Grid node 00+
           k_tmp= kq(ix,iy,iz+1)*ki(3)
           p_mts_yz(j1,iy,iz+1,ptype,iproc)= p_mts_yz(j1,iy,iz+1,ptype,iproc) + k_tmp*vx
           p_mts_yz(j2,iy,iz+1,ptype,iproc)= p_mts_yz(j2,iy,iz+1,ptype,iproc) + k_tmp*vy
           p_mts_yz(j3,iy,iz+1,ptype,iproc)= p_mts_yz(j3,iy,iz+1,ptype,iproc) + k_tmp*vz
           p_mts_yz(Ekg,iy,iz+1,ptype,iproc)= p_mts_yz(Ekg,iy,iz+1,ptype,iproc) + k_tmp*Eki
           
           ! Grid node +0+
           k_tmp= kq(ix,iy+1,iz+1)*ki(4)
           p_mts_yz(j1,iy+1,iz+1,ptype,iproc)= p_mts_yz(j1,iy+1,iz+1,ptype,iproc) + k_tmp*vx
           p_mts_yz(j2,iy+1,iz+1,ptype,iproc)= p_mts_yz(j2,iy+1,iz+1,ptype,iproc) + k_tmp*vy
           p_mts_yz(j3,iy+1,iz+1,ptype,iproc)= p_mts_yz(j3,iy+1,iz+1,ptype,iproc) + k_tmp*vz
           p_mts_yz(Ekg,iy+1,iz+1,ptype,iproc)= p_mts_yz(Ekg,iy+1,iz+1,ptype,iproc) + k_tmp*Eki
           
        endif

     endif

100 enddo

  return
end subroutine part_moments

subroutine calc_avg(n,h,np,p_mts_xy,p_mts_xz,p_mts_yz,data_pavg_xy,&
     data_pavg_xz,data_pavg_yz,sour_xy,sour_xz,sour_yz,ntype,n_mts,&
      ss2D_xy,ss2D_xz,ss2D_yz,nproc,nproc_mpi)
!     ==============================================================
!     VERSION:         0.3
!     LAST MOD:      Dec/23
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:  Calculate average value of particle moments
!     NOTE: 
!     --------------------------------------------------------------
  use mpi
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: ix,iy,iz,ptype,n(3),iproc,nproc,ntype,n_mts,nproc_mpi,ierr
  integer:: np_avg,Tp_avg,j1_avg,j2_avg,j3_avg,dummy,sour_avg,sink_avg
  parameter ( np_avg=1, Tp_avg=2, j1_avg=3, j2_avg=4, j3_avg=5, dummy=6,&
       sour_avg=7, sink_avg=8 ) 
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc), &
       p_mts_xy(n_mts,0:n(1)+2,0:n(2)+2,ntype,nproc),&
       p_mts_xz(n_mts,0:n(1)+2,0:n(3)+2,ntype,nproc),&
       p_mts_yz(n_mts,0:n(2)+2,0:n(3)+2,ntype,nproc),&
       data_pavg_xy(8,0:n(1)+2,0:n(2)+2,ntype),&
       data_pavg_xz(8,0:n(1)+2,0:n(3)+2,ntype),&
       data_pavg_yz(8,0:n(2)+2,0:n(3)+2,ntype),&
       ss2D_xy(2,0:n(1)+2,0:n(2)+2,ntype,nproc),&
       ss2D_xz(2,0:n(1)+2,0:n(3)+2,ntype,nproc),&
       ss2D_yz(2,0:n(2)+2,0:n(3)+2,ntype,nproc)
  real(kind=8):: h(3),k,np_tmp,jx_tmp,jy_tmp,jz_tmp,Tp_tmp
  integer:: sour_xy(0:n(1)+2,0:n(2)+2,ntype,nproc), &
       sour_xz(0:n(1)+2,0:n(3)+2,ntype,nproc),&
       sour_yz(0:n(2)+2,0:n(3)+2,ntype,nproc)
  real(kind=8),allocatable:: data_pavg_tmp_xy(:,:,:,:),data_pavg_tmp_xz(:,:,:,:),&
       data_pavg_tmp_yz(:,:,:,:),data_pavg_tmp(:,:,:,:)

  allocate( data_pavg_tmp_xy(8,0:n(1)+2,0:n(2)+2,ntype), &
       data_pavg_tmp_xz(8,0:n(1)+2,0:n(3)+2,ntype), &
       data_pavg_tmp_yz(8,0:n(2)+2,0:n(3)+2,ntype) )
  
  ! Initialize
  data_pavg_tmp_xy=0.d0
  data_pavg_tmp_xz=0.d0
  data_pavg_tmp_yz=0.d0
           
  !
  ! Reduction
  !
  do iproc=1,nproc
     do ptype=1,ntype

        ! Charge density element
        k=Nm(ptype)/(h(1)*h(2)*h(3))
        
        ! XY plane
        !$OMP PARALLEL
        !$OMP DO
        do iy=0,n(2)+2
           do ix=0,n(1)+2
              data_pavg_tmp_xy(np_avg,ix,iy,ptype)= data_pavg_tmp_xy(np_avg,ix,iy,ptype) + &
                   np(ix,iy,iz_pl,ptype,iproc)                           
              data_pavg_tmp_xy(j1_avg,ix,iy,ptype)= data_pavg_tmp_xy(j1_avg,ix,iy,ptype) + &
                   p_mts_xy(2,ix,iy,ptype,iproc)              
              data_pavg_tmp_xy(j2_avg,ix,iy,ptype)= data_pavg_tmp_xy(j2_avg,ix,iy,ptype) + &
                   p_mts_xy(3,ix,iy,ptype,iproc)
              data_pavg_tmp_xy(j3_avg,ix,iy,ptype)= data_pavg_tmp_xy(j3_avg,ix,iy,ptype) + &
                   p_mts_xy(4,ix,iy,ptype,iproc)
              data_pavg_tmp_xy(Tp_avg,ix,iy,ptype)= data_pavg_tmp_xy(Tp_avg,ix,iy,ptype) + &
                   p_mts_xy(1,ix,iy,ptype,iproc) ! n*<v²>
              if(plt_src.eq.0) then
                 data_pavg_tmp_xy(sour_avg,ix,iy,ptype)= data_pavg_tmp_xy(sour_avg,ix,iy,ptype) + &
                      real(sour_xy(ix,iy,ptype,iproc))*k
              else
                 data_pavg_tmp_xy(sour_avg,ix,iy,ptype)= data_pavg_tmp_xy(sour_avg,ix,iy,ptype) + &
                      ss2D_xy(1,ix,iy,ptype,iproc)
                 data_pavg_tmp_xy(sink_avg,ix,iy,ptype)= data_pavg_tmp_xy(sink_avg,ix,iy,ptype) + &
                      ss2D_xy(2,ix,iy,ptype,iproc)
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
              data_pavg_tmp_xz(np_avg,ix,iz,ptype)= data_pavg_tmp_xz(np_avg,ix,iz,ptype) + &
                   np(ix,n(2)/2+1,iz,ptype,iproc) 
              data_pavg_tmp_xz(j1_avg,ix,iz,ptype)= data_pavg_tmp_xz(j1_avg,ix,iz,ptype) + &
                   p_mts_xz(2,ix,iz,ptype,iproc)
              data_pavg_tmp_xz(j2_avg,ix,iz,ptype)= data_pavg_tmp_xz(j2_avg,ix,iz,ptype) + &
                   p_mts_xz(3,ix,iz,ptype,iproc)
              data_pavg_tmp_xz(j3_avg,ix,iz,ptype)= data_pavg_tmp_xz(j3_avg,ix,iz,ptype) + &
                   p_mts_xz(4,ix,iz,ptype,iproc)
              data_pavg_tmp_xz(Tp_avg,ix,iz,ptype)= data_pavg_tmp_xz(Tp_avg,ix,iz,ptype) + &
                   p_mts_xz(1,ix,iz,ptype,iproc) ! n*<v²>
              if(plt_src.eq.0) then
                 data_pavg_tmp_xz(sour_avg,ix,iz,ptype)= data_pavg_tmp_xz(sour_avg,ix,iz,ptype) + &
                      real(sour_xz(ix,iz,ptype,iproc))*k
              else
                 data_pavg_tmp_xz(sour_avg,ix,iz,ptype)= data_pavg_tmp_xz(sour_avg,ix,iz,ptype) + &
                      ss2D_xz(1,ix,iz,ptype,iproc)
                 data_pavg_tmp_xz(sink_avg,ix,iz,ptype)= data_pavg_tmp_xz(sink_avg,ix,iz,ptype) + &
                      ss2D_xz(2,ix,iz,ptype,iproc)
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
              data_pavg_tmp_yz(np_avg,iy,iz,ptype)= data_pavg_tmp_yz(np_avg,iy,iz,ptype) + &
                   np(ix_pl,iy,iz,ptype,iproc)
              data_pavg_tmp_yz(j1_avg,iy,iz,ptype)= data_pavg_tmp_yz(j1_avg,iy,iz,ptype) + &
                   p_mts_yz(2,iy,iz,ptype,iproc)
              data_pavg_tmp_yz(j2_avg,iy,iz,ptype)= data_pavg_tmp_yz(j2_avg,iy,iz,ptype) + &
                   p_mts_yz(3,iy,iz,ptype,iproc)
              data_pavg_tmp_yz(j3_avg,iy,iz,ptype)= data_pavg_tmp_yz(j3_avg,iy,iz,ptype) + &
                   p_mts_yz(4,iy,iz,ptype,iproc)
              data_pavg_tmp_yz(Tp_avg,iy,iz,ptype)= data_pavg_tmp_yz(Tp_avg,iy,iz,ptype) + &
                   p_mts_yz(1,iy,iz,ptype,iproc) ! n*<v²>
              if(plt_src.eq.0) then
                 data_pavg_tmp_yz(sour_avg,iy,iz,ptype)= data_pavg_tmp_yz(sour_avg,iy,iz,ptype) + &
                      real(sour_yz(iy,iz,ptype,iproc))*k
              else
                 data_pavg_tmp_yz(sour_avg,iy,iz,ptype)= data_pavg_tmp_yz(sour_avg,iy,iz,ptype) + &
                      ss2D_yz(1,iy,iz,ptype,iproc)
                 data_pavg_tmp_yz(sink_avg,iy,iz,ptype)= data_pavg_tmp_yz(sink_avg,iy,iz,ptype) + &
                      ss2D_yz(2,iy,iz,ptype,iproc)
              endif
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL  

     enddo
  enddo
  
  if(nproc_mpi.gt.1) then
     
     ! XY plane
     allocate( data_pavg_tmp(8,0:n(1)+2,0:n(2)+2,ntype) )
     data_pavg_tmp=0.d0
     call MPI_REDUCE(data_pavg_tmp_xy(:,:,:,:), data_pavg_tmp, 8*(n(1)+3)*(n(2)+3)*ntype, &
          MPI_REAL8, MPI_SUM,0,MPI_COMM_WORLD, ierr)
     data_pavg_tmp_xy= data_pavg_tmp
     deallocate( data_pavg_tmp )

     ! XZ plane
     allocate( data_pavg_tmp(8,0:n(1)+2,0:n(3)+2,ntype) )
     data_pavg_tmp=0.d0
     call MPI_REDUCE(data_pavg_tmp_xz(:,:,:,:), data_pavg_tmp, 8*(n(1)+3)*(n(3)+3)*ntype, &
          MPI_REAL8, MPI_SUM,0,MPI_COMM_WORLD, ierr)
     data_pavg_tmp_xz= data_pavg_tmp
     deallocate( data_pavg_tmp )

     ! YZ plane
     allocate( data_pavg_tmp(8,0:n(2)+2,0:n(3)+2,ntype) )
     data_pavg_tmp=0.d0
     call MPI_REDUCE(data_pavg_tmp_yz(:,:,:,:), data_pavg_tmp, 8*(n(2)+3)*(n(3)+3)*ntype, &
          MPI_REAL8, MPI_SUM,0,MPI_COMM_WORLD, ierr)
     data_pavg_tmp_yz= data_pavg_tmp
     deallocate( data_pavg_tmp )

  endif
  
  !
  ! Calculate Tp
  !
  do ptype=1,ntype
        
     ! XY plane
     !$OMP PARALLEL PRIVATE(np_tmp,jx_tmp,jy_tmp,jz_tmp,Tp_tmp)
     !$OMP DO
     do iy=0,n(2)+2
        do ix=0,n(1)+2
           np_tmp= data_pavg_tmp_xy(np_avg,ix,iy,ptype)
           data_pavg_xy(np_avg,ix,iy,ptype)= data_pavg_xy(np_avg,ix,iy,ptype) + np_tmp
           jx_tmp= data_pavg_tmp_xy(j1_avg,ix,iy,ptype)
           data_pavg_xy(j1_avg,ix,iy,ptype)= data_pavg_xy(j1_avg,ix,iy,ptype) + jx_tmp
           jy_tmp= data_pavg_tmp_xy(j2_avg,ix,iy,ptype)
           data_pavg_xy(j2_avg,ix,iy,ptype)= data_pavg_xy(j2_avg,ix,iy,ptype) + jy_tmp
           jz_tmp= data_pavg_tmp_xy(j3_avg,ix,iy,ptype)
           data_pavg_xy(j3_avg,ix,iy,ptype)= data_pavg_xy(j3_avg,ix,iy,ptype) + jz_tmp
           if(np_tmp.gt.0) then
              Tp_tmp= mass(ptype)/qe*(1.d0/3.d0)*( data_pavg_tmp_xy(Tp_avg,ix,iy,ptype)/np_tmp - &
                   (jx_tmp*jx_tmp + jy_tmp*jy_tmp + jz_tmp*jz_tmp)/(np_tmp*np_tmp) )
           else
              Tp_tmp=0.d0
           endif
           ! Issue which may occur in corners or lack of statistics
           if(Tp_tmp.ge.0.d0) data_pavg_xy(Tp_avg,ix,iy,ptype)= data_pavg_xy(Tp_avg,ix,iy,ptype) + Tp_tmp
           data_pavg_xy(sour_avg,ix,iy,ptype)= data_pavg_xy(sour_avg,ix,iy,ptype) + &
                data_pavg_tmp_xy(sour_avg,ix,iy,ptype)
           if(plt_src.eq.1) &
                data_pavg_xy(sink_avg,ix,iy,ptype)= data_pavg_xy(sink_avg,ix,iy,ptype) + &
                data_pavg_tmp_xy(sink_avg,ix,iy,ptype)
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL  
     
     ! XZ plane
     !$OMP PARALLEL PRIVATE(np_tmp,jx_tmp,jy_tmp,jz_tmp,Tp_tmp)
     !$OMP DO
     do iz=0,n(3)+2
        do ix=0,n(1)+2
           np_tmp= data_pavg_tmp_xz(np_avg,ix,iz,ptype)
           data_pavg_xz(np_avg,ix,iz,ptype)= data_pavg_xz(np_avg,ix,iz,ptype) + np_tmp
           jx_tmp= data_pavg_tmp_xz(j1_avg,ix,iz,ptype)
           data_pavg_xz(j1_avg,ix,iz,ptype)= data_pavg_xz(j1_avg,ix,iz,ptype) + jx_tmp
           jy_tmp= data_pavg_tmp_xz(j2_avg,ix,iz,ptype)
           data_pavg_xz(j2_avg,ix,iz,ptype)= data_pavg_xz(j2_avg,ix,iz,ptype) + jy_tmp
           jz_tmp= data_pavg_tmp_xz(j3_avg,ix,iz,ptype)                
           data_pavg_xz(j3_avg,ix,iz,ptype)= data_pavg_xz(j3_avg,ix,iz,ptype) + jz_tmp                
           if(np_tmp.gt.0) then
              Tp_tmp= mass(ptype)/qe*(1.d0/3.d0)*( data_pavg_tmp_xz(Tp_avg,ix,iz,ptype)/np_tmp - &
                   (jx_tmp*jx_tmp + jy_tmp*jy_tmp + jz_tmp*jz_tmp)/(np_tmp*np_tmp) )
           else
              Tp_tmp=0.d0
           endif
           ! Issue which may occur in corners or lack of statistics
           if(Tp_tmp.gt.0) data_pavg_xz(Tp_avg,ix,iz,ptype)= data_pavg_xz(Tp_avg,ix,iz,ptype) + Tp_tmp                
           data_pavg_xz(sour_avg,ix,iz,ptype)= data_pavg_xz(sour_avg,ix,iz,ptype) + &
                data_pavg_tmp_xz(sour_avg,ix,iz,ptype)
           if(plt_src.eq.1) &
                data_pavg_xz(sink_avg,ix,iz,ptype)= data_pavg_xz(sink_avg,ix,iz,ptype) + &
                data_pavg_tmp_xz(sink_avg,ix,iz,ptype)
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL  
     
     ! YZ plane
     !$OMP PARALLEL PRIVATE(np_tmp,jx_tmp,jy_tmp,jz_tmp,Tp_tmp)
     !$OMP DO
     do iz=0,n(3)+2
        do iy=0,n(2)+2
           np_tmp= data_pavg_tmp_yz(np_avg,iy,iz,ptype)
           data_pavg_yz(np_avg,iy,iz,ptype)= data_pavg_yz(np_avg,iy,iz,ptype) + np_tmp
           jx_tmp= data_pavg_tmp_yz(j1_avg,iy,iz,ptype)
           data_pavg_yz(j1_avg,iy,iz,ptype)= data_pavg_yz(j1_avg,iy,iz,ptype) + jx_tmp
           jy_tmp= data_pavg_tmp_yz(j2_avg,iy,iz,ptype)
           data_pavg_yz(j2_avg,iy,iz,ptype)= data_pavg_yz(j2_avg,iy,iz,ptype) + jy_tmp
           jz_tmp= data_pavg_tmp_yz(j3_avg,iy,iz,ptype)
           data_pavg_yz(j3_avg,iy,iz,ptype)= data_pavg_yz(j3_avg,iy,iz,ptype) + jz_tmp                
           if(np_tmp.gt.0) then
              Tp_tmp= mass(ptype)/qe*(1.d0/3.d0)* ( data_pavg_tmp_yz(Tp_avg,iy,iz,ptype)/np_tmp  - &
                   (jx_tmp*jx_tmp + jy_tmp*jy_tmp + jz_tmp*jz_tmp)/(np_tmp*np_tmp) )
           else
              Tp_tmp=0.d0
           endif
           ! Issue which may occur in corners or lack of statistics
           if(Tp_tmp.ge.0.d0) data_pavg_yz(Tp_avg,iy,iz,ptype)= data_pavg_yz(Tp_avg,iy,iz,ptype) + Tp_tmp                
           data_pavg_yz(sour_avg,iy,iz,ptype)= data_pavg_yz(sour_avg,iy,iz,ptype) + &
                data_pavg_tmp_yz(sour_avg,iy,iz,ptype)
           if(plt_src.eq.1) &
                data_pavg_yz(sink_avg,iy,iz,ptype)= data_pavg_yz(sink_avg,iy,iz,ptype) + &
                data_pavg_tmp_yz(sink_avg,iy,iz,ptype)
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL  

  enddo

  ! Deallocate temporary arrays
  deallocate( data_pavg_tmp_xy, data_pavg_tmp_xz, data_pavg_tmp_yz )
  
  return
end subroutine calc_avg

FUNCTION MSTIMER()
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Mar/10
!     MOD AUTHOR:    G. Hagelaar, G. Fubiani
!     COMMENTS:
!     NOTE: 
!     --------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(10):: cdat,ctim,czon
  INTEGER:: time(8),t_new,t_old,MSTIMER
  SAVE t_old

  CALL date_and_time(cdat,ctim,czon,time)
  t_new=((time(5)*60+time(6))*60+time(7))*1000+time(8)
  IF (t_old.EQ.0.d0) t_old=t_new
  MSTIMER=t_new-t_old
  t_old=t_new

  ! This happens at midnight ...
  if(MSTIMER.lt.0) MSTIMER=0

  RETURN

END FUNCTION MSTIMER

subroutine stop_calculation
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Oct/19
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     -------------------------------------------------------------
  use mpi
  implicit none
  integer ierr

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  STOP
  call MPI_Finalize(ierr)

  return
end subroutine stop_calculation

