subroutine part_injection(n,h,bcnd,vxp,sour_xy,sour_xz,sour_yz,ntype,nmax,&
     I_inj,xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow,np_tot,nproc,iseed,&
     nproc_mpi,N_inj,iproc,P_loss,phi,ni0)
!     ==============================================================
!     VERSION:         0.3
!     LAST MOD:      Apr/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:    This subroutine is designed to inject an
!                  electron-positive ion pair.
!     NOTE:          
!     --------------------------------------------------------------
  use mpi
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer nproc_mpi
  integer:: ix,iy,iz,i,k
  integer:: ntype,ptype,nmax,n(3),iproc,nproc,ng_sec
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2),np_tot(ntype,nproc),&
       sour_xy(0:n(1)+2,0:n(2)+2,ntype,nproc), &
       sour_xz(0:n(1)+2,0:n(3)+2,ntype,nproc),&
       sour_yz(0:n(2)+2,0:n(3)+2,ntype,nproc),&
       N_inj(ntype,nproc),N_inj_tmp
  real(kind=8):: h(3),vxp(6,nmax,ntype,nproc),x,y,z,vx,vy,vz,vt,rnd(3),&
       ran2,I_inj,xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow,dN_inj,&
       vz_sav(ntype),dx,dy,dz,xm,ym,zm,P_loss(4,ntype,nproc),Rb,vb,&
       z_tmp,dt_tmp,vmax,fmax,ni0(npart)
    real(kind=8):: phi(0:n(1)+2,0:n(2)+2,0:n(3)+2),phip,ki(8),px,py,pz
  ! Macroscopic parameters 
  integer:: iseed

  ! Initialize particle counter, variables & arrays
  vz_sav=0.d0
  xm= (xl_pow+xr_pow)/2.d0 
  dx= xr_pow-xl_pow
  ym= (yl_pow+yr_pow)/2.d0 
  dy= yr_pow-yl_pow
  zm= (zl_pow+zr_pow)/2.d0 
  dz= zr_pow-zl_pow
  Rb=MIN(dx,dy)/2.d0

  ! Option for an emissive cathode at a Z=cste plane
  ng_sec= 1 ! @ bottom 
  if(dir_sec.eq.-1) ng_sec=2 ! Top

  ! opt=3 : inject electron beam
  ! opt=4 : electron emission off a cathode along +/- Z
  vb=0.d0
  if( ABS(opt_inj).eq.3 .or. ABS(opt_inj).eq.4 ) then 
     vb= dsqrt(2.d0*qe*THm/ABS(mass(1)))
     vt= vt0(1)
     vmax= (vb + dsqrt(vb**2 + 2.d0*vt**2))/2.d0
     fmax= vmax*dexp(-(vmax-vb)**2/vt**2) 
  endif
  
  ! Inject a fixed particle current 
  if(ABS(opt_inj).ne.2) then
     ! Number of particles to inject per time step per OMP and MPI threads
     dN_inj= I_inj*ns_inj*dt/(qe*Nm(1))/real(nproc_mpi)/real(nproc)
     ! Correct for round-off errors
     N_inj_tmp= INT(dN_inj)
     rnd(1)= ran2(iseed)
     if( rnd(1).le.(dN_inj-N_inj_tmp) ) N_inj_tmp= N_inj_tmp + 1
  else ! Inject a number of electrons equal to ion wall losses
     N_inj_tmp= N_inj(2,iproc)
  endif
  
  ! Inject particles
  do i=1,N_inj_tmp
                   
70   rnd(1)=ran2(iseed)
     rnd(2)=ran2(iseed)
     rnd(3)=ran2(iseed)

     if(opt_inj.gt.0) then
        ! Random location inside [x<,x>], [0,ymax], [0,zmax]
        x= -rnd(1)*dx + xr_pow
        y= -rnd(2)*dy + yr_pow
        z= -rnd(3)*dz + zr_pow
     else ! Cosine & flattop distributions
        ! Source term: S(x,y)= S0*cos[pi*(x-xm)/dx]*cos[pi*(y-ym)/dy]
        ! xm= (x1+x2)/2, dx= x2-x1, ym=ymax/2, dy=ymax (y1=0)
        ! rnd= int^x_x1 int^y_y1 S dxdy/ int^x2_x1 int^y2_y1 S dxdy
        x= xm + (dx/pi)*dasin(2.d0*rnd(1)-1.d0)
        if( flag_pbc.eq.0 .or. (flag_pbc.eq.1 .and. flag_pbcz.eq.1) ) then
           y= ym + (dy/pi)*dasin(2.d0*rnd(2)-1.d0)
        else
           y= -rnd(2)*dy + yr_pow
        endif
        if( flag_pbc.eq.0 ) then
           z= zm + (zmax/pi)*dasin(2.d0*rnd(3)-1.d0)
        else
           z= -rnd(3)*dz + zr_pow
        endif
     endif
     if( ABS(opt_inj).eq.3 .or. ABS(opt_inj).eq.4 ) then
        ! Load a disk
        if( ((x-xm)**2 + (y-ym)**2).gt.Rb**2 ) goto 70       
     endif
        
     ! Get particle left grid index
     ix= FLOOR( x/h(1) ) + 1
     iy= FLOOR( y/h(2) ) + 1
     iz= FLOOR( z/h(3) ) + 1
     
     ! Load uniquely inside simulation domain
     if( bcnd(ix,iy,iz).ge.1 .and. bcnd(ix+1,iy,iz).ge.1 .and. &
          bcnd(ix+1,iy+1,iz).ge.1 .and. bcnd(ix,iy+1,iz).ge.1 .and. & 
          bcnd(ix,iy,iz+1).ge.1 .and. bcnd(ix+1,iy,iz+1).ge.1 .and. &
          bcnd(ix+1,iy+1,iz+1).ge.1 .and. bcnd(ix,iy+1,iz+1).ge.1 ) goto 70
     
     ! Loop over ptype particles
     do ptype= 1,ntype

        ! e-beam or electron emission off a cathode along the Z-axis
        if( ABS(opt_inj).eq.3 .or. ABS(opt_inj).eq.4 ) then
           ! Inject only electrons
           if(tag_beam.eq.0) then
              if(ptype.ge.2) goto 80
           else ! pname(tag_beam)=[eb]
              if(ptype.ne.tag_beam) goto 80
           endif

           if( n_cath.eq.2 .and. ABS(opt_inj).eq.4 ) then
              ! Cathode both at top and bottom
              ng_sec= ng_sec + 1
              ! Alternate injection location
              if(ng_sec.gt.2) ng_sec=1
              if(ng_sec.eq.1) then
                 dir_sec=1
              else
                 dir_sec=-1
              endif
           endif

           ! z-coordinate of emissive surface
           z= zg_sec(ng_sec)
        endif

        ! Inject ions 
        rnd(1)= ran2(iseed)
        if( rnd(1).gt.ni0(ptype) ) goto 80
                
        ! Add particle to counter
        np_tot(ptype,iproc)= np_tot(ptype,iproc) + 1
     
        ! Get thermal velocity
        vt= vt0(ptype)
     
        ! Particle index
        k= np_tot(ptype,iproc)
        
        ! Warning
        if(k.gt.nmax) then
           print*, 'k > nmax in part_injection'
           print*, 'Abort calculation ...'
           call stop_calculation
        endif

        ! Electron flux injection
        if(ABS(opt_inj).eq.3 .or. ABS(opt_inj).eq.4) then
           z_tmp= z

           ! e-beam
75         if(ABS(opt_inj).eq.3) then
              ! Shifted Maxwellian flux distribution with u=vB and T=Te
              call shifted_maxwellian_flux(vz,vb,vt,fmax,iseed)     
              rnd(1)= ran2(iseed)
              rnd(2)= ran2(iseed)
              call load_gauss(vx,vy,vt,rnd)
           endif

           ! Emission from a cathode
           if(ABS(opt_inj).eq.4) then
              ! Half Maxwellian flux distribution with T=Tb
              rnd(1)=ran2(iseed)
              vz = real(dir_sec)*vb*dsqrt( -dlog(1-rnd(1)) )     
              rnd(1)= ran2(iseed)
              rnd(2)= ran2(iseed)
              call load_gauss(vx,vy,vb,rnd)
           endif
           
           ! Spread the position over the distance traveled during one time step ns_inj*dt
           rnd(1)= ran2(iseed)
           dt_tmp= rnd(1)*(real(ns_inj)*dt)
           z= z_tmp + vz*dt_tmp

           ! Check bounds
           if(z.lt.0.d0 .or. z.gt.zmax) goto 75
        else
           ! Gaussian distribution
           rnd(1)= ran2(iseed)
           rnd(2)= ran2(iseed)
           call load_gauss(vx,vy,vt,rnd)
           if(vz_sav(ptype).eq.0.d0) then
              rnd(1)= ran2(iseed)
              rnd(2)= ran2(iseed)
              call load_gauss(vx,vy,vt,rnd)
              vz= vx 
              vz_sav(ptype)= vy
           else
              vz= vz_sav(ptype) 
              vz_sav(ptype)= 0.d0
           endif
        endif

        ! Save particle 6d-coordinates
        vxp(1,k,ptype,iproc)= x
        vxp(2,k,ptype,iproc)= y
        vxp(3,k,ptype,iproc)= z        
        vxp(4,k,ptype,iproc)= vx
        vxp(5,k,ptype,iproc)= vy
        vxp(6,k,ptype,iproc)= vz

        ! Source term
        ix= FLOOR( x/h(1) ) + 1
        iy= FLOOR( y/h(2) ) + 1
        iz= FLOOR( z/h(3) ) + 1
        
        if(iz.eq.iz_pl) &
             sour_xy(ix,iy,ptype,iproc)= sour_xy(ix,iy,ptype,iproc) + 1
        if(iy.eq.n(2)/2+1) &
             sour_xz(ix,iz,ptype,iproc)= sour_xz(ix,iz,ptype,iproc) + 1
        if(ix.eq.ix_pl) &
             sour_yz(iy,iz,ptype,iproc)= sour_yz(iy,iz,ptype,iproc) + 1

        ! Save injected power
        P_loss(4,ptype,iproc)= P_loss(4,ptype,iproc) + 0.5d0*Nm(ptype)*mass(ptype)*( &
             vxp(4,k,ptype,iproc)*vxp(4,k,ptype,iproc) + &
             vxp(5,k,ptype,iproc)*vxp(5,k,ptype,iproc) + &
             vxp(6,k,ptype,iproc)*vxp(6,k,ptype,iproc) )

        if(ABS(opt_inj).eq.3) then
           px=( ix*h(1) - x )/h(1)
           py=( iy*h(2) - y )/h(2)
           pz=( iz*h(3) - z )/h(3)
           
           ki(1)= px*py*pz
           ki(2)= (1.d0-px)*py*pz
           ki(3)= (1.d0-px)*(1.d0-py)*pz
           ki(4)= px*(1.d0-py)*pz
           ki(5)= px*py*(1.d0-pz)
           ki(6)= (1.d0-px)*py*(1.d0-pz)
           ki(7)= (1.d0-px)*(1.d0-py)*(1.d0-pz)
           ki(8)= px*(1.d0-py)*(1.d0-pz)
           
           phip= ki(1)*phi(ix,iy,iz) + &
                ki(2)*phi(ix+1,iy,iz) + &
                ki(3)*phi(ix+1,iy+1,iz) + &
                ki(4)*phi(ix,iy+1,iz) + &
                ki(5)*phi(ix,iy,iz+1) + &
                ki(6)*phi(ix+1,iy,iz+1) + &
                ki(7)*phi(ix+1,iy+1,iz+1) + &
                ki(8)*phi(ix,iy+1,iz+1)
           
           P_loss(4,ptype,iproc)= P_loss(4,ptype,iproc) + Nm(ptype)*charge(ptype)*phip
        endif

80      continue
     enddo ! end-loop over ptype particles
  enddo ! end-loop over N_inj

  return  
end subroutine part_injection

subroutine shifted_maxwellian_flux(v,vb,vt,fmax,iseed)  
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Apr/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:    This subroutine is designed to generate a shifted
!                  Maxwellian flux distribution
!     NOTE:          
!     --------------------------------------------------------------
  use mpi
  implicit none
  integer iseed
  real(kind=8):: v,vb,vt,fmax,f,ran2,vm,vp,rnd(2)

  vm= MAX(0.d0,vb-4.d0*vt) ! vmin
  vp= vb+4.d0*vt ! vmax
  
90 rnd(1)= ran2(iseed)
  rnd(2)= ran2(iseed)
  v= vm + rnd(1)*(vp-vm)
  f= v*dexp(-(v-vb)**2/vt**2)
  
  ! acceptance/rejection
  if(rnd(2).gt.f/fmax) goto 90

  return
end subroutine shifted_maxwellian_flux
