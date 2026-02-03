subroutine part_mover(n,h,Ei,Bi,p_mac,P_loss,vxp,&
     bcnd,nmax,ntype,ngrid,flag_dead,nproc,np_tot,&
     iproc,ptype,flag_cex,cnt_cex,cnt_dead,beam_div,sum_q_xz,&
     sum_q_yz,dtype,iseed,n_B,h_B)
!     ==============================================================
!     VERSION:         0.4
!     LAST MOD:      Dec/23
!     MOD AUTHOR
!     NOTE:       Interpolation |------|-------------|
!                              ip      xp           ip+1 
!                              <------><------------>
!                                 1-p         p  
!
!                 Indexes in vxp(): 1== x
!                                   2== y    
!                                   3== z
!                                   4== vx
!                                   5== vy
!                                   6== vz                
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: n(3),nmax,ntype,ngrid,igrid,n_B(3)
  integer:: ix,iy,iz,i,ptype,flag_B_tmp,flag_lost,iproc,nproc
  ! Field arrays
  real(kind=8):: Ei(3,0:n(1)+2,0:n(2)+2,0:n(3)+2),Exp,Eyp,Ezp,h(3), &
       Bi(4,0:n_B(1)+2,0:n_B(2)+2,0:n_B(3)+2),ki(8),h_B(3)
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2),flag_cex(nmax,nproc),&
       cnt_cex(3,nproc)
  integer(kind=1):: flag_dead(nmax,ntype,nproc)
  real(kind=8):: vxp(6,nmax,ntype,nproc),xp_new,yp_new,zp_new,vpx_new,&
       vpy_new,vpz_new,Eki,px,py,pz,cnt_dead(nproc),&
       beam_div(4,nproc)
  real(kind=8):: v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,Bpx,Bpy,Bpz,B_tot
  ! Particle densities
  real(kind=8):: k1,k2,k3,sum_q_xz(0:n(1)+2,0:n(3)+2,ntype,nproc),&
       sum_q_yz(2,0:n(2)+2,0:n(3)+2,ntype,nproc)
  ! Macroscopic parameters
  integer:: np_tot(ntype,nproc),np_lost(ntype,nproc),i_shift,dtype(0:ngrid),&
       d_ind,iseed,n_sec,ip_sec
  real(kind=8):: p_mac(ntype,2,0:ngrid,nproc),P_loss(4,ntype,nproc),theta,&
       ran2,rnd(2),vx_sec,vy_sec,vz_sec,vt

  np_lost(ptype,iproc)=0

  ! Set coefficients for Boris solver
  k1= dt*charge(ptype)/(2.d0*mass(ptype))

  !
  ! Move particles
  !
  do i=1,np_tot(ptype,iproc),1

     !
     ! Save actual velocity & position (avoid calling arrays many times)
     !
     xp_new= vxp(1,i,ptype,iproc)
     yp_new= vxp(2,i,ptype,iproc)
     zp_new= vxp(3,i,ptype,iproc)
     vpx_new= vxp(4,i,ptype,iproc)
     vpy_new= vxp(5,i,ptype,iproc)
     vpz_new= vxp(6,i,ptype,iproc)

     !
     ! Remove particles killed during a collision in the previous step 
     !
     if(flag_dead(i,ptype,iproc).eq.1) then
        ! Add another lost particle
        np_lost(ptype,iproc)= np_lost(ptype,iproc) + 1
        ! Re-initialize 
        flag_dead(i,ptype,iproc)= 0
        ! Counter
        if(ptype.eq.tag_neg) cnt_dead(iproc)= cnt_dead(iproc) + Nm(ptype)*charge(ptype)
        ! Calculate kinetic energy of macroparticle
        Eki= 0.5d0*Nm(ptype)*mass(ptype)*( vpx_new*vpx_new + &
             vpy_new*vpy_new + vpz_new*vpz_new ) 
        ! Total power lost per time step
        P_loss(3,ptype,iproc)= P_loss(3,ptype,iproc) - Eki
        ! Jump to the end of the loop
        goto 110
     endif
     
     !
     ! Get particle left grid index
     !
     ix= INT( xp_new/h(1) ) + 1
     iy= INT( yp_new/h(2) ) + 1
     iz= INT( zp_new/h(3) ) + 1

     !
     ! Calculate fields at position x,y
     !
     px=( ix*h(1) - xp_new )/h(1)
     py=( iy*h(2) - yp_new )/h(2)
     pz=( iz*h(3) - zp_new )/h(3)

     ki(1)= px*py*pz
     ki(2)= (1.d0-px)*py*pz
     ki(3)= (1.d0-px)*(1.d0-py)*pz
     ki(4)= px*(1.d0-py)*pz
     ki(5)= px*py*(1.d0-pz)
     ki(6)= (1.d0-px)*py*(1.d0-pz)
     ki(7)= (1.d0-px)*(1.d0-py)*(1.d0-pz)
     ki(8)= px*(1.d0-py)*(1.d0-pz)

     Exp= ki(1)*Ei(1,ix,iy,iz) + &
          ki(2)*Ei(1,ix+1,iy,iz) + &
          ki(3)*Ei(1,ix+1,iy+1,iz) + &
          ki(4)*Ei(1,ix,iy+1,iz) + &
          ki(5)*Ei(1,ix,iy,iz+1) + &
          ki(6)*Ei(1,ix+1,iy,iz+1) + &
          ki(7)*Ei(1,ix+1,iy+1,iz+1) + &
          ki(8)*Ei(1,ix,iy+1,iz+1)
     
     Eyp= ki(1)*Ei(2,ix,iy,iz) + &
          ki(2)*Ei(2,ix+1,iy,iz) + &
          ki(3)*Ei(2,ix+1,iy+1,iz) + &
          ki(4)*Ei(2,ix,iy+1,iz) + &
          ki(5)*Ei(2,ix,iy,iz+1) + &
          ki(6)*Ei(2,ix+1,iy,iz+1) + &
          ki(7)*Ei(2,ix+1,iy+1,iz+1) + &
          ki(8)*Ei(2,ix,iy+1,iz+1)
     
     Ezp= ki(1)*Ei(3,ix,iy,iz) + &
          ki(2)*Ei(3,ix+1,iy,iz) + &
          ki(3)*Ei(3,ix+1,iy+1,iz) + &
          ki(4)*Ei(3,ix,iy+1,iz) + &
          ki(5)*Ei(3,ix,iy,iz+1) + &
          ki(6)*Ei(3,ix+1,iy,iz+1) + &
          ki(7)*Ei(3,ix+1,iy+1,iz+1) + &
          ki(8)*Ei(3,ix,iy+1,iz+1)

     !
     ! Calculate new velocity
     !
     flag_B_tmp=0
     if(flag_B.eq.1)then

        if(n_B(1).eq.1) then ! Constant B-field

           Bpx= Bi(1,1,1,1)
           Bpy= Bi(2,1,1,1)
           Bpz= Bi(3,1,1,1)
           
           flag_B_tmp=1
        else ! B-field map
           
           if(flag_gridB.eq.1) then
              ! Get particle left grid index on B-field map
              ix= INT( xp_new/h_B(1) ) + 1
              iy= INT( yp_new/h_B(2) ) + 1
              iz= INT( zp_new/h_B(3) ) + 1
           endif
        
           if( Bi(4,ix,iy,iz).gt.1.d-5 ) then  ! > 0.1G
              
              if(flag_gridB.eq.1) then
                 ! Calculate B-fields at position x,y and z on the map
                 px=( ix*h_B(1) - xp_new )/h_B(1)
                 py=( iy*h_B(2) - yp_new )/h_B(2)
                 pz=( iz*h_B(3) - zp_new )/h_B(3)
                 
                 ki(1)= px*py*pz
                 ki(2)= (1.d0-px)*py*pz
                 ki(3)= (1.d0-px)*(1.d0-py)*pz
                 ki(4)= px*(1.d0-py)*pz
                 ki(5)= px*py*(1.d0-pz)
                 ki(6)= (1.d0-px)*py*(1.d0-pz)
                 ki(7)= (1.d0-px)*(1.d0-py)*(1.d0-pz)
                 ki(8)= px*(1.d0-py)*(1.d0-pz)
              endif
              
              ! Bx
              Bpx= ki(1)*Bi(1,ix,iy,iz) + &
                   ki(2)*Bi(1,ix+1,iy,iz) + &
                   ki(3)*Bi(1,ix+1,iy+1,iz) + &
                   ki(4)*Bi(1,ix,iy+1,iz) + &
                   ki(5)*Bi(1,ix,iy,iz+1) + &
                   ki(6)*Bi(1,ix+1,iy,iz+1) + &
                   ki(7)*Bi(1,ix+1,iy+1,iz+1) + &
                   ki(8)*Bi(1,ix,iy+1,iz+1)
              
              ! By
              Bpy= ki(1)*Bi(2,ix,iy,iz) + &
                   ki(2)*Bi(2,ix+1,iy,iz) + &
                   ki(3)*Bi(2,ix+1,iy+1,iz) + &
                   ki(4)*Bi(2,ix,iy+1,iz) + &
                   ki(5)*Bi(2,ix,iy,iz+1) + &
                   ki(6)*Bi(2,ix+1,iy,iz+1) + &
                   ki(7)*Bi(2,ix+1,iy+1,iz+1) + &
                   ki(8)*Bi(2,ix,iy+1,iz+1)
              
              ! Bz
              Bpz= ki(1)*Bi(3,ix,iy,iz) + &
                   ki(2)*Bi(3,ix+1,iy,iz) + &
                   ki(3)*Bi(3,ix+1,iy+1,iz) + &
                   ki(4)*Bi(3,ix,iy+1,iz) + &
                   ki(5)*Bi(3,ix,iy,iz+1) + &
                   ki(6)*Bi(3,ix+1,iy,iz+1) + &
                   ki(7)*Bi(3,ix+1,iy+1,iz+1) + &
                   ki(8)*Bi(3,ix,iy+1,iz+1)
          
              flag_B_tmp=1           
           endif
        endif
     endif

     ! Positive ions not magnetized
     if( flag_B_pos.eq.1 .and. charge(ptype).gt.0.d0 ) flag_B_tmp=0

     if(flag_B_tmp.eq.1) then ! Boris scheme

        B_tot= dsqrt( Bpx*Bpx + Bpy*Bpy + Bpz*Bpz ) 
        k2= k1*( 2.d0/(1.d0 + (k1*B_tot)*(k1*B_tot)) )

        v1x= vpx_new + k1*Exp
        v1y= vpy_new + k1*Eyp
        v1z= vpz_new + k1*Ezp

        v3x= v1x + k1*( v1y*Bpz - v1z*Bpy )
        v3y= v1y + k1*( v1z*Bpx - v1x*Bpz )
        v3z= v1z + k1*( v1x*Bpy - v1y*Bpx )

        v2x= v1x + k2*( v3y*Bpz - v3z*Bpy )
        v2y= v1y + k2*( v3z*Bpx - v3x*Bpz )
        v2z= v1z + k2*( v3x*Bpy - v3y*Bpx ) 

        vpx_new= v2x + k1*Exp
        vpy_new= v2y + k1*Eyp
        vpz_new= v2z + k1*Ezp

     else ! Electrostatic Leap-Frog solver

        vpx_new= vpx_new + 2.d0*k1*Exp
        vpy_new= vpy_new + 2.d0*k1*Eyp
        vpz_new= vpz_new + 2.d0*k1*Ezp
           
     endif

     ! Calculate new position  
     xp_new= xp_new + dt*vpx_new
     yp_new= yp_new + dt*vpy_new
     zp_new= zp_new + dt*vpz_new
     

     ! Get particle new grid location
     ix= FLOOR( xp_new/h(1) ) + 1
     iy= FLOOR( yp_new/h(2) ) + 1
     iz= FLOOR( zp_new/h(3) ) + 1


     !
     ! Correct for particles out of bounds
     !
     if(ix.lt.0) ix=0
     if(ix.gt.n(1)+1) ix=n(1)+1
     if(iy.lt.0) iy=0
     if(iy.gt.n(2)+1) iy=n(2)+1
     if(iz.lt.0) iz=0
     if(iz.gt.n(3)+1) iz=n(3)+1

     !
     ! Check location of particle inside the simulation domain
     !
     flag_lost= 0
     if( bcnd(ix,iy,iz).ge.1 .and. bcnd(ix+1,iy,iz).ge.1 .and. &
          bcnd(ix+1,iy+1,iz).ge.1 .and. bcnd(ix,iy+1,iz).ge.1 .and. & 
          bcnd(ix,iy,iz+1).ge.1 .and. bcnd(ix+1,iy,iz+1).ge.1 .and. &
          bcnd(ix+1,iy+1,iz+1).ge.1 .and. bcnd(ix,iy+1,iz+1).ge.1 ) flag_lost=1


     !
     ! Negative ions and Neumann BC's (no specular reflexion)
     !
     if(ptype.eq.tag_neg) then
        if( xp_new.lt.0.d0 .and. flag_nmn.eq.1 ) flag_lost=2
     endif

     !
     ! Lost particles
     !
     if(flag_lost.ge.1) then

        ! Add another lost particle
        np_lost(ptype,iproc)= np_lost(ptype,iproc) + 1
        
        ! Calculate kinetic energy of macroparticle
        Eki= 0.5d0*Nm(ptype)*mass(ptype)*( vpx_new*vpx_new + &
             vpy_new*vpy_new + vpz_new*vpz_new ) 
        
        ! Total lost power per time step
        P_loss(1,ptype,iproc)= P_loss(1,ptype,iproc) + Eki
        
        ! Grid index
        igrid= bcnd(ix,iy,iz)
        
        ! Neumann BC's and negative ions
        if( flag_lost.eq.2 .and. ptype.eq.tag_neg ) igrid= 0

        ! Dielectric surfaces 
        if( flag_die.eq.1 .and. dtype(igrid).gt.1 ) then

           ! Dielectric combined with periodic BC's along (OZ)
           if(flag_pbc.eq.1) then
              if(zp_new.ge.zmax) then 
                 zp_new= zp_new - zmax
                 iz= FLOOR( zp_new/h(3) ) + 1
              endif
              if(zp_new.le.0.d0 ) then 
                 zp_new= zmax + zp_new
                 iz= FLOOR( zp_new/h(3) ) + 1
              endif
           endif

           ! Linear interpolation
           pz=( iz*h(3) - zp_new )/h(3)
           k3= Nm(ptype)*charge(ptype)

           ! Y=cste planes
           if( dtype(igrid).eq.2 ) then
              px=( ix*h(1) - xp_new )/h(1)
              ki(1)= k3*px*pz
              ki(2)= k3*(1.d0-px)*pz
              ki(3)= k3*(1.d0-px)*(1.d0-pz)
              ki(4)= k3*px*(1.d0-pz)

              sum_q_xz(ix,iz,ptype,iproc)= sum_q_xz(ix,iz,ptype,iproc) + ki(1)
              sum_q_xz(ix+1,iz,ptype,iproc)= sum_q_xz(ix+1,iz,ptype,iproc) + ki(2)
              sum_q_xz(ix+1,iz+1,ptype,iproc)= sum_q_xz(ix+1,iz+1,ptype,iproc) + ki(3)
              sum_q_xz(ix,iz+1,ptype,iproc)= sum_q_xz(ix,iz+1,ptype,iproc) + ki(4)
           endif

           ! X=cste planes #1 and 2
           if( dtype(igrid).eq.3 .or. dtype(igrid).eq.4 ) then 
              py=( iy*h(2) - yp_new )/h(2)
              ki(1)= k3*py*pz
              ki(2)= k3*(1.d0-py)*pz
              ki(3)= k3*(1.d0-py)*(1.d0-pz)
              ki(4)= k3*py*(1.d0-pz)

              d_ind= dtype(igrid)
              sum_q_yz(d_ind-2,iy,iz,ptype,iproc)= sum_q_yz(d_ind-2,iy,iz,ptype,iproc) + ki(1)
              sum_q_yz(d_ind-2,iy+1,iz,ptype,iproc)= sum_q_yz(d_ind-2,iy+1,iz,ptype,iproc) + ki(2)
              sum_q_yz(d_ind-2,iy+1,iz+1,ptype,iproc)= sum_q_yz(d_ind-2,iy+1,iz+1,ptype,iproc) + ki(3)
              sum_q_yz(d_ind-2,iy,iz+1,ptype,iproc)= sum_q_yz(d_ind-2,iy,iz+1,ptype,iproc) + ki(4)
           endif
        endif

        ! Secondary particle emission
        if(charge(ptype).gt.0 .and. gam_sec.gt.0.d0) then
           if(igrid.eq.igrid_sec) then
              rnd(1)= ran2(iseed)
              n_sec= INT(gam_sec)
              if(rnd(1) .le. (gam_sec-n_sec)) n_sec= n_sec+1
              if(n_sec.eq.0) goto 130
              ! Generate secondary particle
              do ip_sec=1,n_sec
                 
                 ! Create a new electron
                 np_tot(1,iproc)= np_tot(1,iproc) + 1
                 i_shift= np_tot(1,iproc)
                 
                 ! Use THm for electron temperature
                 vt= dsqrt(2.d0*qe*ABS(THm)/ABS(mass(1))) 
                 rnd(1)=ran2(iseed)
                 !  dir== -sign(1.d0,vpz_new)
                 vz_sec = -sign(1.d0,vpz_new)*vt*dsqrt( -dlog(1-rnd(1)) )
                 rnd(1)= ran2(iseed)
                 rnd(2)= ran2(iseed)
                 ! Gaussian loading (flux normal to the grid surface)
                 call load_gauss(vx_sec,vy_sec,vt,rnd)
                 
                 ! x and y-locations of impacting ion.
                 vxp(1,i_shift,1,iproc)= xp_new
                 vxp(2,i_shift,1,iproc)= yp_new
                 if(vpz_new.lt.0.d0) then 
                    vxp(3,i_shift,1,iproc)= zg_sec(1)
                 else
                    vxp(3,i_shift,1,iproc)= zg_sec(2)
                 endif
                 vxp(4,i_shift,1,iproc)= vx_sec 
                 vxp(5,i_shift,1,iproc)= vy_sec 
                 vxp(6,i_shift,1,iproc)= vz_sec
                 
                 ! Count power injected into the plasma electrons
                 P_loss(4,1,iproc)= P_loss(4,1,iproc) + 0.5d0*Nm(1)*mass(1)*( &
                      vx_sec*vx_sec + vz_sec*vz_sec + vz_sec*vz_sec )
                 
              enddo
130           continue
           endif
        endif
        
        ! Total number of particle lost at the wall per time step
        p_mac(ptype,np_loss,igrid,iproc)= p_mac(ptype,np_loss,igrid,iproc) + 1.
        
        ! Total lost power at the wall per time step
        p_mac(ptype,P_w,igrid,iproc)= p_mac(ptype,P_w,igrid,iproc) + Eki   

        ! Characteristics of negative ion population impacting num_grd 
        if(ptype.eq.tag_neg) then        
           if(igrid.eq.num_grd) then
              cnt_cex(1,iproc)= cnt_cex(1,iproc) + 1
              if(ABS(flag_cex(i,iproc)).eq.tag_neu) &
                   cnt_cex(2,iproc)= cnt_cex(2,iproc) + 1
              if(flag_cex(i,iproc).eq.-tag_neu) &
                   cnt_cex(3,iproc)= cnt_cex(3,iproc) + 1
              ! RMS beam divergence
              theta= datan(vpy_new/vpx_new)
              beam_div(1,iproc)= beam_div(1,iproc) + theta
              beam_div(2,iproc)= beam_div(2,iproc) + theta*theta
              theta= datan(vpz_new/vpx_new)
              beam_div(3,iproc)= beam_div(3,iproc) + theta
              beam_div(4,iproc)= beam_div(4,iproc) + theta*theta
           endif
        endif
     
     endif

     !
     ! Particles inside the simulation domain
     !
     if(flag_lost.eq.0) then
   
        ! Neumann BCs (LHS only)
        if(flag_nmn.eq.1) then
           if( xp_new.le.0.d0 ) then
              ! Specular reflection
              xp_new= -xp_new
              vpx_new= -vpx_new
           endif
        endif

        ! Periodic BCs
        if(flag_pbc.eq.1) then
           if( yp_new.ge.ymax ) yp_new= yp_new - ymax
           if( yp_new.le.0.d0 ) yp_new= ymax + yp_new
           if( zp_new.ge.zmax ) zp_new= zp_new - zmax
           if( zp_new.le.0.d0 ) zp_new= zmax + zp_new
        endif

        ! Shift location of particle inside array
        i_shift= i-np_lost(ptype,iproc)

        vxp(1,i_shift,ptype,iproc)= xp_new
        vxp(2,i_shift,ptype,iproc)= yp_new
        vxp(3,i_shift,ptype,iproc)= zp_new
        vxp(4,i_shift,ptype,iproc)= vpx_new
        vxp(5,i_shift,ptype,iproc)= vpy_new
        vxp(6,i_shift,ptype,iproc)= vpz_new

        if(ptype.eq.tag_neg) &
             flag_cex(i_shift,iproc)= flag_cex(i,iproc) 
           
     endif

110  continue
  enddo ! end loop over np_tot(ptype) particles

  !
  ! Update particle counter
  !
  np_tot(ptype,iproc)= np_tot(ptype,iproc) - np_lost(ptype,iproc)

  return
 
end subroutine part_mover

subroutine eheating(h,vxp,nmax,ntype,nproc,iseed,np_tot,vt,iproc,P_loss)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      MAR/10
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:   Electrons Maxwellian heating
!                 1== xp
!                 4== vpx
!                 5== vpy
!                 6== vpz
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: nmax,ntype,iseed,np_tot(ntype,nproc),np_tot_tmp,&
       ix,i,ptype,iproc,nproc
  real(kind=8):: h(3),vxp(6,nmax,ntype,nproc),xp_new,vt,vz_sav, &
       ran2,rnd(2),P_loss(4,ntype,nproc),v2old,yp_new,zp_new

  ! Only electrons
  ptype=1

  ! Initialization
  vz_sav=0.d0

  !
  ! Maxwellian heating
  !
  np_tot_tmp= np_tot(ptype,iproc)
  do i=1,np_tot_tmp,1


     ! Get particle left grid index
     xp_new= vxp(1,i,ptype,iproc)
     ix= INT( xp_new/h(1) ) + 1

     ! Add heating by external source
     if( ix.ge.ixl_pow .and. ix.le.ixr_pow ) then 

        if(flag_circxh.eq.1) then
           yp_new= vxp(2,i,ptype,iproc)
           zp_new= vxp(3,i,ptype,iproc)
           if(flag_ahp.eq.0) then 
              if( ((yp_new-ymax/2.d0)**2 + (zp_new-zmax/2.d0)**2).gt.R_ahp**2 ) goto 140
           else
              if( ((yp_new-ymax/2.d0)**2 + (zp_new-zmax/2.d0)**2).lt.R_ahp**2 ) goto 140
           endif              
        endif
        
        ! Use Maxwellian with self-consistent temperature 
        rnd(1)= ran2(iseed)
        if( rnd(1).le.nudt ) then
             
           v2old= vxp(4,i,ptype,iproc)*vxp(4,i,ptype,iproc) + &
                vxp(5,i,ptype,iproc)*vxp(5,i,ptype,iproc) + &
                vxp(6,i,ptype,iproc)*vxp(6,i,ptype,iproc)

           ! Calculate new velocity 
           rnd(1)= ran2(iseed)
           rnd(2)= ran2(iseed)
           call load_gauss(vxp(4,i,ptype,iproc),vxp(5,i,ptype,iproc),vt,rnd)
           if(vz_sav.eq.0.d0) then
              rnd(1)= ran2(iseed)
              rnd(2)= ran2(iseed)
              call load_gauss(vxp(6,i,ptype,iproc),vz_sav,vt,rnd)
           else
              vxp(6,i,ptype,iproc)= vz_sav
              vz_sav= 0.d0
           endif
        
        P_loss(2,ptype,iproc)= P_loss(2,ptype,iproc) + 0.5d0*Nm(ptype)*mass(ptype)*( &
             vxp(4,i,ptype,iproc)*vxp(4,i,ptype,iproc) + &
             vxp(5,i,ptype,iproc)*vxp(5,i,ptype,iproc) + &
             vxp(6,i,ptype,iproc)*vxp(6,i,ptype,iproc) - v2old ) 

        endif ! endif R<= nu*dt

140     continue        
     endif

  enddo ! end loop over np_tot(ptype) particles

  return
 
end subroutine eheating

subroutine charge_deposition(n,h,vxp,nmax,ntype,kq,np,nproc,np_tot,sum_dEk,&
     Nh,iproc,ptype)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Apr/24
!     MOD AUTHOR     G. Fubiani
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: n(3),nmax,ntype,Nh(nproc)
  integer:: ix,iy,iz,i,ptype,iproc,nproc
  real(kind=8):: h(3),ki(8),vxp(6,nmax,ntype,nproc),xp_new,yp_new,zp_new,vpx_new,&
       vpy_new,vpz_new,Eki,px,py,pz,sum_dEk(nproc)
  ! Particle densities
  real(kind=8):: k1,kq(0:n(1)+2,0:n(2)+2,0:n(3)+2),&
       np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)
  ! Macroscopic parameters
  integer:: np_tot(ntype,nproc)

  ! Loop over the particle list
  do i=1,np_tot(ptype,iproc),1

     ! 
     ! Save particle 6D-coordinates
     !
     xp_new= vxp(1,i,ptype,iproc)
     yp_new= vxp(2,i,ptype,iproc)
     zp_new= vxp(3,i,ptype,iproc)
     vpx_new= vxp(4,i,ptype,iproc)
     vpy_new= vxp(5,i,ptype,iproc)
     vpz_new= vxp(6,i,ptype,iproc)

     !
     ! Get particle left grid index
     !
     ix= INT( xp_new/h(1) ) + 1
     iy= INT( yp_new/h(2) ) + 1
     iz= INT( zp_new/h(3) ) + 1
     
     !
     ! Calculate density
     !
     px=( ix*h(1) - xp_new )/h(1)
     py=( iy*h(2) - yp_new )/h(2)
     pz=( iz*h(3) - zp_new )/h(3)
     k1=Nm(ptype)/(h(1)*h(2)*h(3))
     ki(1)= k1*px*py*pz
     ki(2)= k1*(1.d0-px)*py*pz
     ki(3)= k1*(1.d0-px)*(1.d0-py)*pz
     ki(4)= k1*px*(1.d0-py)*pz
     ki(5)= k1*px*py*(1.d0-pz)
     ki(6)= k1*(1.d0-px)*py*(1.d0-pz)
     ki(7)= k1*(1.d0-px)*(1.d0-py)*(1.d0-pz)
     ki(8)= k1*px*(1.d0-py)*(1.d0-pz)
     
     ! Charge assigned to the grid node 000
     np(ix,iy,iz,ptype,iproc)= np(ix,iy,iz,ptype,iproc) + &
          kq(ix,iy,iz)*ki(1)
     
     ! Charge assigned to the grid node +00
     np(ix+1,iy,iz,ptype,iproc)= np(ix+1,iy,iz,ptype,iproc) + &
          kq(ix+1,iy,iz)*ki(2)
     
     ! Charge assigned to the grid node ++0
     np(ix+1,iy+1,iz,ptype,iproc)= np(ix+1,iy+1,iz,ptype,iproc) + &
          kq(ix+1,iy+1,iz)*ki(3)
     
     ! Charge assigned to the grid node 0+0
     np(ix,iy+1,iz,ptype,iproc)= np(ix,iy+1,iz,ptype,iproc) + &
          kq(ix,iy+1,iz)*ki(4)
     
     ! Charge assigned to the grid node 00+
     np(ix,iy,iz+1,ptype,iproc)= np(ix,iy,iz+1,ptype,iproc) + &
          kq(ix,iy,iz+1)*ki(5)
     
     ! Charge assigned to the grid node +0+
     np(ix+1,iy,iz+1,ptype,iproc)= np(ix+1,iy,iz+1,ptype,iproc) + &
          kq(ix+1,iy,iz+1)*ki(6)
     
     ! Charge assigned to the grid node +++
     np(ix+1,iy+1,iz+1,ptype,iproc)= np(ix+1,iy+1,iz+1,ptype,iproc) + &
          kq(ix+1,iy+1,iz+1)*ki(7)
     
     ! Charge assigned to the grid node 0++
     np(ix,iy+1,iz+1,ptype,iproc)= np(ix,iy+1,iz+1,ptype,iproc) + &
          kq(ix,iy+1,iz+1)*ki(8)

     !
     ! Count number of electrons and total kinetic energy in predefined heating region.
     !
     if( ptype.eq.1 .and. Pabs.gt.0.d0 ) then
        if( ix.ge.ixl_pow .and. ix.le.ixr_pow ) then
           
           if(flag_circxh.eq.1) then
              if(flag_ahp.eq.0) then
                 ! Option to deposit power @ r <= Rahp
                 if( ((yp_new-ymax/2.d0)**2 + (zp_new-zmax/2.d0)**2).gt.R_ahp**2 ) goto 100
              else
                 ! Option to deposit power @ r >= Rahp
                 if( ((yp_new-ymax/2.d0)**2 + (zp_new-zmax/2.d0)**2).lt.R_ahp**2 ) goto 100
              endif
           endif
           
           ! Calculate kinetic energy of macroparticle
           Eki= 0.5d0*Nm(ptype)*mass(ptype)*( vpx_new*vpx_new + &
                vpy_new*vpy_new + vpz_new*vpz_new )
           sum_dEk(iproc)= sum_dEk(iproc) + Eki
           Nh(iproc)= Nh(iproc) + 1
           
100        continue
        endif
     endif

  enddo
  
end subroutine charge_deposition
