! subroutine load_part(n,h,bcnd,np,vxp,ntype,nmax,kq,&
!      ni0,np_tot,nproc,iseed,sum_dEk,Nh,mpi_rank,nproc_mpi)
! !     ==============================================================
! !     VERSION:         0.3
! !     LAST MOD:      FEB/15
! !     MOD AUTHOR:    G. Fubiani
! !     COMMENTS:   
! !     NOTE:          
! !     --------------------------------------------------------------
!   use omp_lib
!   use mpi
!   use mod_legacy_particle_globals
!   use mod_constants
!   use mod_utils, only: stop_calculation
!   implicit none
!   include 'particle_info.h'
!   integer ierr,mpi_rank,nproc_mpi
!   integer:: ntype,ptype,nmax,n(3),iproc,nproc,Nh(nproc)
!   real(kind=8):: h(3)
!   ! Particle arrays
!   integer:: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2),np_tot(ntype,nproc), &
!        sum_np_tot_OMP(ntype),sum_np_tot(ntype)
!   real(kind=8):: vxp(6,nmax,ntype,nproc)
!   real(kind=8):: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc),ni0(npart),&
!        sum_dEk(nproc),kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)
!   ! Macroscopic parameters 
!   integer:: iseed(nproc),iseed_OMP
  
!   !
!   ! Loop over OpenMP processes
!   !
  
!   !$OMP PARALLEL PRIVATE(iproc,ptype,iseed_OMP)
!   ! Get processor id (from 0 to nproc-1)
!   iproc= omp_get_thread_num() + 1
!   iseed_OMP= iseed(iproc)
!   call load_part_OMP(n,h,bcnd,np,vxp,ntype,nmax,kq, &
!        ni0,np_tot,nproc,iseed_OMP,sum_dEk,Nh,nproc_mpi,iproc)
!   iseed(iproc)=iseed_OMP
!   !$OMP END PARALLEL

!   sum_np_tot= 0
!   sum_np_tot_OMP= SUM(np_tot,DIM=2)
!   call MPI_ALLREDUCE(sum_np_tot_OMP, sum_np_tot, ntype, MPI_INTEGER, MPI_SUM, &
!        MPI_COMM_WORLD, ierr)
  
!   return
! end subroutine load_part

! subroutine load_part_OMP(n,h,bcnd,np,vxp,ntype,nmax,kq, &
!      ni0,np_tot,nproc,iseed,sum_dEk,Nh,nproc_mpi,iproc)
! !     ==============================================================
! !     VERSION:         0.2
! !     LAST MOD:      NOV/12
! !     MOD AUTHOR:    G. Fubiani
! !     COMMENTS:   Loop order not optimized. This does not matter 
! !                 much in terms of calculation time because this
! !                 subroutine is only called once at the beginning.   
! !     NOTE:          
! !     --------------------------------------------------------------
!   use mpi
!   use mod_utils, only: stop_calculation
!   use mod_rng, only: ran2
!   use mod_constants
!   use mod_legacy_particle_globals
  
!   implicit none
!   include 'particle_info.h'
!   integer:: ix,iy,iz,j,jmax,k,nproc_mpi
!   integer:: ntype,ptype,nmax,n(3),iproc,nproc,Nh(nproc)
!   ! Particle arrays
!   integer:: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2),np_tot(ntype,nproc)
!   real(kind=8):: h(3),vxp(6,nmax,ntype,nproc),x,y,z,vx,vy,vt,dEk,rnd(2),&
!      ki(8),kp,px,py,pz,ni0(npart),sum_dEk(nproc)
!   real(kind=8):: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc), &
!        kq(0:n(1)+2,0:n(2)+2,0:n(3)+2),vz_sav(ntype)
!   ! Macroscopic parameters 
!   integer:: iseed

!   !
!   ! Initialize particle counter, variables & arrays
!   !
!   k=0
!   np_tot(:,iproc)=0
!   vz_sav=0.d0
!   np(:,:,:,:,iproc)=0.d0
!   sum_dEk(iproc)=0.d0
!   Nh(iproc)=0

!   ! Total number of particles
!   jmax= NINT(real(np_cell*n_cell)/real(nproc_mpi)/real(nproc))
!   do j=1,jmax
                 
!      ! random loading
! 70   continue
!      rnd(1)=ran2(iseed)
!      x= rnd(1)*x_load
!      rnd(1)=ran2(iseed)
!      y= rnd(1)*ymax
!      rnd(1)=ran2(iseed)
!      z= rnd(1)*zmax
        
!      ! Get particle left grid index
!      ix= INT( x/h(1) ) + 1
!      iy= INT( y/h(2) ) + 1
!      iz= INT( z/h(3) ) + 1

!      ! Load uniquely inside simulation domain
!      if( bcnd(ix,iy,iz).ge.1 .and. bcnd(ix+1,iy,iz).ge.1 .and. &
!           bcnd(ix+1,iy+1,iz).ge.1 .and. bcnd(ix,iy+1,iz).ge.1 .and. & 
!           bcnd(ix,iy,iz+1).ge.1 .and. bcnd(ix+1,iy,iz+1).ge.1 .and. &
!           bcnd(ix+1,iy+1,iz+1).ge.1 .and. bcnd(ix,iy+1,iz+1).ge.1 ) then
!         goto 70
!      endif

!      do ptype=1,ntype
        
!         ! Check if max. particle per cell is reached for ptype
!         rnd(1)= ran2(iseed)
!         if( rnd(1).gt.ni0(ptype) ) goto 80
        
!         ! Add particle to counter
!         np_tot(ptype,iproc)= np_tot(ptype,iproc) + 1
        
!         ! Get thermal velocity
!         vt= vt0(ptype)
           
!         ! Particle index
!         k= np_tot(ptype,iproc)
        
!         ! Warning
!         if(k.gt.nmax) then
!            print*, 'k > nmax in load_part'
!            print*, 'Abort calculation ...'
!            call stop_calculation
!         endif
        
!         ! Get particle location
!         vxp(1,k,ptype,iproc)= x ! same location for all particles
!         vxp(2,k,ptype,iproc)= y
!         vxp(3,k,ptype,iproc)= z
        
!         rnd(1)= ran2(iseed)
!         rnd(2)= ran2(iseed)
!         call load_gauss(vx,vy,vt,rnd)
!         vxp(4,k,ptype,iproc)= vx
!         vxp(5,k,ptype,iproc)= vy
!         if(vz_sav(ptype).eq.0.d0) then
!            rnd(1)= ran2(iseed)
!            rnd(2)= ran2(iseed)
!            call load_gauss(vx,vy,vt,rnd)
!            vxp(6,k,ptype,iproc)= vx
!            vz_sav(ptype)= vy
!         else
!            vxp(6,k,ptype,iproc)= vz_sav(ptype)
!            vz_sav(ptype)= 0.d0
!         endif
           
!         !
!         ! Calculate initial density
!         !
!         px=( ix*h(1) - x )/h(1)
!         py=( iy*h(2) - y )/h(2)
!         pz=( iz*h(3) - z )/h(3)

!         kp=Nm(ptype)/(h(1)*h(2)*h(3))
!         ki(1)= kp*px*py*pz
!         ki(2)= kp*(1.d0-px)*py*pz
!         ki(3)= kp*(1.d0-px)*(1.d0-py)*pz
!         ki(4)= kp*px*(1.d0-py)*pz
!         ki(5)= kp*px*py*(1.d0-pz)
!         ki(6)= kp*(1.d0-px)*py*(1.d0-pz)
!         ki(7)= kp*(1.d0-px)*(1.d0-py)*(1.d0-pz)
!         ki(8)= kp*px*(1.d0-py)*(1.d0-pz)
        
!         ! Charge assigned to the grid node 000
!         np(ix,iy,iz,ptype,iproc)= np(ix,iy,iz,ptype,iproc) + &
!              kq(ix,iy,iz)*ki(1)
        
!         ! Charge assigned to the grid node +00
!         np(ix+1,iy,iz,ptype,iproc)= np(ix+1,iy,iz,ptype,iproc) + &
!              kq(ix+1,iy,iz)*ki(2)
        
!         ! Charge assigned to the grid node ++0
!         np(ix+1,iy+1,iz,ptype,iproc)= np(ix+1,iy+1,iz,ptype,iproc) + &
!              kq(ix+1,iy+1,iz)*ki(3)
        
!         ! Charge assigned to the grid node 0+0
!         np(ix,iy+1,iz,ptype,iproc)= np(ix,iy+1,iz,ptype,iproc) + &
!              kq(ix,iy+1,iz)*ki(4)
        
!         ! Charge assigned to the grid node 00+
!         np(ix,iy,iz+1,ptype,iproc)= np(ix,iy,iz+1,ptype,iproc) + &
!              kq(ix,iy,iz+1)*ki(5)
        
!         ! Charge assigned to the grid node +0+
!         np(ix+1,iy,iz+1,ptype,iproc)= np(ix+1,iy,iz+1,ptype,iproc) + &
!              kq(ix+1,iy,iz+1)*ki(6)
        
!         ! Charge assigned to the grid node +++
!         np(ix+1,iy+1,iz+1,ptype,iproc)= np(ix+1,iy+1,iz+1,ptype,iproc) + &
!              kq(ix+1,iy+1,iz+1)*ki(7)
        
!         ! Charge assigned to the grid node 0++
!         np(ix,iy+1,iz+1,ptype,iproc)= np(ix,iy+1,iz+1,ptype,iproc) + &
!              kq(ix,iy+1,iz+1)*ki(8)

!         if( ptype.eq.1 .and. Pabs.gt.0.d0 ) then
!            if( ix.ge.ixl_pow .and. ix.le.ixr_pow ) then           
!               ! Calculate kinetic energy of macroparticle
!               dEk= 0.5d0*Nm(ptype)*mass(ptype)*( &
!                    vxp(4,k,ptype,iproc)*vxp(4,k,ptype,iproc) + &
!                    vxp(5,k,ptype,iproc)*vxp(5,k,ptype,iproc) + &
!                    vxp(6,k,ptype,iproc)*vxp(6,k,ptype,iproc) )
              
!               sum_dEk(iproc)= sum_dEk(iproc) + dEk
!               Nh(iproc)= Nh(iproc) + 1           
!            endif
!         endif

! 80      continue
!      enddo
!   enddo

!   return

! end subroutine load_part_OMP

! subroutine load_gauss(vx,vy,vt,rnd)
! !     ==============================================================
! !     VERSION:         0.1
! !     LAST MOD:       AUG/07
! !     MOD AUTHOR:    G. Fubiani
! !     COMMENTS:      Gaussian temperature loading
! !     NOTE:          
! !     --------------------------------------------------------------
!   use mod_legacy_particle_globals
!   use mod_constants
!   implicit none
!   real(kind=8):: theta,vp,vx,vy,vt,rnd(2)

!   ! Gaussian loading
!   vp = vt*dsqrt( -dlog(1-rnd(1)) )
  
!   ! Update transverse normalized momentum
!   theta = 2.d0*pi*rnd(2)
  
!   ! Calculate new normalized velocity
!   vx = vp*dcos(theta)
!   vy = vp*dsin(theta)
  
!   return
! end subroutine load_gauss

