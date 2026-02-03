subroutine part_sorting_OMP(vxp,sorting,nmax,ntype, &
                n,h,Plist,pl_max,nproc,np_tot,iproc,ptype)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      FEB/15
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Sort particles per OMP thread. Particles in the
!                    same cell have neighbor indexes in the array.
!     --------------------------------------------------------------
  implicit none
  integer:: nmax,ntype
  integer:: n(3),i,j,k,ix,iy,iz,ic,ptype,iproc,nproc,pl_max
  ! Particle arrays
  integer:: Plist(0:pl_max,ntype,nproc)
  real(kind=8):: h(3),sorting(6,nmax,nproc),vxp(6,nmax,ntype,nproc),&
       xp_new,yp_new,zp_new
  ! Macroscopic parameters
  integer:: np_tot(ntype,nproc)

  ! Initialization
  Plist(:,ptype,iproc)=0

  !
  ! Count the number of particles in each cell
  !
  do i=1,np_tot(ptype,iproc),1

     ! Particle grid index
     xp_new= vxp(1,i,ptype,iproc)
     yp_new= vxp(2,i,ptype,iproc)
     zp_new= vxp(3,i,ptype,iproc)

     ix= INT( xp_new/h(1) ) + 1
     iy= INT( yp_new/h(2) ) + 1
     iz= INT( zp_new/h(3) ) + 1

     ! Get cell index
     ic= (ix-1) + n(1)*((iy-1) + n(2)*(iz-1)) + 1

     ! Number of particles per cell
     Plist(ic,ptype,iproc)= Plist(ic,ptype,iproc) + 1

     ! Store temporarly
     sorting(1,i,iproc)= xp_new
     sorting(2,i,iproc)= yp_new
     sorting(3,i,iproc)= zp_new
     sorting(4,i,iproc)= vxp(4,i,ptype,iproc)
     sorting(5,i,iproc)= vxp(5,i,ptype,iproc)
     sorting(6,i,iproc)= vxp(6,i,ptype,iproc)

  enddo

  !
  ! Convert Plist into an allocation
  !
  k=0
  do i=1,pl_max     
     j= Plist(i,ptype,iproc)
     Plist(i,ptype,iproc)= k
     k= k + j      
  enddo

  !
  ! Sort array
  !
  do i=1,np_tot(ptype,iproc),1

     ! Particle grid index
     xp_new= sorting(1,i,iproc)
     yp_new= sorting(2,i,iproc)
     zp_new= sorting(3,i,iproc)

     ix= INT( xp_new/h(1) ) + 1
     iy= INT( yp_new/h(2) ) + 1
     iz= INT( zp_new/h(3) ) + 1

     !Get cell index
     ic= (ix-1) + n(1)*((iy-1) + n(2)*(iz-1)) + 1

     j= Plist(ic,ptype,iproc) + 1
     Plist(ic,ptype,iproc)= Plist(ic,ptype,iproc) + 1


    ! Save velocity & position
     vxp(1,j,ptype,iproc)= xp_new
     vxp(2,j,ptype,iproc)= yp_new
     vxp(3,j,ptype,iproc)= zp_new
     vxp(4,j,ptype,iproc)= sorting(4,i,iproc)
     vxp(5,j,ptype,iproc)= sorting(5,i,iproc)
     vxp(6,j,ptype,iproc)= sorting(6,i,iproc)

  enddo

  return
end subroutine part_sorting_OMP

subroutine part_sorting(vxp,sorting,nmax,ntype,n,h,Plist_thd,pl_max, &
     nproc,np_tot,ptype)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      APRIL/17
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Sort particles via its position inside the 
!                    simulation domain instead of per OMP thread. 
!                    vxp() is collapsed.
!     --------------------------------------------------------------
  implicit none
  integer:: nmax,nproc,n(3),i,j,k,ix,iy,iz,ic,iproc,ip,ntype,pl_max, &
       ptype,np_tot(ntype,nproc),sum_np_tot(ntype),np_tot_tmp
  real(kind=8):: xp_new,yp_new,zp_new,vxp_new,vyp_new,vzp_new
  ! Particle arrays
  real(kind=8):: vxp(6,nmax,ntype,nproc),h(3)
  integer :: Plist(0:pl_max),ip_shift(0:nproc), &
       Plist_thd(0:pl_max,ntype,nproc)
  real(kind=8):: sorting(6,nmax*nproc)

  ! Initialization
  Plist(:)= 0 ! Particle list
  ip_shift(0)=0

  ! Sum over OMP proc
  sum_np_tot= SUM(np_tot,DIM=2)
  ! Total number of particles per OMP thread
  np_tot_tmp= INT( dble(sum_np_tot(ptype))/dble(nproc) ) + 1

  !
  ! Count the number of particles in each cell
  !
  !$OMP PARALLEL PRIVATE(iproc,i,ip,xp_new,yp_new,zp_new,&
  !$OMP  vxp_new,vyp_new,vzp_new,ix,iy,iz,ic)
  !$OMP DO
  do iproc=1,nproc
     
     ! Particle list per OMP thread
     Plist_thd(:,ptype,iproc)=0

     ! Shift particle index when vxp() is collapsed
     ip_shift(iproc)= SUM(np_tot(ptype,1:iproc-1))

     do i=1,np_tot(ptype,iproc)
        
        ip= ip_shift(iproc) + i
        
        ! Particle grid index
        xp_new= vxp(1,i,ptype,iproc)
        yp_new= vxp(2,i,ptype,iproc)
        zp_new= vxp(3,i,ptype,iproc)
        vxp_new= vxp(4,i,ptype,iproc)
        vyp_new= vxp(5,i,ptype,iproc)
        vzp_new= vxp(6,i,ptype,iproc)
        
        ix= INT( xp_new/h(1) ) + 1
        iy= INT( yp_new/h(2) ) + 1
        iz= INT( zp_new/h(3) ) + 1
        
        ! Get cell index
        ic= (ix-1) + n(1)*((iy-1) + n(2)*(iz-1)) + 1
        
        ! Number of particles per cell
        !$OMP ATOMIC
        Plist(ic)= Plist(ic) + 1
        
        ! Store temporarly
        sorting(1,ip)= xp_new
        sorting(2,ip)= yp_new
        sorting(3,ip)= zp_new
        sorting(4,ip)= vxp_new
        sorting(5,ip)= vyp_new
        sorting(6,ip)= vzp_new
        
     end do
  enddo
  !$OMP END DO

  !
  ! Convert Plist into an allocation
  !
  !$OMP SINGLE
  k=0
  do i=1,pl_max
     j= Plist(i)
     Plist(i)= k
     k= k + j
  enddo
  !$OMP END SINGLE
  
  !
  ! Sort array
  !
  !$OMP DO PRIVATE(j,k,&
  !$OMP xp_new, yp_new,zp_new,&
  !$OMP vxp_new,vyp_new,vzp_new,&
  !$OMP ix,iy,iz,ic,iproc)
  do i=1,sum_np_tot(ptype)
     
     ! Particle grid index
     xp_new= sorting(1,i)
     yp_new= sorting(2,i)
     zp_new= sorting(3,i)
     vxp_new= sorting(4,i)
     vyp_new= sorting(5,i)
     vzp_new= sorting(6,i)
     
     ix= INT( xp_new/h(1) ) + 1
     iy= INT( yp_new/h(2) ) + 1
     iz= INT( zp_new/h(3) ) + 1
     
     ! Get cell index
     ic= (ix-1) + n(1)*((iy-1) + n(2)*(iz-1)) + 1
     
     !$OMP ATOMIC CAPTURE
     Plist(ic)= Plist(ic) + 1
     j= Plist(ic)
     !$OMP END ATOMIC
     
     k= MOD(j-1,np_tot_tmp) + 1
     iproc= (j-1)/np_tot_tmp + 1
     
     ! Save velocity & position
     vxp(1,k,ptype,iproc)= xp_new
     vxp(2,k,ptype,iproc)= yp_new
     vxp(3,k,ptype,iproc)= zp_new
     vxp(4,k,ptype,iproc)= vxp_new
     vxp(5,k,ptype,iproc)= vyp_new
     vxp(6,k,ptype,iproc)= vzp_new

  enddo
  !$OMP END DO

  !
  ! Particle list per OMP thread
  !
  !$OMP DO PRIVATE(k,iproc)
  do i=1,pl_max
     k= Plist(i)
     if(k.gt.0) then
        iproc= (k-1)/np_tot_tmp+1
        Plist_thd(i,ptype,iproc)= k - (iproc-1)*np_tot_tmp
        if(iproc.gt.1) Plist_thd(i,ptype,iproc-1)= np_tot_tmp
     endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !
  ! Update number of particles per OMP thread
  !
  do iproc=1,nproc
     np_tot(ptype,iproc)= np_tot_tmp
     if( SUM(np_tot(ptype,1:iproc)).gt.sum_np_tot(ptype) ) then
        np_tot(ptype,iproc)= sum_np_tot(ptype) - SUM(np_tot(ptype,1:iproc-1))
        np_tot(ptype,iproc+1:nproc)= 0
        exit
     endif
  enddo

  return
end subroutine part_sorting
