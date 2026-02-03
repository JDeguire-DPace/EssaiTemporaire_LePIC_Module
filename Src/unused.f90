subroutine part_sorting_dom(istep,vxp,nmax,ntype,n,h,pl_max,nproc,np_tot)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      FEB/15
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Sort particles per block. A domain decomposition
!                    is performed. Particles are grouped and sorted 
!                    based on their location inside the simulation 
!                    domain. Each block corresponds to an OMP processor,
!                    Same number of particles per block. 
!     --------------------------------------------------------------
  use omp_lib 
  implicit none
  integer:: nmax,ntype,n(3),i,j,k,ix,iy,iz,ic,ptype,istep,iproc, &
       nproc,pl_max,i_shift,iblock,blocksize,mpi_rank
  ! Particle arrays
  integer:: Plist(pl_max,nproc),Plist_red(pl_max)
  real(kind=8):: h(3),vxp(6,nmax,ntype,nproc),xp_new,yp_new,zp_new,&
       vpx_new,vpy_new,vpz_new
  real(kind=8), allocatable:: sorting(:,:)
  ! Macroscopic parameters
  integer:: np_tot(ntype,nproc),sum_np_tot,np_tot_tmp(0:nproc)
  include 'particle_info.h'

  np_tot_tmp(0)=0

  ! Allocate array
  allocate ( sorting(6,nmax*nproc) )

  do ptype=1,ntype ! Loop over ptype particles

     if(mass(ptype).lt.0) goto 100
     if( SUM(np_tot(ptype,1:nproc)).eq.0 ) goto 100

     ! Initialization
     Plist_red=0
     np_tot_tmp(1:nproc)= np_tot(ptype,1:nproc)

     !$OMP PARALLEL PRIVATE(iproc,xp_new,yp_new,zp_new,ix,iy,iz, &
     !$OMP i,j,k,ic,i_shift)
     ! Get processor id (from 0 to nproc-1)
     iproc= omp_get_thread_num() + 1     

     Plist(:,iproc)= 0

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

        !Get cell index
        ic= (ix-1) + n(1)*((iy-1) + n(2)*(iz-1)) + 1

        ! Number of particles per cell
        Plist(ic,iproc)= Plist(ic,iproc) + 1

        ! Store temporarly
        i_shift= SUM(np_tot_tmp(0:iproc-1)) + i 

        sorting(1,i_shift)= xp_new
        sorting(2,i_shift)= yp_new
        sorting(3,i_shift)= zp_new
        sorting(4,i_shift)= vxp(4,i,ptype,iproc)
        sorting(5,i_shift)= vxp(5,i,ptype,iproc)
        sorting(6,i_shift)= vxp(6,i,ptype,iproc)
        
     enddo

     !
     ! Convert Plist into an allocation
     !
     k=0
     do ic=1,pl_max     
        j= Plist(ic,iproc)
        Plist(ic,iproc)= k
        k= k + j      
     enddo
     !$OMP END PARALLEL

     ! Reduction
     Plist_red(:)=SUM(Plist,DIM=2)

     !
     ! Sort particles per block
     !
     sum_np_tot= SUM(np_tot(ptype,1:nproc))

     ! Balance number of particles
     np_tot(ptype,:)=INT( real(sum_np_tot)/real(nproc) )
     np_tot(ptype,nproc)= np_tot(ptype,nproc) + &
          ( sum_np_tot - SUM(np_tot(ptype,1:nproc)) )

     np_tot_tmp(1:nproc)= np_tot(ptype,1:nproc)
     blocksize= np_tot_tmp(1)

     !$OMP PARALLEL PRIVATE(iproc,xp_new,yp_new,zp_new,ix,iy,iz, &
     !$OMP ic,i,j,iblock)
     ! Get processor id (from 0 to nproc-1)
     iproc= omp_get_thread_num() + 1

     do i=SUM(np_tot_tmp(0:iproc-1))+1,SUM(np_tot_tmp(0:iproc)),1
        
        ! Particle grid index        
        xp_new= sorting(1,i)
        yp_new= sorting(2,i)
        zp_new= sorting(3,i)
        
        ix= INT( xp_new/h(1) ) + 1
        iy= INT( yp_new/h(2) ) + 1
        iz= INT( zp_new/h(3) ) + 1
        
        !Get cell index
        ic= (ix-1) + n(1)*((iy-1) + n(2)*(iz-1)) + 1

        !$OMP critical
        j= Plist_red(ic) + 1
        Plist_red(ic)= Plist_red(ic) + 1
        !$OMP end critical

        iblock= INT( (j-1)/blocksize ) + 1
        if(iblock.gt.nproc) iblock= nproc
        j= j - SUM(np_tot_tmp(0:iblock-1))

        ! Save velocity & position
        vxp(1,j,ptype,iblock)= xp_new
        vxp(2,j,ptype,iblock)= yp_new
        vxp(3,j,ptype,iblock)= zp_new
        vxp(4,j,ptype,iblock)= sorting(4,i)
        vxp(5,j,ptype,iblock)= sorting(5,i)
        vxp(6,j,ptype,iblock)= sorting(6,i)

     enddo
     !$OMP END PARALLEL

100  continue
  enddo ! end-do over ptype particles

  ! Deallocate array
  deallocate ( sorting )

  return
end subroutine part_sorting_dom
