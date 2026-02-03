subroutine calc_rho(n,np,rhs_dom,bcnd,ntype,nproc,nproc_mpi,mpi_rank)
!     ==============================================================
!     VERSION:         0.4
!     LAST MOD:      Oct/16
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  use omp_lib
  use mpi
  implicit none
  integer ierr
  integer:: ix,iy,iz,ptype,n(3),ntype,iproc,nproc,nproc_mpi, &
       mpi_rank,m,k,opt
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc), &
       rhs_dom(0:n(1)+1,0:n(2)+1,0:n(3)/nproc_mpi+1)
  real(kind=8), allocatable:: rhs(:,:,:)
  include 'particle_info.h'

  allocate ( rhs(0:n(1)+1,0:n(2)+1,0:n(3)+1) )

  ! Initialize arrays
  rhs= 0.d0
  rhs_dom=0.d0
  m= n(3)/nproc_mpi

  ! Periodic boundary conditions (performed in parallel)
  if(flag_pbc.eq.1) then

     !$OMP PARALLEL PRIVATE(ix,iy,iz,iproc)
     ! Get processor id (from 0 to nproc-1)
     iproc= omp_get_thread_num() + 1

     ! j=1 & ny+1 planes
     do iz=1,n(3)+1
        do ix=1,n(1)+1
           if( bcnd(ix,1,iz).eq.0 ) then
              np(ix,1,iz,:,iproc)= 0.5d0*( np(ix,1,iz,:,iproc) + &
                   np(ix,n(2)+1,iz,:,iproc) )
              np(ix,0,iz,:,iproc)= np(ix,n(2),iz,:,iproc)
           endif
        
           if( bcnd(ix,n(2)+1,iz).eq.0 ) then
              np(ix,n(2)+2,iz,:,iproc)= np(ix,2,iz,:,iproc)
              np(ix,n(2)+1,iz,:,iproc)= np(ix,1,iz,:,iproc)   
           endif
        enddo
     enddo

     ! k=1 & nz+1 planes
     do iy=1,n(2)+1
        do ix=1,n(1)+1
           if( bcnd(ix,iy,1).eq.0 ) then
              np(ix,iy,1,:,iproc)= 0.5d0*( np(ix,iy,1,:,iproc) + &
                   np(ix,iy,n(3)+1,:,iproc) )
              np(ix,iy,0,:,iproc)= np(ix,iy,n(3),:,iproc)
           endif
        
           if( bcnd(ix,iy,n(3)+1).eq.0 ) then
              np(ix,iy,n(3)+2,:,iproc)= np(ix,iy,2,:,iproc)
              np(ix,iy,n(3)+1,:,iproc)= np(ix,iy,1,:,iproc)   
           endif
        enddo
     enddo
     !$OMP END PARALLEL

  endif
  
  !
  ! Reduction
  !
  !$OMP PARALLEL
  !$OMP DO
  do iz=1,n(3)+1
     do iy=1,n(2)+1
        do ix=1,n(1)+1
           do iproc=1,nproc 
              do ptype=1,ntype
!                 if(charge(ptype).ne.0) then
                 rhs(ix,iy,iz)= rhs(ix,iy,iz) - &
                      charge(ptype)*np(ix,iy,iz,ptype,iproc)
!                 endif
              enddo
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL  

  if(nproc_mpi.gt.1) then
     opt= 2
     if(opt.eq.1) then
        do k=0,nproc_mpi-1
           call MPI_REDUCE(rhs(:,:,k*m:(k+1)*m+1), rhs_dom, (n(1)+2)*(n(2)+2)*(m+2), &
                MPI_REAL8, MPI_SUM, k,MPI_COMM_WORLD, ierr)
        enddo
     else
        call MPI_ALLREDUCE(MPI_IN_PLACE, rhs, (n(1)+2)*(n(2)+2)*(n(3)+2), &
             MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        rhs_dom(:,:,0:m+1)= rhs(:,:,mpi_rank*m:(mpi_rank+1)*m+1)
     endif
  else
     rhs_dom(:,:,0:n(3)+1)= rhs(:,:,0:n(3)+1)
  endif

  deallocate ( rhs )

  return

end subroutine calc_rho

subroutine dens_red(n,np,np_red,bcnd,ntype,nproc,nproc_mpi)
!     ==============================================================
!     VERSION:         0.3
!     LAST MOD:      Oct/16
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  use omp_lib
  use mpi
  implicit none
  integer ierr
  integer:: ix,iy,iz,ptype,n(3),ntype,iproc,nproc,nproc_mpi
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc), &
       np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
  real(kind=8), allocatable:: np_red_tmp(:,:,:,:) 
  include 'particle_info.h'

  ! Initialize arrays
  np_mx(1:ntype)= 0.d0
  np_red=0.d0

  !
  ! Reduction
  !
  !$OMP PARALLEL
  !$OMP DO
  do iz=1,n(3)+1
     do iy=1,n(2)+1
        do ix=1,n(1)+1  
           do iproc=1,nproc 
              do ptype=1,ntype
                 np_red(ix,iy,iz,ptype)= np_red(ix,iy,iz,ptype) + &
                      np(ix,iy,iz,ptype,iproc)
              enddo
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL  

  if(nproc_mpi.gt.1) then
     allocate( np_red_tmp(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype) )
     np_red_tmp= 0.d0
     call MPI_ALLREDUCE(np_red(:,:,:,:), np_red_tmp, (n(1)+3)*(n(2)+3)*(n(3)+3)*ntype, &
          MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     np_red= np_red_tmp
     deallocate( np_red_tmp )
  endif

  do ptype=1,ntype
     !$OMP PARALLEL REDUCTION(MAX:np_mx)
     !$OMP DO
     do iz=1,n(3)+1
        do iy=1,n(2)+1
           do ix=1,n(1)+1
              ! Calculate maximum density
              if( bcnd(ix,iy,iz).eq.-1 ) & ! inside simulation domain
                   np_mx(ptype)= MAX(np_mx(ptype), np_red(ix,iy,iz,ptype))
           enddo
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL  
  enddo

  return

end subroutine dens_red

