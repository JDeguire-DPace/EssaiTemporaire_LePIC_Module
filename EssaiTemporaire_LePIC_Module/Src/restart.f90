subroutine restart(n,h,np,vxp,nmax,ntype,kq,time,nproc,iseed,&
     np_tot,sum_dEk,Nh,flag_cex,mpi_rank,nproc_mpi,np_dup)
!     ==============================================================
!     VERSION:         0.3
!     LAST MOD:      Mar/15
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  use omp_lib
  implicit none
  include 'particle_info.h'
  integer:: i,n(3),nmax,ntype,ptype,iproc,nproc,Nh(nproc),mpi_rank,&
       nproc_mpi,irank,iproc_tmp,np_tot_ext(ntype,100),n1(ntype),n2(ntype)
  ! Particle arrays
  integer:: np_tot(ntype,nproc),iseed(nproc),lgh,flag_cex(nmax,nproc),&
       i_shift,np_lost(ntype,nproc),iseed_OMP,np_tot_tmp
  real(kind=8):: h(3),vxp(6,nmax,ntype,nproc),sum_dEk(nproc),time,&
       np_dup,ran2,rnd
  ! Particle densities
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc),&
       kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)
  character:: pnum*3

  ! Initialize counters and arrays
  np_tot= 0

  !
  ! Set random seed number
  !
  do iproc=1,nproc
     iseed(iproc)= 123456*iproc*(10*mpi_rank+1)
  enddo

  !
  ! Open backup file
  !

  ! Warnings
  if(flag_restart.eq.2) then
     if(mpi_rank_max/nproc_mpi.eq.0) then
        print*, 'mpi_rank_max < # of MPI threads in backup file'
        call stop_calculation
     endif
     if( ((mpi_rank_max/nproc_mpi)-real(mpi_rank_max)/real(nproc_mpi)).ne.0.d0 ) then
        print*, 'mpi_rank_max in backup file / # of MPI threads is not an integer'
        call stop_calculation
     endif
     if(omp_rank_max.gt.100) then
        print*, 'omp_rank_max > 100 in backup file'
        call stop_calculation
     endif
  endif

  ! Read files written with the same number of MPI and OMP threads 
  if(flag_restart.eq.1) then

     write (pnum,'(i3)'),mpi_rank
     if(mpi_rank.le.9) lgh=3
     if(mpi_rank.ge.10 .and. mpi_rank.le.99 ) lgh=2
     if(mpi_rank.ge.100 .and. mpi_rank.le.999 ) lgh=1
     
     open(41+mpi_rank,file='DATA.BAK/particles'//pnum(lgh:3)//'.bak',form='UNFORMATTED')
     read(41+mpi_rank,err=999) time
     read(41+mpi_rank,err=999) np_tot(1:ntype,1:nproc)
     do iproc=1,nproc
        do ptype=1,ntype
           read(41+mpi_rank,err=999) vxp(:,1:np_tot(ptype,iproc),ptype,iproc)
        enddo
        if(tag_neg.gt.0) read(41+mpi_rank,err=999) flag_cex(1:np_tot(tag_neg,iproc),iproc)
     enddo
     close(41+mpi_rank)
     
  else ! Read files written with a different amount of MPI and OMP threads 

     mpi_rank_max= mpi_rank_max/nproc_mpi

     do irank= mpi_rank*mpi_rank_max,(mpi_rank+1)*mpi_rank_max-1 

        write (pnum,'(i3)'),irank
        if(irank.le.9) lgh=3
        if(irank.ge.10 .and. irank.le.99 ) lgh=2
        if(irank.ge.100 .and. irank.le.999 ) lgh=1
        
        open(41+mpi_rank,file='SAV_DATA/particles'//pnum(lgh:3)//'.bak',form='UNFORMATTED')
        print*, 'Reading external file SAV_DATA/particles'//pnum(lgh:3)//'.bak, MPI rank=',mpi_rank

        read(41+mpi_rank,err=999) time
        read(41+mpi_rank,err=999) np_tot_ext(1:ntype,1:omp_rank_max)

        iproc_tmp= 0
        do iproc=1,omp_rank_max

           if(omp_rank_max.eq.nproc) then 
              iproc_tmp= iproc
           else
              iproc_tmp= iproc_tmp + 1
              if(iproc_tmp.gt.nproc) iproc_tmp= 1
           endif

           do ptype=1,ntype
              n1(ptype)= np_tot(ptype,iproc_tmp)
              n2(ptype)= n1(ptype) + np_tot_ext(ptype,iproc)
              ! Warning 
              if(n2(ptype).gt.nmax) then 
                 print*, 'np_tot > nmax in restart(), please correct...'
                 call stop_calculation
              endif

              read(41+mpi_rank,err=999) vxp(:,n1(ptype)+1:n2(ptype),ptype,iproc_tmp)
              np_tot(ptype,iproc_tmp)= n2(ptype)
           enddo

           if(tag_neg.gt.0) read(41+mpi_rank,err=999) flag_cex(n1(tag_neg)+1:n2(tag_neg),iproc_tmp)

        enddo

        close(41+mpi_rank)
     enddo
     
  endif

  !
  ! Duplicate macroparticles
  !
  if(ABS(np_dup).gt.1) then     
     ! Warning
     if(ABS(np_dup)*MAXVAL(np_tot).gt.nmax) then
        print*, 'Warning: |np_dup|*np_tot > nmax in restart(). Please correct...'
        call stop_calculation
     endif

     do iproc=1,nproc
        do ptype=1,ntype
           do i=1,ABS(INT(np_dup))-1
              np_tot_tmp= np_tot(ptype,iproc)
              if(np_tot_tmp.gt.0) then
                 vxp(:,i*np_tot_tmp+1:(i+1)*np_tot_tmp,ptype,iproc)= vxp(:,1:np_tot_tmp,ptype,iproc)
              endif
           enddo
           np_tot(ptype,iproc)= ABS(np_dup)*np_tot(ptype,iproc)
        enddo
     enddo
  endif

  !
  ! Remove a fraction of macroparticles
  !
  if(ABS(np_dup).lt.1) then     
     !$OMP PARALLEL PRIVATE(iproc,iseed_OMP,i_shift,np_tot_tmp,rnd)
     iproc= omp_get_thread_num() + 1
     iseed_OMP= iseed(iproc)     
     
     do ptype=1,ntype
        np_lost(ptype,iproc)= 0
        np_tot_tmp= np_tot(ptype,iproc)
        do i=1,np_tot_tmp
           rnd= ran2(iseed_OMP)
           if( rnd.le.ABS(np_dup)) then
              i_shift= i-np_lost(ptype,iproc)
              vxp(:,i_shift,ptype,iproc)= vxp(:,i,ptype,iproc)
           else
              np_lost(ptype,iproc)= np_lost(ptype,iproc) + 1
           endif           
        enddo
        np_tot(ptype,iproc)= np_tot(ptype,iproc) - np_lost(ptype,iproc)
        vxp(:,np_tot(ptype,iproc)+1:np_tot_tmp,ptype,iproc)= 0.d0
     enddo

     iseed(iproc)= iseed_OMP  
     !$OMP END PARALLEL     
  endif

  
  !
  ! Loop over OpenMP processes
  !  
  !$OMP PARALLEL PRIVATE(iproc,ptype)
  iproc= omp_get_thread_num() + 1
  do ptype=1,ntype
     call restart_OMP(n,h,np,vxp,nmax,ptype,ntype,kq,iproc,nproc,&
          np_tot,sum_dEk,Nh)
  enddo
  !$OMP END PARALLEL

  return
  
  !
  ! Files not read correctly
  !
999 if(mpi_rank.eq.0) print*, 'Backup files were not read correctly!'  
  call stop_calculation
  
end subroutine restart

subroutine restart_OMP(n,h,np,vxp,nmax,ptype,ntype,kq,&
     iproc,nproc,np_tot,sum_dEk,Nh)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Mar/15
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: ip,ix,iy,iz,n(3),nmax,ntype,ptype,iproc,nproc,Nh(nproc)
  ! Particle arrays
  integer:: np_tot(ntype,nproc)
  real(kind=8):: h(3),vxp(6,nmax,ntype,nproc),sum_dEk(nproc),x,y,z,vx,vy,&
       vz,dEk,kp,ki(8),px,py,pz
  ! Particle densities
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc),&
       kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)

  !
  ! Initialize counters and arrays
  !
  np(:,:,:,ptype,iproc)=0.d0
  sum_dEk(iproc)=0.d0
  Nh(iproc)=0

  !
  ! Start iteration over particles
  !
  
  do ip=1,np_tot(ptype,iproc),1
     
     x= vxp(1,ip,ptype,iproc)
     y= vxp(2,ip,ptype,iproc)
     z= vxp(3,ip,ptype,iproc)
     vx= vxp(4,ip,ptype,iproc)
     vy= vxp(5,ip,ptype,iproc)
     vz= vxp(6,ip,ptype,iproc)
     
     !
     ! Calculate density
     !
     ix= INT( x/h(1) ) + 1
     iy= INT( y/h(2) ) + 1
     iz= INT( z/h(3) ) + 1
           
     px=( ix*h(1) - x )/h(1)
     py=( iy*h(2) - y )/h(2)
     pz=( iz*h(3) - z )/h(3)
           
     kp=Nm(ptype)/(h(1)*h(2)*h(3))
     ki(1)= kp*px*py*pz
     ki(2)= kp*(1.d0-px)*py*pz
     ki(3)= kp*(1.d0-px)*(1.d0-py)*pz
     ki(4)= kp*px*(1.d0-py)*pz
     ki(5)= kp*px*py*(1.d0-pz)
     ki(6)= kp*(1.d0-px)*py*(1.d0-pz)
     ki(7)= kp*(1.d0-px)*(1.d0-py)*(1.d0-pz)
     ki(8)= kp*px*(1.d0-py)*(1.d0-pz)
     
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
        
     if( ptype.eq.1 .and. Pabs.gt.0.d0 ) then
        if( ix.ge.ixl_pow .and. ix.le.ixr_pow ) then
           ! Calculate kinetic energy of macroparticle
           dEk= 0.5d0*Nm(ptype)*mass(ptype)*( vx*vx + vy*vy + vz*vz ) 
           sum_dEk(iproc)= sum_dEk(iproc) + dEk
           Nh(iproc)= Nh(iproc) + 1
        endif
     endif
        
  enddo ! end loop over particles
 
  return 
end subroutine restart_OMP


