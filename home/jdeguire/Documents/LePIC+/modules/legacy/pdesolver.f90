subroutine pdesolver(u,b,bcnd,h,n,ncycl,eps,omega,k,ktot,res,ng,rank,nproc)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      DEC/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Solve diffusion equation using a V-shaped
!                    multi-grid method. SOR algorithm is used for 
!                    relaxation.
!     NOTE:          u(x,y,t) is defined as u(0:nx+2,0:ny+2) where
!                    1 and n+1 are for the boundary conditions.
!     --------------------------------------------------------------
  use mod_constants, only: eps0
  use mod_utils, only: stop_calculation
  implicit none
  integer:: k,n(3),ng,ncycl,rank,nproc,n_mg(3,12), &
       bcnd(0:n(1)+2,0:n(2)+2,0:n(3)/nproc+2),flag
  real(kind=8):: h(3),h_mg(3,12),u(0:n(1)+2,0:n(2)+2,-1:n(3)/nproc+2), &
       b(0:n(1)+1,0:n(2)+1,0:n(3)/nproc+1),nr
  real(kind=8), allocatable:: e2(:,:,:),e4(:,:,:),e8(:,:,:),e16(:,:,:), &
       e32(:,:,:),e64(:,:,:),e128(:,:,:),e256(:,:,:),e512(:,:,:), &
       e1024(:,:,:),e2048(:,:,:),r2(:,:,:),r4(:,:,:),r8(:,:,:),r16(:,:,:), &
       r32(:,:,:),r64(:,:,:),r128(:,:,:),r256(:,:,:),r512(:,:,:), &
       r1024(:,:,:),r2048(:,:,:)
  integer, allocatable::  bcnd2(:,:,:),bcnd4(:,:,:),bcnd8(:,:,:), &
       bcnd16(:,:,:),bcnd32(:,:,:),bcnd64(:,:,:),bcnd128(:,:,:), &
       bcnd256(:,:,:),bcnd512(:,:,:),bcnd1024(:,:,:),bcnd2048(:,:,:)
  real(kind=8):: eps,omega,ktot,u0,alpha
  real(kind=8):: res,res_s ! residuals
  logical:: converged
  include 'mg.h'

  ! Parameters for convergence test
  u0= 1.d0 ! 1V
  alpha= ABS((n(1)+1)*(n(2)+1)*(n(3)+1)*ac)*u0

  ! Initialize variables and arrays
  ktot=0
  res_s=0.d0

  ! Define array length and mesh size in MG sub-levels
  flag=0
  do k=1,12
     n_mg(1,k)=n(1)/2**(k-1)
     n_mg(2,k)=n(2)/2**(k-1)
     n_mg(3,k)=n(3)/2**(k-1)/nproc
     nr=real(n(3))/real(2**(k-1))/real(nproc)
     if( (real(n_mg(3,k))-nr).eq.0.d0 .and. (real(n_mg(3,k)/2)-nr/2.d0).lt.0 ) flag=1
     if(flag.eq.1) n_mg(3,k)=n(3)/2**(k-1)

     h_mg(:,k)= h*2**(k-1)
  enddo

  allocate( e2(0:n_mg(1,2)+2,0:n_mg(2,2)+2,-1:n_mg(3,2)+2), &
       e4(0:n_mg(1,3)+2,0:n_mg(2,3)+2,-1:n_mg(3,3)+2), &
       e8(0:n_mg(1,4)+2,0:n_mg(2,4)+2,-1:n_mg(3,4)+2), &
       e16(0:n_mg(1,5)+2,0:n_mg(2,5)+2,-1:n_mg(3,5)+2), &
       e32(0:n_mg(1,6)+2,0:n_mg(2,6)+2,-1:n_mg(3,6)+2), &
       e64(0:n_mg(1,7)+2,0:n_mg(2,7)+2,-1:n_mg(3,7)+2), &
       e128(0:n_mg(1,8)+2,0:n_mg(2,8)+2,-1:n_mg(3,8)+2), &
       e256(0:n_mg(1,9)+2,0:n_mg(2,9)+2,-1:n_mg(3,9)+2), &
       e512(0:n_mg(1,10)+2,0:n_mg(2,10)+2,-1:n_mg(3,10)+2), &
       e1024(0:n_mg(1,11)+2,0:n_mg(2,11)+2,-1:n_mg(3,11)+2), &
       e2048(0:n_mg(1,12)+2,0:n_mg(2,12)+2,-1:n_mg(3,12)+2) )
  allocate( r2(0:n_mg(1,2)+1,0:n_mg(2,2)+1,0:n_mg(3,2)+1), &
       r4(0:n_mg(1,3)+1,0:n_mg(2,3)+1,0:n_mg(3,3)+1), &
       r8(0:n_mg(1,4)+1,0:n_mg(2,4)+1,0:n_mg(3,4)+1), &
       r16(0:n_mg(1,5)+1,0:n_mg(2,5)+1,0:n_mg(3,5)+1), &
       r32(0:n_mg(1,6)+1,0:n_mg(2,6)+1,0:n_mg(3,6)+1), &
       r64(0:n_mg(1,7)+1,0:n_mg(2,7)+1,0:n_mg(3,7)+1), &
       r128(0:n_mg(1,8)+1,0:n_mg(2,8)+1,0:n_mg(3,8)+1), &
       r256(0:n_mg(1,9)+1,0:n_mg(2,9)+1,0:n_mg(3,9)+1), &
       r512(0:n_mg(1,10)+1,0:n_mg(2,10)+1,0:n_mg(3,10)+1), &
       r1024(0:n_mg(1,11)+1,0:n_mg(2,11)+1,0:n_mg(3,11)+1), &
       r2048(0:n_mg(1,12)+1,0:n_mg(2,12)+1,0:n_mg(3,12)+1) )
  allocate( bcnd2(0:n_mg(1,2)+2,0:n_mg(2,2)+2,0:n_mg(3,2)+2), &
       bcnd4(0:n_mg(1,3)+2,0:n_mg(2,3)+2,0:n_mg(3,3)+2), &
       bcnd8(0:n_mg(1,4)+2,0:n_mg(2,4)+2,0:n_mg(3,4)+2), &
       bcnd16(0:n_mg(1,5)+2,0:n_mg(2,5)+2,0:n_mg(3,5)+2), &
       bcnd32(0:n_mg(1,6)+2,0:n_mg(2,6)+2,0:n_mg(3,6)+2), &
       bcnd64(0:n_mg(1,7)+2,0:n_mg(2,7)+2,0:n_mg(3,7)+2), &
       bcnd128(0:n_mg(1,8)+2,0:n_mg(2,8)+2,0:n_mg(3,8)+2), &
       bcnd256(0:n_mg(1,9)+2,0:n_mg(2,9)+2,0:n_mg(3,9)+2), &
       bcnd512(0:n_mg(1,10)+2,0:n_mg(2,10)+2,0:n_mg(3,10)+2), &
       bcnd1024(0:n_mg(1,11)+2,0:n_mg(2,11)+2,0:n_mg(3,11)+2), &
       bcnd2048(0:n_mg(1,12)+2,0:n_mg(2,12)+2,0:n_mg(3,12)+2) )

  ! Calculate Solution
  k=1
  converged=.FALSE.
  do while (.NOT.converged)
     
     call mg(u,b,bcnd,h_mg,res,n_mg,omega,ng,alpha*eps,k,ktot,rank,nproc,&
          e2,e4,e8,e16,e32,e64,e128,e256,e512,e1024,e2048,&
          r2,r4,r8,r16,r32,r64,r128,r256,r512,r1024,r2048,&
          bcnd2,bcnd4,bcnd8,bcnd16,bcnd32,bcnd64,bcnd128, &
          bcnd256,bcnd512,bcnd1024,bcnd2048) 
     res= res/alpha

     ! Some warnings 
     if(k.eq.1) then 
        res_s=res
     else
        if(ABS(res/res_s).gt.10.d0) then
           if(rank.eq.0) then
              print*, ' '
              print*, 'Warning: PDE solver is diverging'
              print*, 'Abort calculation ...'
              print*, 'Final sol: k=',k,'ratio=',res/res_s
           endif
           call stop_calculation
        endif
     endif

     ! Test for convergence
     if( res.lt.eps .or. k.eq.ncycl ) converged=.TRUE.

     ! Print info. on screen
     if(MOD(k,1000).eq.0) then
        if(rank.eq.0) then
           print*, 'k=',k,'res=',res
           print*, 'Equivalent SOR iterations:',nint(ktot)
        endif
     endif

     k= k+1

  enddo ! end dowhile

  deallocate( e2,e4,e8,e16,e32,e64,e128,e256,e512,e1024,e2048 )
  deallocate( r2,r4,r8,r16,r32,r64,r128,r256,r512,r1024,r2048 )
  deallocate( bcnd2,bcnd4,bcnd8,bcnd16,bcnd32,bcnd64,bcnd128, &
       bcnd256,bcnd512,bcnd1024,bcnd2048 )

  return
end subroutine pdesolver
