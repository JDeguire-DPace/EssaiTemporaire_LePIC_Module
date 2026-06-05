subroutine mg(u,b,bcnd,h,res,n,omega,ng,eps,k,ktot,rank,nproc,&
          e2,e4,e8,e16,e32,e64,e128,e256,e512,e1024,e2048,&
          r2,r4,r8,r16,r32,r64,r128,r256,r512,r1024,r2048,&
          bcnd2,bcnd4,bcnd8,bcnd16,bcnd32,bcnd64,bcnd128, &
          bcnd256,bcnd512,bcnd1024,bcnd2048)
   !     ==============================================================
   !     VERSION:         0.1
   !     LAST MOD:      DEC/14
   !     MOD AUTHOR:    G. Fubiani
   !     COMMENTS:      Solve poisson equation using a V-shaped
   !                    multi-grid method. SOR algorithm is used for 
   !                    relaxation.
   !     NOTE:          1) u(x,y,t) is defined as u(0:nx+2,0:ny+2) where
   !                    1 and n+1 are for the boundary conditions.
   !                    2) Graphical representation of MG scheme:
   !                     
   !                             Fine grid (n=8)
   !                    |---|---|---|---|---|---|---|---|
   !                    1   2   3   4   5   6   7   8   9
   !                    1(BC)                          n+1(BC)  
   !                             
   !                          First level coarsed grid
   !                    |-------|-------|-------|-------|
   !                    1       2       3       4       5
   !                    1(BC)                         n/2+1(BC)  
   !
   !     --------------------------------------------------------------
  use mod_constants, only: eps0
  implicit none
  integer:: i,k,ig,n(3,12),ng,ksor,iter(ng),il,rank,nproc, &
       bcnd(0:n(1,1)+2,0:n(2,1)+2,0:n(3,1)+2)
  real(kind=8):: h(3,12),u(0:n(1,1)+2,0:n(2,1)+2,-1:n(3,1)+2),&
       b(0:n(1,1)+1,0:n(2,1)+1,0:n(3,1)+1)
  ! Define multigrid arrays
  real(kind=8):: e2(0:n(1,2)+2,0:n(2,2)+2,-1:n(3,2)+2), &
       e4(0:n(1,3)+2,0:n(2,3)+2,-1:n(3,3)+2), &
       e8(0:n(1,4)+2,0:n(2,4)+2,-1:n(3,4)+2), &
       e16(0:n(1,5)+2,0:n(2,5)+2,-1:n(3,5)+2), &
       e32(0:n(1,6)+2,0:n(2,6)+2,-1:n(3,6)+2), &
       e64(0:n(1,7)+2,0:n(2,7)+2,-1:n(3,7)+2), &
       e128(0:n(1,8)+2,0:n(2,8)+2,-1:n(3,8)+2), &
       e256(0:n(1,9)+2,0:n(2,9)+2,-1:n(3,9)+2), &
       e512(0:n(1,10)+2,0:n(2,10)+2,-1:n(3,10)+2), &
       e1024(0:n(1,11)+2,0:n(2,11)+2,-1:n(3,11)+2), &
       e2048(0:n(1,12)+2,0:n(2,12)+2,-1:n(3,12)+2)
  real(kind=8):: r2(0:n(1,2)+1,0:n(2,2)+1,0:n(3,2)+1), &
       r4(0:n(1,3)+1,0:n(2,3)+1,0:n(3,3)+1), &
       r8(0:n(1,4)+1,0:n(2,4)+1,0:n(3,4)+1), &
       r16(0:n(1,5)+1,0:n(2,5)+1,0:n(3,5)+1), &
       r32(0:n(1,6)+1,0:n(2,6)+1,0:n(3,6)+1), &
       r64(0:n(1,7)+1,0:n(2,7)+1,0:n(3,7)+1), &
       r128(0:n(1,8)+1,0:n(2,8)+1,0:n(3,8)+1), &
       r256(0:n(1,9)+1,0:n(2,9)+1,0:n(3,9)+1), &
       r512(0:n(1,10)+1,0:n(2,10)+1,0:n(3,10)+1), &
       r1024(0:n(1,11)+1,0:n(2,11)+1,0:n(3,11)+1), &
       r2048(0:n(1,12)+1,0:n(2,12)+1,0:n(3,12)+1)
  integer:: bcnd2(0:n(1,2)+2,0:n(2,2)+2,0:n(3,2)+2), &
       bcnd4(0:n(1,3)+2,0:n(2,3)+2,0:n(3,3)+2), &
       bcnd8(0:n(1,4)+2,0:n(2,4)+2,0:n(3,4)+2), &
       bcnd16(0:n(1,5)+2,0:n(2,5)+2,0:n(3,5)+2), &
       bcnd32(0:n(1,6)+2,0:n(2,6)+2,0:n(3,6)+2), &
       bcnd64(0:n(1,7)+2,0:n(2,7)+2,0:n(3,7)+2), &
       bcnd128(0:n(1,8)+2,0:n(2,8)+2,0:n(3,8)+2), &
       bcnd256(0:n(1,9)+2,0:n(2,9)+2,0:n(3,9)+2), &
       bcnd512(0:n(1,10)+2,0:n(2,10)+2,0:n(3,10)+2), &
       bcnd1024(0:n(1,11)+2,0:n(2,11)+2,0:n(3,11)+2), &
       bcnd2048(0:n(1,12)+2,0:n(2,12)+2,0:n(3,12)+2)
  integer:: n1,n2,n4,n8,n16,n32,n64,n128,n256,n512,n1024,n2048
  parameter ( n1=1, n2=2, n4=3, n8=4, n16=5, n32=6, n64=7, &
       n128=8, n256=9, n512=10, n1024=11, n2048=12 )
 integer:: n3_1,n3_2,n3_4,n3_8,n3_16,n3_32,n3_64,n3_128,n3_256, &
       n3_512,n3_1024,n3_2048
  real(kind=8):: ktot,eps,omega
  real(kind=8):: res ! residual

  !
  ! Initialization
  !
  do i=1,ng
     iter(i)=0
  enddo

  n3_1=n(3,n1)*nproc
  n3_2=n3_1/2
  n3_4=n3_1/4
  n3_8=n3_1/8
  n3_16=n3_1/16
  n3_32=n3_1/32
  n3_64=n3_1/64
  n3_128=n3_1/128
  n3_256=n3_1/256
  n3_512=n3_1/512
  n3_1024=n3_1/1024
  n3_2048=n3_1/2048

  !
  ! Start iteration 
  !
  do ig=1,ng*2-1 ! 2*ng-1 for V shape iteration
     
     ! Level #1
     il=1
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(u,b,h(:,n1),bcnd,res,n(:,n1),n3_1,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(u,b,h(:,n1),bcnd,bcnd2,r2,e2,n(:,n1),n3_1,rank,nproc)
     endif

     ! Level #2
     il=2
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then
        call sor_rb(e2,r2,h(:,n2),bcnd2,res,n(:,n2),n3_2,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e2,r2,h(:,n2),bcnd2,bcnd4,r4,e4,n(:,n2),n3_2,rank,nproc)
        if( ig.gt.ng ) call prolongation(u,e2,bcnd,n(:,n1),n3_1,nproc)
     endif

     ! Level #3
     il=3
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then
        call sor_rb(e4,r4,h(:,n4),bcnd4,res,n(:,n4),n3_4,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e4,r4,h(:,n4),bcnd4,bcnd8,r8,e8,n(:,n4),n3_4,rank,nproc)
        if( ig.gt.ng ) call prolongation(e2,e4,bcnd2,n(:,n2),n3_2,nproc)
     endif

     ! Level #4
     il=4
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e8,r8,h(:,n8),bcnd8,res,n(:,n8),n3_8,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e8,r8,h(:,n8),bcnd8,bcnd16,r16,e16,n(:,n8),n3_8,rank,nproc)
        if( ig.gt.ng ) call prolongation(e4,e8,bcnd4,n(:,n4),n3_4,nproc)
     endif

     ! Level #5
     il=5
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e16,r16,h(:,n16),bcnd16,res,n(:,n16),n3_16,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e16,r16,h(:,n16),bcnd16,bcnd32,r32,e32,n(:,n16),n3_16,rank,nproc)
        if( ig.gt.ng ) call prolongation(e8,e16,bcnd8,n(:,n8),n3_8,nproc)
     endif

     ! Level #6
     il=6
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e32,r32,h(:,n32),bcnd32,res,n(:,n32),n3_32,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e32,r32,h(:,n32),bcnd32,bcnd64,r64,e64,n(:,n32),n3_32,rank,nproc)
        if( ig.gt.ng ) call prolongation(e16,e32,bcnd16,n(:,n16),n3_16,nproc)
     endif

     ! Level #7
     il=7
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e64,r64,h(:,n64),bcnd64,res,n(:,n64),n3_64,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e64,r64,h(:,n64),bcnd64,bcnd128,r128,e128,n(:,n64),n3_64,rank,nproc)
        if( ig.gt.ng ) call prolongation(e32,e64,bcnd32,n(:,n32),n3_32,nproc)
     endif

     ! Level #8
     il=8
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e128,r128,h(:,n128),bcnd128,res,n(:,n128),n3_128,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e128,r128,h(:,n128),bcnd128,bcnd256,r256,e256,n(:,n128),n3_128,rank,nproc)
        if( ig.gt.ng ) call prolongation(e64,e128,bcnd64,n(:,n64),n3_64,nproc)
     endif

     ! Level #9
     il=9
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e256,r256,h(:,n256),bcnd256,res,n(:,n256),n3_256,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e256,r256,h(:,n256),bcnd256,bcnd512,r512,e512,n(:,n256),n3_256,rank,nproc)
        if( ig.gt.ng ) call prolongation(e128,e256,bcnd128,n(:,n128),n3_128,nproc)
     endif

     ! Level #10
     il=10
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e512,r512,h(:,n512),bcnd512,res,n(:,n512),n3_512,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e512,r512,h(:,n512),bcnd512,bcnd1024,r1024,e1024,n(:,n512),n3_512,rank,nproc)
        if( ig.gt.ng ) call prolongation(e256,e512,bcnd256,n(:,n256),n3_256,nproc)
     endif

     ! Level #11
     il=11
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e1024,r1024,h(:,n1024),bcnd1024,res,n(:,n1024),n3_1024,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e1024,r1024,h(:,n1024),bcnd1024,bcnd2048,r2048,e2048,n(:,n1024),n3_1024,rank,nproc)
        if( ig.gt.ng ) call prolongation(e512,e1024,bcnd512,n(:,n512),n3_512,nproc)
     endif

     ! Level max.
     if( ig.eq.ng ) then
        if(ng.eq.1) call sor_rb(u,b,h(:,n1),bcnd,res,n(:,n1),n3_1,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.2) call sor_rb(e2,r2,h(:,n2),bcnd2,res,n(:,n2),n3_2,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.3) call sor_rb(e4,r4,h(:,n4),bcnd4,res,n(:,n4),n3_4,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.4) call sor_rb(e8,r8,h(:,n8),bcnd8,res,n(:,n8),n3_8,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.5) call sor_rb(e16,r16,h(:,n16),bcnd16,res,n(:,n16),n3_16,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.6) call sor_rb(e32,r32,h(:,n32),bcnd32,res,n(:,n32),n3_32,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.7) call sor_rb(e64,r64,h(:,n64),bcnd64,res,n(:,n64),n3_64,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.8) call sor_rb(e128,r128,h(:,n128),bcnd128,res,n(:,n128),n3_128,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.9) call sor_rb(e256,r256,h(:,n256),bcnd256,res,n(:,n256),n3_256,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.10) call sor_rb(e512,r512,h(:,n512),bcnd512,res,n(:,n512),n3_512,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.11) call sor_rb(e1024,r1024,h(:,n1024),bcnd1024,res,n(:,n1024),n3_1024,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.12) call sor_rb(e2048,r2048,h(:,n2048),bcnd2048,res,n(:,n2048),n3_2048,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(ig)= iter(ig) + ksor
           
        if(ng.eq.2) call prolongation(u,e2,bcnd,n(:,n1),n3_1,nproc)
        if(ng.eq.3) call prolongation(e2,e4,bcnd2,n(:,n2),n3_2,nproc)
        if(ng.eq.4) call prolongation(e4,e8,bcnd4,n(:,n4),n3_4,nproc)
        if(ng.eq.5) call prolongation(e8,e16,bcnd8,n(:,n8),n3_8,nproc)
        if(ng.eq.6) call prolongation(e16,e32,bcnd16,n(:,n16),n3_16,nproc)
        if(ng.eq.7) call prolongation(e32,e64,bcnd32,n(:,n32),n3_32,nproc)
        if(ng.eq.8) call prolongation(e64,e128,bcnd64,n(:,n64),n3_64,nproc)
        if(ng.eq.9) call prolongation(e128,e256,bcnd128,n(:,n128),n3_128,nproc)
        if(ng.eq.10) call prolongation(e256,e512,bcnd256,n(:,n256),n3_256,nproc)
        if(ng.eq.11) call prolongation(e512,e1024,bcnd512,n(:,n512),n3_512,nproc)
        if(ng.eq.12) call prolongation(e1024,e2048,bcnd1024,n(:,n1024),n3_1024,nproc)
     endif
     
  enddo
  
  do i=1,ng
     ktot= ktot + real(iter(i))/2.d0**(2*(i-1))
  enddo
  
  return
end subroutine mg


subroutine restriction(u,b,h,bcnd,bcnd2h,r2h,e2h,n,n3,rank,nproc)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      DEC/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      h -> 2h grid mapping using a 27 point trilinear 
!                    interpolation
!     NOTE:          1) We do not need to compute 1 and n/2+1
!                    2) For the stencyl [a() array] restriction, see. 
!                       A Multigrid Tutorial, 2nd Ed., SIAM, p. 130-131
!                       Note that their domain range from 0 to n 
!                       (we use 1 to n+1) consequently we must at 1 to 
!                       each indexes in the book in order to properly 
!                       convert to our conventions. 
!
!                          W   E   W   E   W   E   W 
!                      W   E   W   E   W   E   W   E
!                    |-*-|-*-|-*-|-*-|-*-|-*-|-*-|-*-|
!                    1   2   3   4   5   6   7   8   9
!                    1(BC)      Fine grid           n+1(BC)  
!                             
!                        W       E       W       E
!                    |---*---|---*---|---*---|---*---|
!                    1       2       3       4       5
!                    1(BC)     Coarsed grid       n/2+1(BC)  
! 
!                    Legend= W==aw, E==ae, ac= -( aw + ae )
!                            aw==a^i-1/2, ae==a^i+1/2
!
!     ---------------------------------------------------------------
  use mpi
  use mod_constants, only: eps0
  use mod_part_info, only: flag_pbc, flag_nmn
  implicit none
  integer:: i,j,k,ir,jr,kr,km,kp,n(3),n3,m,flag
  integer rank,nproc,ierr,status(MPI_STATUS_SIZE)
  real(kind=8):: h(3),u(0:n(1)+2,0:n(2)+2,-1:n(3)+2), &
       b(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
       e2h(0:n(1)/2+2,0:n(2)/2+2,-1:n(3)/2+2), &
       r2h(0:n(1)/2+1,0:n(2)/2+1,0:n(3)/2+1), &
       r_pbc(n(1)+1,n(2)+1),mr
  integer:: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2), &
       bcnd2h(0:n(1)/2+2,0:n(2)/2+2,0:n(3)/2+2)
  real(kind=8):: re,r000,rpp0,rp00,rpm0,r0m0,rmm0,rm00,rmp0,r0p0, &
       r00p,rppp,rp0p,rpmp,r0mp,rmmp,rm0p,rmpp,r0pp, &
       r00m,rppm,rp0m,rpmm,r0mm,rmmm,rm0m,rmpm,r0pm
  include 'mg.h'

  !
  ! Compute size of local block
  !
  m= n3/nproc
  mr= real(n3)/real(nproc)

  flag=0
  if((real(m/2)-(mr/2.d0)).eq.0.d0) flag=1

  km= 1
  if(flag.eq.1) then 
     ! At least 1 node per proc  
     if( (rank+1).eq.nproc ) then    
        kp= m/2+1
     else
        kp= m/2
     endif
  else ! Less than 1 node per proc
     m=n3
     kp= n3/2+1
  endif

  !
  ! Calculate redidual @ k=n(3) for periodic BCs
  !
  if(flag_pbc.eq.0) goto 90

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ir,jr,re)
  !$OMP DO
  do j=1,n(2)/2+1
     do i=1,n(1)/2+1
        ir=2*i-1
        jr=2*j-1

        if(bcnd(ir,jr,m+1).le.0) then

           if( ir.ne.1 .and. jr.ne.1 ) then 
              call getres(u,b,h,ir-1,jr-1,m,n,re,flag_nmn)
              r_pbc(ir-1,jr-1)= re
           endif

           if( ir.ne.1 ) then
              call getres(u,b,h,ir-1,jr,m,n,re,flag_nmn)
              r_pbc(ir-1,jr)= re
           endif

           if( ir.ne.1 .and. jr.ne.(n(2)+1) ) then
              call getres(u,b,h,ir-1,jr+1,m,n,re,flag_nmn)
              r_pbc(ir-1,jr+1)= re
           endif

           if( ir.ne.1 .and. jr.ne.1 ) then
              call getres(u,b,h,ir,jr-1,m,n,re,flag_nmn)
              r_pbc(ir,jr-1)= re
           endif

           call getres(u,b,h,ir,jr,m,n,re,flag_nmn)
           r_pbc(ir,jr)= re

           if( jr.ne.(n(2)+1) ) then
              call getres(u,b,h,ir,jr+1,m,n,re,flag_nmn)
              r_pbc(ir,jr+1)= re
           endif

           if( ir.ne.(n(1)+1) .and. jr.ne.1 ) then
              call getres(u,b,h,ir+1,jr-1,m,n,re,flag_nmn)
              r_pbc(ir+1,jr-1)= re
           endif

           if( ir.ne.(n(1)+1) ) then
              call getres(u,b,h,ir+1,jr,m,n,re,flag_nmn)
              r_pbc(ir+1,jr)= re
           endif

           if( ir.ne.(n(1)+1) .and. jr.ne.(n(2)+1) ) then
              call getres(u,b,h,ir+1,jr+1,m,n,re,flag_nmn)
              r_pbc(ir+1,jr+1)= re
           endif

        else
           
           if( ir.ne.1 .and. jr.ne.1 ) r_pbc(ir-1,jr-1)= 0.d0
           if(ir.ne.1) r_pbc(ir-1,jr)= 0.d0
           if( ir.ne.1 .and. jr.ne.n(2)+1 ) r_pbc(ir-1,jr+1)= 0.d0                      
           if(jr.ne.1) r_pbc(ir,jr-1)= 0.d0           
           r_pbc(ir,jr)= 0.d0           
           if(jr.ne.n(2)+1) r_pbc(ir,jr+1)= 0.d0          
           if( ir.ne.n(1)+1 .and. jr.ne.1 ) r_pbc(ir+1,jr-1)= 0.d0           
           if(ir.ne.n(1)+1) r_pbc(ir+1,jr)= 0.d0           
           if( ir.ne.n(1)+1 .and. jr.ne.n(2)+1 ) r_pbc(ir+1,jr+1)= 0.d0
           
        endif
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  if( nproc.gt.1 .and. flag.eq.1 ) then
     if(rank.eq.(nproc-1)) then
        ! Send periodic BC r(:,:,m) to rank #0
        call MPI_SEND(r_pbc,(n(1)+1)*(n(2)+1), MPI_REAL8, 0, &
             2, MPI_COMM_WORLD, ierr)
     endif
     if(rank.eq.0) then
        ! Received periodic value from nproc-1 : r(:,:,0)=r(:,:,m)
        call MPI_RECV(r_pbc,(n(1)+1)*(n(2)+1), MPI_REAL8, nproc-1, &
             2, MPI_COMM_WORLD, status, ierr)
     endif
  endif

90 continue

  !
  ! Trilinear 27 point interpolation:
  !

  ! Sweep the coarse grid

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,ir,jr,kr,re, &
  !$OMP       r000,rpp0,rp00,rpm0,r0m0,rmm0,rm00,rmp0,r0p0, &
  !$OMP       r00p,rppp,rp0p,rpmp,r0mp,rmmp,rm0p,rmpp,r0pp, &
  !$OMP       r00m,rppm,rp0m,rpmm,r0mm,rmmm,rm0m,rmpm,r0pm)
  !$OMP DO COLLAPSE(2)
  do k=km,kp
     do j=1,n(2)/2+1
        do i=1,n(1)/2+1
        
           ! Address of the fine grid point
           ir=2*i-1
           jr=2*j-1
           kr=2*k-1

           ! Update bcnd() array for higher grid level
           bcnd2h(i,j,k)=bcnd(ir,jr,kr)

           if(ir.eq.1) bcnd2h(0,j,k)=bcnd(0,jr,kr)
           if(ir.eq.n(1)+1) bcnd2h(n(1)/2+2,j,k)=bcnd(n(1)+2,jr,kr)
           if(jr.eq.1) bcnd2h(i,0,k)=bcnd(ir,0,kr)
           if(jr.eq.n(2)+1) bcnd2h(i,n(2)/2+2,k)=bcnd(ir,n(2)+2,kr)
           if(kr.eq.1) bcnd2h(i,j,0)=bcnd(ir,jr,0)
           if(kr.eq.m+1) bcnd2h(i,j,m/2+2)=bcnd(ir,jr,m+2)


           ! Computes the residuals around fine grid point
           if(bcnd2h(i,j,k).le.0) then

              ! Initialization
              r000=0.d0 
              rpp0=0.d0
              rp00=0.d0
              rpm0=0.d0
              r0m0=0.d0
              rmm0=0.d0
              rm00=0.d0
              rmp0=0.d0
              r0p0=0.d0
              r00p=0.d0
              rppp=0.d0
              rp0p=0.d0
              rpmp=0.d0
              r0mp=0.d0
              rmmp=0.d0
              rm0p=0.d0
              rmpp=0.d0
              r0pp=0.d0
              r00m=0.d0
              rppm=0.d0
              rp0m=0.d0
              rpmm=0.d0
              r0mm=0.d0
              rmmm=0.d0
              rm0m=0.d0
              rmpm=0.d0
              r0pm=0.d0

              ! 1: 000
              call getres(u,b,h,ir,jr,kr,n,r000,flag_nmn)
              
              ! 2: ++0
              if( ir.ne.(n(1)+1) .and. jr.ne.(n(2)+1) ) then
                 if(bcnd(ir+1,jr+1,kr).le.0) &
                      call getres(u,b,h,ir+1,jr+1,kr,n,rpp0,flag_nmn)
              endif

              ! 3: +00
              if( ir.ne.(n(1)+1) ) then
                 if(bcnd(ir+1,jr,kr).le.0) &
                      call getres(u,b,h,ir+1,jr,kr,n,rp00,flag_nmn)
              endif

              ! 4: +-0
              if( bcnd(ir+1,jr-1,kr).le.0 ) then
                 if( ir.ne.(n(1)+1) .and. jr.ne.1 ) then
                    call getres(u,b,h,ir+1,jr-1,kr,n,rpm0,flag_nmn)
                 else
                    if( jr.eq.1 ) &
                         call getres(u,b,h,ir+1,n(2),kr,n,rpm0,flag_nmn)
                 endif
              endif

              ! 5: 0-0
              if( bcnd(ir,jr-1,kr).le.0 ) then
                 if( jr.eq.1 ) then
                    call getres(u,b,h,ir,n(2),kr,n,r0m0,flag_nmn)
                 else
                    call getres(u,b,h,ir,jr-1,kr,n,r0m0,flag_nmn)
                 endif
              endif

              ! 6: --0
              if( bcnd(ir-1,jr-1,kr).le.0 .and. ir.ne.1 ) then
                 if(jr.ne.1) then
                    call getres(u,b,h,ir-1,jr-1,kr,n,rmm0,flag_nmn)
                 else
                    call getres(u,b,h,ir-1,n(2),kr,n,rmm0,flag_nmn)
                 endif
              endif
              
              ! 7: -00
              if( ir.ne.1 ) then
                 if(bcnd(ir-1,jr,kr).le.0) &
                      call getres(u,b,h,ir-1,jr,kr,n,rm00,flag_nmn)
              endif
              
              ! 8: -+0
              if( ir.ne.1 .and. jr.ne.(n(2)+1) ) then
                 if(bcnd(ir-1,jr+1,kr).le.0) &
                      call getres(u,b,h,ir-1,jr+1,kr,n,rmp0,flag_nmn)
              endif

              ! 9: 0+0
              if( jr.ne.(n(2)+1) ) then
                 if(bcnd(ir,jr+1,kr).le.0) &
                      call getres(u,b,h,ir,jr+1,kr,n,r0p0,flag_nmn)
              endif
              
              ! 10: 00+
              if( kr.ne.(m+1) ) then
                 if(bcnd(ir,jr,kr+1).le.0) &
                      call getres(u,b,h,ir,jr,kr+1,n,r00p,flag_nmn)
              endif
              
              ! 11: +++
              if( ir.ne.(n(1)+1) .and. jr.ne.(n(2)+1) .and. &
                   kr.ne.(m+1) ) then
                 if(bcnd(ir+1,jr+1,kr+1).le.0) &
                      call getres(u,b,h,ir+1,jr+1,kr+1,n,rppp,flag_nmn)
              endif

              ! 12: +0+
              if( ir.ne.(n(1)+1) .and. kr.ne.(m+1) ) then
                 if(bcnd(ir+1,jr,kr+1).le.0) &
                      call getres(u,b,h,ir+1,jr,kr+1,n,rp0p,flag_nmn)
              endif

              ! 13: +-+
              if(bcnd(ir+1,jr-1,kr+1).le.0)  then
                 if( ir.ne.(n(1)+1) .and. jr.ne.1 .and. &
                      kr.ne.(m+1) ) then
                    call getres(u,b,h,ir+1,jr-1,kr+1,n,rpmp,flag_nmn)
                 else
                    if( jr.eq.1 .and. kr.ne.(m+1) ) &
                         call getres(u,b,h,ir+1,n(2),kr+1,n,rpmp,flag_nmn)
                 endif   
              endif

              ! 14: 0-+
              if(bcnd(ir,jr-1,kr+1).le.0)  then
                 if( jr.ne.1  .and. kr.ne.(m+1) ) then
                    call getres(u,b,h,ir,jr-1,kr+1,n,r0mp,flag_nmn)
                 else
                    if( jr.eq.1 .and. kr.ne.(m+1) ) &
                         call getres(u,b,h,ir,n(2),kr+1,n,r0mp,flag_nmn)
                 endif
              endif

              ! 15: --+
              if( bcnd(ir-1,jr-1,kr+1).le.0 .and. ir.ne.1 ) then 
                 if( jr.ne.1 .and. kr.ne.(m+1) ) then
                    call getres(u,b,h,ir-1,jr-1,kr+1,n,rmmp,flag_nmn)
                 else
                    if( jr.eq.1 .and. kr.ne.(m+1) ) &
                         call getres(u,b,h,ir-1,n(2),kr+1,n,rmmp,flag_nmn)
                 endif
              endif
              
              ! 16: -0+
              if( ir.ne.1 .and. kr.ne.(m+1)  ) then
                 if(bcnd(ir-1,jr,kr+1).le.0) &
                      call getres(u,b,h,ir-1,jr,kr+1,n,rm0p,flag_nmn)
              endif

              ! 17: -++
              if( ir.ne.1 .and. jr.ne.(n(2)+1) .and. &
                   kr.ne.(m+1)  ) then
                 if(bcnd(ir-1,jr+1,kr+1).le.0) &
                      call getres(u,b,h,ir-1,jr+1,kr+1,n,rmpp,flag_nmn)
              endif

              ! 18: 0++
              if( jr.ne.(n(2)+1) .and. kr.ne.(m+1)  ) then
                 if(bcnd(ir,jr+1,kr+1).le.0) &
                      call getres(u,b,h,ir,jr+1,kr+1,n,r0pp,flag_nmn)
              endif

              ! 19: 00-
              if(bcnd(ir,jr,kr-1).le.0) then
                 if( kr.eq.1 .and. (rank.eq.0 .or. flag.eq.0) ) then
                    r00m= r_pbc(ir,jr)
                 else
                    call getres(u,b,h,ir,jr,kr-1,n,r00m,flag_nmn)
                 endif
              endif
              
              ! 20: ++-
              if(bcnd(ir+1,jr+1,kr-1).le.0)  then
                 if( kr.eq.1 .and. (rank.eq.0 .or. flag.eq.0) ) then
                    if( ir.ne.(n(1)+1) .and. jr.ne.(n(2)+1) ) &
                         rppm= r_pbc(ir+1,jr+1)
                 else
                    if( ir.ne.(n(1)+1) .and. jr.ne.(n(2)+1) ) &
                         call getres(u,b,h,ir+1,jr+1,kr-1,n,rppm,flag_nmn)                      
                 endif
              endif

              ! 21: +0-
              if(bcnd(ir+1,jr,kr-1).le.0)  then
                 if( kr.eq.1 .and. (rank.eq.0 .or. flag.eq.0) ) then
                   if(ir.ne.(n(1)+1)) rp0m= r_pbc(ir+1,jr)
                 else
                   if(ir.ne.(n(1)+1)) call getres(u,b,h,ir+1,jr,kr-1,n,rp0m,flag_nmn)
                 endif
              endif

              ! 22: +--
              if(bcnd(ir+1,jr-1,kr-1).le.0) then
                 if( kr.eq.1 .and. (rank.eq.0 .or. flag.eq.0) ) then
                    if( ir.ne.(n(1)+1) .and. jr.ne.1 ) rpmm= r_pbc(ir+1,jr-1)
                 else
                    if( ir.ne.(n(1)+1) .and. jr.ne.1 ) &
                         call getres(u,b,h,ir+1,jr-1,kr-1,n,rpmm,flag_nmn)
                    if( jr.eq.1 ) &
                         call getres(u,b,h,ir+1,n(2),kr-1,n,rpmm,flag_nmn)
                 endif
              endif

              ! 23: 0--
              if(bcnd(ir,jr-1,kr-1).le.0) then
                 if( kr.eq.1 .and. (rank.eq.0 .or. flag.eq.0) ) then
                    if(jr.ne.1) r0mm= r_pbc(ir,jr-1)
                 else
                    if(jr.ne.1) &
                         call getres(u,b,h,ir,jr-1,kr-1,n,r0mm,flag_nmn)
                    if(jr.eq.1) &
                         call getres(u,b,h,ir,n(2),kr-1,n,r0mm,flag_nmn)
                 endif
              endif

              ! 24: ---
              if( bcnd(ir-1,jr-1,kr-1).le.0 .and. ir.ne.1 ) then
                 if( kr.eq.1 .and. (rank.eq.0 .or. flag.eq.0) ) then
                    if(jr.ne.1) rmmm= r_pbc(ir-1,jr-1)
                 else
                    if(jr.ne.1) &
                         call getres(u,b,h,ir-1,jr-1,kr-1,n,rmmm,flag_nmn)
                    if(jr.eq.1) &
                         call getres(u,b,h,ir-1,n(2),kr-1,n,rmmm,flag_nmn)
                 endif
              endif
              
              ! 25: -0-
              if( bcnd(ir-1,jr,kr-1).le.0 .and. ir.ne.1 ) then
                 if( kr.eq.1 .and. (rank.eq.0 .or. flag.eq.0) ) then
                    rm0m= r_pbc(ir-1,jr)
                 else
                    call getres(u,b,h,ir-1,jr,kr-1,n,rm0m,flag_nmn)
                 endif
              endif

              ! 26: -+-
              if(bcnd(ir-1,jr+1,kr-1).le.0 .and. ir.ne.1 ) then
                 if( kr.eq.1 .and. (rank.eq.0 .or. flag.eq.0) ) then
                    if(jr.ne.(n(2)+1)) &
                         rmpm= r_pbc(ir-1,jr+1)
                 else
                    if(jr.ne.(n(2)+1)) &
                         call getres(u,b,h,ir-1,jr+1,kr-1,n,rmpm,flag_nmn)
                 endif
              endif

              ! 27: 0+-
              if(bcnd(ir,jr+1,kr-1).le.0) then
                 if( kr.eq.1 .and. (rank.eq.0 .or. flag.eq.0) ) then
                    if(jr.ne.(n(2)+1)) r0pm= r_pbc(ir,jr+1)      
                 else
                    if(jr.ne.(n(2)+1)) call getres(u,b,h,ir,jr+1,kr-1,n,r0pm,flag_nmn)
                 endif
              endif

              ! Compute residual 
              re=( 8.d0*r000 + &
                   4.d0*( rp00 + rm00 + r0p0 + r0m0 + r00p + r00m ) + &
                   2.d0*( rpp0 + rmp0 + rpm0 + rmm0 + &
                          rp0p + rm0p + r0pp + r0mp + rp0m + &
                          rm0m + r0pm + r0mm ) + &
                          rppp + rmpp + rpmp + rmmp + &
                          rppm + rmpm + rpmm + rmmm )/64.0d0
              r2h(i,j,k)=re

           else

              r2h(i,j,k)=0.d0
           
           endif
           
           ! Zero initial solution on coarse grid
           e2h(i,j,k)=0.0d0

        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  ! Zero the remaining unknowns on the coarse grid
  if(flag.eq.1) then
     e2h(0,:,0:m/2+2)=0.d0 ! i=0 plane
     e2h(n(1)/2+2,:,0:m/2+2)=0.d0 ! i= nx/2+2 plane
     
     e2h(:,0,0:m/2+2)=0.d0  ! j=0 plane
     e2h(:,n(2)/2+2,0:m/2+2)=0.d0! j= ny/2+2 plane
     
     e2h(:,:,0)=0.d0 ! k=0 plane
     e2h(:,:,m/2+1:m/2+2)=0.d0 ! k= nz/2+2 plane
  else
     e2h(0,:,:)=0.d0 ! i=0 plane
     e2h(n(1)/2+2,:,:)=0.d0 ! i= nx/2+2 plane
     
     e2h(:,0,:)=0.d0  ! j=0 plane
     e2h(:,n(2)/2+2,:)=0.d0! j= ny/2+2 plane
     
     e2h(:,:,0)=0.d0 ! k=0 plane
     e2h(:,:,n(3)/2+2)=0.d0 ! k= nz/2+2 plane
  endif
 
  ! Skip inter-node communications for m<1
  if( nproc.gt.1 .and. flag.eq.0 ) goto 100 

  !
  ! Send ghost nodes for r2h to processes of higher rank 
  !
  if(MOD(rank,2).eq.1) then ! rank is odd
     ! Received ghost value from rank-1 stored in r2h(:,:,k=0)
     call MPI_RECV(r2h(0,0,0),(n(1)/2+2)*(n(2)/2+2), MPI_REAL8, rank-1, &
          0, MPI_COMM_WORLD, status, ierr)
     if(rank.lt.nproc-1) then
        ! Send ghost value r2h(:,:,k=m) to rank+1
        call MPI_SEND(r2h(0,0,m/2),(n(1)/2+2)*(n(2)/2+2), MPI_REAL8, rank+1, &
             0, MPI_COMM_WORLD, ierr)
     endif
  else ! Rank is even 
     if(rank.gt.0) then
        call MPI_RECV(r2h(0,0,0),(n(1)/2+2)*(n(2)/2+2), MPI_REAL8, rank-1, &
             0, MPI_COMM_WORLD, status, ierr)
     endif
     if(rank.lt.nproc-1) then
        call MPI_SEND(r2h(0,0,m/2),(n(1)/2+2)*(n(2)/2+2), MPI_REAL8, rank+1, &
                0, MPI_COMM_WORLD, ierr)
     endif
  endif

  !
  ! Send/receive ghost nodes for bcnd2h to/from neighbour processes
  !
  if(MOD(rank,2).eq.1) then ! rank is odd
     ! Send ghost value bcnd2h(:,j=1) to rank-1
     call MPI_SEND(bcnd2h(0,0,1),(n(1)/2+3)*(n(2)/2+3), MPI_INTEGER, rank-1, &
          0, MPI_COMM_WORLD, ierr)
     ! Received ghost value from rank-1 stored in bcnd2h(:,-1:0)
     call MPI_RECV(bcnd2h(0,0,0),(n(1)/2+3)*(n(2)/2+3), MPI_INTEGER, rank-1, &
          0, MPI_COMM_WORLD, status, ierr)
     
     if(rank.lt.nproc-1) then
        ! Send ghost value bcnd2h(:,m-1:m) to rank+1
        call MPI_SEND(bcnd2h(0,0,m/2),(n(1)/2+3)*(n(2)/2+3), MPI_INTEGER, rank+1, &
             0, MPI_COMM_WORLD, ierr)
        ! Received ghost value from rank+1 stored in bcnd2h(:,j=m+1)
        call MPI_RECV(bcnd2h(0,0,m/2+1),(n(1)/2+3)*(n(2)/2+3), MPI_INTEGER, rank+1, &
             0, MPI_COMM_WORLD, status, ierr)
     endif
  else ! Rank is even 
     if(rank.gt.0) then
        call MPI_RECV(bcnd2h(0,0,0),(n(1)/2+3)*(n(2)/2+3), MPI_INTEGER, rank-1, & 
                0, MPI_COMM_WORLD, status, ierr)
        call MPI_SEND(bcnd2h(0,0,1),(n(1)/2+3)*(n(2)/2+3), MPI_INTEGER, rank-1, &
             0, MPI_COMM_WORLD, ierr)
     endif
     
        if(rank.lt.nproc-1) then
           call MPI_RECV(bcnd2h(0,0,m/2+1),(n(1)/2+3)*(n(2)/2+3), MPI_INTEGER, rank+1, &
                0, MPI_COMM_WORLD, status, ierr)
           call MPI_SEND(bcnd2h(0,0,m/2),(n(1)/2+3)*(n(2)/2+3), MPI_INTEGER, rank+1, &
                0, MPI_COMM_WORLD, ierr)
        endif
     endif

100 continue

  return
end subroutine restriction


subroutine getres(uorg,rhsorg,h,ir,jr,kr,n,r,flag)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      DEC/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Calculate residual
!     NOTE:          
!     --------------------------------------------------------------
  use mod_constants, only: eps0
  implicit none
  integer:: ir,jr,kr,n(3),flag
  real(kind=8):: h(3), &
       uorg(0:n(1)+2,0:n(2)+2,-1:n(3)+2), &
       rhsorg(0:n(1)+1,0:n(2)+1,0:n(3)+1)
  real(kind=8):: uij,ue,uw,un,us,ut,ub,urhs,r
  include 'mg.h'

  ! Data points used in the calculation of the residual
  uij=uorg(ir,jr,kr)
  un=uorg(ir,jr+1,kr)
  us=uorg(ir,jr-1,kr)
  ue=uorg(ir+1,jr,kr)
  ut=uorg(ir,jr,kr+1)
  ub=uorg(ir,jr,kr-1)
  if( ir.eq.1 .and. flag.eq.1 ) then
     ! Neumann BC's (LHS)
     uw=0.d0
     aw=0.d0
     ae=2.d0*ae
  else
     uw=uorg(ir-1,jr,kr)
  endif
  urhs=rhsorg(ir,jr,kr)
  
  ! ...and the residual
  r= urhs - ( aw*uw + ae*ue + an*un + as*us + ac*uij &
       + at*ut + ab*ub )

  return
end subroutine getres

subroutine prolongation(u,e2h,bcnd,n,n3,nproc)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      DEC/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      2h -> h grid mapping using a 3D tricubic 
!                    (Full Weighted) interpolation      
!     NOTE:          We do not need to compute 0 and n/2+1  
!     --------------------------------------------------------------
  implicit none
  integer:: i,j,k,ir,jr,kr,km,kp,n(3),n3,m,nproc
  real(kind=8):: u(0:n(1)+2,0:n(2)+2,-1:n(3)+2), &
       e2h(0:n(1)/2+2,0:n(2)/2+2,-1:n(3)/2+2), &
       r000,rm00,rmm0,r0m0,r00m,r0mm,rmmm,rm0m,mr
  integer:: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)

  !
  ! Compute size of local block
  !
  m=(n3/2)/nproc
  mr= (real(n3)/2.d0)/real(nproc)

  km= 1
  if( (real(m)-mr).eq.0.d0 ) then         
     kp= m+1
  else ! Less than 1 node per proc
     kp= n3/2+1
  endif

  !
  ! Sweep the coarse grid
  !

  !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(i,j,k,ir,jr,kr,r000,rm00, &
  !$OMP    rmm0,r0m0,r00m,r0mm,rmmm,rm0m)
  !$OMP DO COLLAPSE(2)
  do k=km,kp
     do j=1,n(2)/2+1
        do i=1,n(1)/2+1
        
           ! Load all the required data points
           r000=e2h(i,j,k)
           r0m0=e2h(i,j-1,k)
           rmm0=e2h(i-1,j-1,k)
           rm00=e2h(i-1,j,k)
           r00m=e2h(i,j,k-1)
           r0mm=e2h(i,j-1,k-1)
           rmmm=e2h(i-1,j-1,k-1)
           rm0m=e2h(i-1,j,k-1)

           ! Compute fine grid address
           ir=2*i-1
           jr=2*j-1
           kr=2*k-1

           ! Prolong
           if(bcnd(ir,jr,kr).le.0) u(ir,jr,kr)= u(ir,jr,kr) + r000
           
           if( bcnd(ir,jr-1,kr).le.0 ) &
                u(ir,jr-1,kr)= u(ir,jr-1,kr) + (r000+r0m0)/2.0d0        

           if( bcnd(ir-1,jr-1,kr).le.0 ) &
                u(ir-1,jr-1,kr)= u(ir-1,jr-1,kr) + &
                (r000+r0m0+rmm0+rm00)/4.0d0
           
           if(bcnd(ir-1,jr,kr).le.0 ) &
                u(ir-1,jr,kr)= u(ir-1,jr,kr) + (r000+rm00)/2.0d0 

           if(bcnd(ir,jr,kr-1).le.0 ) &
                u(ir,jr,kr-1)= u(ir,jr,kr-1) + (r000+r00m)/2.0d0 

           if(bcnd(ir,jr-1,kr-1).le.0 ) &
                u(ir,jr-1,kr-1)= u(ir,jr-1,kr-1) + &
                (r000+r00m+r0mm+r0m0)/4.0d0 

           if(bcnd(ir-1,jr,kr-1).le.0 ) &
                u(ir-1,jr,kr-1)= u(ir-1,jr,kr-1) + &
                (r000+r00m+rm0m+rm00)/4.0d0 
 
           if(bcnd(ir-1,jr-1,kr-1).le.0 ) &
                u(ir-1,jr-1,kr-1)= u(ir-1,jr-1,kr-1) + &
                (r000+r00m+rm0m+rm00+r0m0+r0mm+rmmm+rmm0)/8.0d0 
           
        enddo ! repeat over the whole coarse grid
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  return
end subroutine prolongation
