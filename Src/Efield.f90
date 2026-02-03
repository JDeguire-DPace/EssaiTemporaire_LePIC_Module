subroutine calc_Efield(n,h,phi,Ei,bcnd)
!     ===================================================================
!     VERSION:         0.2
!     LAST MOD:      Sep/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: Scalars are evaluated at grid points, vectors in between
!     NOTE:     What about BCs and corner points?
!     -------------------------------------------------------------------
  implicit none
  integer:: ix,iy,iz,n(3),bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
  real(kind=8):: h(3),phi(0:n(1)+2,0:n(2)+2,0:n(3)+2), &
       Ei(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
  include 'particle_info.h'

  !
  ! Interior points
  !
  !$OMP PARALLEL
  !$OMP DO
  do iz=1,n(3)+1
     do iy=1,n(2)+1
        do ix=1,n(1)+1

           ! Interior points
           if( bcnd(ix,iy,iz).le.0 ) then

              ! Second order correct in dx
              Ei(1,ix,iy,iz)= -( phi(ix+1,iy,iz)-phi(ix-1,iy,iz) )/(2.d0*h(1))

              ! Second order correct in dy
              Ei(2,ix,iy,iz)= -( phi(ix,iy+1,iz)-phi(ix,iy-1,iz) )/(2.d0*h(2))

              ! Second order correct in dz
              Ei(3,ix,iy,iz)= -( phi(ix,iy,iz+1)-phi(ix,iy,iz-1) )/(2.d0*h(3))

           endif

        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  !
  ! Boundary conditions
  !
  !$OMP PARALLEL
  !$OMP DO
  do iz=1,n(3)+1
     do iy=1,n(2)+1
        do ix=1,n(1)+1

           ! Skip everything 
           if( bcnd(ix,iy,iz).eq.-1 ) goto 10

           !
           ! Neumann BCs
           !
           if( bcnd(ix,iy,iz).eq.-2 ) then ! YZ plane, LHS only
           
              ! Calculate Ex
              Ei(1,1,iy,iz)= 0.d0

              if( iy.gt.1 .and. iy.lt.n(2)+1 .and. &
                   iz.gt.1 .and. iz.lt.n(3)+1 ) then
                 ! Calculate Ey
                 Ei(2,1,iy,iz)= -( phi(1,iy+1,iz)-phi(1,iy-1,iz) )/(2.d0*h(2))
                 ! Calculate Ez
                 Ei(3,1,iy,iz)= -( phi(1,iy,iz+1)-phi(1,iy,iz-1) )/(2.d0*h(3))
                 goto 10
              endif
              
           endif


           !
           ! Periodic BCs
           !
           if( bcnd(ix,iy,iz).eq.0 .or. bcnd(ix,iy,iz).eq.-2 ) then

              ! iy= 0,1,ny+1 and ny+2 simulation planes
              if( iy.eq.1 ) then
                 ! Calculate Ey
                 Ei(2,ix,1,iz)= -( phi(ix,2,iz)-phi(ix,n(2),iz) )/(2.d0*h(2))
                 Ei(2,ix,0,iz)= Ei(2,ix,n(2),iz)
                 Ei(2,ix,n(2)+1,iz)= Ei(2,ix,1,iz)
                 Ei(2,ix,n(2)+2,iz)= Ei(2,ix,2,iz)
                 ! Calculate Ex
                 Ei(1,ix,0,iz)= Ei(1,ix,n(2),iz)
                 Ei(1,ix,n(2)+2,iz)= Ei(1,ix,2,iz)
                 ! Calculate Ez
                 Ei(3,ix,0,iz)= Ei(3,ix,n(2),iz)
                 Ei(3,ix,n(2)+2,iz)= Ei(3,ix,2,iz)
              endif

              ! iz= 0,1,nz+1 and nz+2 planes
              if( iz.eq.1 ) then
                 ! Calculate Ez
                 Ei(3,ix,iy,1)= -( phi(ix,iy,2)-phi(ix,iy,n(3)) )/(2.d0*h(3))
                 Ei(3,ix,iy,0)= Ei(3,ix,iy,n(3))
                 Ei(3,ix,iy,n(3)+1)= Ei(3,ix,iy,1)
                 Ei(3,ix,iy,n(3)+2)= Ei(3,ix,iy,2)
                 ! Calculate Ex
                 Ei(1,ix,iy,0)= Ei(1,ix,iy,n(3))
                 Ei(1,ix,iy,n(3)+2)= Ei(1,ix,iy,2)
                 ! Calculate Ey
                 Ei(2,ix,iy,0)= Ei(2,ix,iy,n(3))
                 Ei(2,ix,iy,n(3)+2)= Ei(2,ix,iy,2)
              endif
              
              goto 10
           endif
           
           !
           ! Walls
           !
           if( bcnd(ix,iy,iz).ge.1 ) then
              
              ! Wall on the west side
              if( bcnd(ix+1,iy,iz).le.0 ) then          
                 ! Ex
                 Ei(1,ix,iy,iz)= 2.d0*Ei(1,ix+1,iy,iz) - Ei(1,ix+2,iy,iz)
                 ! Ey
                 Ei(2,ix,iy,iz)= -( phi(ix,iy+1,iz)-phi(ix,iy-1,iz) )/(2.d0*h(2))
                 ! Ez
                 Ei(3,ix,iy,iz)= -( phi(ix,iy,iz+1)-phi(ix,iy,iz-1) )/(2.d0*h(3))
              endif

              ! Wall on the east side
              if( bcnd(ix-1,iy,iz).le.0 ) then          
                 ! Ex
                 Ei(1,ix,iy,iz)= 2.d0*Ei(1,ix-1,iy,iz) - Ei(1,ix-2,iy,iz)
                 ! Ey
                 Ei(2,ix,iy,iz)= -( phi(ix,iy+1,iz)-phi(ix,iy-1,iz) )/(2.d0*h(2))
                 ! Ez
                 Ei(3,ix,iy,iz)= -( phi(ix,iy,iz+1)-phi(ix,iy,iz-1) )/(2.d0*h(3))
              endif

              ! Wall at the south side
              if( bcnd(ix,iy+1,iz).le.0 ) then
                 ! Ex
                 Ei(1,ix,iy,iz)= -( phi(ix+1,iy,iz)-phi(ix-1,iy,iz) )/(2.d0*h(1))
                 ! Ey
                 Ei(2,ix,iy,iz)= 2.d0*Ei(2,ix,iy+1,iz) - Ei(2,ix,iy+2,iz)
                 ! Ez
                 Ei(3,ix,iy,iz)= -( phi(ix,iy,iz+1)-phi(ix,iy,iz-1) )/(2.d0*h(3))
              endif

              ! Wall at the north side
              if( bcnd(ix,iy-1,iz).le.0 ) then
                 ! Ex
                 Ei(1,ix,iy,iz)= -( phi(ix+1,iy,iz)-phi(ix-1,iy,iz) )/(2.d0*h(1))
                 ! Ey
                 Ei(2,ix,iy,iz)= 2.d0*Ei(2,ix,iy-1,iz) - Ei(2,ix,iy-2,iz)
                 ! Ez
                 Ei(3,ix,iy,iz)= -( phi(ix,iy,iz+1)-phi(ix,iy,iz-1) )/(2.d0*h(3))
              endif

              ! Wall at the bottom
              if( bcnd(ix,iy,iz+1).le.0 ) then
                 ! Ex
                 Ei(1,ix,iy,iz)= -( phi(ix+1,iy,iz)-phi(ix-1,iy,iz) )/(2.d0*h(1))
                 ! Ey
                 Ei(2,ix,iy,iz)= -( phi(ix,iy+1,iz)-phi(ix,iy-1,iz) )/(2.d0*h(2))                 
                 ! Ez
                 Ei(3,ix,iy,iz)= 2.d0*Ei(3,ix,iy,iz+1) - Ei(3,ix,iy,iz+2)   
              endif

              ! Wall at the top
              if( bcnd(ix,iy,iz-1).le.0 ) then
                 ! Ex
                 Ei(1,ix,iy,iz)= -( phi(ix+1,iy,iz)-phi(ix-1,iy,iz) )/(2.d0*h(1))
                 ! Ey
                 Ei(2,ix,iy,iz)= -( phi(ix,iy+1,iz)-phi(ix,iy-1,iz) )/(2.d0*h(2))   
                 ! Ez
                 Ei(3,ix,iy,iz)= 2.d0*Ei(3,ix,iy,iz-1) - Ei(3,ix,iy,iz-2) 
              endif
              
           endif

10         continue       

        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  return
  
end subroutine calc_Efield
