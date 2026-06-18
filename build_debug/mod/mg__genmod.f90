        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 16 22:31:33 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MG__genmod
          INTERFACE 
            SUBROUTINE MG(U,B,BCND,H,RES,N,OMEGA,NG,EPS,K,KTOT,RANK,    &
     &NPROC,E2,E4,E8,E16,E32,E64,E128,E256,E512,E1024,E2048,R2,R4,R8,R16&
     &,R32,R64,R128,R256,R512,R1024,R2048,BCND2,BCND4,BCND8,BCND16,     &
     &BCND32,BCND64,BCND128,BCND256,BCND512,BCND1024,BCND2048)
              INTEGER(KIND=4) :: NG
              INTEGER(KIND=4) :: N(3,12)
              REAL(KIND=8) :: U(0:N((1,1))+2,0:N((2,1))+2,-1:N((3,1))+2)
              REAL(KIND=8) :: B(0:N((1,1))+1,0:N((2,1))+1,0:N((3,1))+1)
              INTEGER(KIND=4) :: BCND(0:N((1,1))+2,0:N((2,1))+2,0:N((3,1&
     &))+2)
              REAL(KIND=8) :: H(3,12)
              REAL(KIND=8) :: RES
              REAL(KIND=8) :: OMEGA
              REAL(KIND=8) :: EPS
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: KTOT
              INTEGER(KIND=4) :: RANK
              INTEGER(KIND=4) :: NPROC
              REAL(KIND=8) :: E2(0:N((1,2))+2,0:N((2,2))+2,-1:N((3,2))+2&
     &)
              REAL(KIND=8) :: E4(0:N((1,3))+2,0:N((2,3))+2,-1:N((3,3))+2&
     &)
              REAL(KIND=8) :: E8(0:N((1,4))+2,0:N((2,4))+2,-1:N((3,4))+2&
     &)
              REAL(KIND=8) :: E16(0:N((1,5))+2,0:N((2,5))+2,-1:N((3,5))+&
     &2)
              REAL(KIND=8) :: E32(0:N((1,6))+2,0:N((2,6))+2,-1:N((3,6))+&
     &2)
              REAL(KIND=8) :: E64(0:N((1,7))+2,0:N((2,7))+2,-1:N((3,7))+&
     &2)
              REAL(KIND=8) :: E128(0:N((1,8))+2,0:N((2,8))+2,-1:N((3,8))&
     &+2)
              REAL(KIND=8) :: E256(0:N((1,9))+2,0:N((2,9))+2,-1:N((3,9))&
     &+2)
              REAL(KIND=8) :: E512(0:N((1,10))+2,0:N((2,10))+2,-1:N((3, &
     &10))+2)
              REAL(KIND=8) :: E1024(0:N((1,11))+2,0:N((2,11))+2,-1:N((3,&
     &11))+2)
              REAL(KIND=8) :: E2048(0:N((1,12))+2,0:N((2,12))+2,-1:N((3,&
     &12))+2)
              REAL(KIND=8) :: R2(0:N((1,2))+1,0:N((2,2))+1,0:N((3,2))+1)
              REAL(KIND=8) :: R4(0:N((1,3))+1,0:N((2,3))+1,0:N((3,3))+1)
              REAL(KIND=8) :: R8(0:N((1,4))+1,0:N((2,4))+1,0:N((3,4))+1)
              REAL(KIND=8) :: R16(0:N((1,5))+1,0:N((2,5))+1,0:N((3,5))+1&
     &)
              REAL(KIND=8) :: R32(0:N((1,6))+1,0:N((2,6))+1,0:N((3,6))+1&
     &)
              REAL(KIND=8) :: R64(0:N((1,7))+1,0:N((2,7))+1,0:N((3,7))+1&
     &)
              REAL(KIND=8) :: R128(0:N((1,8))+1,0:N((2,8))+1,0:N((3,8))+&
     &1)
              REAL(KIND=8) :: R256(0:N((1,9))+1,0:N((2,9))+1,0:N((3,9))+&
     &1)
              REAL(KIND=8) :: R512(0:N((1,10))+1,0:N((2,10))+1,0:N((3,10&
     &))+1)
              REAL(KIND=8) :: R1024(0:N((1,11))+1,0:N((2,11))+1,0:N((3, &
     &11))+1)
              REAL(KIND=8) :: R2048(0:N((1,12))+1,0:N((2,12))+1,0:N((3, &
     &12))+1)
              INTEGER(KIND=4) :: BCND2(0:N((1,2))+2,0:N((2,2))+2,0:N((3,&
     &2))+2)
              INTEGER(KIND=4) :: BCND4(0:N((1,3))+2,0:N((2,3))+2,0:N((3,&
     &3))+2)
              INTEGER(KIND=4) :: BCND8(0:N((1,4))+2,0:N((2,4))+2,0:N((3,&
     &4))+2)
              INTEGER(KIND=4) :: BCND16(0:N((1,5))+2,0:N((2,5))+2,0:N((3&
     &,5))+2)
              INTEGER(KIND=4) :: BCND32(0:N((1,6))+2,0:N((2,6))+2,0:N((3&
     &,6))+2)
              INTEGER(KIND=4) :: BCND64(0:N((1,7))+2,0:N((2,7))+2,0:N((3&
     &,7))+2)
              INTEGER(KIND=4) :: BCND128(0:N((1,8))+2,0:N((2,8))+2,0:N((&
     &3,8))+2)
              INTEGER(KIND=4) :: BCND256(0:N((1,9))+2,0:N((2,9))+2,0:N((&
     &3,9))+2)
              INTEGER(KIND=4) :: BCND512(0:N((1,10))+2,0:N((2,10))+2,0:N&
     &((3,10))+2)
              INTEGER(KIND=4) :: BCND1024(0:N((1,11))+2,0:N((2,11))+2,0:&
     &N((3,11))+2)
              INTEGER(KIND=4) :: BCND2048(0:N((1,12))+2,0:N((2,12))+2,0:&
     &N((3,12))+2)
            END SUBROUTINE MG
          END INTERFACE 
        END MODULE MG__genmod
