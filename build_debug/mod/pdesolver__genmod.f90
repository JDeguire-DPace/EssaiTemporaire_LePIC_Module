        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 16 22:31:33 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PDESOLVER__genmod
          INTERFACE 
            SUBROUTINE PDESOLVER(U,B,BCND,H,N,NCYCL,EPS,OMEGA,K,KTOT,RES&
     &,NG,RANK,NPROC)
              INTEGER(KIND=4) :: NPROC
              INTEGER(KIND=4) :: N(3)
              REAL(KIND=8) :: U(0:N((1))+2,0:N((2))+2,-1:N((3))/NPROC+2)
              REAL(KIND=8) :: B(0:N((1))+1,0:N((2))+1,0:N((3))/NPROC+1)
              INTEGER(KIND=4) :: BCND(0:N((1))+2,0:N((2))+2,0:N((3))/   &
     &NPROC+2)
              REAL(KIND=8) :: H(3)
              INTEGER(KIND=4) :: NCYCL
              REAL(KIND=8) :: EPS
              REAL(KIND=8) :: OMEGA
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: KTOT
              REAL(KIND=8) :: RES
              INTEGER(KIND=4) :: NG
              INTEGER(KIND=4) :: RANK
            END SUBROUTINE PDESOLVER
          END INTERFACE 
        END MODULE PDESOLVER__genmod
