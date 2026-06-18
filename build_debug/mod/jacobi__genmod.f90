        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 16 22:31:33 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE JACOBI__genmod
          INTERFACE 
            SUBROUTINE JACOBI(U,B,H,BCND,RES,N,N3,EPS,KSOR,KMG,DIG,RANK,&
     &NPROC)
              INTEGER(KIND=4) :: NPROC
              INTEGER(KIND=4) :: N(3)
              REAL(KIND=8) :: U(0:N((1))+2,0:N((2))+2,-1:N((3))+2)
              REAL(KIND=8) :: B(0:N((1))+1,0:N((2))+1,0:N((3))+1)
              REAL(KIND=8) :: H(3)
              INTEGER(KIND=4) :: BCND(0:N((1))+2,0:N((2))+2,0:N((3))+2)
              REAL(KIND=8) :: RES
              INTEGER(KIND=4) :: N3
              REAL(KIND=8) :: EPS
              INTEGER(KIND=4) :: KSOR
              INTEGER(KIND=4) :: KMG
              INTEGER(KIND=4) :: DIG
              INTEGER(KIND=4) :: RANK
            END SUBROUTINE JACOBI
          END INTERFACE 
        END MODULE JACOBI__genmod
