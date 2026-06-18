        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 16 22:31:33 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RESTRICTION__genmod
          INTERFACE 
            SUBROUTINE RESTRICTION(U,B,H,BCND,BCND2H,R2H,E2H,N,N3,RANK, &
     &NPROC)
              INTEGER(KIND=4) :: N(3)
              REAL(KIND=8) :: U(0:N((1))+2,0:N((2))+2,-1:N((3))+2)
              REAL(KIND=8) :: B(0:N((1))+1,0:N((2))+1,0:N((3))+1)
              REAL(KIND=8) :: H(3)
              INTEGER(KIND=4) :: BCND(0:N((1))+2,0:N((2))+2,0:N((3))+2)
              INTEGER(KIND=4) :: BCND2H(0:N((1))/2+2,0:N((2))/2+2,0:N((3&
     &))/2+2)
              REAL(KIND=8) :: R2H(0:N((1))/2+1,0:N((2))/2+1,0:N((3))/2+1&
     &)
              REAL(KIND=8) :: E2H(0:N((1))/2+2,0:N((2))/2+2,-1:N((3))/2+&
     &2)
              INTEGER(KIND=4) :: N3
              INTEGER(KIND=4) :: RANK
              INTEGER(KIND=4) :: NPROC
            END SUBROUTINE RESTRICTION
          END INTERFACE 
        END MODULE RESTRICTION__genmod
