        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 16 22:31:33 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETRES__genmod
          INTERFACE 
            SUBROUTINE GETRES(UORG,RHSORG,H,IR,JR,KR,N,R,FLAG)
              INTEGER(KIND=4) :: N(3)
              REAL(KIND=8) :: UORG(0:N((1))+2,0:N((2))+2,-1:N((3))+2)
              REAL(KIND=8) :: RHSORG(0:N((1))+1,0:N((2))+1,0:N((3))+1)
              REAL(KIND=8) :: H(3)
              INTEGER(KIND=4) :: IR
              INTEGER(KIND=4) :: JR
              INTEGER(KIND=4) :: KR
              REAL(KIND=8) :: R
              INTEGER(KIND=4) :: FLAG
            END SUBROUTINE GETRES
          END INTERFACE 
        END MODULE GETRES__genmod
