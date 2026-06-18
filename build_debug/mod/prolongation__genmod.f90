        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 16 22:31:33 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROLONGATION__genmod
          INTERFACE 
            SUBROUTINE PROLONGATION(U,E2H,BCND,N,N3,NPROC)
              INTEGER(KIND=4) :: N(3)
              REAL(KIND=8) :: U(0:N((1))+2,0:N((2))+2,-1:N((3))+2)
              REAL(KIND=8) :: E2H(0:N((1))/2+2,0:N((2))/2+2,-1:N((3))/2+&
     &2)
              INTEGER(KIND=4) :: BCND(0:N((1))+2,0:N((2))+2,0:N((3))+2)
              INTEGER(KIND=4) :: N3
              INTEGER(KIND=4) :: NPROC
            END SUBROUTINE PROLONGATION
          END INTERFACE 
        END MODULE PROLONGATION__genmod
