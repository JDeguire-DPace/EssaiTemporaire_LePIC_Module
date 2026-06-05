module mod_utils
    use iso_fortran_env, only: real64, int32
    implicit none
    private
    public :: stop_calculation, indexx, ran2
contains
    subroutine stop_calculation
        use mpi
        implicit none
        integer ierr

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        STOP
        call MPI_Finalize(ierr)

        return
    end subroutine stop_calculation

    function ran2(irand) result(r)
        use iso_fortran_env, only: int32, real64
        implicit none
        integer(int32), intent(inout) :: irand
        integer(int32), parameter :: ia=16807_int32, im=2147483647_int32, iq=127773_int32, ir=2836_int32
        real(real64),   parameter :: am=1.0_real64/real(im,real64)
        real(real64) :: r
        integer(int32) :: k

        k = irand / iq
        irand = ia*(irand - k*iq) - ir*k
        if (irand < 0_int32) irand = irand + im
        r = am * real(irand, real64)
    end function ran2



    SUBROUTINE indexx(n,arr,indx)
    !     ==============================================================
    !     VERSION:          /
    !     LAST MOD:      Years 1986-92
    !     MOD AUTHOR:    Numerical Recipes Software
    !     COMMENTS:         /
    !     --------------------------------------------------------------
    INTEGER:: n,indx(n),M,NSTACK
    REAL(KIND=8):: arr(n)
    PARAMETER (M=7,NSTACK=50)
    INTEGER:: i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
    REAL(KIND=8) a

    do  j=1,n
        indx(j)=j
    enddo
    jstack=0
    l=1
    ir=n
    1 if(ir-l.lt.M)then
        do j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
            enddo
            i=0
    2       indx(i+1)=indxt
        enddo
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
    else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
            itemp=indx(l+1)
            indx(l+1)=indx(l)
            indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
    3    continue
        i=i+1
        if(arr(indx(i)).lt.a)goto 3
    4    continue
        j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
    5    indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then 
            print*, 'NSTACK too small in indexx'
            STOP
        endif
        if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
        else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
        endif
    endif
    goto 1
    END SUBROUTINE indexx

end module mod_utils