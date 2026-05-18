function ran2(irand)
  implicit none
  integer,parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
  real(kind=8),parameter :: am=1.0d0/im
  real(kind=8):: ran2
  integer :: k,irand

  k=irand/iq
  irand=ia*(irand-k*iq)-ir*k
  if (irand.lt.0) irand=irand+im
  ran2=am*irand
  
  return

end function ran2