       module parameters
!       implicit none
!       integer :: mxci,lci,msqrd16,kci,irwnd,mxvcr,nfl       
!       integer :: noo,non,maxel,mxsml,irmx,maxes,iblck
!       integer :: irwndp,irwndm,mxvcrp,mxvcrm
      integer, parameter ::mxci=280,lci=4*mxci,msqrd16=lci*lci
      integer, parameter ::kci=4300000000,irwnd=40,mxvcr=kci*irwnd,nfl=400
      integer, parameter ::noo=20,non=20,maxel=6225000000,mxsml=500,irmx=800
      integer, parameter ::maxes=35000000,iblck=5000
      integer, parameter ::irwndp=20,irwndm=20,mxvcrp=kci*irwndp
      integer, parameter ::mxvcrm=kci*irwndm
       end module parameters
