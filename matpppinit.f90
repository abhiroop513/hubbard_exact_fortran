! This program generates the Hamiltonian Matrix. Both diagonal and 
! off-diagonal parts are generated together. The input files are
! output of gendet which is detgen.out and magmol.inp. The output
! is written to matgen.out.
       use parameters
       use interfaces
      implicit real(8) (a-h,o-z)
      implicit integer(4) (i-n)
!     include "param2.cmn"
      complex(8) :: hmat,htemp,imag,phaset
      character(len=11),allocatable :: dstnms(:)
      allocatable :: idgseq(:),hmat(:),jhcol(:),irocom(:)
      dimension :: vbond(40)
      dimension :: iorb(20),ibond(2,40),jorb(20),nto2(30),jbond(2,40)
      dimension :: trnsfr(40),htemp(600),jhtemp(600),v(40,40),nzocc(20),uhub(20)
      dimension :: esite(20),nup(20),ndn(20),nocc(20),uhubi(20),coord(3,30)
      dimension :: inbond(2,40),trnsfr2(40),xp(20),yp(20),zp(20),vatt(20)
      complex(8) :: spn, spz
      real(8), allocatable :: amg(:)
!     common /srch/ idgseq,istates
!     common /all/ nto2
!***********************************************
!     interface
!      subroutine binsrch(id1,ifound,ibin,nto2,idet,islt)
!      implicit real(8)(a-h,o-z)
!      implicit integer(4)(i-n)
!      dimension :: idet(:),nto2(:)
!     end subroutine binsrch
!     end interface
!****************************************
!     write(*,*)'entered matppp'
      open(unit=1,file='initial.inp')
      open(unit=2,file='emfield.inp')
      open(unit=8,file='gen.out',form='unformatted')
!     open(unit=7,file='matpppq.out',form='unformatted')
      open(unit=7,file='matppp.out')
!     open(unit=7,file='matpppq.out')
      open(unit=9,file='dst.inp')
      open(unit=19,file='phase.inp')
  55  format(1x,a8)

      read(1,*) norb,nlctrns,amltplt,min,max,field,ipot
      read(1,*) accu,nbonds,nchrg,vrp
      read(1,*) (trnsfr(i),i=1,nbonds)
      read(1,*) (ibond(1,i),ibond(2,i),i=1,nbonds)
      read(1,*) (uhub(i),i=1,norb)
      read(1,*) (esite(i),i=1,norb)
      read(1,*) (nzocc(i),i=1,norb)
      read(1,*) (vatt(i),i=1,norb-1)

      allocate(amg(norb))
      read(1,*) (amg(i),i=1,norb)
      write(*,*) (amg(i),i=1,norb)
!     stop
!   For Next Nearest Neighbours
!     read(1,*)nnbonds
!     read(1,*) (inbond(1,i),inbond(2,i),i=1,nnbonds)
!     read(1,*) (trnsfr2(i),i=1,nnbonds)

!     read(2,*)idire1,idire2,idirm
      read(2,*)elfx1,elfy1,elfz1
!     read(2,*)elfx2,elfy2,elfz2
!     read(2,*)amgfx,amgfy,amgfz

!     read(12,*)isteps,deltat
!      read(17,*)field
!      time=isteps*deltat

       do i=1,norb
       read(9,*)xp(i),yp(i),zp(i)
!      write(*,*) xp(i),yp(i),zp(i)
       enddo

      if(nbonds.gt.40) write(*,*)'error - nbonds exceeds 40'
      if(norb.gt.20) write(*,*)'error - norb exceeds 20'
       imag=dcmplx(1.0d0,0.0d0)
      read(8) nsite,istates
      allocate(idgseq(istates))
      read(8) (idgseq(i),i=1,istates)
!     read(19,*)phaset
      phaset = (0.0,0.0)
      nto2(1)=1
      do i=2,30
      nto2(i)=nto2(i-1)*2
      if(nto2(i).le.istates) cycle
      ibin=i
      end do

!------------------dipolar interaction----------------------------
       v = 0.0
       do i =1,norb-1
       do j=i+1,norb
         if(i.ne.j) then
        r=(j-i)
        v1=vrp/r**3
         if(dabs(v1).le.1.0d-3)then
         v(i,j)=0.0d0
         else
         v(i,j)=v1
       endif
       endif
       enddo
       enddo

      do i=1,norb
      do j=1,norb
      v(j,i)=v(i,j)
      enddo
      enddo

!----------------------------------------------------------------------------------

!      if(ipot.eq.1) then 
!      do i=1,norb
!      read(9,*) xp(i),yp(i),zp(i) 
!      write(*,*) xp(i),yp(i),zp(i)
!      enddo
!      do i=1,norb
!         do j=1,norb
!            if(i.ne.j) then
!              rijsq=(xp(i)-xp(j))**2+(yp(i)-yp(j))**2+(zp(i)-zp(j))**2 
!              if((uhub(i)+uhub(j)).gt.1.0d-6) then
!              aux=((28.794d0/(uhub(i)+uhub(j)))**2+rijsq)
!              v(i,j)=14.397d0/dsqrt(aux)
!              else
!              endif
!            endif
!         enddo
!      enddo
!      endif
!------------------------------------------------
      allocate(irocom(kci))      
      ielem=1
!     write(*,*)'istates',istates
       allocate(hmat(maxel))
       allocate(jhcol(maxel))
      do i=1,istates   ! MAIN DO LOOP
         itot=1
         irocom(i)=ielem
         id=idgseq(i)
!        write(*,*)'id',id
         call decomp(id,norb,iorb)
         temp=0.0
         do k=1,norb
            nocc(k)=0
            nup(k)=0
            ndn(k)=0
            if(iorb(k).eq.3) then
               nocc(k)=2
               nup(k)=1
               ndn(k)=1
            endif
            if(iorb(k).eq.2) then
               nocc(k)=1
               nup(k)=1
            endif
            if(iorb(k).eq.1) then
               nocc(k)=1
               ndn(k)=1
            endif
         enddo
         do k=1,norb
            temp=temp+esite(k)*nocc(k)+uhub(k)*nup(k)*ndn(k)
            if(k.ne.norb) then
!             temp=temp+vatt(k)*(nocc(k)-nzocc(k))*(nocc(k+1)-nzocc(k+1))
             temp=temp+vatt(k)*(nocc(k))*(nocc(k+1))
            end if

!     1+nocc(k)*cos(0.0d0)
!**         kp1=k+1
!**         do j=kp1,norb
!**            temp=temp+v(k,j)*(nocc(k)-nzocc(k))*(nocc(j)-nzocc(j))
!**         enddo
         enddo
         htemp(itot)=temp*imag
!        write(*,*) 'i,htemp',i,htemp(itot)
         jhtemp(itot)=i
         itot=itot+1

!CCCC  electric field operator begins
!        write(*,*)'idirm,amg'
!        write(*,*)idirm,amg
!        do k=1,norb
!        htemp(itot)=htemp(itot)+(elfx1*xp(k)*nocc(k)+elfy1*yp(k)*nocc(k)+elfz1*zp(k)*nocc(k))*sin(fomegat)
!        enddo

!        write(*,113) itot,htemp(itot)

!        jhtemp(itot)=i
!        itot=itot+1

!      Electric field operator ends
!==========================================================================================
!     Magnetic field operator begins
      spz=(0.5,0.0)
      do ij=1, norb
         if((iorb(ij).eq.0).or.(iorb(ij).eq.3)) cycle

         if(iorb(ij).eq.1) spn=-spz
         if(iorb(ij).eq.2) spn=spz

         htemp(itot)=amg(ij)*spn
         jhtemp(itot)=i
         itot=itot+1
      end do

!     Magnetic field operator ends

         do 20 inb=1,nbonds
            ihop=0
 22         call trnsfrdn(iorb,inb,jorb,ires,norb,ibond,ihop,phs)
            if(ires.eq.0) go to 30
            call convrt(jorb,norb,idup)
            call binsrch(idup,ifound,ibin,nto2,idgseq,istates)
            jhtemp(itot)=ifound
            if(ihop.eq.0)htemp(itot)=trnsfr(inb)*phs*exp(phaset)
            if(ihop.eq.1)htemp(itot)=trnsfr(inb)*phs*exp(-phaset)

!            htemp(itot)=trnsfr(inb)*phs*imag
            itot=itot+1
  30        call trnsfrup(iorb,inb,jorb,ires,norb,ibond,ihop,phs)
            if(ires.eq.0) go to 21
            call convrt(jorb,norb,iddn)
            call binsrch(iddn,ifound,ibin,nto2,idgseq,istates)
            jhtemp(itot)=ifound
            if(ihop.eq.0)htemp(itot)=trnsfr(inb)*phs*exp(phaset)
            if(ihop.eq.1)htemp(itot)=trnsfr(inb)*phs*exp(-phaset)
!            htemp(itot)=trnsfr(inb)*phs*imag
            itot=itot+1
  21        ihop=ihop+1
            if(ihop.gt.1) go to 20
            go to 22
  20     continue
          go to 79
!cccc   NEXT-NEAREST NEIGHBOUR
!@          do 120 inb=1,nnbonds
!@          ihop=0
!@122       call trnsfrdn(iorb,inb,jorb,ires,norb,inbond,ihop,phs)
!@          if(ires.eq.0) go to 130
!@          call convrt(jorb,norb,idup)
!@          call binsrch(idup,ifound,ibin,nto2,idgseq,istates)
!@          jhtemp(itot)=ifound
!@          if(ihop.eq.0)htemp(itot)=trnsfr2(inb)*phs*exp(phaset)
!@          if(ihop.eq.1)htemp(itot)=trnsfr2(inb)*phs*exp(-phaset)

!c            htemp(itot)=trnsfr(inb)*phs*imag
!@          itot=itot+1
!@130        call trnsfrup(iorb,inb,jorb,ires,norb,inbond,ihop,phs)
!@          if(ires.eq.0) go to 121
!@          call convrt(jorb,norb,iddn)
!@          call binsrch(iddn,ifound,ibin,nto2,idgseq,istates)
!@          jhtemp(itot)=ifound
!@          if(ihop.eq.0)htemp(itot)=trnsfr2(inb)*phs*exp(phaset)
!@          if(ihop.eq.1)htemp(itot)=trnsfr2(inb)*phs*exp(-phaset)
!c            htemp(itot)=trnsfr(inb)*phs*imag
!@          itot=itot+1
!@121        ihop=ihop+1
!@          if(ihop.gt.1) go to 120
!@          go to 122
!@ 120     continue
!CC    NNN Hopping ends
  79     itot=itot-1
         if(itot.gt.500) write(*,*) 'error itot exceeds 500'
         call order(htemp,jhtemp,hmat,jhcol,itot,ielem)
      end do
      ielem=ielem-1
      write(*,*)'ielem,istates',ielem,istates
!     write(7) istates,ielem
!     write(7) (irocom(k),k=1,istates)
!     write(7) (hmat(k),k=1,ielem)
!     write(7) (jhcol(k),k=1,ielem)
      write(7,*) istates,ielem
      write(7,*) (irocom(k),k=1,istates)
      write(7,*) (hmat(k),k=1,ielem)
      write(7,*) (jhcol(k),k=1,ielem)
!     write(*,991) (hmat(k),k=1,ielem)
!     write(*,992) (irocom(k),k=1,istates)
!     write(*,992) (jhcol(k),k=1,ielem)
 991  format(1x,5e16.8)
 992  format(1x,10i7)
      stop
      end
      subroutine decomp(id,norb,iorb)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension iorb(:)
      i2=2
      do i=1,norb
      idp=ishft(id,-i2)
      idp=ishft(idp,i2)
      iorb(norb-i+1)=id-idp
      id=ishft(idp,-i2)
      end do
!     return
      end
      subroutine trnsfrdn(iorb,inb,jorb,ires,norb,ibond,ihop,phs)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension :: iorb(:),ibond(:,:),jorb(:),korb(20)
!***************************************************
      interface
      subroutine annphase(j,jorb,jspin,jphase)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension jorb(:)
      end subroutine annphase
      end interface
!************************************************
      interface
      subroutine crphase(j,jorb,jspin,jphase)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension jorb(:)
      end subroutine crphase
      end interface
!************************************************
      phs=1.0
      idon=ibond(1,inb)
      iacc=ibond(2,inb)
      if(ihop.eq.1) then
      idon=ibond(2,inb)
      iacc=ibond(1,inb)
      endif
      ires=0
      do i=1,norb
      jorb(i)=iorb(i)
      korb(i)=iorb(i)
      end do
      if(jorb(iacc).eq.3.or.jorb(iacc).eq.1) go to 20
      if(jorb(idon).eq.0.or.jorb(idon).eq.2) go to 20
      if(jorb(iacc).ne.0) go to 21
      jorb(iacc)=1
      if(jorb(idon).eq.3) jorb(idon)=2
      if(jorb(idon).eq.1) jorb(idon)=0
      go to 25
  21  jorb(iacc)=3
      if(jorb(idon).eq.3) jorb(idon)=2
      if(jorb(idon).eq.1) jorb(idon)=0
  25  ires=1
      jphase=2
      jspin=-1
      call annphase(idon,iorb,jspin,jphase)
!     write(*,*)'jphase after annphase',jphase
      korb(idon)=iorb(idon)-1
      call crphase(iacc,korb,jspin,jphase)
!     write(*,*)'jphase after crphase',jphase
      if((jphase/2)*2.ne.jphase) phs=-1.0
  20  end        
      subroutine trnsfrup(iorb,inb,jorb,ires,norb,ibond,ihop,phs)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension :: iorb(:),ibond(:,:),jorb(:),korb(20)
!***************************************************
      interface
      subroutine annphase(j,jorb,jspin,jphase)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension jorb(:)
      end subroutine annphase
      end interface
!************************************************
      interface
      subroutine crphase(j,jorb,jspin,jphase)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension jorb(:)
      end subroutine crphase
      end interface
!************************************************
      phs=1.0
      idon=ibond(1,inb)
      iacc=ibond(2,inb)
      if(ihop.eq.1) then
      idon=ibond(2,inb)
      iacc=ibond(1,inb)
      endif
      ires=0
      do i=1,norb
      jorb(i)=iorb(i)
      korb(i)=iorb(i)
      end do
      if(jorb(iacc).eq.3.or.jorb(iacc).eq.2) go to 20
      if(jorb(idon).eq.0.or.jorb(idon).eq.1) go to 20
      if(jorb(iacc).ne.0) go to 21
      jorb(iacc)=2
      if(jorb(idon).eq.3) jorb(idon)=1
      if(jorb(idon).eq.2) jorb(idon)=0
      go to 25
  21  jorb(iacc)=3
      if(jorb(idon).eq.3) jorb(idon)=1
      if(jorb(idon).eq.2) jorb(idon)=0
  25  ires=1
      jphase=2
      jspin=1
      call annphase(idon,iorb,jspin,jphase)
      korb(idon)=iorb(idon)-2
      call crphase(iacc,korb,jspin,jphase)
!     write(*,*)'jphase',jphase
      if((jphase/2)*2.ne.jphase) phs=-1.0
  20  end        
      subroutine convrt(jorb,norb,id)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension jorb(:)
      id=0
      do 10 i=1,norb
      if(jorb(i).eq.0) go to 10
      j=norb-i+1
      js=jorb(i)
      go to (2,3,4),js
   2  jl=2*j-2
      id=ibset(id,jl)
      go to 10
   3  jh=2*j-1
      id=ibset(id,jh)
      go to 10
   4  jl=2*j-2
      jh=jl+1
      id=ibset(id,jl)
      id=ibset(id,jh)
  10  continue
      end
      subroutine order(cotemp,itdg,conew,itnew,itot,ielem)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      complex(8) :: atrns,cotemp,conew
!     include 'param.cmn'
      dimension :: itdg(:),cotemp(:),conew(:),itnew(:)
      do 10 i=1,itot
      if(i.eq.itot) go to 10
      ip1=i+1
      do 15 j=ip1,itot
      if(itdg(j).eq.0) go to 15
      if(itdg(i).lt.itdg(j)) go to 15
      if(itdg(i).eq.itdg(j)) go to 20
      itrns=itdg(i)
      itdg(i)=itdg(j)
      itdg(j)=itrns
      atrns=cotemp(i)
      cotemp(i)=cotemp(j)
      cotemp(j)=atrns
      go to 15
  20  itdg(j)=0
      cotemp(i)=cotemp(i)+cotemp(j)
      cotemp(j)=0.0
  15  continue
  10  continue
      j=0
      do 30 i=1,itot
      if(itdg(i).eq.0) go to 30
!      if(dabs(cotemp(i)).lt.1.0e-8) go to 30
         if((abs(dreal(cotemp(i))).lt.1.0d-5).and.(abs(dimag(cotemp(i))).lt.1.0d-5)) go to 30
      itnew(ielem+j)=itdg(i)
      conew(ielem+j)=cotemp(i)
      j=j+1
  30  continue
      ielem=ielem+j
      end
      subroutine binsrch(id1,ifound,ibin,nto2,idet,islt)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension :: idet(:),nto2(:)
!     common /srch/idet,islt
      ibeg=ibin-1
      itst=nto2(ibeg)
  35  if(itst.gt.islt) go to 50
      if(id1.eq.idet(itst)) go to 25
      if(id1.gt.idet(itst)) go to 15
  50  ibeg=ibeg-1
      if(ibeg.eq.0) go to 60
      itst=itst-nto2(ibeg)
!      write(*,*) 'itst,ibeg',itst,ibeg
      go to 35
  15  ibeg=ibeg-1
      if(ibeg.eq.0) go to 60
      itst=itst+nto2(ibeg)
!      write(*,*) 'itst,ibeg',itst,ibeg
      go to 35
  60  write(*,*)'error final state not found in 2nd diagram search'
      go to 40
  25  ifound=itst
!      write(*,*) ifound
  40  end
      subroutine annphase(j,jorb,jspin,jphase)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension jorb(:)
      jm1=j-1
      if(j.eq.1) go to 35
      do 30 iphs=1,jm1
      if(jorb(iphs).eq.3.or.jorb(iphs).eq.0) go to 30
      jphase=jphase+1
  30  continue
  35  if(jspin.eq.-1.and.jorb(j).eq.3) jphase=jphase+1
      return
      end
      subroutine crphase(j,jorb,jspin,jphase)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension jorb(:)
!     write(*,*)'j in crphs',j
      jm1=j-1
      if(j.eq.1) go to 35
      do 30 iphs=1,jm1
      if(jorb(iphs).eq.3.or.jorb(iphs).eq.0) go to 30
      jphase=jphase+1
  30  continue
  35  if(jspin.eq.-1.and.jorb(j).eq.2) jphase=jphase+1
      end
