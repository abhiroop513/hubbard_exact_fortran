! This program generates the integers corresponding to the Slater
! determinants for sites=nsite, electrons=nlctrn and amltplt is the
! desired Ms value (not 2S+1 but any S 0, 0.5 etc).
! output is written to detgen.out and the input file is magmol.inp
! itrtngen.inp is a file to keep track of the number of sites later.
      use parameters
      implicit real(8)  (a-h,o-z)
      implicit integer(4)  (i-n)
      dimension :: ksite(30),lbit(30)
      allocatable :: idiag(:),n_op(:,:),spin(:,:),n_up(:,:),n_dn(:,:)
      complex(8), allocatable:: evec(:)
      allocate(idiag(kci))
      allocate(spin(kci,20))
      allocate(n_op(kci,30))
      allocate(n_up(kci,30))
      allocate(n_dn(kci,30))
      write(*,*)'entered gendet'
      open(unit=1,file='initial.inp')
      open(unit=7,file='gen.out',form='unformatted')
!     open(unit=17,file='gen2.out')
      open(unit=18,file='sztot.out',form='unformatted')
      open(unit=28,file='chtot.out',form='unformatted')
      open(unit=29,file='nup.out',form='unformatted')
      open(unit=30,file='ndown.out',form='unformatted')
!     open(unit=18,file='sztot.out')
!     open(unit=2,file='itrtngen.inp')
      read(1,*) nsite,nlctrn,amltplt
      nmltplt=4
      nbitsps=2
      write(2,*) nsite,1
! nsite=# of sites, nmltplt=multiplicity of site spin and amltplt
! is the multiplicity of the ground state of the chain.
!     write(*,*)'nsite,nlctrn,amltplt'
!     write(*,*) nsite,nlctrn,amltplt
      izero=0
      imin=0
      imax=2**(nbitsps*nsite)
      nstspn=nmltplt-1
      lbit(1)=(nsite-1)*nbitsps
      do i=2,nsite
       lbit(i)=lbit(i-1)-nbitsps
      end do
!     write(*,*)'imin,imax',imin,imax
!     write(*,*)'lbit(i),i=1,nsite'
!     write(*,*)(lbit(i),i=1,nsite)
      istates=0
      do 10 i=imin,imax
      do 16 k=1,nsite
      ksite(k)=0
  16  continue
      check=0
      lctrn=0
      do 15 j=1,nsite
      temp=0.0
      ibeg=lbit(j)
      call mvbits(i,ibeg,nbitsps,ksite(j),izero)
      if(ksite(j).gt.nstspn) go to 10
      lctrn=lctrn+1
      if(ksite(j).eq.0) lctrn=lctrn-1
      if(ksite(j).eq.3) lctrn=lctrn+1
      if(ksite(j).eq.1) then
       temp = -0.5
       check=check-0.5
      end if
      if(ksite(j).eq.2) then 
       temp = 0.5
       check=check+0.5
      end if
      spin(istates+1,j) = temp

      if(ksite(j).eq.0) then
        n_up(istates+1,j) = 0
        n_dn(istates+1,j) = 0
        n_op(istates+1,j) = 0
      end if  
      if(ksite(j).eq.1) then
        n_up(istates+1,j) = 0
        n_dn(istates+1,j) = 1
        n_op(istates+1,j) = 1
      end if
      if(ksite(j).eq.2) then
        n_up(istates+1,j) = 1
        n_dn(istates+1,j) = 0
        n_op(istates+1,j) = 1
      end if
      if(ksite(j).eq.3) then
        n_up(istates+1,j) = 1
        n_dn(istates+1,j) = 1
        n_op(istates+1,j) = 2
      end if

  15  continue
      if(dabs(check-amltplt).gt.1.0d-8) go to 10
      if(lctrn.ne.nlctrn) go to 10
      istates=istates+1
      idiag(istates)=i
      
!     write(*,*)'i,istates,ksite(i),i=1,nsite'
!     write(*,33)i,istates,(ksite(k),k=1,nsite)
  10  continue
      write(*,*)'idiag',istates
      if(istates.gt.kci) write(*,*)'error istates gt kci'
      write(7)nsite,istates
      write(7)(idiag(i),i=1,istates)
      write(18)((spin(i,j),j=1,nsite),i=1,istates)
      write(28)((n_op(i,j),j=1,nsite),i=1,istates)
      write(29)((n_up(i,j),j=1,nsite),i=1,istates)
      write(30)((n_dn(i,j),j=1,nsite),i=1,istates)


!=========================================================
!            Create a Doublon State
!=========================================================
      allocate(evec(istates))
      evec = (0.0,0.0)
      nsiteby2 = nsite/2
      write(*,*)nsiteby2
      do i = 1, istates
       icount = 0
       do j = 1, nsite, 2
         if((n_op(i,j).eq.2).and.(n_op(i,j+1).eq.0))  icount = icount+1
       end do
        if(icount.eq.nsiteby2) then
          evec(i) = (1.0, 0.0)
           write(*,*)(n_op(i,j), j = 1, nsite)
        end if
      end do
      open(3, file='unsymdet.out')
      write(3,*)ajunk,istates
      write(3,*)(evec(i),i=1,istates)
!      write(*,*)(ksite(j), j = 1, nsite)

!     do i = 1,istates
!       write(17,*)idiag(i)
!       write(40,*)(n_op(i,j),j=1,nsite)
!     end do
!     write(*,*)nsite,nlctrn,amltplt,istates
!      write(*,33)(idiag(i),i=1,istates)
  33  format(1x,10i7)
      end
