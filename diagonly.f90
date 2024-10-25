      implicit none
      real(8),allocatable ::hlanc_d(:),hlanc_od(:),W(:),Z(:,:),work(:)
      real(8),allocatable :: hd(:),hod(:)
!     real(8),allocatable ::hreal(:,:),himg(:,:)
      real(8) :: a,b,norm,norm_1,deltat,time,prob1,perp,perp1,dcabs,temp2
      integer :: nsite,istates,ielem,ksteps,isteps,junk,junk2,nbonds
      integer :: info,i,j,k,i2,fact,t,tmax,j1,isite,jsite,ksite,nsiteby2
      integer,allocatable :: row_ptr(:),jhcol(:),iwsp(:),idg(:),nocc(:,:),jcorr(:)
      complex(8)::temp1,sum1,prob_c,prd3,scorr,afcorr,bndcorr
      real(8) :: prob2,sum2,sigma,summ,beta1,tol,q,ajunk, ipr
      complex(8),allocatable ::evec(:),evec_1(:),hmat(:),phinp1(:),zlanc(:,:)
      complex(8),allocatable ::temp(:,:),tzlanc(:,:),hdcomplex(:),prd(:,:)
      complex(8),allocatable :: U(:,:),psi_0(:),gs(:),psizero(:),psit(:),evec_init(:)
      character(len=1)::jobz
      integer :: ideg,lwsp,liwsp,iexp,iflag,ns,nav, fl
      complex(8),allocatable ::evec_2(:),wspc(:),amat(:,:),bndord(:),bndcorr2(:)
      real(8) :: entropy,szcorr,chcorr,temp3,temp4,chfluc,szfluc
      real(8),allocatable ::spin(:,:),szz(:),chrd(:),szcorr2(:),chcorr2(:),docc(:)
      character(len=20), external :: str
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      interface
       subroutine sprmul(istates,ielem,hmat,col_ind,row_ptr,b,c)
        implicit none
        complex(8) ::c(:),hmat(:),b(:)
        integer :: row_ptr(:),col_ind(:),istates,ielem
        integer :: n,j,ind1,ind2,diff
      end subroutine sprmul
      end interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      interface
      subroutine conduct(istates,ielem,hmat,jhcol,row_ptr,ksteps,hlanc1,gs,sigma)   
      implicit none
      integer :: i,j,i1,j1,nbonds,nbnd,norb,nlctrn,ifound
      integer :: id,idup,jphase,istates,junk,ires,ielem
      integer :: ksteps,k,inb,ihop,jcount,jhcol(:),row_ptr(:)
      integer,allocatable :: ibond(:,:),iorb(:),jorb(:),idgseq(:)
      complex(8),allocatable :: opr(:),R(:),d1(:),phinp1(:),evec(:),evec_1(:)
      real(8) :: a,b,perp,norm
      complex(8)::gs(:),konst,temp1,temp2,green,opr_sum,jc,hmat(:)
      real :: ajunk,aphase,crphase,w,eps,sgn
      real(8) :: sigma,hlanc1,trnsfr,phs
      real(8),allocatable :: e1(:),hd(:),hod(:)
      end subroutine conduct
      end interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      interface
      subroutine bond_order(norb,istates,coef,nbonds)
      implicit none
      integer :: i,j,i1,j1,nbonds,nbnd,norb,nlctrn,ifound,junk2
      integer :: id,idup,jphase,istates,junk,ires,ielem,info
      integer :: jnb,ihop1,ihop2,iddnup,iddndn,idupdn,idupup,indcor,ires1,ires2
      integer :: ksteps,k,inb,ihop,jcount,isite,ksite
      integer,allocatable ::ibond(:,:),iorb(:),jorb(:),korb(:), idgseq(:),jcorr(:)
      complex(8),allocatable ::R(:),d1(:),phinp1(:),evec(:),evec_1(:),evec2(:),bndtemp(:,:)
      real(8) :: a,b,perp,norm
      complex(8)::coef(:),konst,temp1,temp2,green,opr_sum,jc,opr,bndcorr
      real :: ajunk,aphase,crphase,w,eps,sgn
      real(8) :: sigma,hlanc1,trnsfr,phs,phs1,phs2,cormat
      real(8),allocatable :: e1(:),hd(:),hod(:),Z(:,:),work(:)
      character(len=1)::jobz
      real, allocatable :: bndcor(:)
      end subroutine bond_order
      end interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     interface
!      subroutine entanglement(nst,c,entropy)
!      implicit none
!      integer :: nsite,nst,nby2,nstby2,i,j,k,l
!      integer,allocatable :: indexx(:,:)
!      real(8),allocatable :: rho_r(:,:),rho_i(:,:),zr(:,:),zi(:,:),eig(:)
!      real(8),allocatable :: fv1(:),fv2(:),fm1(:,:)
!      complex(8) :: tempc,c(:)
!      real(8) :: entropy
!      integer :: matz,ierr
!      end subroutine entanglement
!     end interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      interface
      subroutine fecorr(nsite,istates,idg,evec,scorr,afcorr)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension :: ksite(16),idg(:)
      complex(8) :: evec(:)
      complex(8) :: scorr,afcorr
      end subroutine fecorr
      end interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      open(23,file='gen.out',form='unformatted')
      open(24,file = 'sztot.out',form='unformatted')
      open(34,file = 'chtot.out',form='unformatted')
!     open(unit=3,file='matpppq.out',form='unformatted')
      open(unit=3,file='matppp.out')
!     open(unit=12,file='unsymdet.out',form='unformatted')
      open(unit=11,file='initial.inp')
      open(unit=12,file='unsymdet.out')
      open(unit=17,file='evec_init.out')
      ksteps = 100
      

      read(3,*)istates,ielem
      allocate(evec(istates))
      read(12,*)ajunk,junk2
      read(12,*)(evec(i),i=1,istates)
      close (12)
      read(11,*) junk
      read(11,*) ajunk,nbonds
      close(11)

!     call density_index(istates)
      allocate(hmat(ielem))
      allocate(jhcol(ielem))
      allocate(row_ptr(istates+1))
      allocate(hlanc_d(ksteps))
      allocate(hlanc_od(ksteps-1))
      allocate(idg(istates))

      read(23)nsite,junk
      read(23)(idg(i),i=1,istates)
      close(23)
      allocate(spin(istates,nsite))
      allocate(nocc(istates,nsite))
      allocate(szz(nsite))
      allocate(docc(nsite))
      allocate(chrd(nsite))
      allocate(szcorr2(nsite))
      allocate(chcorr2(nsite))
      allocate(bndord(nbonds))
      allocate(bndcorr2(nbonds))
      allocate(jcorr(nsite))
      read(24)((spin(i,isite),isite = 1,nsite),i=1,istates)
      read(34)((nocc(i,isite),isite = 1,nsite),i=1,istates)
      close(24)
      close(34)
      fact = 1
!     do i = 2,ksteps
!      fact = fact*i
!     end do
!      fact = fact*(1.0d-3) 
      read(3,*) (row_ptr(k),k=1,istates)
      read(3,*)(hmat(k),k=1,ielem)
      read(3,*) (jhcol(k),k=1,ielem)
      row_ptr(istates+1)=ielem+1
      row_ptr=row_ptr-1

      jobz = "V"                 ! Lapack diagonal routine requires
      open(71,file='maxsteps.inp')
      read(71,*)tmax,deltat
      close(71)

!    Create an initial state vector
     allocate(evec_1(istates))
     allocate(evec_init(istates))  !array for storing initial vector
!      write(*,*)"HELLO4"
       allocate(zlanc(istates,ksteps))
!     open(44,file='time_steps')
      open(54,file='prob.dat',position='append')
      open(57,file='szcorr.dat',position='append')
      open(58,file='chcorr.dat',position='append')
      open(75,file='fecorr.dat',position='append')
      open(77,file='afcorr.dat')
      open(76,file='spinden.dat',position='append')
      open(86,file='charden.dat',position='append')
      allocate(psit(ksteps))
      allocate(psizero(ksteps))
      psizero = (0.0,0.0)
      psizero(1) =(1.0,0.0) 
      allocate(hd(ksteps))
      allocate(hod(ksteps-1))
      allocate(gs(istates))
      open(15,file='eigvals.out')
       tol = 1.0d-7
       time = 0.0
!      evec_init = evec    ! Store the initial vector in this array
!      read(17,*)ajunk,junk2
!      read(17,*)(evec_init(i),i=1,istates)

!       deltat = 0.01
       do t = 1, tmax       !  Loop for Time Steps
!      if(mod(t,10).eq.0) write(*,*)t
       norm_1 = 1.0
       evec_1 = (0.0,0.0)
       b = 0.0
!--------MAKE THE TRIDIAGONAL MATRIX --------------
       do,k=1,istates
        zlanc(k,1)=evec(k)
       end do
      beta1 = 1.0
      allocate(phinp1(istates))
      do j = 1, ksteps  !! (TILL CONVERGENCE)
       norm = 0.0
       a = 0.0
        call sprmul(istates,ielem,hmat,jhcol,row_ptr,evec,phinp1)
        do i=1,istates
         a = a + dconjg(evec(i))*phinp1(i)
         norm = norm + dconjg(evec(i))*evec(i)
        end do
         write(45,*)"kstep=",j,"a=",a,"norm=",norm
!        a = a/norm
!       b = norm/norm_1
!       hreal(j,j)=a
        hlanc_d(j)=a
        if(j.eq.ksteps) exit
        perp = 0.0
        do k= 1,istates
         phinp1(k)= phinp1(k) - a*evec(k) - b*evec_1(k)
         perp = perp + (dconjg(phinp1(k)))*phinp1(k)
        end do
!~~~~~~~~~Reorthogonalization~~~~~~~~~~~~~
!       do j1 = 1,j
!        q = 0.0
!        do k=1,istates
!         q = q + dconjg(zlanc(k,j1))*phinp1(k)
!        end do
!        phinp1 = (phinp1 - q*zlanc(:,j1))
!        phinp1 = phinp1/(1-q**2)
!       end do
!==========================================
!        b = perp/norm
         b = sqrt(perp)
         if(b.le.tol) exit
!        phinp1 = (sqrt(b))*phinp1
         phinp1 = phinp1/(dcmplx(b,0.0))
!       hlanc_od(j)=sqrt(b)
        hlanc_od(j)=b
        beta1 = beta1*hlanc_od(j)
        do k=1,istates
        zlanc(k,j+1) = phinp1(k)   ! Storing the Krylov vectors
        end do
        evec_1 = evec
        evec = phinp1
         perp = 0.0
        norm_1= norm
      end do
!-------- Tridiagonal matrix constructed ------------

      deallocate(phinp1)
!     deallocate(hmat)
      allocate(Z(ksteps,ksteps))
      allocate(work(20*ksteps))
      hd = hlanc_d
      hod = hlanc_od
!     write(*,*)"ksteps =",ksteps
!     CALL RST(ksteps, ksteps, hlanc_d, hlanc_od, MATZ, Z, IERR)
!     call ch(ksteps,ksteps,hreal,himg,W,MATZ,zr,zi,fv1,fv2,fm1,IERR)
       CALL dstev(jobz,ksteps,hlanc_d,hlanc_od,Z,ksteps,work,info)
       deallocate(work)
!      write(*,*)"INFO =",info
!     if(t.lt.338) go to 55
!     do i = 1,2
!     write(*,*)hlanc_d(i)
!     end do
!     write(15,*)
!******* Calculate Spin density and Correlation *********** 
!     read(24)((spin(i,isite),isite = 1,nsite),i=1,idim)
!      write(37,*)((spin(i,isite),isite = 1,nsite),i=1,idim)
!     do isite = 1,nsite
                !     szz(isite)=0.0
!     do i = 1, idim
!       call szcalc()
!      szz(isite) = szz(isite) + conjg(evec(i))*spin(i,isite)*evec(i)
!      scorr = scorr + conjg(evec(i))*spin(i,isite)*spin(i,isite)*evec(i)
!     end do
!     end do
!

!     stop
!---------------------------------------
  55  if(t.gt.1) go to 25
      gs = (0.0,0.0)
      summ = 0.0
      do i = 1, istates
       do j = 1, ksteps
        gs(i) = gs(i) + zlanc(i,j)*Z(j,1)
       end do
       summ = summ + (gs(i))*(dconjg(gs(i)))
      end do
       gs = gs/(sqrt(summ))

      ! deltat = ((1.0d-3)*fact)/beta1
      ! deltat = deltat**(1.0/real(ksteps))
      ! write(44,*) "DELTA = ",deltat
!     open(43,file='ground_st.inp')
      open(unit=12,file='unsymdet.out')
      write(12,*)hlanc_d(1),junk2
      write(12,*)(gs(i),i=1,istates)
      write(*,*)"<<<<<<<<<<<", hlanc_d(1),">>>>>>>>>>>>>>>"
!     close(43)
      stop
        szcorr = 0.0
        chcorr = 0.0
        szfluc = 0.0
        chfluc = 0.0
        szcorr2 = 0.0
        chcorr2 = 0.0
        ipr = 0.0
       do isite = 1,nsite
        szz(isite)=0.0
        chrd(isite)=0.0
        do i = 1, istates
         ipr = ipr + zabs(gs(i)*gs(i)*gs(i)*gs(i))
         szz(isite) = szz(isite) + conjg(gs(i))*spin(i,isite)*gs(i)
         chrd(isite) = chrd(isite) + conjg(gs(i))*nocc(i,isite)*gs(i)
          !docc(isite) = docc(isite) +conjg(gs(i))*nocc(i,isite)*gs(i)
!        jsite = mod(isite+1,nsite) !  For closed ring (PBC)
         do ksite = 1,nsite
          temp3 = conjg(gs(i))*spin(i,isite)*spin(i,ksite)*gs(i)
          szcorr2(abs(isite-ksite)+1)=szcorr2(abs(isite-ksite)+1) + temp3
          temp4 =conjg(gs(i))*(nocc(i,isite)-nav)*(nocc(i,ksite)-nav)*gs(i)
          chcorr2(abs(isite-ksite)+1)=chcorr2(abs(isite-ksite)+1) + temp4
          if((isite.le.nsiteby2).and.(ksite.le.nsiteby2)) then
           szfluc = szfluc + conjg(gs(i))*spin(i,isite)*spin(i,ksite)*gs(i)
           chfluc = chfluc + conjg(gs(i))*nocc(i,isite)*nocc(i,ksite)*gs(i)
          end if
         end do
         jsite = isite+1
         if(isite.eq.nsite) jsite=1
         !if(isite.eq.nsite) cycle
         szcorr = szcorr +conjg(gs(i))*spin(i,isite)*spin(i,jsite)*gs(i)
         chcorr = chcorr +conjg(gs(i))*nocc(i,isite)*nocc(i,jsite)*gs(i)
        end do
        szfluc = szfluc - (sum(szz))**2
        chfluc = chfluc - (sum(chrd))**2

        ipr = 1.0/ipr
          write(76,*)isite,szz(isite)
          write(86,*)isite,chrd(isite)
       end do
       jcorr = 0.0
       do i = 1, nsite
        do j = 1, nsite
         jcorr(abs(i-j)+1)=jcorr(abs(i-j)+1)+1
        end do
       end do
       write(*,*)"HELLO"
       call bond_order(nsite,istates,gs,nbonds)

        do i = 0,nsite-1
         write(96,*)i,szcorr2(i+1)/real(jcorr(i+1))
         write(97,*)i,chcorr2(i+1)/real(jcorr(i+1))
        end do
        do i = 0,nbonds-1
         write(99,*)i+1,abs(bndord(i+1))
         write(98,*)i,abs(bndcorr2(i+1))
        end do

      stop


! 25  allocate(amat(ksteps,ksteps))
!=================================================
! Time Evolution
!=================================================
  25  allocate(hdcomplex(ksteps))
      do i = 1,ksteps
       temp1 = (-1.0)*deltat*(dcmplx(0.0,hlanc_d(i)))
       hdcomplex(i) = exp(temp1)
!      temp2 = deltat*(hlanc_d(i))
!      hdcomplex(i) = dcmplx(cos(temp2),(-1)*sin(temp2))
      end do           

!       allocate(prd(ksteps,ksteps))
!      allocate(prd2(ksteps,ksteps))
      do j = 1,ksteps
       prd3 = (0.0,0.0)
       do k = 1,ksteps
        prd3 = prd3 + Z(j,k)*hdcomplex(k)*Z(1,k)
       end do
       do i = 1,istates
        zlanc(i,j) = prd3*zlanc(i,j)
       end do
      end do

      prob_c = (0.0,0.0)
      perp = 0.0
      time = time + deltat
      evec = (0.0,0.0)

       do i = 1,istates
        do k = 1,ksteps
         evec(i) = evec(i) + zlanc(i,k)
        end do
         perp = perp + (dconjg(evec(i)))*evec(i)
         prob_c = prob_c + dconjg(evec_init(i))*evec(i)   ! Overlap with Initial State
       end do
      
       evec = evec/(sqrt(perp))
       prob1 =(dconjg(prob_c)*prob_c)/perp
!******* Calculate Spin density and Correlation *********** 
!     read(24)((spin(i,isite),isite = 1,nsite),i=1,idim)
!      write(37,*)((spin(i,isite),isite = 1,nsite),i=1,idim)
!     if(mod(tmax,500).eq.0) then
        
        szcorr = 0.0
        chcorr = 0.0
       do isite = 1,nsite
        szz(isite)=0.0
        chrd(isite)=0.0
        do i = 1, istates
         szz(isite) = szz(isite) + conjg(evec(i))*spin(i,isite)*evec(i)
         chrd(isite) = chrd(isite) + conjg(evec(i))*nocc(i,isite)*evec(i)
         jsite = mod(isite+1,nsite)
         if(isite.eq.nsite-1) jsite=nsite
         szcorr = szcorr + conjg(evec(i))*spin(i,isite)*spin(i,jsite)* evec(i)
         chcorr = chcorr + conjg(evec(i))*nocc(i,isite)*nocc(i,jsite)* evec(i)
        end do
!         write(76,*)isite,szz(isite)
!         write(86,*)isite,chrd(isite)
       end do
!       write(76,*)
!     end if
       chcorr = chcorr/real(nsite)
       szcorr = szcorr/real(nsite)
!       call fecorr(nsite,istates,idg,evec,scorr,afcorr)
!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
!    Entanglement Entropy (von Neumann)
!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
!      call entanglement(istates,evec,entropy)
!      call entanglement(norb,istates,idgseq,evec,entropy)
!      prob1 =(dconjg(prob_c)*prob_c)
      write(54,*)time,prob1
      write(57,*)time,szcorr
      write(58,*)time,chcorr
!     write(*,*)time,prob1
!     write(64,*)time,entropy
      write(75,*)time,realpart(scorr)
      write(77,*)time,realpart(afcorr)
      open(unit=12,file='unsymdet.out')
      write(12,*)ajunk,junk2
      write(12,*)(evec(i),i=1,istates)
      close(12)
       psizero = psit
      deallocate(Z)
      deallocate(hdcomplex)
!     deallocate(prd)
!     deallocate(amat)
!     deallocate(U)
!      open(89, file = 'time_steps')
!       write(89, *)t, deltat
!      close(89)
      end do      ! Time Steps
!     stop
!     call conduct(istates,ielem,hmat,jhcol,row_ptr,ksteps,hlanc_d(1),gs,sigma)
      end
!,,,,,,,,,,,,,,,,,,,,,,,
! Sparse Multiplication
!'''''''''''''''''''''''
      subroutine sprmul(istates,ielem,hmat,col_ind,row_ptr,b,c)
!                sprmul(istates,ielem,hmat,jhcol,row_ptr,evec,phinp1)
      implicit none
       complex(8) ::c(:),hmat(:),b(:)
       integer :: row_ptr(:),col_ind(:),istates,ielem
       integer :: n,j,ind1,ind2,diff
!      j = j+1
       c = (0.0,0.0)
       do n = 1,istates
        j =1
        diff = row_ptr(n+1)-row_ptr(n)
        do while (j.le.diff)
!       do while (j.le.row_ptr(n+1))
         ind1 = row_ptr(n)+j
         ind2 = col_ind(ind1)
         c(n) = c(n) + hmat(ind1)*b(ind2)
         j = j+1
        end do
       end do
      end subroutine sprmul
!************************************************************
       character(len=20) function str(k)
       integer, intent(in) :: k
       write (str, 30) k
       str = adjustl(str)
   30  format(I2.2)
       end function str


