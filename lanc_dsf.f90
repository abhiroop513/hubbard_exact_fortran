      implicit none
      real(8),allocatable ::hlanc_d(:),hlanc_od(:),W(:),Z(:,:),work(:)
      real(8),allocatable :: hd(:),hod(:)
!     real(8),allocatable ::hreal(:,:),himg(:,:)
      real(8) :: a,b,norm,norm_1,deltat,time,prob1,perp,perp1,dcabs
      integer :: nsite,istates,ielem,ksteps,isteps,junk,junk2,nbonds
      integer :: info,i,j,k,i2,fact,t,tmax,j1,isite,jsite,ksite,nsiteby2
      integer,allocatable :: row_ptr(:),jhcol(:),iwsp(:),idg(:),nocc(:,:),jcorr(:)
      complex(8)::temp1,sum1,prob_c,prd3,scorr,afcorr,bndcorr
      real(8) :: prob2,sum2,sigma,summ,beta1,tol,q,ajunk, ipr
      complex(8),allocatable ::evec(:),evec_1(:),hmat(:),phinp1(:),zlanc(:,:),d1(:),R(:)
      complex(8),allocatable ::temp(:,:),tzlanc(:,:),hdcomplex(:),prd(:,:)
      complex(8),allocatable :: U(:,:),psi_0(:),gs(:),psizero(:),psit(:),evec_init(:)
      character(len=1)::jobz
      integer :: ideg,lwsp,liwsp,iexp,iflag,ns,nav, fl, io, jcount
      complex(8),allocatable ::evec_2(:),wspc(:),amat(:,:),bndord(:),bndcorr2(:), sq_gs(:)
      real(8) :: entropy,szcorr,chcorr,temp3,temp4,chfluc,szfluc
      real(8),allocatable ::spin(:,:),szz(:),chrd(:),docc(:), e1(:), rvec(:,:), qvec(:)
      character(len=20), external :: str
      real(8):: omega, eps, gs_energy
      complex(8):: konst,temp2,green1, green2, opr_sum,jc,opr, icom
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
      ksteps = 50
      

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
      allocate(jcorr(nsite))
      allocate(sq_gs(istates))
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
      !open(54,file='prob.dat',position='append')
      !open(57,file='szcorr.dat',position='append')
      !open(58,file='chcorr.dat',position='append')
      !open(75,file='fecorr.dat',position='append')
      !open(77,file='afcorr.dat')
      !open(76,file='spinden.dat',position='append')
      !open(86,file='charden.dat',position='append')
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
!       do t = 1, tmax       !  Loop for Time Steps
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

!---------------------------------------
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
      gs_energy = hlanc_d(1)
!     close(43)
      write(*,*) "Diagonalization Done"      
      open(57,file='coords.inp')
      allocate(rvec(nsite,2))

      ! Read arrays from file
      do isite = 1, nsite
       read(57, *) rvec(isite, :)
      end do
      
      allocate(qvec(2))       
      !qvec = 3.141592653589*(/0.66667, 0.0/)
      qvec = 3.141592653589*(/0.5, 0.0/)

!-------- Initialize vector for Dynamical sructure factor ----------------
! ----   S_q |Ψ0> = Σ exp(-iq.r)Sz_i |Ψ0> ---------------------
      icom = (0.0, 1.0)
      sq_gs = 0.0
      do i = 1,istates
        do isite = 1,nsite
         sq_gs(i) = sq_gs(i) +  exp(icom*dot_product(qvec,rvec(isite,:)))*gs(i)*spin(i,isite)
        end do
      end do  
      !sq_gs = sq_gs/norm2(sq_gs)

!------- Second Lanczos Step ---------------------------------------------      
       norm_1 = 1.0
       evec_1 = (0.0,0.0)

         
       evec = sq_gs/sqrt(abs(dot_product(sq_gs, sq_gs)))
       b = 0.0
!--------MAKE THE TRIDIAGONAL MATRIX --------------
      zlanc(:,1)=evec
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
        do j1 = 1,j
         q = 0.0
         do k=1,istates
          q = q + dconjg(zlanc(k,j1))*phinp1(k)
         end do
         phinp1 = (phinp1 - q*zlanc(:,j1))
         phinp1 = phinp1/(1-q**2)
        end do
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
      !allocate(Z(ks3.141592653589teps,ksteps))
      !allocate(work(20*ksteps))
      hd = hlanc_d
      hod = hlanc_od
!----------------------------------------------------
!      gs = (0.0,0.0)
!      summ = 0.0
!      do i = 1, istates
!       do j = 1, ksteps
!        gs(i) = gs(i) + zlanc(i,j)*Z(j,1)
!       end do
!       summ = summ + (gs(i))*(dconjg(gs(i)))
!      end do
!      gs = gs/(sqrt(summ))
!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
!    DYNAMICAL STRUCTURE FACTOR
!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

!========  Perfrom the Continued Fraction =================
      allocate(R(ksteps))
      allocate(e1(ksteps))
      allocate(d1(ksteps))
      omega = -8.0
      eps = 0.05
      jcount = 0
      open(38,file='dsf_vs_omega.dat')
      do k = 1,2400
        konst = dcmplx(omega,0.0) + dcmplx(0.0,eps) + dcmplx(gs_energy,0.0)
        R(ksteps) = (0.0,0.0)
        e1(ksteps) = 0.0
        d1(ksteps) = konst - dcmplx(hd(ksteps),0.0)
        do i=ksteps-1,1,-1
          d1(i) = konst - dcmplx(hd(i),0)  ! hlanc1 = hlanc_d(1)
          e1(i) = (hod(i))**2
          temp2 =d1(i+1)-R(i+1)
          R(i) = e1(i)/(temp2)
        end do
  
        temp1 = d1(1) - R(1)
        !write(*,*)"Shape TEMP1", shape(temp1)
        
        green1 = (dot_product(sq_gs, sq_gs))*(dcmplx(1.0,0.0))/temp1
  
        !sigma = imagpart(green*opr_sum)
        !sigma = sigma/(omega*(3.141592653589))
  !     write(38,*)w,sigma
  !     write(38,*)w,(imagpart(green)*(-1.0/3.141592653589))
      !end do
      !write(*,*)'jcount=',jcount
      !write(*,*)"COMPLETE"

! ========  Perform the Second Continued Fraction ===============
      !jcount = 0
      !open(38,file='dsf.dat')
      
      !do k = 1,800
        konst = dcmplx(omega,0.0) + dcmplx(0.0,eps) - dcmplx(gs_energy,0.0)
        R(ksteps) = (0.0,0.0)
        e1(ksteps) = 0.0
        d1(ksteps) = konst + dcmplx(hd(ksteps),0.0)
        do i=ksteps-1,1,-1
          d1(i) = konst + dcmplx(hd(i),0)  ! hlanc1 = hlanc_d(1)
          e1(i) = (hod(i))**2
          temp2 =d1(i+1)+R(i+1)
          R(i) = e1(i)/(temp2)
        end do

        temp1 = d1(1) + R(1)
        !green = (dcmplx(1.0,0.0))/temp1
        green2 = (dot_product(sq_gs, sq_gs))*(dcmplx(1.0,0.0))/temp1

  !     write(38,*)w,sigma
  !     write(38,*)w,(imagpart(green)*(-1.0/3.141592653589))
  !     write(38,*) omega, imagpart(green1 - green2)/real(nsite)
        write(38,*) omega, -imagpart(green1)/real(nsite)
        omega = omega + 0.02
        jcount = jcount + 1
      end do
      write(*,*)'jcount=',jcount
      write(*,*)"COMPLETE both continued fractions"

!       call fecorr(nsite,istates,idg,evec,scorr,afcorr)
!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
!    Entanglement Entropy (von Neumann)
!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
!      call entanglement(istates,evec,entropy)
!      call entanglement(norb,istates,idgseq,evec,entropy)
!      prob1 =(dconjg(prob_c)*prob_c)

!      write(54,*)time,prob1
!      write(57,*)time,szcorr
!      write(58,*)time,chcorr
!!     write(*,*)time,prob1
!!     write(64,*)time,entropy
!      write(75,*)time,realpart(scorr)
!      write(77,*)time,realpart(afcorr)
!      open(unit=12,file='unsymdet.out')
!      write(12,*)ajunk,junk2
!      write(12,*)(evec(i),i=1,istates)
!      close(12)
!       psizero = psit
!      deallocate(Z)
!      deallocate(hdcomplex)
!!     deallocate(prd)
!!     deallocate(amat)
!!     deallocate(U)
!!      open(89, file = 'time_steps')
!!       write(89, *)t, deltat
!!      close(89)
!      end do      ! Time Steps
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


