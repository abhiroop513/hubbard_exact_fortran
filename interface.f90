      module interfaces
      interface
      subroutine decomp(id,norb,iorb)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension iorb(:)
      end subroutine decomp
      end interface
!***********************************************
      interface
      subroutine trnsfrdn(iorb,inb,jorb,ires,norb,ibond,ihop,phs)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension :: iorb(:),ibond(:,:),jorb(:),korb(20)
      end subroutine trnsfrdn
      end interface
!******************************************************
      interface
      subroutine trnsfrup(iorb,inb,jorb,ires,norb,ibond,ihop,phs)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension :: iorb(:),ibond(:,:),jorb(:),korb(20)
      end subroutine trnsfrup
      end interface
!********************************************************
      interface
      subroutine convrt(jorb,norb,id)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension jorb(:)
      end subroutine convrt
      end interface
!******************************************************
      interface
      subroutine binsrch(id1,ifound,ibin,nto2,idet,islt)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      dimension :: idet(:),nto2(:)
      end subroutine binsrch
      end interface
!***************************************************
      interface
      subroutine order(cotemp,itdg,conew,itnew,itot,ielem)
      implicit real(8)(a-h,o-z)
      implicit integer(4)(i-n)
      complex(8) :: atrns,cotemp,conew
      dimension :: itdg(:),cotemp(:),conew(:),itnew(:)
      end subroutine order
      end interface
!*************************************************

      end module interfaces
