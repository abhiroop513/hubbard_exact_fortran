      implicit none
      integer :: nsite,istates,i,n,a
      complex(8),allocatable :: idg(:)
      real(8) :: norm,ran,ran2,junk
      open(2,file='gen.out',form='unformatted')
      read(2)nsite,istates
      junk = 50.0
      allocate(idg(istates))
      idg = (0.0,0.0)
      norm = 0.0
!     n = nsite/2
      do i = 1,istates
       call random_number(ran)
       call random_number(ran2)
       idg(i) = ran+ran2
       norm = norm + (idg(i))**2
      end do
      idg = idg/(sqrt(norm))
!     open(1,file='unsymdet.out',form='unformatted')
      open(1,file='unsymdet.out')
      write(1,*)junk,istates
      write(1,*)(idg(i),i=1,istates)
      write(*,*)nsite,istates
      end
         
