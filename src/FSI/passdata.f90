      subroutine passdata(A,B,C,D,E)
      implicit none

      integer:: A
      real(8) :: B(A),C(A),D(A),E(A)

      integer :: i
      print *, "NumNodes :", A

!      do i =1, A
!       print *, B(i), C(i), D(i)
!      enddo

      E=B

      return
      end subroutine passdata
      subroutine passpara(A, B)

      implicit none
       integer :: A
       real(8) :: B

       print *, "Number A and B are ", A, B

       end subroutine passpara

      subroutine passintdata(A,B)

      implicit none

      integer:: A
      integer:: B(A,3)

      integer:: i

      do i=1,A
       print *, B(i,1:3)
      enddo
      end subroutine passintdata
