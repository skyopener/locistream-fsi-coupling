! // read grid and comainNumNodesectivitiy inforamtion and return them
      subroutine readgrid
      use maindata

      implicit none
!      integer :: mainNumNodes, mainNumElems
!      real(8),allocatable:: mainXgl(:),mainYgl(:),mainZgl(:)
!      integer, allocatable:: mainConnect(:,:)

     !Local variables
      integer:: io,np,node,eletype,num

      open(1020,file='trimesh.dat',action='read')

      np=0
      do
       read(1020,*,iostat=io)
       if (io < 0) exit
       np = np + 1
      enddo

      rewind(1020)
      mainNumNodes = np

      allocate(mainXgl(mainNumNodes))
      allocate(mainYgl(mainNumNodes))
      allocate(mainZgl(mainNumNodes))
!      allocate(uf(mainNumNodes,3))
 
      do num=1,mainNumNodes
      read(1020,*) node,mainXgl(num),mainYgl(num),mainZgl(num)
!        uf(num,1) = mainXgl(num)
!        uf(num,2) = mainYgl(num)
!        uf(num,3) = mainZgl(num)
      enddo

      open(1021,file='connect.dat',action='read')
      np=0

      do
       read(1021,*,iostat=io)
       if (io < 0) exit
       np = np + 1
      enddo

      rewind(1021)
      mainNumElems = np
      print *,'Number of nodes in the FE mesh is:',mainNumNodes
      print *,'Number of elements in the FE mesh is:',mainNumElems

!      read(*,*)
      allocate(mainConnect(mainNumElems,3));

      do num=1,mainNumElems
       read(1021,*)node,eletype,mainConnect(num,1),mainConnect(num,2),mainConnect(num,3)
!       print *, mainConnect(num,1),mainConnect(num,2),mainConnect(num,3)
      enddo
      rewind(1021)

!       read(*,*)
!       stop
      close(1020)
      close(1021)

      end subroutine readgrid

