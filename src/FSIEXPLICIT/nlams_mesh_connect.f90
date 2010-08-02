!Read triangle mesh and connectivity from two different files trimesh.dat and connect.dat.

module mesh_connect_vars
 real(8),save,dimension(:),allocatable::Xgl,Ygl,Zgl
! real(8),save,dimension(:,:),allocatable::UndefGrid,CoordGridsUnd
 real(8),save,dimension(:,:),allocatable::UndefGrid
 integer,save,dimension(:,:),allocatable::Connect
 integer,save ::NumNodes
 integer,save ::NumElems
 integer,save ::NumNodesPerEl,NumDofPerNode,NumElDof    !Number of nodes per element, number of DOF per node
end module mesh_connect_vars

subroutine nlams_mesh_connect
 use mesh_connect_vars
 use LociSTREAMData, only:LociSTREAMNumNodes,LociSTREAMXgl,LociSTREAMYgl,LociSTREAMZgl,LociSTREAMNumElems,LociSTREAMConnect

 implicit none
 integer:: i

 NumNodes = LociSTREAMNumNodes

! print *, 'NumNodes :',NumNodes

 allocate(Xgl(NumNodes)); Xgl = 0.0 ;
 allocate(Ygl(NumNodes)); Ygl = 0.0 ;
 allocate(Zgl(NumNodes)); Zgl = 0.0 ;
 allocate(UndefGrid(NumNodes,3)); UndefGrid = 0.0 ;
! allocate(CoordGridsUnd(NumNodes,3)); CoordGridsUnd = 0.0 ; ! might not need
 
! do i=1,NumNodes
!   print *, 'LociSTREAMXgl: ',LociSTREAMXgl(i)
! enddo


 
 Xgl = LociSTREAMXgl
 Ygl = LociSTREAMYgl
 Zgl = LociSTREAMZgl
 
 do i=1,NumNodes
  UndefGrid(i,1) = Xgl(i)
  UndefGrid(i,2) = Ygl(i)
  UndefGrid(i,3) = Zgl(i)
 enddo

 NumElems = LociSTREAMNumElems

 allocate(Connect(NumElems,3)); Connect = 0
 
 do i=1,NumElems
  Connect(i,1) = LociSTREAMConnect(i,1)
  Connect(i,2) = LociSTREAMConnect(i,2)
  Connect(i,3) = LociSTREAMConnect(i,3)
 enddo

 NumNodesPerEl = 3
 NumDofPerNode = 6
 NumEldof      = NumNodesPerEl*NumDofPerNode


!EleFlag = 0 !Triangle element
! Input element coordinates
!  open(1020,file='trimesh.dat',action='read')
!  open(1021,file='connect.dat',action='read')
!  open(23033,file='temp.dat',action='write')
! Grid coordinates
!  np=0
!  do
!    read(1020,*,iostat=io)
!      if (io < 0) exit
!      np = np + 1
!  enddo
!  rewind(1020)
!  NumNodes = np
!  allocate(Xgl(NumNodes));
!  allocate(Ygl(NumNodes));
!  allocate(Zgl(NumNodes));
!  allocate(UndefGrid(NumNodes,3)); 
!  allocate(CoordGridsUnd(NumNodes,3));
!trimesh.dat is the grid file
!  do num=1,NumNodes
!    read(1020,*)node1,Xgl(num),Ygl(num),Zgl(num)
!    print *, num,Xgl(num),Ygl(num)
!Reading the undeformed nodal coordinates in a single array
!    UndefGrid(num,1) = Xgl(num)
!    UndefGrid(num,2) = Ygl(num)
!    UndefGrid(num,3) = Zgl(num)
!  enddo
!For coupling with STREAM
! CoordGridsUnd = UndefGrid
! Connectivity
!  np=0
!  do
!    read(1021,*,iostat=io)
!      if (io < 0) exit
!      np = np + 1
!  enddo
!
!  rewind(1021)
!  NumElems = np
!  print *,'Number of nodes in the FE mesh is:',NumNodes
!  print *,'Number of elements in the FE mesh is:',NumElems
!
!  allocate(Connect(NumElems,3));
!
!!Reading connectivity from connect.dat
!  do num=1,NumElems
!    read(1021,*)node,eletype,Connect(num,1),Connect(num,2),Connect(num,3)
!!    print *, num, Connect(num,1),Connect(num,2),Connect(num,3)
!  enddo
!  rewind(1021)
!  close(1020)
!  close(1021)
!  
!!Element type - EleFlag indicates whether a triangle element or a rectangular element is being used. For triangle, EleFlag is 0.
!  !Only triangle is available right now.
!  if (EleFlag == 0) then
!     NumNodesPerEl = 3
!     NumDofPerNode = 6
!     NumEldof      = NumNodesPerEl*NumDofPerNode
!  else
!     NumNodesPerEl = 4
!     NumDofPerNode = 6
!     NumEldof      = NumNodesPerEl*NumDofPerNode
!  endif

end subroutine nlams_mesh_connect

